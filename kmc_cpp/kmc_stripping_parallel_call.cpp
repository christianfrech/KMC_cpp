//partition down 3d array
#include <mpi.h>
#include "kmc_stripping_parallel.hpp"

int main() {
    std::cout << "first line \n";
    int nprocs = 4;//argv[0];
    int xprocs;
    int yprocs;
    int zprocs;
    int size; 
    int rank;
    int gcf;    
    std::vector<int> dims = {0,0,0};

    std::ostringstream ss;
    std::string output_filename;
    std::string alltimes_filename;
    std::string times_filename;
    std::string count_filename;
    std::ofstream out_file;
    std::fstream in_file;
    std::vector<std::string> lines;
    std::string line;
    int read_idx = 0;
    std::vector< std::string > infile_name = {"geo_small.txt"};
    std::vector<double> vertex_rates = {7.65954e11};
    std::vector<double> edge_rates = {1.23537e6};
    std::vector<std::string> folders = {"parallel_output"};
    std::vector< std::string > catalog_file = {"ratecatalog_nogb.txt"};
    std::string region_infile;
    std::vector< std::vector< std::vector<double> > > reg_rates = {{{8.28403e10, 8.28403e10}}}; //settings of kb = 8.6173e-5, T = 300, E = 0.05, prefactor = 5e12
    std::tuple< std::vector< std::vector< std::vector<int> > >, std::vector<int>, std::vector<double>, std::vector<double> > return_tuple;
    int iterations = 1;
    double time = 10;
    int waitint = 0;

    std::vector<Matrix<int>*> vacancies; std::vector<int> counts; std::vector<double> times; std::vector<double> all_times;
    

    std::cout << "variables initialized \n";

    for (int h=0; h<(int)folders.size(); h++) {
        for (int i=0; i<iterations; i++) {
            for (int k=0; k<(int)edge_rates.size(); k++) {

                in_file.open(infile_name[h]);
                std::cout << "opened file \n";

                if (in_file.is_open()) {
                    while ( getline (in_file,line) )
                    {
                        lines.push_back(line);
                    }
                    in_file.close();
                }
                // parsing first line to grab dimensions of lattice ###
                std::string dims_str = lines[read_idx]; //getting dimension line
                read_idx ++;
                std::vector<std::string> dims_toks = tokenizer(dims_str," "); 
                
                print_1Dvector_string(dims_toks);

                for (int idx=0; idx<3; idx++) {

                    dims[idx] = std::stoi(dims_toks[idx+1]); 
                }                

                if (nprocs == 0) {
                    std::cout << "ERROR: no command line arguments provided\n";
                    exit(0); 
                }
                else {
                    //gcf = 2;
                    gcf = find_gcf(nprocs, dims[0]); 
                }

                xprocs = gcf;
                yprocs = (int)(nprocs / gcf);

                std::cout << "starting mpi stuff \n";
                // specify number of procs using slurm
                MPI_Init(NULL, NULL);
                MPI_Comm_size(MPI_COMM_WORLD, &size);
                MPI_Comm_rank(MPI_COMM_WORLD, &rank);

                std::cout << "rank: " << rank << "\n";
                std::cout << "mpi initialized \n";
                // chunkify all needed variables
                int x_idx = rank % xprocs;
                int y_idx = floor(rank / xprocs);

                int x_chunk_start = (int)( dims[0]/xprocs * x_idx );
                int x_chunk_end = (int)( dims[0]/xprocs * (x_idx + 1) );
                int y_chunk_start = (int)( dims[1]/yprocs * y_idx );
                int y_chunk_end = (int)( dims[1]/yprocs * (y_idx + 1) );
                std::vector<int> procs = {xprocs, yprocs};

                std::vector< std::vector<int> > chunk_bounds = {{x_chunk_start, x_chunk_end},{y_chunk_start, y_chunk_end}, {0,dims[2]}};
                std::cout << "rank: " << rank <<"\n";
                print_2Dvector(chunk_bounds);
                //exit(0);
                std::chrono::system_clock::time_point start;

                std::cout << "-----------------------------------------\n";
                std::cout << "i: " << i << "\n"; 
                
                char hostname[256];
                gethostname(hostname, sizeof(hostname));
                fflush(stdout); 
                /*
                if (rank == 0) {
                    printf("PID %d on %s ready for attach for rank %d\n", getpid(), hostname, rank);
                    while (waitint==0) {
                        sleep(1);
                    }
                }
                */
                
                MPI_Barrier(MPI_COMM_WORLD);
                printf("past barrier");
                fflush(stdout);

                Lattice* lattice = populate_lattice(infile_name[h], catalog_file[h], region_infile, vertex_rates[k], edge_rates[k], reg_rates[h], dims, chunk_bounds, rank, procs);
                std::cout << "after populate_lattice " << std::endl;
                start = std::chrono::system_clock::now();
                std::cout << "time got \n";
                MPI_Barrier(MPI_COMM_WORLD);
                lattice_return_struct return_tuple = lattice->new_kmc_iterator(time, start, folders[h], i);
                delete lattice;
                
                vacancies = return_tuple.get_all_vacancies(); counts = return_tuple.get_move_counts(); times = return_tuple.get_time_count(); all_times =return_tuple.get_all_times();

                MPI_Barrier(MPI_COMM_WORLD); // all the sub ranks/processes waits here
                /* process 0 will aggregate the results*/
                int size1 = (int)vacancies.size();
                int size2 = (int)vacancies[0]->rows();
                int size3 = (int)vacancies[0]->cols();
                int size = (int)(vacancies.size() * vacancies[0]->rows() * vacancies[0]->cols());
                std::vector<int> output(size);
                int idx = 0;
                Matrix<int>* one_vac;

                for (int frames=0; frames<(int)vacancies.size(); frames++) {
                    for (int vacs=0; vacs<(int)vacancies[frames]->rows(); vacs++) {
                        for (int coord=0; coord<(int)vacancies[frames]->cols(); coord++) {
                            idx = (int)frames * size3 * size2 + vacs * size3 + coord;
                            one_vac = vacancies[frames];
                            output[idx] = (*one_vac)[vacs][coord];
                        }
                    }
                }

                int size_in;
                std::vector<int> vacs_in;
                std::vector< std::vector<int> > all_vacs_flattened;
                std::vector< std::vector< std::vector< std::vector<int> > > > all_vacs;
                FourDBoolArr vacancies_out(2, dims[0], dims[1], dims[2]);
                int i1; int i2; int i3; int i4;
                MPI_Status status;

                if(rank==0) {
                    for (int i=1; i<size; i++) {
                        x_idx = i % dims[0];
                        y_idx = floor(i / xprocs);
                        x_chunk_start = (int)( dims[0]/xprocs * x_idx );
                        y_chunk_start = (int)( dims[1]/yprocs * y_idx );
                    
                        //int MPI_Recv(void *buf, int count, MPI_Datatype datatype, int source, int tag,
                        //MPI_Comm comm, MPI_Status *status)
                        
                        MPI_Recv(
                                &size_in,      //Address of the message we are receiving.
                                1,                  //Number of elements handled by that address.
                                MPI_INT,            //MPI_TYPE of the message we are sending.
                                i,                  //Rank of sending process
                                1,                  //Message Tag
                                MPI_COMM_WORLD,      //MPI Communicator
                                &status
                            );
                        MPI_Recv(
                                &vacs_in,      //Address of the message we are receiving.
                                size_in,                  //Number of elements handled by that address.
                                MPI_INT,            //MPI_TYPE of the message we are sending.
                                i,                  //Rank of sending process
                                1,                  //Message Tag
                                MPI_COMM_WORLD,      //MPI Communicator
                                &status
                            );
                        //MPI_Irecv(&max_rate, 1, MPI_DOUBLE, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &request);
                        
                        printf("Process %d has received \n", i);
                        
                        for (int j=0; j<size1; j++) {
                            for (int k=0; k<size2; k++) {
                                for (int l=0; l<size3; l++) {
                                    idx = (int)(j*size2*size3 + k*size3 + l);
                                    if (l==0) i1 = vacs_in[idx];
                                    else if (l==1) i2 = vacs_in[idx] + x_chunk_start;
                                    else if (l==2) i3 = vacs_in[idx] + y_chunk_start;
                                    else if (l==3) i4 = vacs_in[idx];
                                    vacancies_out(i1,i2,i3,i4) = 1;                                 
                                }
                            }
                        } 
                    }
                }
                else {
                    MPI_Send(
                            &size,      //Address of the message we are receiving.
                            1,                  //Number of elements handled by that address.
                            MPI_INT,            //MPI_TYPE of the message we are sending.
                            rank,                  //Rank of sending process
                            1,                  //Message Tag
                            MPI_COMM_WORLD      //MPI Communicator
                        );
                    MPI_Send(
                            &output,      //Address of the message we are receiving.
                            size,                  //Number of elements handled by that address.
                            MPI_INT,            //MPI_TYPE of the message we are sending.
                            rank,                  //Rank of sending process
                            1,                  //Message Tag
                            MPI_COMM_WORLD      //MPI Communicator
                        );
                    printf("Process %d has sent \n", i);

                }
                      
                printf("\n");
                MPI_Finalize();

                Matrix<int>* all_vacancies = vacancies_out.nonzero();

                for (int l=0; l<all_vacancies->rows(); l++) {
                    ss << folders[h] << "/vacs/vacancies_output_" << (i) << "_" << l << "_" << k << "rate.txt";
                    output_filename = ss.str();
                    write_to_file(output_filename, vacancies[l]);
                    ss.str("");
                    ss.clear();
                }

                for (int l=0; l<all_times.size(); l++) {
                    ss << folders[h] << "/all_times/alltimes_output_" << (i) << "_" << l << "_" << k << "rate.txt";
                    std::cout << ss.str() << "\n";
                    alltimes_filename = ss.str();
                    out_file.open(alltimes_filename);
                    out_file << all_times[l] << "\n";
                    ss.str("");
                    ss.clear();
                    out_file.close();
                } 

                for (int l=0; l<(int)counts.size(); l++) {
                    ss << folders[h] << "/counts/counts_output_" << (i) << "_" << l << "_" << k << "rate.txt";
                    count_filename = ss.str();
                    out_file.open(count_filename);
                    out_file << counts[l] << "\n";
                    ss.str("");
                    ss.clear();
                    out_file.close();

                    ss << folders[h] << "/times/times_output_" << (i) << "_" << l << "_" << k << "rate.txt";
                    times_filename = ss.str();
                    out_file.open(times_filename);
                    out_file << times[l] << "\n";
                    ss.str("");
                    ss.clear();
                    out_file.close();
                }

            }
        }
    }

    for (int i=0; i < (int)vacancies.size(); i++) { delete vacancies[i]; }
    
    return 0;
}
