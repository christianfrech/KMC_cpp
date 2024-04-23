//partition down 3d array
#include <mpi.h>
#include "kmc_stripping_parallel.hpp"

int main(int argc, char *argv[]) {
    int nprocs = 4;//argv[0];
    int xprocs;
    int yprocs;
    int zprocs;
    int size; 
    int rank;
    int gcf;    
    int dims[3] = {0,0,0};

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
    Lattice* lattice;
    
    std::vector< std::string > infile_name = {"geo_in.txt"};
    std::vector<double> vertex_rates = {7.65954e11};
    std::vector<double> edge_rates = {1.23537e6};
    std::vector<std::string> folders = {"parallel_output"};
    std::vector< std::string > catalog_file = {"ratecatalog_nogb.txt"};
    std::string region_infile;
    std::vector< std::vector< std::vector<double> > > reg_rates = {{{8.28403e10, 5e12}}}; //settings of kb = 8.6173e-5, T = 300, E = 0.05, prefactor = 5e12
    std::tuple< std::vector< std::vector< std::vector<int> > >, std::vector<int>, std::vector<double>, std::vector<double> > return_tuple;
    int iterations = 1;
    double time = 5e9;

    std::vector< std::vector< std::vector<int> > > vacancies; std::vector<int> counts; std::vector<double> times; std::vector<double> all_times;
    

    for (int h=0; h<(int)folders.size(); h++) {
        for (int i=0; i<iterations; i++) {
            for (int k=0; k<(int)edge_rates.size(); k++) {

                in_file.open(infile_name[h]);

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

                for (int idx=0; i<3; idx++) {
                    dims[idx] = std::stoi(dims_toks[idx+1]); 
                }

                if (nprocs == 0) {
                    std::cout << "ERROR: no command line arguments provided\n";
                    exit(0); 
                }
                else {
                    gcf = find_gcf(nprocs, dims[0]); 
                }

                xprocs = gcf;
                yprocs = (int)(nprocs / gcf);

                // specify number of procs using slurm
                MPI_Init(NULL, NULL);
                MPI_Comm_size(MPI_COMM_WORLD, &size);
                MPI_Comm_rank(MPI_COMM_WORLD, &rank);
                

                // chunkify all needed variables
                int x_idx = rank % dims[0];
                int y_idx = floor(rank / xprocs);

                int x_chunk_start = (int)( dims[0]/xprocs * x_idx );
                int x_chunk_end = (int)( dims[0]/xprocs * (x_idx + 1) );
                int y_chunk_start = (int)( dims[1]/yprocs * y_idx );
                int y_chunk_end = (int)( dims[1]/yprocs * (y_idx + 1) );
                std::vector<int> procs = {xprocs, yprocs};

                std::vector< std::vector<int> > chunk_bounds = {{x_chunk_start, x_chunk_end},{y_chunk_start, y_chunk_end}, {0,dims[2]}};
                std::chrono::system_clock::time_point start;

                std::cout << "-----------------------------------------\n";
                std::cout << "i: " << i << "\n"; 
                lattice = populate_lattice(infile_name[h], catalog_file[h], region_infile, vertex_rates[k], edge_rates[k], reg_rates[h], chunk_bounds, rank, procs);
                start = std::chrono::system_clock::now();
                lattice_return_struct return_tuple = lattice->new_kmc_iterator(time, start, folders[h], i);
                delete lattice;
                
                vacancies = return_tuple.get_all_vacancies(); counts = return_tuple.get_move_counts(); times = return_tuple.get_time_count(); all_times =return_tuple.get_all_times();

                MPI_Barrier(MPI_COMM_WORLD); // all the sub ranks/processes waits here
                /* process 0 will aggregate the results*/ 
                int size1 = (int)vacancies.size();
                int size2 = (int)vacancies[0].size();
                int size3 = (int)vacancies[0][0].size();
                int size = (int)(vacancies.size() * vacancies[0].size() * vacancies[0][0].size());
                std::vector<int> output(size);
                int idx = 0;

                for (int frames=0; frames<(int)vacancies.size(); frames++) {
                    for (int vacs=0; vacs<(int)vacancies[0].size(); vacs++) {
                        for (int coord=0; coord<(int)vacancies[0][0].size(); coord++) {
                            idx = (int)frames*vacancies[0].size()*vacancies[0][0].size() + (int)vacs*vacancies[0].size() + coord;
                            output[idx] = vacancies[frames][vacs][coord];
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

                std::vector< std::vector<int> > all_vacancies = vacancies_out.nonzero();

                for (int l=0; l<all_vacancies.size(); l++) {
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
    
    
    return 0;
}
