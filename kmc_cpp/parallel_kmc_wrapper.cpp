//partition down 3d array
#include <mpi.h>
#include "kmc_stripping_parallel.hpp"

int main() {
    std::cout << "first line \n";
    int xprocs;
    int yprocs;
    int zprocs;
    int nprocs; 
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
    std::vector< std::string > infile_name = {"hexagon_corrected_10k.txt"};     
    std::vector<double> vertex_rates = {7.65954e11};
    std::vector<double> edge_rates = {1.23537e6};
    std::vector<std::string> folders = {"parallel_output"};
    std::vector< std::string > catalog_file = {"ratecatalog_nogb.txt"};
    std::string region_infile;
    std::vector< std::vector< std::vector<double> > > reg_rates = {{{4.56e7, 5e12}}}; //{{{8.28403e10, 8.28403e10}}}; //settings of kb = 8.6173e-5, T = 300, E = 0.05, prefactor = 5e12
    std::tuple< std::vector< std::vector< std::vector<int> > >, std::vector<int>, std::vector<double>, std::vector<double> > return_tuple;
    int iterations = 1;
    double time = 100;
    int waitint = 0;

    std::vector<int> counts; std::vector<double> times; std::vector<double> all_times;    
    
    std::cout << "variables initialized \n";

    bool restart = false;
    int last_move_tick = 0;
    double last_time = 0;
    std::tuple<int, double> tup_out;
                
    for (int h=0; h<(int)folders.size(); h++) {
        for (int i=0; i<iterations; i++) {
            for (int k=0; k<(int)edge_rates.size(); k++) {
                tup_out = get_last_iter_time(folders[h]);

                last_move_tick = std::get<0>(tup_out);
                last_time = std::get<1>(tup_out);
                std::cout << "last_time: " << last_time << "\n";
                std::cout << "last_move_tick: " << last_move_tick << "\n";

                if (last_move_tick!=0) restart = true;

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
                
                print_1Dvector(dims_toks);

                for (int idx=0; idx<3; idx++) {
                    dims[idx] = std::stoi(dims_toks[idx+1]); 
                }                

                std::cout << "starting mpi stuff \n";

                // specify number of procs using slurm
                MPI_Init(NULL, NULL);
                MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
                MPI_Comm_rank(MPI_COMM_WORLD, &rank);

                std::cout << "nprocs: " << nprocs << "\n";
                if (nprocs == 0) {
                    std::cout << "ERROR: no command line arguments provided\n";
                    exit(0); 
                }
                else {
                    gcf = create_partition(nprocs, dims[0]); 
                }
                std::cout << "gcf: " << gcf << "\n";
                xprocs = gcf;
                yprocs = (int)(nprocs / gcf);
                xprocs = 2;//nprocs;
                yprocs = 2;//1;
                std::cout << "xprocs: " << xprocs << " yprocs: " << yprocs << "\n";

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

                MPI_Barrier(MPI_COMM_WORLD);
                if (rank == 0) {
                    std::cout << "rank: " << rank <<"\n";
                    print_2Dvector(chunk_bounds);
                }
                MPI_Barrier(MPI_COMM_WORLD);
                if (rank == 1) {
                    std::cout << "rank: " << rank <<"\n";
                    print_2Dvector(chunk_bounds);
                }
                MPI_Barrier(MPI_COMM_WORLD);
                if (rank == 2) {
                    std::cout << "rank: " << rank <<"\n";
                    print_2Dvector(chunk_bounds);
                }
                MPI_Barrier(MPI_COMM_WORLD);
                if (rank == 3) {
                    std::cout << "rank: " << rank <<"\n";
                    print_2Dvector(chunk_bounds);
                }
                MPI_Barrier(MPI_COMM_WORLD);

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
                lattice_return_struct return_tuple = lattice->new_kmc_iterator(time, start, folders[h], i, k, last_move_tick, last_time);
                delete lattice;
                
                counts = return_tuple.get_move_counts(); times = return_tuple.get_time_count(); all_times =return_tuple.get_all_times();

                MPI_Barrier(MPI_COMM_WORLD); // all the sub ranks/processes waits here

                std::cout << "rank: " << rank << " post main call barrier \n";
                /* process 0 will aggregate the results*/
                /*
                int size = 0;
                int size1 = (int)vacancies.size();

                for (int frame=0; frame < (int)vacancies.size(); frame++) {
                    size += (int)( vacancies[frame]->rows() * vacancies[frame]->cols());
                }

                std::vector<int> output(size);
                int idx = 0;
                Matrix<int>* one_vac;
                int size2 = 0;
                int size3 = 0;
                int size_in;
                std::vector<int> vacs_in;
                Matrix<int>* proc_vacancies;
                std::vector< Matrix<int>* > vacancies_allframes;
                FourDBoolArr vacancies_out(2, dims[0], dims[1], dims[2]);
                int i1; int i2; int i3; int i4;
                MPI_Status status;
                MPI_Request request;

                std::vector<int> nums_and_proc(2);
                std::vector<int> num_proc_buffer(2);

                for (int frames=0; frames<(int)vacancies.size(); frames++) {
                    vacancies_out.zero();
                    int num_elems = (int)(vacancies[frames]->rows() * 4);
                    nums_and_proc[0] = num_elems;
                    nums_and_proc[1] = rank;
                    std::cout << "rank: " << rank << " nums_and_proc: \n";
                    print_1Dvector(nums_and_proc);
                    std::cout << "rank: " << rank << " vacancies[frames]: \n";
                    vacancies[frames]->print();

                    MPI_Isend(
                        nums_and_proc.data(),      //Address of the message we are receiving.
                        2,                  //Number of elements handled by that address.
                        MPI_INT,            //MPI_TYPE of the message we are sending.
                        0,                  //Rank of receiving process
                        1,                  //Message Tag
                        MPI_COMM_WORLD,      //MPI Communicator
                        &request 
                    );
                    MPI_Isend(
                        vacancies[frames]->data(),      //Address of the message we are receiving.
                        (int)(vacancies[frames]->rows() * 4),        //Number of elements handled by that address.
                        MPI_INT,            //MPI_TYPE of the message we are sending.
                        0,                  //Rank of receiving process
                        2,                  //Message Tag
                        MPI_COMM_WORLD,     //MPI Communicator
                        &request 
                    );

                    printf("Process %d has sent \n", rank);
             

                    MPI_Barrier(MPI_COMM_WORLD); // all the sub ranks/processes waits here

                    if (rank == 0) {
                        for (int rec_proc=0; rec_proc<nprocs; rec_proc++) {
                                    
                            std::cout << "rank: " << rank << " rec_proc " << rec_proc << "\n";
                            MPI_Recv(
                                num_proc_buffer.data(),      //Address of the message we are receiving.
                                2,                  //Number of elements handled by that address.
                                MPI_INT,            //MPI_TYPE of the message we are sending.
                                rec_proc,           //Rank of sending process
                                1,                  //Message Tag
                                MPI_COMM_WORLD,      //MPI Communicator
                                &status 
                            ); 
                           
                            vacs_in.resize(num_proc_buffer[0]);
                            std::cout << "num_proc_buffer[1]: " << num_proc_buffer[1] << " num_proc_buffer[0]: " << num_proc_buffer[0] << "\n";
                            std::cout << "num_proc_buffer[1]: " << num_proc_buffer[1] << " vacs_in.size(): " << vacs_in.size() << "\n";
                            
                            MPI_Recv(
                                    vacs_in.data(),      //Address of the message we are receiving.
                                    num_proc_buffer[0],                  //Number of elements handled by that address.
                                    MPI_INT,            //MPI_TYPE of the message we are sending.
                                    rec_proc,                  //Rank of sending process
                                    2,                  //Message Tag
                                    MPI_COMM_WORLD,      //MPI Communicator
                                    &status 
                            );

                            if (num_proc_buffer[0] != 0) {
                                                        
                                std::cout << "rank: " << rank << " vacs_in: \n";
                                print_1Dvector(vacs_in);
                                x_idx = num_proc_buffer[1] % dims[0];
                                y_idx = floor(num_proc_buffer[1] / xprocs);
                                x_chunk_start = (int)( dims[0]/xprocs * x_idx );
                                y_chunk_start = (int)( dims[1]/yprocs * y_idx );
                                std::cout << "rank: " << rank << " dims[0] " << dims[0] << " dims[1] " << dims[1] << " xprocs " << xprocs << " yprocs " << yprocs << " nums_and_proc[0] " << nums_and_proc[0] << " nums_and_proc[1] " << nums_and_proc[1] << "\n";
                                std::cout << "rank: " << rank << " x_idx " << x_idx << " y_idx " << y_idx << " x_chunk_start " << x_chunk_start << " y_chunk_start " << y_chunk_start << "\n\n";
                                size2 = (int)(vacs_in.size()/4);
                                size3 = 4;
                                for (int k=0; k<size2; k++) {
                                    idx = (int)(k*size3);
                                    i1 = vacs_in[idx];
                                    i2 = vacs_in[idx+1] + x_chunk_start;
                                    i3 = vacs_in[idx+2] + y_chunk_start;
                                    i4 = vacs_in[idx+3];
                                    vacancies_out(i1,i2,i3,i4) = 1;
                                    std::cout << "rank: " << rank << " input:   [ " << i1 << " " << i2 << " " << i3 << " " << i4 << " ]\n";
                                }
                            }                  
                        }
                        proc_vacancies = vacancies_out.nonzero();
                        vacancies_allframes.push_back(proc_vacancies);
                    }

                    MPI_Barrier(MPI_COMM_WORLD); // all the sub ranks/processes waits here
                }
                */
                
                printf("MPI Finalize\n");
                MPI_Finalize();
                /*
                for (int l=0; l<vacancies_allframes.size(); l++) {
                    ss << folders[h] << "/vacs/vacancies_output_" << (i) << "_" << l << "_" << k << "rate.txt";
                    output_filename = ss.str();
                    write_to_file(output_filename, vacancies_allframes[l]);
                    ss.str("");
                    ss.clear();
                }
                */
            }
        }
    }


    return 0;
}
