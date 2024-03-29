//partition down 3d array
#include <mpi.h>


int main(int argc, char *argv[]) {
    nprocs;
    dims = [x,y,z];
    int xprocs;
    int yprocs;
    int zprocs;
    int size; int rank;
    int gcf;

    if (argc == 0) {
        std::cout << "ERROR: no command line arguments provided\n"
        exit(0);
    }
    else {
        gcf = find_gcf(argv[0], dims[0]); 
    }

    xprocs = gcf;
    yprocs = (int)(nprocs / gcf); 

    std::fstream in_file;
    in_file.open(infile_name);
    std::vector<std::string> lines;
    std::string line;
    int read_idx = 0;

    if (in_file.is_open()) {
        while ( getline (in_file,line) )
        {
            lines.push_back(line);
        }
        in_file.close();
    }

    // parsing first line to grab dimensions of lattice ###
    std::string dims = lines[read_idx]; //getting dimension line
    read_idx ++;
    std::vector<std::string> dims_str = tokenizer(dims," "); 
    std::vector<int> dims(3);

    for (int i=0; i<(int)dims.size(); i++) {
        dims[i] = std::stoi(dims_str[i+1]); 
    }
    file.close();

    // specify number of procs using slurm
    MPI_Init();
    MPI_comm_size(MPI_COMM_WORLD, &size);
    MPI_comm_rank(MPI_COMM_WORLD, &rank);
    
    // send one row to each core 

    // chunkify all needed variables
    int x_idx = rank % dims[0];
    int y_idx = floor(rank / x_procs);

    int x_chunk_start = (int)( dims[0]/x_procs * x_idx );
    int x_chunk_end = (int)( dims[0]/x_procs * (x_idx + 1) );
    int y_chunk_start = (int)( dims[1]/y_procs * y_idx );
    int y_chunk_end = (int)( dims[1]/y_procs * (y_idx + 1) );
    std::vector<int> chunk_bounds = {{x_chunk_start, x_chunk_end},{y_chunk_start, y_chunk_end}, {dims[2][0],dims[2][0]}};

    for (int h=0; h<folders.size(); h++) {
        for (int i=0; i<iterations; i++) {
            for (int k=0; k<(int)edge_rates.size(); k++) {
                std::cout << "-----------------------------------------\n";
                std::cout << "i: " << i << "\n"; 
                std::cout << "temp: " << temperatures[k] << "\n"; 
                lattice = populate_lattice(in_file[h], catalog_file[h], vertex_rates[k], edge_rates[k], reg_rates[h], chunk_bounds, rank, nprocs);
                start = std::chrono::system_clock::now();
                return_tuple = lattice->new_kmc_iterator(time, start, folders[h], i);
                delete lattice;
            }
        }
    }
    
    if(process_Rank==0) {
        for(i=1;i<size_Of_Comm;i++)
        {            
            // Sending chunkified variables to each processor
            if(2*i<=N)
            {
                printf("scattering data %f\n", *((distro_Array)+2*i));  
                MPI_Send(
                    distro_Array+2*(i-1),      //Address of the message we are sending.
                    2,                  //Number of elements handled by that address.
                    MPI_DOUBLE,            //MPI_TYPE of the message we are sending.
                    i,                  //Rank of receiving process
                    1,                  //Message Tag
                    MPI_COMM_WORLD      //MPI Communicator
                );
            }
        }
    }
    else {

        // Recieving chunkified variables at each processor

        printf("waiting for data by %d\n", process_Rank);
        if(2*process_Rank<=N)
        {
            MPI_Recv(
                    &scattered_Data,      //Address of the message we are receiving.
                    2,                  //Number of elements handled by that address.
                    MPI_DOUBLE,            //MPI_TYPE of the message we are sending.
                    0,                  //Rank of sending process
                    1,                  //Message Tag
                    MPI_COMM_WORLD,      //MPI Communicator
                    MPI_STATUS_IGNORE   //MPI Status Object
                );

            printf("Process %d has received:  ",process_Rank);
            
            if ((curr_move_num + (num_interface_sites - vacs_on_interface)) >= ((int)moves_shifts.rows() - 20)) {
                // resizing data structures to accommodate all moves 

                int newsize = 2 * moves_shifts.rows();
                rate_cumsum.resize(newsize);
                moves_coords.reshape(newsize, 4);
                moves_shifts.reshape(newsize, 3);
                moves_lattice.reshape(newsize, 1);
                moves_vacs.reshape(newsize, 1);
            }
            /*
            // finding all moves along the {111} family of vectors
            for (int s=0; s < (int)diag_directions.size(); s++) {
                if ((l == 0) && (i == 0) && (diag_directions[s][2] == 1)) {/* checking for leftmost non-periodic boundary along z-axis*/}
                
            /*
                else if ((l == (int)(lattice_dim[2]-1)) && (i == 1) && (diag_directions[s][2] == 1)) {/* checking for rightmost non-periodic boundary along z-axis*///}
            /*
                else {
                    moves_coords[curr_move_num][0] = i;

                    if ((i == 0) && (vacancies(1, (((j - diag_directions[s][0]) % lattice_dim[0] + lattice_dim[0]) % lattice_dim[0]), (((k - diag_directions[s][1]) % lattice_dim[1] + lattice_dim[1]) % lattice_dim[1]), (((l - diag_directions[s][2]) % lattice_dim[2] + lattice_dim[2]) % lattice_dim[2])) == 0)) {
                        // checking that vertex site -> bc site move has new site occupied by atom

                        moves_coords[curr_move_num][1] = (((j - diag_directions[s][0]) % lattice_dim[0] + lattice_dim[0]) % lattice_dim[0]);
                        moves_coords[curr_move_num][2] = (((k - diag_directions[s][1]) % lattice_dim[1] + lattice_dim[1]) % lattice_dim[1]);
                        moves_coords[curr_move_num][3] = (((l - diag_directions[s][2]) % lattice_dim[2] + lattice_dim[2]) % lattice_dim[2]);                                        
                        moves_shifts[curr_move_num][0] = - diag_directions[s][0];
                        moves_shifts[curr_move_num][1] = - diag_directions[s][1];
                        moves_shifts[curr_move_num][2] = - diag_directions[s][2]; 
                        moves_lattice[curr_move_num][0] = 0;
                        moves_vacs[curr_move_num][0] = idx;
                        
                        // getting rate corresponding to move
                        rate = new_get_rateconstants(moves_coords[curr_move_num], moves_shifts[curr_move_num], moves_lattice[curr_move_num][0]);
                        
                        if (rate == -1) {curr_move_num --;}
                        else {
                            if (curr_move_num == 0) {rate_cumsum[curr_move_num] = rate;}
                            else {rate_cumsum[curr_move_num] = rate + rate_cumsum[curr_move_num-1];}
                        }
                        curr_move_num ++;
                        moves[0] ++;
                    }
                    
                    else if ((i == 1) && (vacancies(0, (((j + diag_directions[s][0]) % lattice_dim[0] + lattice_dim[0]) % lattice_dim[0]), (((k + diag_directions[s][1]) % lattice_dim[1] + lattice_dim[1]) % lattice_dim[1]), (((l + diag_directions[s][2]) % lattice_dim[2] + lattice_dim[2]) % lattice_dim[2])) == 0)) {
                        // checking that bc site -> vertex site move has new site occupied by atom
                        moves_coords[curr_move_num][1] = (((j + diag_directions[s][0]) % lattice_dim[0] + lattice_dim[0]) % lattice_dim[0]);
                        moves_coords[curr_move_num][2] = (((k + diag_directions[s][1]) % lattice_dim[1] + lattice_dim[1]) % lattice_dim[1]);
                        moves_coords[curr_move_num][3] = (((l + diag_directions[s][2]) % lattice_dim[2] + lattice_dim[2]) % lattice_dim[2]);
                        moves_shifts[curr_move_num][0] = diag_directions[s][0];
                        moves_shifts[curr_move_num][1] = diag_directions[s][1];
                        moves_shifts[curr_move_num][2] = diag_directions[s][2];
                        moves_lattice[curr_move_num][0] = 1;
                        moves_vacs[curr_move_num][0] = idx; 
                        
                        // getting rate corresponding to move
                        rate = new_get_rateconstants(moves_coords[curr_move_num], moves_shifts[curr_move_num], moves_lattice[curr_move_num][0]);
                        
                        if (rate == -1) {curr_move_num --;}
                        else {
                            if (curr_move_num == 0) {rate_cumsum[curr_move_num] = rate;}
                            else {rate_cumsum[curr_move_num] = rate + rate_cumsum[curr_move_num-1];}
                        }
                        curr_move_num ++;
                        moves[0] ++;
                    }
                }
            }
            */
            /*
            // finding all moves along the {100} family of vectors
            for (int s=0; s < (int)edge_directions.size(); s++) {
                
                if ((l == 0) && (edge_directions[s][2] == -1)) {}  
                else if ((l == (int)(lattice_dim[2]-1)) && (edge_directions[s][2] == 1)) {}
                else if (vacancies(i, (((j + edge_directions[s][0]) % lattice_dim[0] + lattice_dim[0]) % lattice_dim[0]), (((k + edge_directions[s][1]) % lattice_dim[1] + lattice_dim[1]) % lattice_dim[1]), (((l + edge_directions[s][2]) % lattice_dim[2] + lattice_dim[2]) % lattice_dim[2])) == 0) {
                    // checking that vertex site -> vertex site or bc site -> bc site move has new site occupied by atom

                    moves_shifts[curr_move_num][0] = edge_directions[s][0];
                    moves_shifts[curr_move_num][1] = edge_directions[s][1];
                    moves_shifts[curr_move_num][2] = edge_directions[s][2];
                    moves_coords[curr_move_num][0] = i;
                    moves_coords[curr_move_num][1] = (((j + edge_directions[s][0]) % lattice_dim[0] + lattice_dim[0]) % lattice_dim[0]);
                    moves_coords[curr_move_num][2] = (((k + edge_directions[s][1]) % lattice_dim[1] + lattice_dim[1]) % lattice_dim[1]);
                    moves_coords[curr_move_num][3] = (((l + edge_directions[s][2]) % lattice_dim[2] + lattice_dim[2]) % lattice_dim[2]);
                    moves_vacs[curr_move_num][0] = idx; 
                    
                    if (i == 0) {
                        moves_lattice[curr_move_num][0] = 2;
                    }
                    else if (i == 1) {
                        moves_lattice[curr_move_num][0] = 3;
                    }
                    
                    // getting rate corresponding to move
                    rate = new_get_rateconstants(moves_coords[curr_move_num], moves_shifts[curr_move_num], moves_lattice[curr_move_num][0]);
                    
                    if (rate == -1) {curr_move_num --;}
                    else {
                        if (curr_move_num == 0) {rate_cumsum[curr_move_num] = rate;}
                        else {rate_cumsum[curr_move_num] = rate + rate_cumsum[curr_move_num-1];}
                    }
                    curr_move_num ++;
                    moves[1] ++;
                }
            }
        */
        }

        // UPDATING SIZE OF DATA STRUCTURES CONTIANING COORDINATES AND RATES OF MOVES
        num_of_moves = curr_move_num;
        rate_cumsum.resize(num_of_moves);
        moves_vacs.reshape(num_of_moves, 1);
        moves_coords.reshape(num_of_moves, 4);
        moves_shifts.reshape(num_of_moves, 3);
        moves_lattice.reshape(num_of_moves, 1);
        
        MPI_Send(
            &sum,      //Address of the message we are sending.
            1,                  //Number of elements handled by that address.
            MPI_DOUBLE,            //MPI_TYPE of the message we are sending.
            0,                  //Rank of receiving process
            1,                  //Message Tag
            MPI_COMM_WORLD      //MPI Communicator
            );
        printf("\n");

    
    printf("Process %d has received:  ",process_Rank);
    double sum=0;
    
    
    MPI_Barrier(MPI_COMM_WORLD); // all the sub ranks/processes waits here
   /* process 0 will aggregate the results*/
    if(process_Rank==0)
    {
        for(i=1;i<size_Of_Comm;i++)
        {
            double sum=0;
            MPI_Recv(
                    &sum,      //Address of the message we are receiving.
                    1,                  //Number of elements handled by that address.
                    MPI_DOUBLE,            //MPI_TYPE of the message we are sending.
                    i,                  //Rank of sending process
                    1,                  //Message Tag
                    MPI_COMM_WORLD,      //MPI Communicator
                    MPI_STATUS_IGNORE   //MPI Status Object
                );
            printf("Process %d has sent: %f \n", i, sum);
        }
    }
            
    printf("\n");
    MPI_Finalize(); 
    return 0;
}
