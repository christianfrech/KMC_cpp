#include <iostream>
#include <sstream>
#include <stdio.h>
#include <string>
#include <vector>
#include <cstdlib>
#include <fstream>

#include "str_func.hpp"

/*!
 * \brief Prints the contents of a 3D vector of doubles.
 * \param vec The 3D vector to be printed.
 */
void print_3Dvector_float(std::vector<std::vector<std::vector<double>>> vec) {
    int M = (int)vec.size();
    std::cout << "[ ";
    for (int m = 0; m < M; m++) {
        std::cout << "[ ";
        for (int n = 0; n < (int)vec[m].size(); n++) {
            std::cout << "[ ";
            for (int o = 0; o < (int)vec[m][n].size(); o++) {
                std::cout << vec[m][n][o] << " ";
            }
            std::cout << "] ";
        }
        std::cout << "] ";
        std::cout << "\n  ";
    }
    std::cout << "] \n\n";
}

/*!
 * \brief Prints the contents of a 3D vector of integers.
 * \param vec The 3D vector to be printed.
 */
void print_3Dvector(std::vector<std::vector<std::vector<int>>> vec) {
    int M = (int)vec.size();
    std::cout << "[ ";
    for (int m = 0; m < M; m++) {
        std::cout << "[ ";
        for (int n = 0; n < (int)vec[m].size(); n++) {
            std::cout << "[ ";
            for (int o = 0; o < (int)vec[m][n].size(); o++) {
                std::cout << vec[m][n][o] << " ";
            }
            std::cout << "] ";
        }
        std::cout << "] ";
        std::cout << "\n  ";
    }
    std::cout << "] \n\n";
}

/*!
 * \brief Prints the contents of a 2D vector of doubles.
 * \param vec The 2D vector to be printed.
 */
void print_2Dvector(std::vector<std::vector<double>> vec) {
    int M = (int)vec.size();
    std::cout << "[ ";
    for (int m = 0; m < M; m++) {
        std::cout << "[ ";
        for (int n = 0; n < (int)vec[m].size(); n++) {
            std::cout << vec[m][n] << " ";
        }
        std::cout << "] ";
        std::cout << "\n  ";
    }
    std::cout << "] \n\n";
}

/*!
 * \brief Prints the contents of a 2D vector of integers.
 * \param vec The 2D vector to be printed.
 */
void print_2Dvector(std::vector<std::vector<int>> vec) {
    int M = (int)vec.size();
    std::cout << "[ ";
    for (int m = 0; m < M; m++) {
        std::cout << "[ ";
        for (int n = 0; n < (int)vec[m].size(); n++) {
            std::cout << vec[m][n] << " ";
        }
        std::cout << "] ";
        std::cout << "\n  ";
    }
    std::cout << "] \n\n";
}

/*!
 * \brief Prints the contents of a 2D vector of size_t.
 * \param vec The 2D vector to be printed.
 */
void print_2Dvector(std::vector<std::vector<size_t>> vec) {
    int M = (int)vec.size();
    std::cout << "[ ";
    for (int m = 0; m < M; m++) {
        std::cout << "[ ";
        for (int n = 0; n < (int)vec[m].size(); n++) {
            std::cout << vec[m][n] << " ";
        }
        std::cout << "] ";
        std::cout << "\n  ";
    }
    std::cout << "] \n\n";
}

/*!
 * \brief Prints the contents of a 1D vector of integers.
 * \param vec The 1D vector to be printed.
 */
void print_1Dvector(std::vector<int> vec) {
    int N = (int)vec.size();
    std::cout << "[ ";
    for (int n = 0; n < N; n++) {
        std::cout << vec[n] << " ";
    }
    std::cout << "] \n\n";
}

/*!
 * \brief Prints the contents of a 1D array of integers.
 * \param arr Pointer to the integer array.
 * \param size The size of the array.
 */
void print_1Darr(int* arr, int size) {
    std::cout << "[ ";
    for (int n = 0; n < size; n++) {
        std::cout << arr[n] << " ";
    }
    std::cout << "] \n\n";
}

/*!
 * \brief Prints the contents of a 1D vector of doubles.
 * \param vec The 1D vector to be printed.
 */
void print_1Dvector(std::vector<double> vec) {
    int N = (int)vec.size();
    std::cout << "[ ";
    for (int n = 0; n < N; n++) {
        std::cout << vec[n] << " ";
    }
    std::cout << "] \n\n";
}

/*!
 * \brief Prints the contents of a 1D vector of strings.
 * \param vec The 1D vector to be printed.
 */
void print_1Dvector(std::vector<std::string> vec) {
    int N = (int)vec.size();
    std::cout << "[ ";
    for (int n = 0; n < N; n++) {
        std::cout << vec[n] << " ";
    }
    std::cout << "] \n\n";
}

/*!
 * \brief Prints the contents of a 1D vector of size_t values.
 * \param vec Reference to the 1D vector to be printed.
 */
void print_1Dvector(std::vector<size_t>& vec) {
    int N = (int)vec.size();
    std::cout << "[ ";
    for (int n = 0; n < N; n++) {
        std::cout << vec[n] << " ";
    }
    std::cout << "] \n\n";
}

/*!
 * \brief Prints the contents of a 1D vector of boolean values.
 * \param vec Reference to the 1D vector to be printed.
 */
void print_1Dvector(std::vector<bool>& vec) {
    int N = (int)vec.size();
    std::cout << "[ ";
    for (int n = 0; n < N; n++) {
        std::cout << vec[n] << " ";
    }
    std::cout << "] \n\n";
}


/*!
 * \brief Writes output data in a parallel manner using MPI.
 * 
 * This function handles the parallel writing of output data across multiple processes.
 * It communicates with all the processes in the MPI world to gather data and write the output files.
 * The data is gathered and processed in chunks by different MPI ranks and written to files by rank 0.
 * 
 * \param vacancies The matrix of vacancy data to be written.
 * \param dims A vector containing the dimensions of the grid.
 * \param folder The folder where the output will be written.
 * \param xprocs The number of processes in the x direction.
 * \param yprocs The number of processes in the y direction.
 * \param nprocs The total number of processes.
 * \param rank The rank of the current process.
 * \param i The index of the current time step.
 * \param l The index of the current iteration.
 * \param k The current move tick count.
 * \param t The time at the current step.
 * \param get_rank A flag indicating whether to include the rank in the output data (default is false).
 */
/*
void write_output_parallel(const Matrix<int>& vacancies, std::vector<int> dims, std::string folder, int xprocs, int yprocs, int nprocs, int rank, int i, int l, int k, double t, int curr_num_vacs, bool get_rank = false) {
    
    // Initialize variables
    int size = (int)( vacancies.rows() * vacancies.cols() );
    std::vector<int> output(size);
    int idx = 0;
    int size2 = 0; int size3 = 0; int size_in;
    int i1, i2, i3, i4, x_idx, y_idx, x_chunk_start, y_chunk_start;
    std::vector<int> vacs_in;
    MPI_Status status;
    MPI_Request request;
    std::ostringstream ss;
    
    std::vector<int> nums_and_proc(2);
    std::vector<int> num_proc_buffer(2);
    std::vector<int> receive_counts(nprocs);
    std::vector<int> receive_displacements(nprocs, 0);
    std::vector<int> sum_vec(nprocs);

    int num_elems = (int)(vacancies.rows() * 4);
    nums_and_proc[0] = num_elems;
    nums_and_proc[1] = rank;

    int sum_of_elems = 0;
    
    // Send the number of elements and process rank to rank 0
    if (rank != 0) {
        MPI_Isend(
            nums_and_proc.data(),  // Address of the message we are sending
            2,                     // Number of elements handled by the address
            MPI_INT,               // MPI type of the message
            0,                     // Rank of receiving process
            4,                     // Message tag
            MPI_COMM_WORLD,        // MPI communicator
            &request 
        );
    }

    // Get current number of vacancies in the system

    MPI_Barrier(MPI_COMM_WORLD);  // Synchronize all processes

    if (rank == 0) {
        receive_counts[0] = num_elems;
        
        // Receive the number of elements and ranks from other processes
        for (int rec_proc = 0; rec_proc < nprocs; rec_proc++) {
            if (rec_proc != rank) {
                MPI_Recv(
                    num_proc_buffer.data(),  // Address of the message we are receiving
                    2,                       // Number of elements handled by the address
                    MPI_INT,                 // MPI type of the message
                    rec_proc,                // Rank of sending process
                    4,                       // Message tag
                    MPI_COMM_WORLD,          // MPI communicator
                    &status 
                ); 

                receive_counts[num_proc_buffer[1]] = num_proc_buffer[0];
            }              
        }
        
        // Calculate the displacements and sums for each process
        for (int sum_i = 0; sum_i < (int)receive_counts.size(); sum_i++) {
            if (sum_i == 0) {
                receive_displacements[sum_i] = 0;
                sum_vec[sum_i] = receive_counts[sum_i]; 
            } else {
                receive_displacements[sum_i] = receive_counts[sum_i - 1] + receive_displacements[sum_i - 1];
                sum_vec[sum_i] = receive_counts[sum_i] + sum_vec[sum_i - 1];
            }
        }

        for (auto& n : receive_counts) 
            sum_of_elems += n;
    }

    vacs_in.resize(sum_of_elems);
    MPI_Barrier(MPI_COMM_WORLD);

    // Broadcast displacements and counts to all processes
    MPI_Bcast(receive_displacements.data(), (int)receive_displacements.size(), MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(receive_counts.data(), (int)receive_counts.size(), MPI_INT, 0, MPI_COMM_WORLD);

    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Gatherv(vacancies.data(), (vacancies.rows() * 4), MPI_INT, vacs_in.data(), receive_counts.data(), receive_displacements.data(), MPI_INT, 0, MPI_COMM_WORLD);

    MPI_Barrier(MPI_COMM_WORLD);

    int shift_idx = 0;
    int coords_size = (int)(vacs_in.size() / 4);
    
    if (rank == 0) {
        std::cout << "writing move_ticks: " << k << "\n";

        // Initialize output matrices for vacancies and process ranks
        FourDBoolArr vacancies_out((size_t)2, (size_t)dims[0], (size_t)dims[1], (size_t)dims[2]);
        FourDArr proc_rank((size_t)2, (size_t)dims[0], (size_t)dims[1], (size_t)dims[2]);
        vacancies_out.zero();
        proc_rank.zero();
        
        // Print the sum vector
        print_1Dvector(sum_vec);

        // Process data and assign to output arrays
        for (int idx = 0; idx < coords_size; idx++) {
            if ((shift_idx < nprocs) && (idx * 4 == sum_vec[shift_idx])) { 
                shift_idx++; 
            }
        
            x_idx = shift_idx % xprocs;
            y_idx = floor(shift_idx / xprocs);
            
            x_chunk_start = (int)(dims[0] / xprocs * x_idx);
            y_chunk_start = (int)(dims[1] / yprocs * y_idx); 

            i1 = vacs_in[4 * idx];
            i2 = vacs_in[4 * idx + 1] + x_chunk_start;
            i3 = vacs_in[4 * idx + 2] + y_chunk_start;
            i4 = vacs_in[4 * idx + 3];
            
            if (!get_rank) {
                vacancies_out(i1, i2, i3, i4) = 1;
            } else {
                proc_rank(i1, i2, i3, i4) = shift_idx + 1;
            }
        }

        std::string output_filename = ss.str();
        ss << folder << "/vacs/vacancies_output_" << i << "_" << l << "_" << k << "_" << t << "_moves.txt";
        std::cout << "filename: " << ss.str() << "\n";

        // Write the vacancies or vacancies and rank to file
        if (!get_rank) {
            Matrix<int> all_vacancies = vacancies_out.nonzero();
            write_to_file(ss.str(), all_vacancies);
            std::cout << "all_vacancies.rows(): " << all_vacancies.rows() << " curr_num_vacs: " << curr_num_vacs << "\n"; 
            
            if (all_vacancies.rows() != curr_num_vacs) {
                std::cout << "Change in total number of vacancies in simulation -- all_vacancies.rows(): " << all_vacancies.rows() << " curr_num_vacs: " << curr_num_vacs << "\n"; 
                if (k != 0) {
                    exit(0);
                }
            }
        } 
        else { 
            Matrix<int> vacancies_and_rank = proc_rank.nonzero_elems(); 
            write_to_file(ss.str(), vacancies_and_rank);
            std::cout << "vacancies_and_rank.rows(): " << vacancies_and_rank.rows() << " curr_num_vacs: " << curr_num_vacs << "\n"; 
            
            if (vacancies_and_rank.rows() != curr_num_vacs) {
                std::cout << "Change in total number of vacancies in simulation -- vacancies_and_rank.rows(): " << vacancies_and_rank.rows() << " curr_num_vacs: " << curr_num_vacs << "\n"; 
                if (k != 0) {
                    exit(0);
                }
            }
        }
        
        if (reconstruct == true) {
            reconstruct_ghost_sites(&all_vacancies, dims, rank, chunk_bounds, xprocs, yprocs, nprocs)
        }
        

        ss.str("");
        ss.clear();

        
    }

    MPI_Barrier(MPI_COMM_WORLD);  // Synchronize all processes
}
*/

/*!
 * \brief Sums up the vacancies across all processes and returns the total number of vacancies.
 * 
 * This function calculates the total number of vacancies across all processes in an MPI parallel environment.
 * Each process sends its local vacancy data to the root process (rank 0), which gathers and processes the data
 * from all processes, summing the number of vacancies and returning the result.
 * 
 * \param vacancies The matrix of sites of vacancies to be summed over.
 * \param dims A vector containing the dimensions of the lattices (used to calculate chunks per process).
 * \param xprocs The number of processes in the x direction.
 * \param yprocs The number of processes in the y direction.
 * \param nprocs The total number of processes.
 * \param rank The rank of the current process.
 * 
 * \return The total number of vacancies across all processes.
 */
/*
int sum_vacs_allprocs(const Matrix<int>& vacancies, std::vector<int> dims, int xprocs, int yprocs, int nprocs, int rank) {
    // Initialize variables
    int size = (int)( vacancies.rows() * vacancies.cols() );  ///< Size of the vacancy matrix
    std::vector<int> output(size);
    int idx = 0;
    int size2 = 0, size3 = 0, size_in;
    int i1, i2, i3, i4, x_idx, y_idx, x_chunk_start, y_chunk_start;
    std::vector<int> vacs_in;

    MPI_Status status;
    MPI_Request request;
    std::ostringstream ss;

    // MPI communication buffers
    std::vector<int> nums_and_proc(2);  ///< Number of elements and process rank
    std::vector<int> num_proc_buffer(2);  ///< Buffer for receiving data from other processes
    std::vector<int> receive_counts(nprocs);  ///< Number of elements to receive from each process
    std::vector<int> receive_displacements(nprocs, 0);  ///< Displacements of received data
    std::vector<int> sum_vec(nprocs);  ///< Running sum of received data sizes

    int num_elems = (int)(vacancies.rows() * 4);  ///< Total number of elements in the vacancy matrix
    nums_and_proc[0] = num_elems;
    nums_and_proc[1] = rank;

    int sum_of_elems = 0;

    // Send the number of elements and rank to the root process (rank 0)
    if (rank != 0) {
        MPI_Isend(
            nums_and_proc.data(),  // Address of the message we are sending
            2,                     // Number of elements
            MPI_INT,               // MPI data type
            0,                     // Rank of receiving process
            4,                     // Message tag
            MPI_COMM_WORLD,        // MPI communicator
            &request 
        );
    }

    MPI_Barrier(MPI_COMM_WORLD);  ///< Synchronize all processes

    // Root process (rank 0) handles gathering of data
    if (rank == 0) {
        receive_counts[0] = num_elems;

        // Receive the number of elements from each process
        for (int rec_proc = 0; rec_proc < nprocs; rec_proc++) {
            if (rec_proc != rank) {
                MPI_Recv(
                    num_proc_buffer.data(),  // Address of the message we are receiving
                    2,                       // Number of elements
                    MPI_INT,                 // MPI data type
                    rec_proc,                // Rank of sending process
                    4,                       // Message tag
                    MPI_COMM_WORLD,          // MPI communicator
                    &status 
                ); 

                receive_counts[num_proc_buffer[1]] = num_proc_buffer[0];
            }              
        }
        
        // Compute the displacements and total sum of elements from each process
        for (int sum_i = 0; sum_i < (int)receive_counts.size(); sum_i++) {
            if (sum_i == 0) {
                receive_displacements[sum_i] = 0;
                sum_vec[sum_i] = receive_counts[sum_i]; 
            } else {
                receive_displacements[sum_i] = receive_counts[sum_i - 1] + receive_displacements[sum_i - 1];
                sum_vec[sum_i] = receive_counts[sum_i] + sum_vec[sum_i - 1];
            }
        }

        // Sum up the total number of elements across all processes
        for (auto& n : receive_counts) 
            sum_of_elems += n;
    }

    vacs_in.resize(sum_of_elems);  ///< Resize the vacancy data buffer

    MPI_Barrier(MPI_COMM_WORLD);  ///< Synchronize all processes

    // Broadcast the receive displacements and counts to all processes
    MPI_Bcast(receive_displacements.data(), (int)receive_displacements.size(), MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(receive_counts.data(), (int)receive_counts.size(), MPI_INT, 0, MPI_COMM_WORLD);

    MPI_Barrier(MPI_COMM_WORLD);  ///< Synchronize all processes

    // Gather the vacancy data from all processes to rank 0
    MPI_Gatherv(vacancies.data(), (vacancies.rows() * 4), MPI_INT, vacs_in.data(), receive_counts.data(), receive_displacements.data(), MPI_INT, 0, MPI_COMM_WORLD);

    MPI_Barrier(MPI_COMM_WORLD);  ///< Synchronize all processes

    int shift_idx = 0;
    int coords_size = (int)(vacs_in.size() / 4);

    // Root process (rank 0) processes the gathered data
    if (rank == 0) {
        FourDBoolArr vacancies_out((size_t)2, (size_t)dims[0], (size_t)dims[1], (size_t)dims[2]);  ///< Output array for vacancies
        vacancies_out.zero();  ///< Initialize the output array

        // Process each vacancy and update the output array
        for (int idx = 0; idx < coords_size; idx++) {
            if ((shift_idx < nprocs) && (idx * 4 == sum_vec[shift_idx])) { 
                shift_idx++;  ///< Move to the next process
            }

            x_idx = shift_idx % xprocs;
            y_idx = floor(shift_idx / xprocs);

            // Calculate chunk starting indices for x and y
            x_chunk_start = (int)(dims[0] / xprocs * x_idx);
            y_chunk_start = (int)(dims[1] / yprocs * y_idx);

            i1 = vacs_in[4 * idx];
            i2 = vacs_in[4 * idx + 1] + x_chunk_start;
            i3 = vacs_in[4 * idx + 2] + y_chunk_start;
            i4 = vacs_in[4 * idx + 3];

            vacancies_out(i1, i2, i3, i4) = 1;  ///< Set vacancy in the output array
        }

        // Return the total number of non-zero vacancies
        Matrix<int> all_vacancies = vacancies_out.nonzero();
        return (int)all_vacancies.rows();
    }

    // In case of other ranks, no output is generated.
}


void reconstruct_ghost_sites(const Matrix<int>& vacancies, std::vector<int> dims, int rank, std::vector<std::vector<int>> chunk_bounds, int xprocs, int yprocs, int nprocs) {
    int lattice_pos; int x; int y; int z;
    int x_idx; int y_idx; int z_idx;

    int xlo_edge = (((chunk_bounds[0][0]-1) % dims[0] + dims[0]) % dims[0]);
    int xhi_edge = (((chunk_bounds[0][1]) % dims[0] + dims[0]) % dims[0]);
    int ylo_edge = (((chunk_bounds[1][0]-1) % dims[1] + dims[1]) % dims[1]);
    int yhi_edge = (((chunk_bounds[1][1]) % dims[1] + dims[1]) % dims[1]);

    for (int vac_idx=0; vac_idx<(int)vacancies_pos.size(); vac_idx++) {
        
        lattice_pos = vacancies[vac_idx][0];
        x = vacancies[vac_idx][1];
        y = vacancies[vac_idx][2];
        z = vacancies[vac_idx][3];

        if (((procs[0] != 1) || (procs[1] != 1)) && (atomtype == 0)) {

            // first layer ghost sites for (111) and (100) moves
            if ((y >= (chunk_bounds[1][0] - 1)) && (y < (chunk_bounds[1][1]) )) {
                if ((x == xlo_edge) && (lattice_pos == 1)) {
                    
                    y_idx = mod_with_bounds((y - chunk_bounds[1][0] + 1), dims[1]);
                    (*temp_proc_neg_x_neighbors)(0,1,y_idx,z_idx) = 1;
                    
                    //filling in other corner to preserve corner periodic boundary conditions
                    if ((y_idx == x_dims[2]-1) && (temp_proc_neighbors(rank,2) == rank)) { 
                        (*temp_proc_neg_x_neighbors)(0,1,0,z_idx) = 1; }

                    if ((y == ylo_edge)) {
                        x_idx = mod_with_bounds((x - chunk_bounds[0][0] + 1), dims[0]);
                        (*temp_proc_neg_y_neighbors)(0,1,x_idx,z_idx) = 1;
                    }   
                }
            }
            if ((x >= (chunk_bounds[0][0] - 1)) && (x < (chunk_bounds[0][1]) )) {
                if ((y == ylo_edge)) { 
                    if (lattice_pos == 1) {
                        x_idx = mod_with_bounds((x - chunk_bounds[0][0] + 1), dims[0]);
                        (*temp_proc_neg_y_neighbors)(0,1,x_idx,z_idx) = 1;

                        //filling in other corner to preserve corner periodic boundary conditions
                        if ((x_idx == y_dims[2]-1) && (temp_proc_neighbors(rank,0) == rank)) { 
                            (*temp_proc_neg_y_neighbors)(0,1,0,z_idx) = 1; }
                    }
                }
            }

            // second layer ghost sites for (100) moves
            if ((y >= (chunk_bounds[1][0])) && (y < (chunk_bounds[1][1] + 1) )) {
                if ((x == xhi_edge) && (lattice_pos == 1)) {
                    y_idx = mod_with_bounds((y - chunk_bounds[1][0]), dims[1]);
                    (*temp_proc_pos_x_neighbors)(0,1,y_idx,z_idx) = 1;

                    //filling in other corner to preserve corner periodic boundary conditions
                    if ((y_idx == x_dims_pos[2]-1) && (temp_proc_neighbors(rank,2) == rank)) { 
                        (*temp_proc_pos_x_neighbors)(0,1,0,z_idx) = 1; }
                }
            }
            if ((x >= (chunk_bounds[0][0])) && (x < (chunk_bounds[0][1] + 1) )) {
                if ((y == yhi_edge) && (lattice_pos == 1)) {
                    
                    x_idx = mod_with_bounds((x - chunk_bounds[0][0]), dims[0]);
                    (*temp_proc_pos_y_neighbors)(0,1,x_idx,z_idx) = 1;

                    //filling in other corner to preserve corner periodic boundary conditions
                    if ((x_idx == 0) && (temp_proc_neighbors(rank,0) == rank)) { 
                        (*temp_proc_pos_y_neighbors)(0,1,y_dims_pos[2]-1,z_idx) = 1; }
                }
            }

            // first layer ghost sites for (111) and (100) moves
            if ((y >= (chunk_bounds[1][0])) && (y < (chunk_bounds[1][1] + 1) )) {
                if ((x == xhi_edge) && (lattice_pos == 0)) {
                    
                    y_idx = mod_with_bounds((y - chunk_bounds[1][0]), dims[1]);
                    (*temp_proc_pos_x_neighbors)(0,0,y_idx,z_idx) = 1;

                    //filling in other corner to preserve corner periodic boundary conditions
                    if ((y_idx == 0) && (temp_proc_neighbors(rank,2) == rank)) { 
                        (*temp_proc_pos_x_neighbors)(0,0,(size_t)(x_dims_pos[2]-1),z_idx) = 1; 
                    }                    
                    if ((y == yhi_edge)) {
                        x_idx = mod_with_bounds((x - chunk_bounds[0][0]), dims[0]);
                        (*temp_proc_pos_y_neighbors)(0,0,x_idx,z_idx) = 1;
                    }
                }
            }
            if ((x >= (chunk_bounds[0][0])) && (x < (chunk_bounds[0][1] + 1) )) {
                if ((y == yhi_edge)) {
                    x_idx = mod_with_bounds((x - chunk_bounds[0][0]), dims[0]);    
                    (*temp_proc_pos_y_neighbors)(0,0,x_idx,z_idx) = 1; 

                    //filling in other corner to preserve corner periodic boundary conditions
                    if ((x_idx == 0) && (temp_proc_neighbors(rank,0) == rank)) { 
                        (*temp_proc_pos_y_neighbors)(0,0,(size_t)(y_dims_pos[2]-1),z_idx) = 1; 
                        }
                }
            }

            // second layer ghost sites for (100) moves 
            if ((y >= (chunk_bounds[1][0] - 1)) && (y < (chunk_bounds[1][1]) )) {
                if ((x == xlo_edge) && (lattice_pos == 0)) {

                    y_idx = mod_with_bounds((y - chunk_bounds[1][0] + 1), dims[1]);
                    (*temp_proc_neg_x_neighbors)(0,0,y_idx,z_idx) = 1;

                    //filling in other corner to preserve corner periodic boundary conditions
                    //if (y_idx == x_dims[2]-1) { (*temp_proc_neg_x_neighbors)(0,0,0,z_idx) = 1; }
                }
            }
            if ((x >= (chunk_bounds[0][0] - 1)) && (x < (chunk_bounds[0][1]) )) {
                if ((y == ylo_edge) && (lattice_pos == 0)) {
                    
                    x_idx = mod_with_bounds(((x - 0.5) - chunk_bounds[0][0] + 1), dims[0]);
                    (*temp_proc_neg_y_neighbors)(0,0,x_idx,z_idx) = 1;
                }
            }
        }
    }
}
*/