#include <iostream>
#include <sstream>
#include <stdio.h>
#include <string>
#include <vector>
#include <cstdlib>
#include <fstream>

#include "str_func.hpp"


void print_3Dvector_float(std::vector< std::vector< std::vector<double> > > vec) {
    int M = (int)vec.size();;
    std::cout << "[ ";
    for (int m=0; m<M; m++) {
        std::cout << "[ ";
        for (int n=0; n<(int)vec[m].size(); n++) {
            std::cout << "[ ";
            for (int o=0; o<(int)vec[m][n].size(); o++) {
                std::cout << vec[m][n][o] << " ";
            }
            std::cout << "] ";
        }
        std::cout << "] ";
        std::cout << "\n  ";
    }
    std::cout << "] \n\n";
}

void print_3Dvector(std::vector< std::vector< std::vector<int> > > vec) {
    int M = (int)vec.size();;
    std::cout << "[ ";
    for (int m=0; m<M; m++) {
        std::cout << "[ ";
        for (int n=0; n<(int)vec[m].size(); n++) {
            std::cout << "[ ";
            for (int o=0; o<(int)vec[m][n].size(); o++) {
                std::cout << vec[m][n][o] << " ";
            }
            std::cout << "] ";
        }
        std::cout << "] ";
        std::cout << "\n  ";
    }
    std::cout << "] \n\n";
}

void print_2Dvector(std::vector< std::vector<double> > vec) {
    int M = (int)vec.size();;
    std::cout << "[ ";
    for (int m=0; m<M; m++) {
        std::cout << "[ ";
        for (int n=0; n<(int)vec[m].size(); n++) {
            std::cout << vec[m][n] << " ";
        }
        std::cout << "] ";
        std::cout << "\n  ";
    }
    std::cout << "] \n\n";
}

void print_2Dvector(std::vector< std::vector<int> > vec) {
    int M = (int)vec.size();;
    std::cout << "[ ";
    for (int m=0; m<M; m++) {
        std::cout << "[ ";
        for (int n=0; n<(int)vec[m].size(); n++) {
            std::cout << vec[m][n] << " ";
        }
        std::cout << "] ";
        std::cout << "\n  ";
    }
    std::cout << "] \n\n";
}

void print_2Dvector(std::vector< std::vector<size_t> > vec) {
    int M = (int)vec.size();;
    std::cout << "[ ";
    for (int m=0; m<M; m++) {
        std::cout << "[ ";
        for (int n=0; n<(int)vec[m].size(); n++) {
            std::cout << vec[m][n] << " ";
        }
        std::cout << "] ";
        std::cout << "\n  ";
    }
    std::cout << "] \n\n";
}

void print_1Dvector(std::vector<int> vec) {
    int N = (int)vec.size();
    std::cout << "[ ";
    for (int n=0; n<N; n++) {
        std::cout << vec[n] << " ";
    }
    std::cout << "] \n\n";
}

void print_1Darr(int* arr, int size) {
    std::cout << "[ ";
    for (int n=0; n<size; n++) {
        std::cout << arr[n] << " ";
    }
    std::cout << "] \n\n";
}

void print_1Dvector(std::vector<double> vec) {
    int N = (int)vec.size();
    std::cout << "[ ";
    for (int n=0; n<N; n++) {
        std::cout << vec[n] << " ";
    }
    std::cout << "] \n\n";
}

void print_1Dvector(std::vector<std::string> vec) {
    int N = (int)vec.size();
    std::cout << "[ ";
    for (int n=0; n<N; n++) {
        std::cout << vec[n] << " ";
    }
    std::cout << "] \n\n";
}

void print_1Dvector(std::vector<size_t> &vec) {
    int N = (int)vec.size();
    std::cout << "[ ";
    for (int n=0; n<N; n++) {
        std::cout << vec[n] << " ";
    }
    std::cout << "] \n\n";
}

void print_1Dvector(std::vector<bool> &vec) {
    int N = (int)vec.size();
    std::cout << "[ ";
    for (int n=0; n<N; n++) {
        std::cout << vec[n] << " ";
    }
    std::cout << "] \n\n";
}

void write_output_parallel(const Matrix<int>& vacancies, std::vector<int> dims, std::string folder, int xprocs, int yprocs, int nprocs, int rank, int i, int l, int k, double t) {
    //std::cout << "rank: " << rank << " start write parallel \n";
    int size = (int)( vacancies.rows() * vacancies.cols());

    std::vector<int> output(size);
    int idx = 0;
    int size2 = 0; int size3 = 0; int size_in;
    int i1; int i2; int i3; int i4; int x_idx; int y_idx; int x_chunk_start; int y_chunk_start;
    std::vector<int> vacs_in;

    MPI_Status status;
    MPI_Request request;
    std::ostringstream ss;

    std::vector<int> nums_and_proc(2);
    std::vector<int> num_proc_buffer(2);
    std::vector<int> receive_counts(nprocs);
    std::vector<int> receive_displacements(nprocs,0);
    std::vector<int> sum_vec(nprocs);

    int num_elems = (int)(vacancies.rows() * 4);
    nums_and_proc[0] = num_elems;
    nums_and_proc[1] = rank;

    int sum_of_elems = 0;
    
    if (rank != 0) {
        MPI_Isend(
            nums_and_proc.data(),      //Address of the message we are receiving.
            2,                  //Number of elements handled by that address.
            MPI_INT,            //MPI_TYPE of the message we are sending.
            0,                  //Rank of receiving process
            4,                  //Message Tag
            MPI_COMM_WORLD,      //MPI Communicator
            &request 
        );
    }

    MPI_Barrier(MPI_COMM_WORLD); // all the sub ranks/processes waits here

    if (rank == 0) {

        receive_counts[0] = num_elems;
        // Caught signal 11 (Segmentation fault: Sent by the kernel at address (nil)) at rank = 0
        for (int rec_proc=0; rec_proc<nprocs; rec_proc++) {
            if (rec_proc != rank) {

                MPI_Recv(
                    num_proc_buffer.data(),      //Address of the message we are receiving.
                    2,                  //Number of elements handled by that address.
                    MPI_INT,            //MPI_TYPE of the message we are sending.
                    rec_proc,           //Rank of sending process
                    4,                  //Message Tag
                    MPI_COMM_WORLD,      //MPI Communicator
                    &status 
                ); 

                receive_counts[num_proc_buffer[1]] = num_proc_buffer[0];
            }              
        }
        
        for (int sum_i=0; sum_i<(int)receive_counts.size(); sum_i++) {
            if (sum_i==0) {
                receive_displacements[sum_i] = 0;
                sum_vec[sum_i] = receive_counts[sum_i]; 
            }
            else {
                receive_displacements[sum_i] = receive_counts[sum_i-1] + receive_displacements[sum_i-1];
                sum_vec[sum_i] = receive_counts[sum_i] + sum_vec[sum_i-1];
            }
        }

        for (auto& n : receive_counts) 
            sum_of_elems += n;

    }
    vacs_in.resize(sum_of_elems);

    MPI_Barrier(MPI_COMM_WORLD);

    MPI_Bcast(receive_displacements.data(), (int)receive_displacements.size(), MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(receive_counts.data(), (int)receive_counts.size(), MPI_INT, 0, MPI_COMM_WORLD);

    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Gatherv(vacancies.data(), (vacancies.rows()*4), MPI_INT, vacs_in.data(), receive_counts.data(), receive_displacements.data(), MPI_INT, 0, MPI_COMM_WORLD);
    //print_1Dvector(vacs_in);
    MPI_Barrier(MPI_COMM_WORLD);

    int shift_idx = 0;
    int coords_size = (int)(vacs_in.size()/4);
    
    if (rank == 0) {
        std::cout << "writing move_ticks: " << k << "\n";
        // Caught signal 11 (Segmentation fault: Sent by the kernel at address (nil)) when making FourDBoolArr
        FourDBoolArr vacancies_out((size_t)2, (size_t)dims[0], (size_t)dims[1], (size_t)dims[2]);
        vacancies_out.zero();
        
        for (int idx=0; idx<coords_size; idx++) {
            
            //std::cout << "4*idx: "  << 4*idx << " shift_idx: " << shift_idx << "\n";
            if ((shift_idx < nprocs) && (idx*4 == sum_vec[shift_idx])) { shift_idx++; }
        
            x_idx = shift_idx % xprocs;
            y_idx = floor(shift_idx / xprocs);
            
            x_chunk_start = (int)( dims[0]/xprocs * x_idx );
            y_chunk_start = (int)( dims[1]/yprocs * y_idx ); 
            

            i1 = vacs_in[4*idx];
            i2 = vacs_in[4*idx+1] + x_chunk_start;
            i3 = vacs_in[4*idx+2] + y_chunk_start;
            i4 = vacs_in[4*idx+3];
            //std::cout << "i1: "  << i1 << " i2: "  << i2 << " i3: "  << i3 << " i4: "  << i4 << "\n";
            vacancies_out(i1,i2,i3,i4) = 1;
        }

        Matrix<int> all_vacancies = vacancies_out.nonzero();
        std::string output_filename = ss.str();
        ss << folder << "/vacs/vacancies_output_" << i << "_" << l << "_" << k << "_" << t << "_moves.txt";
        std::cout << "filename: " << ss.str() << "\n";
        write_to_file(ss.str(), all_vacancies);
        ss.str("");
        ss.clear();
    }
    
    MPI_Barrier(MPI_COMM_WORLD); // all the sub ranks/processes waits here
}

int sum_vacs_allprocs(const Matrix<int>& vacancies, std::vector<int> dims, int xprocs, int yprocs, int nprocs, int rank) {
    //std::cout << "rank: " << rank << " start write parallel \n";
    int size = (int)( vacancies.rows() * vacancies.cols());

    std::vector<int> output(size);
    int idx = 0;
    int size2 = 0; int size3 = 0; int size_in;
    int i1; int i2; int i3; int i4; int x_idx; int y_idx; int x_chunk_start; int y_chunk_start;
    std::vector<int> vacs_in;

    MPI_Status status;
    MPI_Request request;
    std::ostringstream ss;

    std::vector<int> nums_and_proc(2);
    std::vector<int> num_proc_buffer(2);
    std::vector<int> receive_counts(nprocs);
    std::vector<int> receive_displacements(nprocs,0);
    std::vector<int> sum_vec(nprocs);

    int num_elems = (int)(vacancies.rows() * 4);
    nums_and_proc[0] = num_elems;
    nums_and_proc[1] = rank;

    int sum_of_elems = 0;
    
    if (rank != 0) {
        MPI_Isend(
            nums_and_proc.data(),      //Address of the message we are receiving.
            2,                  //Number of elements handled by that address.
            MPI_INT,            //MPI_TYPE of the message we are sending.
            0,                  //Rank of receiving process
            4,                  //Message Tag
            MPI_COMM_WORLD,      //MPI Communicator
            &request 
        );
    }

    MPI_Barrier(MPI_COMM_WORLD); // all the sub ranks/processes waits here

    if (rank == 0) {

        receive_counts[0] = num_elems;
        // Caught signal 11 (Segmentation fault: Sent by the kernel at address (nil)) at rank = 0
        for (int rec_proc=0; rec_proc<nprocs; rec_proc++) {
            if (rec_proc != rank) {

                MPI_Recv(
                    num_proc_buffer.data(),      //Address of the message we are receiving.
                    2,                  //Number of elements handled by that address.
                    MPI_INT,            //MPI_TYPE of the message we are sending.
                    rec_proc,           //Rank of sending process
                    4,                  //Message Tag
                    MPI_COMM_WORLD,      //MPI Communicator
                    &status 
                ); 

                receive_counts[num_proc_buffer[1]] = num_proc_buffer[0];
            }              
        }
        
        for (int sum_i=0; sum_i<(int)receive_counts.size(); sum_i++) {
            if (sum_i==0) {
                receive_displacements[sum_i] = 0;
                sum_vec[sum_i] = receive_counts[sum_i]; 
            }
            else {
                receive_displacements[sum_i] = receive_counts[sum_i-1] + receive_displacements[sum_i-1];
                sum_vec[sum_i] = receive_counts[sum_i] + sum_vec[sum_i-1];
            }
        }

        for (auto& n : receive_counts) 
            sum_of_elems += n;

    }
    vacs_in.resize(sum_of_elems);

    MPI_Barrier(MPI_COMM_WORLD);

    MPI_Bcast(receive_displacements.data(), (int)receive_displacements.size(), MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(receive_counts.data(), (int)receive_counts.size(), MPI_INT, 0, MPI_COMM_WORLD);

    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Gatherv(vacancies.data(), (vacancies.rows()*4), MPI_INT, vacs_in.data(), receive_counts.data(), receive_displacements.data(), MPI_INT, 0, MPI_COMM_WORLD);
    //print_1Dvector(vacs_in);
    MPI_Barrier(MPI_COMM_WORLD);

    int shift_idx = 0;
    int coords_size = (int)(vacs_in.size()/4);
    
    if (rank == 0) {
        // Caught signal 11 (Segmentation fault: Sent by the kernel at address (nil)) when making FourDBoolArr
        FourDBoolArr vacancies_out((size_t)2, (size_t)dims[0], (size_t)dims[1], (size_t)dims[2]);
        vacancies_out.zero();
        
        for (int idx=0; idx<coords_size; idx++) {
            
            //std::cout << "4*idx: "  << 4*idx << " shift_idx: " << shift_idx << "\n";
            if ((shift_idx < nprocs) && (idx*4 == sum_vec[shift_idx])) { shift_idx++; }
        
            x_idx = shift_idx % xprocs;
            y_idx = floor(shift_idx / xprocs);
            
            x_chunk_start = (int)( dims[0]/xprocs * x_idx );
            y_chunk_start = (int)( dims[1]/yprocs * y_idx ); 
            

            i1 = vacs_in[4*idx];
            i2 = vacs_in[4*idx+1] + x_chunk_start;
            i3 = vacs_in[4*idx+2] + y_chunk_start;
            i4 = vacs_in[4*idx+3];
            //std::cout << "i1: "  << i1 << " i2: "  << i2 << " i3: "  << i3 << " i4: "  << i4 << "\n";
            vacancies_out(i1,i2,i3,i4) = 1;
        }

        Matrix<int> all_vacancies = vacancies_out.nonzero();
        return (int)all_vacancies.rows();
    }
    
}