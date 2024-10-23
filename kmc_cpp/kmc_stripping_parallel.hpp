#include <iostream>
#include <sstream>
#include <stdio.h>
#include <assert.h> 
#include <array> 
#include <random>
#include <string>
#include <vector>
#include <cstdlib>
#include <tuple>        
#include <numeric>  
#include <algorithm>   
#include <fstream>
#include <map>
#include <cmath>
#include <math.h>
#include <chrono>
#include <execution>
#include <cstdint>
#include <mpi.h>

#include "hpp_files/math_func.hpp"
#include "hpp_files/vec_func.hpp"



/*------------------------------------------------------------------------------------*/
 /*! \brief A class for storing a simulation cell of atoms and propogating moves around the 
 lattice */
class Lattice {

    std::vector< std::vector<int> > diag_directions;
    std::vector< std::vector<int> > edge_directions;

    std::map<int, double> rate_typedict;

    std::vector<bool> periodic;

    std::vector<int> ver_tuples_cumsum;
    std::vector<int> bc_tuples_cumsum;
    std::vector<int> ver_edge_tuples_cumsum; 
    std::vector<int> bc_edge_tuples_cumsum;

    public:

        std::vector<Region*> regions;
        std::vector<int> sublattice_dim;
        std::vector<int> total_dims;
        std::map<int, std::string> a_types;
        FourDArr vertex_sites;
        FourDBoolArr vacancies;
        FourDArr bc_sites;
        FourDArr region_sites;
        int rank;
        int num_procs;
        double t;
        FourDBoolArr proc_neg_x_neighbors;
        FourDBoolArr proc_neg_y_neighbors;
        FourDBoolArr proc_pos_x_neighbors;
        FourDBoolArr proc_pos_y_neighbors;
        Matrix<int> proc_neighbors;
        std::vector<int> proc_dims;
        std::vector<double> probs;
        std::vector<double> rates;
        std::vector<double> rate_cumsum;
        std::vector< std::vector< std::vector< std::vector<double> > > > region_energies;
        Matrix<int> moves_coords; 
        Matrix<int> moves_shifts;
        Matrix<int> moves_lattice;
        Matrix<int> moves_vacs;
        Matrix<double> ratecatalog_111; 
        Matrix<double> ratecatalog_100;
        Matrix<double> regionrates_111_L;
        Matrix<double> regionrates_111_R;
        Matrix<double> regionrates_100_L;
        Matrix<double> regionrates_100_R;
        Matrix<int> configs_111;
        Matrix<int> configs_100;
        Matrix<int> vacancies_pos;
        int num_of_moves;
        std::mt19937 mt_obj;
        std::uniform_int_distribution<std::mt19937::result_type> x_rand;
        std::uniform_int_distribution<std::mt19937::result_type> y_rand;
        double temp;
        int watch_var;
        int num_of_vacs;
        std::vector< std::vector<int> > chunk_bounds;
        std::vector<double> prev_times;
        std::vector< std::vector<int> > prev_moves;
        std::vector< std::vector<int> > prev_newlocs;
        std::vector< std::vector<int> > prev_oldlocs;
        std::vector<int> prev_lattice;
        std::vector<int> prev_idxs;
        std::vector< std::vector<int> > par_prev_moves;
        std::vector<int> par_move_ticks;
        std::vector<int> par_prev_idx;
        std::vector<int> prev_move_type;
        std::vector<int> prev_move_type_ticks;
        int ghost_done_tag;
        int par_done_tag;
        int conflict_done_flag;
        int void_threshold;
        double void_barrier;
            

        Lattice(int xdim, int ydim, int zdim, int num_vacancies, int num_regions, int number_procs, int x_neigh, int y_neigh, int xtot, int ytot, int ztot):
            vertex_sites((size_t)1, (size_t)xdim, (size_t)ydim, (size_t)zdim),
            vacancies((size_t)2, (size_t)xdim, (size_t)ydim, (size_t)zdim),
            bc_sites((size_t)1, (size_t)xdim, (size_t)ydim, (size_t)zdim),
            region_sites((size_t)2, (size_t)xdim, (size_t)ydim, (size_t)zdim),
            num_procs((size_t)number_procs),
            proc_neg_x_neighbors((size_t)1, (size_t)2, (size_t)y_neigh, (size_t)zdim), 
            proc_neg_y_neighbors((size_t)1, (size_t)2, (size_t)x_neigh, (size_t)zdim), 
            proc_pos_x_neighbors((size_t)1, (size_t)2, (size_t)y_neigh, (size_t)zdim), 
            proc_pos_y_neighbors((size_t)1, (size_t)2, (size_t)x_neigh, (size_t)zdim),  
            proc_neighbors((size_t)number_procs, (size_t)8),
            moves_coords((size_t)(14 * num_vacancies), 4),
            moves_shifts((size_t)(14 * num_vacancies), 3),
            moves_lattice((size_t)(14 * num_vacancies), 1),
            moves_vacs((size_t)(14 * num_vacancies), 1),
            ratecatalog_111(2, (size_t)exp_int(2,8)),
            ratecatalog_100(2, (size_t)exp_int(2,14)),
            regionrates_111_L((size_t)num_regions, (size_t)exp_int(2,8)),
            regionrates_111_R((size_t)num_regions, (size_t)exp_int(2,8)),
            regionrates_100_L((size_t)num_regions, (size_t)exp_int(2,14)),
            regionrates_100_R((size_t)num_regions, (size_t)exp_int(2,14)),
            configs_111(1, (size_t)exp_int(2,8)),
            configs_100(1, (size_t)exp_int(2,14)),
            vacancies_pos((size_t)num_vacancies, 4),
            mt_obj((unsigned int)(std::chrono::high_resolution_clock::now().time_since_epoch().count())),
            //mt_obj((unsigned int)725863834569),
            x_rand(0, (size_t)xdim),
            y_rand(0, (size_t)ydim)

            {
                diag_directions = {{0,0,0}, {1,0,0}, {0,1,0}, {1,1,0}, {0,0,1}, {1,0,1}, {0,1,1}, {1,1,1}};
                edge_directions = {{0,0,1}, {-1,0,0}, {0,-1,0}, {0,1,0}, {1,0,0}, {0,0,-1}};
                sublattice_dim = {xdim,ydim,zdim};
                total_dims = {xtot,ytot,ztot};
                t = 0;
                num_of_moves = 0;
                num_of_vacs = num_vacancies;
                ghost_done_tag = 10;
                par_done_tag = 11;
                conflict_done_flag = 12;
            }

        /**
        * @brief Subroutine for adding move information to matrices stored as attributes of Lattice structure.
        * 
        * This function calculates the new coordinates, shifts, and lattice for a move, then stores
        * this information in the relevant arrays. It also computes the rate for the move using a rate 
        * constant function and updates the cumulative rate.
        * 
        * @param i            Initial i coordinate.
        * @param j            Initial j coordinate.
        * @param k            Initial k coordinate.
        * @param l            Initial l coordinate.
        * @param curr_move_num Current move number to be processed.
        * @param direc_sign   Direction sign for the move.
        * @param s            Index of direction in diag_directions/edge_directions.
        * @param idx          Index of vacancy in vacancies_pos Matrix.
        * @param lattice      Integer corresponding to move type .
        * @return int         Updated move number after processing.
        */
        int add_move(int i, int j, int k, int l, int curr_move_num, int direc_sign, int s, int idx, int lattice) {
            double rate = 0; 

            if ((lattice == 2) || (lattice == 3)) {
                moves_coords(curr_move_num,0) = i;
                moves_coords(curr_move_num,1) = (((j + edge_directions[s][0]) % sublattice_dim[0] + sublattice_dim[0]) % sublattice_dim[0]);
                moves_coords(curr_move_num,2) = (((k + edge_directions[s][1]) % sublattice_dim[1] + sublattice_dim[1]) % sublattice_dim[1]);
                moves_coords(curr_move_num,3) = (((l + edge_directions[s][2]) % sublattice_dim[2] + sublattice_dim[2]) % sublattice_dim[2]);                           
                moves_shifts(curr_move_num,0) = edge_directions[s][0];
                moves_shifts(curr_move_num,1) = edge_directions[s][1];
                moves_shifts(curr_move_num,2) = edge_directions[s][2];
                moves_lattice(curr_move_num,0) = lattice;
                moves_vacs(curr_move_num,0) = idx;

            } else if ((lattice == 0) || (lattice == 1)) {
                if (lattice == 0) { moves_coords(curr_move_num,0) = 1; }
                else if (lattice == 1) { moves_coords(curr_move_num,0) = 0; }
                moves_coords(curr_move_num,1) = (((j + direc_sign * diag_directions[s][0]) % sublattice_dim[0] + sublattice_dim[0]) % sublattice_dim[0]);
                moves_coords(curr_move_num,2) = (((k + direc_sign * diag_directions[s][1]) % sublattice_dim[1] + sublattice_dim[1]) % sublattice_dim[1]);
                moves_coords(curr_move_num,3) = (((l + direc_sign * diag_directions[s][2]) % sublattice_dim[2] + sublattice_dim[2]) % sublattice_dim[2]);                                    
                moves_shifts(curr_move_num,0) = direc_sign * diag_directions[s][0];
                moves_shifts(curr_move_num,1) = direc_sign * diag_directions[s][1];
                moves_shifts(curr_move_num,2) = direc_sign * diag_directions[s][2]; 
                moves_lattice(curr_move_num,0) = lattice;
                moves_vacs(curr_move_num,0) = idx;
            }

            // getting rate corresponding to move
            rate = new_get_rateconstants(moves_coords[curr_move_num], moves_shifts[curr_move_num], moves_lattice(curr_move_num,0));

            if (rate == -1) {curr_move_num --;} // case for when move not available in rate catalog
            else {
                if (curr_move_num == 0) {rate_cumsum[curr_move_num] = rate;}
                else {rate_cumsum[curr_move_num] = rate + rate_cumsum[curr_move_num-1];}
            }

            curr_move_num++;

            return curr_move_num;
        }

        /**
        * @brief Routine to check if an adjacent site is unoccupied for a move.
        * 
        * This function checks whether the target position in the lattice is free by verifying
        * if the corresponding site in the vacancy matrix is occupied.
        * 
        * @param i           Initial i coordinate.
        * @param j           Initial j coordinate.
        * @param k           Initial k coordinate.
        * @param l           Initial l coordinate.
        * @param direc_sign  Direction sign for the move.
        * @param s           Index of direction in diag_directions/edge_directions.
        * @param lattice     Integer corresponding to move type (0, 1, 2, or 3).
        * @return bool       True if the site is has no vacancy, false if occupied with vacancy.
        */
        bool check_move_free(int i, int j, int k, int l, int direc_sign, int s, int lattice) {
            int i1; int i2; int i3; int i4;
            if ((lattice == 2) || (lattice == 3)) {
                i1 = i;
                i2 = (((j + edge_directions[s][0]) % sublattice_dim[0] + sublattice_dim[0]) % sublattice_dim[0]);
                i3 = (((k + edge_directions[s][1]) % sublattice_dim[1] + sublattice_dim[1]) % sublattice_dim[1]);
                i4 = (((l + edge_directions[s][2]) % sublattice_dim[2] + sublattice_dim[2]) % sublattice_dim[2]); 
            } else if ((lattice == 0) || (lattice == 1)) {
                if (lattice == 0) { i1 = 1; }
                else if (lattice == 1) { i1 = 0; }
                i2 = (((j + direc_sign * diag_directions[s][0]) % sublattice_dim[0] + sublattice_dim[0]) % sublattice_dim[0]);
                i3 = (((k + direc_sign * diag_directions[s][1]) % sublattice_dim[1] + sublattice_dim[1]) % sublattice_dim[1]);
                i4 = (((l + direc_sign * diag_directions[s][2]) % sublattice_dim[2] + sublattice_dim[2]) % sublattice_dim[2]);
            }
            if (vacancies(i1,i2,i3,i4)) return false;

            return true;
        }

        
        /** 
        * @brief Routine to count the number of vacancies in the nearest neighbor (NN) shell.
        * 
        * This function calculates and returns the number of neighboring vacancies for a given
        * lattice site, taking into account the periodic boundary conditions along all directions.
        * It checks various diagonal and edge directions and communicates with neighboring processors
        * if the vacancy is located at the boundary.
        * 
        * @param idx The index of the vacancy whose neighbors are to be counted.
        * @return int The number of nearest neighbor vacancies.
        */
        int get_NN_count(int idx) {
            int i = vacancies_pos(idx,0); int j = vacancies_pos(idx,1); int k = vacancies_pos(idx,2); int l = vacancies_pos(idx,3);
            int new_j=0; int new_k=0; int new_l=0;
            int NN_count = 0;

            std::vector<size_t> x_dims = proc_pos_x_neighbors.size_vec;
            std::vector<size_t> y_dims = proc_pos_y_neighbors.size_vec;

            for (int s=0; s < (int)diag_directions.size(); s++) {

                if (i == 0) {
                    new_j = (((j - diag_directions[s][0]) % sublattice_dim[0] + sublattice_dim[0]) % sublattice_dim[0]);
                    new_k = (((k - diag_directions[s][1]) % sublattice_dim[1] + sublattice_dim[1]) % sublattice_dim[1]);
                    new_l = (((l - diag_directions[s][2]) % sublattice_dim[2] + sublattice_dim[2]) % sublattice_dim[2]);
                }
                else if (i == 1) {
                    new_j = (((j + diag_directions[s][0]) % sublattice_dim[0] + sublattice_dim[0]) % sublattice_dim[0]);
                    new_k = (((k + diag_directions[s][1]) % sublattice_dim[1] + sublattice_dim[1]) % sublattice_dim[1]);
                    new_l = (((l + diag_directions[s][2]) % sublattice_dim[2] + sublattice_dim[2]) % sublattice_dim[2]);
                }

                if ((l == 0) && (i == 0) && (diag_directions[s][2] == 1)) {/* checking for leftmost non-periodic boundary along z-axis*/}
                
                else if ((l == (int)(sublattice_dim[2]-1)) && (i == 1) && (diag_directions[s][2] == 1)) {/* checking for rightmost non-periodic boundary along z-axis*/}
                
                else if ((i == 0) && (j == 0) && (diag_directions[s][0] == 1)) {/*communicate with proc to -x direction*/
                    
                    if ((k == 0) && (diag_directions[s][1] == 1)) {
                        /*check neighbor -x,-y array */
                        if ((proc_neighbors(rank,5) == rank) && (!check_move_free(i,j,k,l,-1,s,0))) {NN_count ++;}
                        else if ((proc_neighbors(rank,5) != rank) && (proc_neg_x_neighbors(0, 1, 0, (size_t)(new_l)))) {NN_count ++;}
                    }
                    else {
                        /*check neighbor -x array */
                        if ((proc_neighbors(rank,4) == rank) && (!check_move_free(i,j,k,l,-1,s,0))) {NN_count ++;}
                        else if ((proc_neighbors(rank,4) != rank) && (proc_neg_x_neighbors(0, 1, (size_t)(new_k+1), (size_t)(new_l)))) {NN_count ++;}
                    }
                }
                
                else if ((i == 1) && (j == (sublattice_dim[0] - 1)) && (diag_directions[s][0] == 1)) {/*communicate with proc to +x direction*/ 
                    
                    if ((k == ((sublattice_dim[1] - 1))) && (diag_directions[s][1] == 1)) {
                        /*check neighbor +x,+y array */
                        if (((proc_neighbors(rank,1) == rank)) && (!check_move_free(i,j,k,l,1,s,1))) {NN_count ++;}
                        else if (((proc_neighbors(rank,1) != rank)) && (proc_pos_x_neighbors(0, 0, (size_t)(x_dims[2]-1), (size_t)(new_l)))) {NN_count ++;}
                    }
                    else {
                        /*check neighbor +x array */
                        if (((proc_neighbors(rank,0) == rank)) && (!check_move_free(i,j,k,l,1,s,1))) {NN_count ++;}
                        else if (((proc_neighbors(rank,0) != rank)) && (proc_pos_x_neighbors(0, 0, (size_t)(new_k), (size_t)(new_l)))) {NN_count ++;}
                    }
                }

                else if ((i == 0) && (k == 0) && (diag_directions[s][1] == 1)) {/*communicate with proc to -y direction*/
                    
                    if ((j == 0) && (diag_directions[s][0] == 1)) {
                        /*check neighbor -x,-y array */
                        if ((proc_neighbors(rank,5) == rank) && (!check_move_free(i,j,k,l,-1,s,0))) {NN_count ++;}
                        else if ((proc_neighbors(rank,5) != rank) && (proc_neg_y_neighbors(0, 1, 0, (size_t)(new_l)))) {NN_count ++;}
                    }
                    else if (((proc_neighbors(rank,6) == rank)) && (!check_move_free(i,j,k,l,-1,s,0))) {NN_count ++;}
                    else if (((proc_neighbors(rank,6) != rank)) && (proc_neg_y_neighbors(0, 1, (size_t)(new_j+1), (size_t)(new_l)))) {NN_count ++;}
                }

                else if ((i == 1) && (k == (sublattice_dim[1] - 1)) && (diag_directions[s][1] == 1)) {/*communicate with proc to +y direction*/
                    
                    if ((j == ((sublattice_dim[0] - 1))) && (diag_directions[s][0] == 1)) {
                        /*check neighbor +x,+y array */
                        if ((proc_neighbors(rank,1) == rank) && (!check_move_free(i,j,k,l,1,s,1))) {NN_count ++;}
                        else if (((proc_neighbors(rank,1) != rank)) && (proc_pos_y_neighbors(0, 0, (size_t)(y_dims[2]-1), (size_t)(new_l)))) {NN_count ++;}
                    }
                    else {
                        /*check neighbor +y array */
                        if ((proc_neighbors(rank,2) == rank) && (!check_move_free(i,j,k,l,1,s,1))) {NN_count ++;}
                        else if ((proc_neighbors(rank,2) != rank) && (proc_pos_y_neighbors(0, 0, (size_t)(new_j), (size_t)(new_l)))) {NN_count ++;}
                    }
                }

                else {
                    if ((i == 0) && (vacancies(1, (((j - diag_directions[s][0]) % sublattice_dim[0] + sublattice_dim[0]) % sublattice_dim[0]), (((k - diag_directions[s][1]) % sublattice_dim[1] + sublattice_dim[1]) % sublattice_dim[1]), (((l - diag_directions[s][2]) % sublattice_dim[2] + sublattice_dim[2]) % sublattice_dim[2])) == 0)) {
                        // checking that vertex site -> bc site move has new site occupied by atom
                        NN_count ++;
                    }
                    else if ((i == 1) && (vacancies(0, (((j + diag_directions[s][0]) % sublattice_dim[0] + sublattice_dim[0]) % sublattice_dim[0]), (((k + diag_directions[s][1]) % sublattice_dim[1] + sublattice_dim[1]) % sublattice_dim[1]), (((l + diag_directions[s][2]) % sublattice_dim[2] + sublattice_dim[2]) % sublattice_dim[2])) == 0)) {
                        // checking that bc site -> vertex site move has new site occupied by atom
                        NN_count ++;
                    }
                }

            }

            return NN_count;
        }


        /** 
        * @brief Find actions in the subdomain assigned to the current processor.
        * 
        * This function loops over all vacancies in the system and identifies possible moves
        * based on the nearest neighbor and edge directions. It performs communication between
        * processors when necessary, handles boundary conditions, and resizes data structures 
        * when required to store new moves.
        * 
        * The function also updates the cumulative sum of move rates.
        */
        void parallel_get_actions() {
            //std::cout << "rank: " << rank << " enter parallel get actions \n";
            int curr_move_num = 0; // total number of moves at this current timestep  
            int vacs_on_interface = 0; // vacancies at last z-index of lattice (used to calculate rate for stripping)
            int num_interface_sites = sublattice_dim[0] * sublattice_dim[1]; // total number of sites at last z-index 
            std::fill(rate_cumsum.begin(), rate_cumsum.end(), 0);


            rate_cumsum.resize((int)moves_coords.rows());
            int i=0; int j=0; int k=0; int l=0; 
            int new_j=0; int new_k=0; int new_l=0;

            std::vector<size_t> x_dims = proc_pos_x_neighbors.size_vec;
            std::vector<size_t> y_dims = proc_pos_y_neighbors.size_vec;

            // looping over all vacancies in system

            //std::cout << "rank: " << rank << " pre loop \n";
            for (int idx=0; idx < (int)vacancies_pos.rows(); idx++) {
                // position in lattice of vacancy 
                i = vacancies_pos(idx,0);
                j = vacancies_pos(idx,1);
                k = vacancies_pos(idx,2);
                l = vacancies_pos(idx,3);

                if ((curr_move_num  >= ((int)moves_shifts.rows() - 20))) {
                    // resizing data structures to accommodate all moves 

                    int newsize = 2 * ((int)moves_shifts.rows() + (int)vacancies_pos.rows() * 14);
                    //std::cout << "newsize: " << newsize << "\n";
                    rate_cumsum.resize(newsize);
                    moves_coords.reshape(newsize, 4, rank);
                    moves_shifts.reshape(newsize, 3, rank);
                    moves_lattice.reshape(newsize, 1, rank);
                    moves_vacs.reshape(newsize, 1, rank);

                    //std::cout << "rank: " << rank << " moves_coords.rows(): " << moves_coords.rows() << " moves_shifts.rows(): " << moves_shifts.rows() << " moves_lattice.rows(): " << moves_lattice.rows() << " moves_vacs.rows(): " << moves_vacs.rows() << "\n";
                    //std::cout << "vacancies_pos.rows(): " << vacancies_pos.rows() << "\n";
                }

                //NN_vac = get_NN_count(vacancies_pos[idx], i);

                // finding all moves along the {111} family of vectors
                for (int s=0; s < (int)diag_directions.size(); s++) {
                    if (i == 0) {
                        new_j = (((j - diag_directions[s][0]) % sublattice_dim[0] + sublattice_dim[0]) % sublattice_dim[0]);
                        new_k = (((k - diag_directions[s][1]) % sublattice_dim[1] + sublattice_dim[1]) % sublattice_dim[1]);
                        new_l = (((l - diag_directions[s][2]) % sublattice_dim[2] + sublattice_dim[2]) % sublattice_dim[2]);
                    }
                    
                    else if (i == 1) {
                        new_j = (((j + diag_directions[s][0]) % sublattice_dim[0] + sublattice_dim[0]) % sublattice_dim[0]);
                        new_k = (((k + diag_directions[s][1]) % sublattice_dim[1] + sublattice_dim[1]) % sublattice_dim[1]);
                        new_l = (((l + diag_directions[s][2]) % sublattice_dim[2] + sublattice_dim[2]) % sublattice_dim[2]);
                    }

                    if ((l == 0) && (i == 0) && (diag_directions[s][2] == 1)) {/* checking for leftmost non-periodic boundary along z-axis*/}
                    
                    else if ((l == (int)(sublattice_dim[2]-1)) && (i == 1) && (diag_directions[s][2] == 1)) {/* checking for rightmost non-periodic boundary along z-axis*/}
                    
                    else if ((i == 0) && (j == 0) && (diag_directions[s][0] == 1)) {/*communicate with proc to -x direction*/
                        
                        if ((k == 0) && (diag_directions[s][1] == 1)) {
                            /*check neighbor -x,-y array */
                            if ((proc_neighbors(rank,5) == rank) && check_move_free(i,j,k,l,-1,s,0)) {curr_move_num = add_move(i,j,k,l,curr_move_num,-1,s,idx,0);}
                            else if ((proc_neighbors(rank,5) != rank) && !(proc_neg_x_neighbors(0, 1, 0, (size_t)(new_l)))) {curr_move_num = add_move(i,j,k,l,curr_move_num,-1,s,idx,0);}
                        }
                        else {
                            /*check neighbor -x array */
                            if ((proc_neighbors(rank,4) == rank) && check_move_free(i,j,k,l,-1,s,0)) {curr_move_num = add_move(i,j,k,l,curr_move_num,-1,s,idx,0);}
                            else if ((proc_neighbors(rank,4) != rank) && !(proc_neg_x_neighbors(0, 1, (size_t)(new_k+1), (size_t)(new_l)))) {curr_move_num = add_move(i,j,k,l,curr_move_num,-1,s,idx,0);}
                        }
                    }
                    
                    else if ((i == 1) && (j == (sublattice_dim[0] - 1)) && (diag_directions[s][0] == 1)) {/*communicate with proc to +x direction*/ 
                        
                        if ((k == ((sublattice_dim[1] - 1))) && (diag_directions[s][1] == 1)) {
                            /*check neighbor +x,+y array */
                            if (((proc_neighbors(rank,1) == rank)) && check_move_free(i,j,k,l,1,s,1)) {curr_move_num = add_move(i,j,k,l,curr_move_num,1,s,idx,1);}
                            else if (((proc_neighbors(rank,1) != rank)) && !(proc_pos_x_neighbors(0, 0, (size_t)(x_dims[2]-1), (size_t)(new_l)))) {curr_move_num = add_move(i,j,k,l,curr_move_num,1,s,idx,1);}
                        }
                        else {
                            /*check neighbor +x array */
                            if (((proc_neighbors(rank,0) == rank)) && check_move_free(i,j,k,l,1,s,1)) {curr_move_num = add_move(i,j,k,l,curr_move_num,1,s,idx,1);}
                            else if (((proc_neighbors(rank,0) != rank)) && !(proc_pos_x_neighbors(0, 0, (size_t)(new_k), (size_t)(new_l)))) {curr_move_num = add_move(i,j,k,l,curr_move_num,1,s,idx,1);}
                        }
                    }

                    else if ((i == 0) && (k == 0) && (diag_directions[s][1] == 1)) {/*communicate with proc to -y direction*/
                        
                        if ((j == 0) && (diag_directions[s][0] == 1)) {
                            /*check neighbor -x,-y array */
                            if ((proc_neighbors(rank,5) == rank) && check_move_free(i,j,k,l,-1,s,0)) {curr_move_num = add_move(i,j,k,l,curr_move_num,-1,s,idx,0);}
                            else if ((proc_neighbors(rank,5) != rank) && !(proc_neg_y_neighbors(0, 1, 0, (size_t)(new_l)))) {curr_move_num = add_move(i,j,k,l,curr_move_num,-1,s,idx,0);}
                        }
                        else if (((proc_neighbors(rank,6) == rank)) && check_move_free(i,j,k,l,-1,s,0)) {curr_move_num = add_move(i,j,k,l,curr_move_num,-1,s,idx,0);}
                        else if (((proc_neighbors(rank,6) != rank)) && !(proc_neg_y_neighbors(0, 1, (size_t)(new_j+1), (size_t)(new_l)))) {curr_move_num = add_move(i,j,k,l,curr_move_num,-1,s,idx,0);}
                    }

                    else if ((i == 1) && (k == (sublattice_dim[1] - 1)) && (diag_directions[s][1] == 1)) {/*communicate with proc to +y direction*/
                        
                        if ((j == ((sublattice_dim[0] - 1))) && (diag_directions[s][0] == 1)) {
                            /*check neighbor +x,+y array */
                            if ((proc_neighbors(rank,1) == rank) && check_move_free(i,j,k,l,1,s,1)) {curr_move_num = add_move(i,j,k,l,curr_move_num,1,s,idx,1);}
                            else if (((proc_neighbors(rank,1) != rank)) && (!(proc_pos_y_neighbors(0, 0, (size_t)(y_dims[2]-1), (size_t)(new_l))))) {curr_move_num = add_move(i,j,k,l,curr_move_num,1,s,idx,1);}
                        }
                        else {
                            /*check neighbor +y array */
                            if ((proc_neighbors(rank,2) == rank) && check_move_free(i,j,k,l,1,s,1)) {curr_move_num = add_move(i,j,k,l,curr_move_num,1,s,idx,1);}
                            else if ((proc_neighbors(rank,2) != rank) && !(proc_pos_y_neighbors(0, 0, (size_t)(new_j), (size_t)(new_l)))) {curr_move_num = add_move(i,j,k,l,curr_move_num,1,s,idx,1);}
                        }
                    }

                    else {
                        if ((i == 0) && (vacancies(1, (((j - diag_directions[s][0]) % sublattice_dim[0] + sublattice_dim[0]) % sublattice_dim[0]), (((k - diag_directions[s][1]) % sublattice_dim[1] + sublattice_dim[1]) % sublattice_dim[1]), (((l - diag_directions[s][2]) % sublattice_dim[2] + sublattice_dim[2]) % sublattice_dim[2])) == 0)) {
                            // checking that vertex site -> bc site move has new site occupied by atom
                            curr_move_num = add_move(1,j,k,l,curr_move_num,-1,s,idx,0);
                        }
                        else if ((i == 1) && (vacancies(0, (((j + diag_directions[s][0]) % sublattice_dim[0] + sublattice_dim[0]) % sublattice_dim[0]), (((k + diag_directions[s][1]) % sublattice_dim[1] + sublattice_dim[1]) % sublattice_dim[1]), (((l + diag_directions[s][2]) % sublattice_dim[2] + sublattice_dim[2]) % sublattice_dim[2])) == 0)) {
                            // checking that bc site -> vertex site move has new site occupied by atom
                            curr_move_num = add_move(0,j,k,l,curr_move_num,1,s,idx,1);
                        }
                    }
                }

                // finding all moves along the {100} family of vectors
                for (int s=0; s < (int)edge_directions.size(); s++) {

                    new_j = (((j + edge_directions[s][0]) % sublattice_dim[0] + sublattice_dim[0]) % sublattice_dim[0]);
                    new_k = (((k + edge_directions[s][1]) % sublattice_dim[1] + sublattice_dim[1]) % sublattice_dim[1]);
                    new_l = (((l + edge_directions[s][2]) % sublattice_dim[2] + sublattice_dim[2]) % sublattice_dim[2]); 

                    if ((l == 0) && (edge_directions[s][2] == -1)) {}  

                    else if ((l == (int)(sublattice_dim[2]-1)) && (edge_directions[s][2] == 1)) {}

                    else if ((j == 0) && (edge_directions[s][0] == -1)) {/*communicate with proc to -x direction*/
                        if ((proc_neighbors(rank,4) == rank) && check_move_free(i,j,k,l,1,s,(i+2))) {curr_move_num = add_move(i,j,k,l,curr_move_num,1,s,idx,(i+2));}       
                        else if ((proc_neighbors(rank,4) != rank) && !(proc_neg_x_neighbors(0, i, (size_t)(new_k+1),(size_t)(new_l)))) {curr_move_num = add_move(i,j,k,l,curr_move_num,1,s,idx,(i+2));}                    
                    }
                    else if ((j == (sublattice_dim[0] - 1)) && (edge_directions[s][0] == 1))  {/*communicate with proc to +x direction*/
                        if ((proc_neighbors(rank,0) == rank) && check_move_free(i,j,k,l,1,s,(i+2))) {curr_move_num = add_move(i,j,k,l,curr_move_num,1,s,idx,(i+2));}
                        else if ((proc_neighbors(rank,0) != rank) && !(proc_pos_x_neighbors(0, i, (size_t)(new_k),(size_t)(new_l)))) {curr_move_num = add_move(i,j,k,l,curr_move_num,1,s,idx,(i+2));}                    
                    }
                    else if ((k == 0) && (edge_directions[s][1] == -1)) {/*communicate with proc to -y direction*/
                        if ((proc_neighbors(rank,6) == rank) && check_move_free(i,j,k,l,1,s,(i+2))) {curr_move_num = add_move(i,j,k,l,curr_move_num,1,s,idx,(i+2));}       
                        if ((proc_neighbors(rank,6) != rank) && !(proc_neg_y_neighbors(0, i, (size_t)(new_j+1),(size_t)(new_l)))) {curr_move_num = add_move(i,j,k,l,curr_move_num,1,s,idx,(i+2));}                    
                    }
                    else if ((k == (sublattice_dim[1] - 1)) && (edge_directions[s][1] == 1)) {/*communicate with proc to +y direction*/
                        if ((proc_neighbors(rank,2) == rank) && check_move_free(i,j,k,l,1,s,(i+2))) {curr_move_num = add_move(i,j,k,l,curr_move_num,1,s,idx,(i+2));}     
                        if ((proc_neighbors(rank,2) != rank) && !(proc_pos_y_neighbors(0, i, (size_t)(new_j),(size_t)(new_l)))) {curr_move_num = add_move(i,j,k,l,curr_move_num,1,s,idx,(i+2));}                    
                    }
                    else if (vacancies(i, (((j + edge_directions[s][0]) % sublattice_dim[0] + sublattice_dim[0]) % sublattice_dim[0]), 
                        (((k + edge_directions[s][1]) % sublattice_dim[1] + sublattice_dim[1]) % sublattice_dim[1]), 
                        (((l + edge_directions[s][2]) % sublattice_dim[2] + sublattice_dim[2]) % sublattice_dim[2])) == 0) {
                        // checking that vertex site -> vertex site or bc site -> bc site move has new site occupied by atom
                        curr_move_num = add_move(i,j,k,l,curr_move_num,1,s,idx,(i+2));
                    }
                }
            }

            // UPDATING SIZE OF DATA STRUCTURES CONTIANING COORDINATES AND RATES OF MOVES
            //std::cout << "rank: " << rank << " curr_move_num: " << curr_move_num << " moves_coords.rows(): " << moves_coords.rows() << " moves_shifts.rows(): " << moves_shifts.rows() << " moves_lattice.rows(): " << moves_lattice.rows() << " moves_vacs.rows(): " << moves_vacs.rows() << "\n";
            
            num_of_moves = curr_move_num;
            rate_cumsum.resize(num_of_moves);
            moves_vacs.reshape(num_of_moves, 1, rank);
            moves_coords.reshape(num_of_moves, 4, rank);
            moves_shifts.reshape(num_of_moves, 3, rank);
            moves_lattice.reshape(num_of_moves, 1, rank);

            //std::cout << "rank: " << rank << " end of get actions\n";
        }


        /**
        * @brief Determines the rate constant corresponding to the nearest neighbor (NN) encoding 
        *        and move type of a vacancy in the lattice.
        * 
        * This function computes the rate constant for a vacancy move based on the configuration 
        * of the lattice, the coordinates of the vacancy, and the direction of the move. It handles 
        * both vertex and boundary condition (bc) sites and checks for special regions with 
        * predefined rate constants.
        * 
        * @param coord A pointer to an array of integers representing the coordinates of the vacancy 
        *              in the lattice. The array has 4 elements: coord[0] to coord[3] are the lattice coordinates.
        * @param shift A pointer to an array of integers representing the direction of the move.
        *              The array has 3 elements corresponding to the movement along each axis.
        * @param lattice An integer representing the type of move:
        *                - 0: Moving vacancy from vertex site to bc site
        *                - 1: Moving vacancy from bc site to vertex site
        *                - 2 or 3: Moving vacancy from vertex site to vertex site, or bc site to bc site
        * 
        * @return A double representing the rate constant for the vacancy move. If the move is not allowed, 
        *         the function returns -1.
        * 
        * @throws std::exception If the lattice type is invalid, an exception is thrown.
        */
        double new_get_rateconstants(int* coord, int* shift, int lattice) {  
            
            double rate = -1; 
            int LR_idx;  // Index corresponding to the direction of movement in the lattice (left/right)
            int idx = 0;

            // DETERMINING DIRECTION OF MOVE //

            // Case 1: Moving vacancy from vertex site to bc site
            if (lattice == 0) {
                if (shift[2] == 0) {
                    LR_idx = 0;
                } else {
                    LR_idx = 1;
                }
            }

            // Case 2: Moving vacancy from bc site to vertex site
            else if (lattice == 1) {
                if (shift[2] == 1) {
                    LR_idx = 0;
                } else {
                    LR_idx = 1;
                }
            }

            // Case 3: Moving vacancy from vertex site to vertex site OR bc site to bc site
            else if ((lattice == 2) || (lattice == 3)) {
                if (shift[2] == 1) {
                    LR_idx = 0;
                } else {
                    LR_idx = 1;
                }
            }

            // Return -1 if move not allowed
            if (idx == -1) {
                return -1;
            }

            // Determine region ID of the site
            int reg_id = region_sites(coord[0], coord[1], coord[2], coord[3]);

            if (reg_id != 0) {
                // Case where move occurs in a region with specially defined rate constants
                if ((lattice == 0) || (lattice == 1)) {
                    if (LR_idx == 1) {
                        rate = regionrates_111_L[(reg_id-1)][idx];
                    } else if (LR_idx == 0) {
                        rate = regionrates_111_R[(reg_id-1)][idx];
                    }
                }
                else if ((lattice == 2) || (lattice == 3)) {
                    if (LR_idx == 1) {
                        rate = regionrates_100_L[(reg_id-1)][idx];
                    } else if (LR_idx == 0) {
                        rate = regionrates_100_R[(reg_id-1)][idx];
                    }
                }
                else {
                    std::string str_output = "Error: invalid lattice type in search_catalog()";
                    printf("%s", str_output.c_str());
                    throw std::exception();
                }

            } else {
                // Case where no pre-defined region exists, use bulk rate constants
                if ((lattice == 1) || (lattice == 0)) {
                    rate = ratecatalog_111[LR_idx][idx];
                } else if ((lattice == 2) || (lattice == 3)) {
                    rate = ratecatalog_100[LR_idx][idx];
                }
            }

            return rate;
        }


        /**
        * @brief Generates the encoding corresponding to the configuration of nearest neighbors
        *        for a vacancy in the bulk of the processor domain.
        * 
        * This function computes an integer corresponding to a binary encoding for configuration of 
        * nearest neighbors around the vacancy's position
        * 
        * @param vac A pointer to an array of integers representing the vacancy's position in the lattice. 
        *            The array has 4 elements: vac[0] is the type of site (0 or 1), and vac[1], vac[2], vac[3]
        *            represent the coordinates of the vacancy in the lattice.
        * @param lattice An integer representing the type of move:
        *                - 3: Moving vacancy from bc site to bc site
        *                - 2: Moving vacancy from vertex site to vertex site
        *                - 1: Moving vacancy from bc site to vertex site
        *                - 0: Moving vacancy from vertex site to bc site
        * 
        * @return An integer representing the sum of the encoding for nearest neighbor interactions.
        */
        int new_get_neighbors(int* vac, int lattice) {
            // Ensure that the input coordinates are valid
            assert((vac[0] == 0) || (vac[0] == 1));
            assert( ((vac[1] >= 0) && (vac[1] < sublattice_dim[0])) );
            assert( ((vac[2] >= 0) && (vac[2] < sublattice_dim[1])) );
            assert( ((vac[3] >= 0) && (vac[3] < sublattice_dim[2])) );

            int sum = 0;
            int m = (int)a_types.size(); // Number of atom types in the system

            // Case 1: Moving vacancy from bc site to bc site
            if (lattice == 3) {
                for (int i=0; i<(int)diag_directions.size(); i++) {
                    sum += exp_int(m,i) * vertex_sites(0, 
                        (((vac[0] - diag_directions[i][0]) % sublattice_dim[0] + sublattice_dim[0]) % sublattice_dim[0]),
                        (((vac[1] - diag_directions[i][1]) % sublattice_dim[1] + sublattice_dim[1]) % sublattice_dim[1]),
                        (((vac[2] - diag_directions[i][2]) % sublattice_dim[2] + sublattice_dim[2]) % sublattice_dim[2]));
                }
                for (int i=0; i<(int)edge_directions.size(); i++) {
                    sum += exp_int(m,i) * bc_sites(0, 
                        (((vac[0] + edge_directions[i][0]) % sublattice_dim[0] + sublattice_dim[0]) % sublattice_dim[0]),
                        (((vac[1] + edge_directions[i][1]) % sublattice_dim[1] + sublattice_dim[1]) % sublattice_dim[1]),
                        (((vac[2] + edge_directions[i][2]) % sublattice_dim[2] + sublattice_dim[2]) % sublattice_dim[2]));
                }
            }
            
            // Case 2: Moving vacancy from vertex site to vertex site
            else if (lattice == 2) {
                for (int i=0; i<(int)diag_directions.size(); i++) {
                    sum += exp_int(m,i) * bc_sites(0, 
                        (((vac[0] + diag_directions[i][0]) % sublattice_dim[0] + sublattice_dim[0]) % sublattice_dim[0]),
                        (((vac[1] + diag_directions[i][1]) % sublattice_dim[1] + sublattice_dim[1]) % sublattice_dim[1]),
                        (((vac[2] + diag_directions[i][2]) % sublattice_dim[2] + sublattice_dim[2]) % sublattice_dim[2]));
                }
                for (int i=0; i<(int)edge_directions.size(); i++) {
                    sum += exp_int(m,i) * vertex_sites(0, 
                        (((vac[0] + edge_directions[i][0]) % sublattice_dim[0] + sublattice_dim[0]) % sublattice_dim[0]),
                        (((vac[1] + edge_directions[i][1]) % sublattice_dim[1] + sublattice_dim[1]) % sublattice_dim[1]),
                        (((vac[2] + edge_directions[i][2]) % sublattice_dim[2] + sublattice_dim[2]) % sublattice_dim[2]));
                }
            }
            
            // Case 3: Moving vacancy from bc site to vertex site
            else if (lattice == 1) {
                for (int i=0; i<(int)diag_directions.size(); i++) {
                    sum += exp_int(m,i) * vertex_sites(0, 
                        (((vac[0] - diag_directions[i][0]) % sublattice_dim[0] + sublattice_dim[0]) % sublattice_dim[0]),
                        (((vac[1] - diag_directions[i][1]) % sublattice_dim[1] + sublattice_dim[1]) % sublattice_dim[1]),
                        (((vac[2] - diag_directions[i][2]) % sublattice_dim[2] + sublattice_dim[2]) % sublattice_dim[2]));
                }
            }
            
            // Case 4: Moving vacancy from vertex site to bc site
            else if (lattice == 0) {
                for (int i=0; i<(int)diag_directions.size(); i++) {
                    sum += exp_int(m,i) * bc_sites(0, 
                        (((vac[0] + diag_directions[i][0]) % sublattice_dim[0] + sublattice_dim[0]) % sublattice_dim[0]),
                        (((vac[1] + diag_directions[i][1]) % sublattice_dim[1] + sublattice_dim[1]) % sublattice_dim[1]),
                        (((vac[2] + diag_directions[i][2]) % sublattice_dim[2] + sublattice_dim[2]) % sublattice_dim[2]));
                }
            }

            return sum;
        }


        /**
        * @brief Checks if a move exceeds the bounds of the processor domain and communicates the 
        *        information to adjacent processors if necessary.
        *
        * This function verifies if a vacancy move crosses the boundary of the current processor domain. 
        * If the move exceeds the domain, the function sends the updated information to adjacent processors 
        * using MPI communication. The function handles different boundary conditions and directions of the move.
        *
        * @param i_old The old lattice layer index of the vacancy before the move.
        * @param j_old The old x-coordinate of the vacancy before the move.
        * @param k_old The old y-coordinate of the vacancy before the move.
        * @param l_old The old z-coordinate of the vacancy before the move.
        * @param move_idx The index of the move information .
        * @param new_loc A constant reference to a vector of integers representing the new coordinates 
        *                of the vacancy after the move. The vector contains the lattice layer indices, x, y, and z.
        * 
        * @return True if the move crosses the boundary of the processor domain and information is sent 
        *         to adjacent processors, otherwise False.
        * 
        * @note The function uses MPI non-blocking communication (MPI_Isend) to send data to neighboring 
        *       processors. It handles multiple ghost regions and different data structures for various 
        *       lattice configurations.
        */
        bool parallel_processes_check(int i_old, int j_old, int k_old, int l_old, int move_idx, const std::vector<int>& new_loc) {
            //
            // Assertions to ensure the validity of the old and new coordinates
            assert(((i_old == 0) || (i_old == 1)) && ((new_loc[0] == 0) || (new_loc[0] == 1)));
            assert(((j_old >= 0) && (j_old < sublattice_dim[0])) && ((new_loc[1] >= 0) && (new_loc[1] < sublattice_dim[0])));
            assert(((k_old >= 0) && (k_old < sublattice_dim[1])) && ((new_loc[2] >= 0) && (new_loc[2] < sublattice_dim[1])));
            assert(((l_old >= 0) && (l_old < sublattice_dim[2])) && ((new_loc[3] >= 0) && (new_loc[3] < sublattice_dim[2])));
            //

            int new_proc;
            int bufferlen = 9;

            std::vector<int> loc_buffer1(bufferlen); // Buffer for sending new vacancy location
            std::vector<int> loc_buffer2(bufferlen); // Additional buffer for multi-processor communication
            std::vector<int> loc_buffer3(bufferlen); // Additional buffer for multi-processor communication

            MPI_Request request1;
            MPI_Request request2;
            MPI_Request request3;

            int i = new_loc[0];
            int j = new_loc[1];
            int k = new_loc[2];
            int l = new_loc[3];

            // Printing the new location for debugging
            print_1Dvector(new_loc);

            // Dimensions of neighboring processors in x and y directions
            std::vector<size_t> x_dims = proc_pos_x_neighbors.size_vec;
            std::vector<size_t> y_dims = proc_pos_y_neighbors.size_vec;

            // (111) moves handling for boundary conditions
            if ((moves_lattice(move_idx, 0) == 0) || (moves_lattice(move_idx, 0) == 1)) {
                if ((j_old == 0) && (moves_shifts(move_idx, 0) == -1)) {
                    // Communicate with the processor in the -x direction
                    if ((k_old == 0) && (moves_shifts(move_idx, 1) == -1)) {
                        // Communicate with the processor in the (-x, -y) direction
                        if ((proc_neighbors(rank, 6) == proc_neighbors(rank, 5)) && (proc_neighbors(rank, 4) == proc_neighbors(rank, 5))) {
                            // No communication needed
                        } else if ((proc_neighbors(rank, 6) == proc_neighbors(rank, 5)) && (proc_neighbors(rank, 4) == rank)) {
                            new_proc = proc_neighbors(rank, 6);
                            for (int idx = 0; idx < new_loc.size(); idx++) {
                                loc_buffer1[idx] = new_loc[idx];
                            }
                            loc_buffer1[4] = rank;
                            loc_buffer1[5] = i_old;
                            loc_buffer1[6] = j_old;
                            loc_buffer1[7] = k_old;
                            loc_buffer1[8] = l_old;

                            // Send data to the neighboring processor
                            MPI_Isend(loc_buffer1.data(), bufferlen, MPI_INT, new_proc, 1, MPI_COMM_WORLD, &request1);
                        }
                        else if ((proc_neighbors(rank, 4) == proc_neighbors(rank, 5)) && (proc_neighbors(rank, 6) == rank)) {
                            new_proc = proc_neighbors(rank, 4);
                            for (int idx = 0; idx < new_loc.size(); idx++) {
                                loc_buffer1[idx] = new_loc[idx];
                            }
                            loc_buffer1[4] = rank;
                            loc_buffer1[5] = i_old;
                            loc_buffer1[6] = j_old;
                            loc_buffer1[7] = k_old;
                            loc_buffer1[8] = l_old;

                            // Send data to the neighboring processor
                            MPI_Isend(loc_buffer1.data(), bufferlen, MPI_INT, new_proc, 1, MPI_COMM_WORLD, &request1);
                        }
                        else {
                            // Send data to multiple neighboring processors
                            new_proc = proc_neighbors(rank, 6);
                            for (int idx = 0; idx < new_loc.size(); idx++) {
                                loc_buffer1[idx] = new_loc[idx];
                            }
                            loc_buffer1[4] = rank;
                            loc_buffer1[5] = i_old;
                            loc_buffer1[6] = j_old;
                            loc_buffer1[7] = k_old;
                            loc_buffer1[8] = l_old;

                            MPI_Isend(loc_buffer1.data(), bufferlen, MPI_INT, new_proc, 1, MPI_COMM_WORLD, &request1);

                            new_proc = proc_neighbors(rank, 5);
                            for (int idx = 0; idx < new_loc.size(); idx++) {
                                loc_buffer2[idx] = new_loc[idx];
                            }
                            loc_buffer2[4] = rank;
                            loc_buffer2[5] = i_old;
                            loc_buffer2[6] = j_old;
                            loc_buffer2[7] = k_old;
                            loc_buffer2[8] = l_old;

                            MPI_Isend(loc_buffer2.data(), bufferlen, MPI_INT, new_proc, 1, MPI_COMM_WORLD, &request2);

                            new_proc = proc_neighbors(rank, 4);
                            for (int idx = 0; idx < new_loc.size(); idx++) {
                                loc_buffer3[idx] = new_loc[idx];
                            }
                            loc_buffer3[4] = rank;
                            loc_buffer3[5] = i_old;
                            loc_buffer3[6] = j_old;
                            loc_buffer3[7] = k_old;
                            loc_buffer3[8] = l_old;

                            MPI_Isend(loc_buffer3.data(), bufferlen, MPI_INT, new_proc, 1, MPI_COMM_WORLD, &request3);
                        }
                        proc_neg_y_neighbors(0, i, (size_t)0, l) = 1;
                        proc_neg_x_neighbors(0, i, (size_t)0, l) = 1;
                    }
                    else {
                        // Communicate with the processor in the -x direction
                        new_proc = proc_neighbors(rank, 4);
                        if (new_proc == rank) {
                            return false;
                        } else {
                            for (int idx = 0; idx < new_loc.size(); idx++) {
                                loc_buffer1[idx] = new_loc[idx];
                            }
                            loc_buffer1[4] = rank;
                            loc_buffer1[5] = i_old;
                            loc_buffer1[6] = j_old;
                            loc_buffer1[7] = k_old;
                            loc_buffer1[8] = l_old;

                            MPI_Isend(loc_buffer1.data(), bufferlen, MPI_INT, new_proc, 1, MPI_COMM_WORLD, &request1);
                        }
                        proc_neg_x_neighbors(0, i, k+1, l) = 1;
                    }

                    return true;
                }
                else if ((j_old == (sublattice_dim[0] - 1)) && (moves_shifts(move_idx,0) == 1)) {/*communicate with proc to +x direction*/
                    if ((k_old == (sublattice_dim[1] - 1))  && (moves_shifts(move_idx,1) == 1)) {
                        //std::cout << "rank: " << rank << " communicate with proc to (+x,+y) direction \n\n";
                        /*communicate with proc to (+x,+y) direction*/

                        if ((proc_neighbors(rank,2) == proc_neighbors(rank,1)) && (proc_neighbors(rank,2) == proc_neighbors(rank,0))) {return false;}
                        
                        else if ((proc_neighbors(rank,2) == proc_neighbors(rank,1)) && (proc_neighbors(rank,0) == rank)) {
                            new_proc = proc_neighbors(rank,2);
                            for (int idx=0; idx<new_loc.size(); idx++) {loc_buffer1[idx] = new_loc[idx];}
                            loc_buffer1[4] = rank;
                            loc_buffer1[5] = i_old;
                            loc_buffer1[6] = j_old;
                            loc_buffer1[7] = k_old;
                            loc_buffer1[8] = l_old;
                            //std::cout << "rank: " << rank << " sending to proc: " << new_proc <<  "\n";
                            MPI_Isend(
                                loc_buffer1.data(),        //Address of the message we are sending.
                                bufferlen,                  //Number of elements handled by that address.
                                MPI_INT,            //MPI_TYPE of the message we are sending.
                                new_proc,           //Rank of receiving process
                                1,                  //Message Tag
                                MPI_COMM_WORLD,      //MPI Communicator
                                &request1 ); 
                            ////std::cout << "rank: " << rank << " sent to proc: " << new_proc <<  "\n";
                        }
                        else if ((proc_neighbors(rank,0) == proc_neighbors(rank,1)) && (proc_neighbors(rank,2) == rank)) {
                            new_proc = proc_neighbors(rank,0);
                            for (int idx=0; idx<new_loc.size(); idx++) {loc_buffer1[idx] = new_loc[idx];}
                            loc_buffer1[4] = rank;
                            loc_buffer1[5] = i_old;
                            loc_buffer1[6] = j_old;
                            loc_buffer1[7] = k_old;
                            loc_buffer1[8] = l_old;
                            //std::cout << "rank: " << rank << " sending to proc: " << new_proc <<  "\n";
                            MPI_Isend(
                                loc_buffer1.data(),        //Address of the message we are sending.
                                bufferlen,                  //Number of elements handled by that address.
                                MPI_INT,            //MPI_TYPE of the message we are sending.
                                new_proc,           //Rank of receiving process
                                1,                  //Message Tag
                                MPI_COMM_WORLD,      //MPI Communicator
                                &request1 ); 
                            //std::cout << "rank: " << rank << " sent to proc: " << new_proc <<  "\n";
                        }
                        else {
                            new_proc = proc_neighbors(rank,2);
                            for (int idx=0; idx<new_loc.size(); idx++) {loc_buffer1[idx] = new_loc[idx];}
                            loc_buffer1[4] = rank;
                            loc_buffer1[5] = i_old;
                            loc_buffer1[6] = j_old;
                            loc_buffer1[7] = k_old;
                            loc_buffer1[8] = l_old;
                            //std::cout << "rank: " << rank << " sending to proc: " << new_proc <<  "\n";
                            MPI_Isend(
                                loc_buffer1.data(),        //Address of the message we are sending.
                                bufferlen,                  //Number of elements handled by that address.
                                MPI_INT,            //MPI_TYPE of the message we are sending.
                                new_proc,           //Rank of receiving process
                                1,                  //Message Tag
                                MPI_COMM_WORLD,      //MPI Communicator
                                &request1 ); 
                            //std::cout << "rank: " << rank << " sent to proc: " << new_proc <<  "\n";
                            new_proc = proc_neighbors(rank,1);
                            for (int idx=0; idx<new_loc.size(); idx++) {loc_buffer2[idx] = new_loc[idx];}
                            loc_buffer2[4] = rank;
                            loc_buffer2[5] = i_old;
                            loc_buffer2[6] = j_old;
                            loc_buffer2[7] = k_old;
                            loc_buffer2[8] = l_old;
                            //std::cout << "rank: " << rank << " sending to proc: " << new_proc <<  "\n";
                            MPI_Isend(loc_buffer2.data(), bufferlen, MPI_INT, new_proc, 1, MPI_COMM_WORLD, &request2);
                            //std::cout << "rank: " << rank << " sent to proc: " << new_proc <<  "\n";
                            new_proc = proc_neighbors(rank,0);
                            for (int idx=0; idx<new_loc.size(); idx++) {loc_buffer3[idx] = new_loc[idx];}
                            loc_buffer3[4] = rank;
                            loc_buffer3[5] = i_old;
                            loc_buffer3[6] = j_old;
                            loc_buffer3[7] = k_old;
                            loc_buffer3[8] = l_old;
                            //std::cout << "rank: " << rank << " sending to proc: " << new_proc <<  "\n";
                            MPI_Isend(loc_buffer3.data(), bufferlen, MPI_INT, new_proc, 1, MPI_COMM_WORLD, &request3);
                            //std::cout << "rank: " << rank << " sent to proc: " << new_proc <<  "\n";
                        }

                        proc_pos_x_neighbors(0, i, (size_t)(x_dims[2]-1), (size_t)(l)) = 1;
                        proc_pos_y_neighbors(0, i, (size_t)(y_dims[2]-1), (size_t)(l)) = 1;
                    }
                    else {
                        //std::cout << "rank: " << rank << " communicate with proc to (+x) direction \n\n";
                        /*communicate with proc to (+x) direction*/
                        
                        new_proc = proc_neighbors(rank,0); // proc_neighbors indices follow: 0:+x, 1:(+x,+y), 2:+y, 3:(-x,+y), 4:-x, 5:(-x,-y), 6:-y, 7:(+x,-y) 
                        if (new_proc == rank) {return false;}
                        else {
                            for (int idx=0; idx<new_loc.size(); idx++) {loc_buffer1[idx] = new_loc[idx];}
                            loc_buffer1[4] = rank;
                            loc_buffer1[5] = i_old;
                            loc_buffer1[6] = j_old;
                            loc_buffer1[7] = k_old;
                            loc_buffer1[8] = l_old;
                            //std::cout << "rank: " << rank << " sending to proc: " << new_proc <<  "\n";
                            MPI_Isend(loc_buffer1.data(), bufferlen, MPI_INT, new_proc, 1, MPI_COMM_WORLD, &request1);
                            //std::cout << "rank: " << rank << " sent to proc: " << new_proc <<  "\n";
                        }
                        proc_pos_x_neighbors(0, i, (size_t)(k), (size_t)(l)) = 1;
                    }

                    return true;
                }
                else if ((k_old == 0) && (moves_shifts(move_idx,1) == -1)) {
                    //std::cout << "rank: " << rank << " communicate with proc to (-y) direction \n\n";
                    /*communicate with proc to -y direction*/ 
                    new_proc = proc_neighbors(rank,6); // proc_neighbors indices follow: 0:+x, 1:(+x,+y), 2:+y, 3:(-x,+y), 4:-x, 5:(-x,-y), 6:-y, 7:(+x,-y) 

                    if (new_proc == rank) { return false; }
                    else {
                        for (int idx=0; idx<new_loc.size(); idx++) {loc_buffer1[idx] = new_loc[idx];}
                        loc_buffer1[4] = rank;
                        loc_buffer1[5] = i_old;
                        loc_buffer1[6] = j_old;
                        loc_buffer1[7] = k_old;
                        loc_buffer1[8] = l_old; 
                        //std::cout << "rank: " << rank << " sending to proc: " << new_proc <<  "\n";
                        MPI_Isend(loc_buffer1.data(), bufferlen, MPI_INT, new_proc, 1, MPI_COMM_WORLD, &request1);
                        //std::cout << "rank: " << rank << " sent to proc: " << new_proc <<  "\n";
                    }
                    proc_neg_y_neighbors(0, i, (size_t)(j+1), (size_t)(l)) = 1;
                    return true;
                }
                else if ((k_old == (sublattice_dim[1] - 1)) && (moves_shifts(move_idx,1) == 1)) {
                    //std::cout << "rank: " << rank << " communicate with proc to (+y) direction \n\n";
                    /*communicate with proc to +y direction*/
                    new_proc = proc_neighbors(rank,2); // proc_neighbors indices follow: 0:+x, 1:(+x,+y), 2:+y, 3:(-x,+y), 4:-x, 5:(-x,-y), 6:-y, 7:(+x,-y) 

                    if (new_proc == rank) { return false; }
                    else {
                        for (int idx=0; idx<new_loc.size(); idx++) {loc_buffer1[idx] = new_loc[idx];}
                        loc_buffer1[4] = rank;
                        loc_buffer1[5] = i_old;
                        loc_buffer1[6] = j_old;
                        loc_buffer1[7] = k_old;
                        loc_buffer1[8] = l_old;
                        //std::cout << "rank: " << rank << " sending to proc: " << new_proc <<  "\n";
                        MPI_Isend(loc_buffer1.data(), bufferlen, MPI_INT, new_proc, 1, MPI_COMM_WORLD, &request1);
                        //std::cout << "rank: " << rank << " sent to proc: " << new_proc <<  "\n";
                    }

                    proc_pos_y_neighbors(0, i, (size_t)(j), (size_t)(l)) = 1;

                    return true;
                }
            }
            if ((moves_lattice(move_idx,0) == 2) || (moves_lattice(move_idx,0) == 3)) {
                if ((moves_shifts(move_idx,0) == 1) && (j_old == (sublattice_dim[0] - 1))) {
                    if (i_old == 0) proc_pos_x_neighbors(0, i, (size_t)(k+1), (size_t)(l)) = 1;
                    else proc_pos_x_neighbors(0, i, (size_t)(k+1), (size_t)(l)) = 1;

                    new_proc = proc_neighbors(rank,0);
                    if (new_proc == rank) { return false; }

                    for (int idx=0; idx<new_loc.size(); idx++) {loc_buffer1[idx] = new_loc[idx];}
                    loc_buffer1[4] = rank;
                    loc_buffer1[5] = i_old;
                    loc_buffer1[6] = j_old;
                    loc_buffer1[7] = k_old;
                    loc_buffer1[8] = l_old; 
                    //std::cout << "rank: " << rank << " sending to proc: " << new_proc <<  "\n";
                    MPI_Isend(loc_buffer1.data(), bufferlen, MPI_INT, new_proc, 1, MPI_COMM_WORLD, &request1);
                    //std::cout << "rank: " << rank << " sent to proc: " << new_proc <<  "\n";

                    return true;
                }
                else if ((moves_shifts(move_idx,0) == -1) && (j_old == 0)) {
                    if (i_old == 1) proc_neg_x_neighbors(0, i, (size_t)(k+1), (size_t)(l)) = 1;
                    else proc_neg_x_neighbors(0, i, (size_t)(k+1), (size_t)(l)) = 1;

                    new_proc = proc_neighbors(rank,4);
                    if (new_proc == rank) { return false; }

                    for (int idx=0; idx<new_loc.size(); idx++) {loc_buffer1[idx] = new_loc[idx];}
                    loc_buffer1[4] = rank;
                    loc_buffer1[5] = i_old;
                    loc_buffer1[6] = j_old;
                    loc_buffer1[7] = k_old;
                    loc_buffer1[8] = l_old; 
                    //std::cout << "rank: " << rank << " sending to proc: " << new_proc <<  "\n";
                    MPI_Isend(loc_buffer1.data(), bufferlen, MPI_INT, new_proc, 1, MPI_COMM_WORLD, &request1);
                    //std::cout << "rank: " << rank << " sent to proc: " << new_proc <<  "\n";

                    return true;
                }
                else if ((moves_shifts(move_idx,1) == 1) && (k_old == (sublattice_dim[1] - 1))) {
                    if (i_old == 1) proc_pos_y_neighbors(0, i, (size_t)(j+1), (size_t)(l)) = 1;
                    else proc_pos_y_neighbors(0, i, (size_t)(j+1), (size_t)(l)) = 1;

                    new_proc = proc_neighbors(rank,2);
                    if (new_proc == rank) { return false; }

                    for (int idx=0; idx<new_loc.size(); idx++) {loc_buffer1[idx] = new_loc[idx];}
                    loc_buffer1[4] = rank;
                    loc_buffer1[5] = i_old;
                    loc_buffer1[6] = j_old;
                    loc_buffer1[7] = k_old;
                    loc_buffer1[8] = l_old; 
                    //std::cout << "rank: " << rank << " sending to proc: " << new_proc <<  "\n";
                    MPI_Isend(loc_buffer1.data(), bufferlen, MPI_INT, new_proc, 1, MPI_COMM_WORLD, &request1);
                    //std::cout << "rank: " << rank << " sent to proc: " << new_proc <<  "\n";

                    return true;
                }
                else if ((moves_shifts(move_idx,1) == -1) && (k_old == 0)) {
                    if (i_old == 1) proc_neg_y_neighbors(0, i, (size_t)(j+1), (size_t)(l)) = 1;
                    else proc_neg_y_neighbors(0, i, (size_t)(j+1), (size_t)(l)) = 1;
                    
                    new_proc = proc_neighbors(rank,6);
                    if (new_proc == rank) { return false; }

                    for (int idx=0; idx<new_loc.size(); idx++) {loc_buffer1[idx] = new_loc[idx];}
                    loc_buffer1[4] = rank;
                    loc_buffer1[5] = i_old;
                    loc_buffer1[6] = j_old;
                    loc_buffer1[7] = k_old;
                    loc_buffer1[8] = l_old; 
                    //std::cout << "rank: " << rank << " sending to proc: " << new_proc <<  "\n";
                    MPI_Isend(loc_buffer1.data(), bufferlen, MPI_INT, new_proc, 1, MPI_COMM_WORLD, &request1);
                    //std::cout << "rank: " << rank << " sent to proc: " << new_proc <<  "\n";

                    return true;
                }
            }
            
            return false;
        }


        /**
        * @brief Sends information about ghost sites to other processors.
        *
        * This method communicates the status of ghost sites to neighboring processors.
        * It performs checks to ensure that the old and new locations are within valid ranges.
        * The data includes the indices of the ghost sites and the rank of the sending processor.
        * Depending on the values of the input parameters, it sends the data to the appropriate neighboring processor
        * using MPI_Isend for parallel transfer.
        *
        * @param i_old Old lattice coordinate (vertex/bc site) of site in ghost zone.
        * @param j_old Old x-coordinate of the ghost site.
        * @param k_old Old y-coordinate of the ghost site.
        * @param l_old Old z-coordinate of the ghost site.
        * @param idx Index of the ghost site (used for data management).
        * @param new_loc Vector containing the new locations of the ghost sites.
        * @param parallel_transfer Boolean flag indicating if the transfer is parallel.
        */
        void ghost_site_send(int i_old, int j_old, int k_old, int l_old, int idx, const std::vector<int>& new_loc, bool parallel_transfer) {         
            //
            assert( ((i_old == 0) || (i_old == 1)) && ((new_loc[0] == 0) || (new_loc[0] == 1)) );
            assert( ((j_old >= 0) && (j_old < sublattice_dim[0])) && ((new_loc[1] >= 0) && (new_loc[1] < sublattice_dim[0])) );
            assert( ((k_old >= 0) && (k_old < sublattice_dim[1])) && ((new_loc[2] >= 0) && (new_loc[2] < sublattice_dim[1])) );
            assert( ((l_old >= 0) && (l_old < sublattice_dim[2])) && ((new_loc[3] >= 0) && (new_loc[3] < sublattice_dim[2])) );
            //

            int new_proc;
            int bufferlen = 9;
            int side_idx = -1;
            
            std::vector<int> loc_buffer1(bufferlen); // new location of vacancy
            std::vector<int> loc_buffer2(bufferlen); // new location of vacancy
            std::vector<int> loc_buffer3(bufferlen); // new location of vacancy

            MPI_Request request1;
            MPI_Request request2;
            MPI_Request request3;

            int i = new_loc[0];
            int j = new_loc[1];
            int k = new_loc[2];
            int l = new_loc[3];
            //std::cout << "rank: " << rank <<  " ghost_site_send()\n";
            //std::cout << "rank: " << rank <<  " i: " << i <<  " j: " << j <<  " k: " << k <<  " l: " << l << "\n";
            //std::cout << "rank: " << rank << " new_loc\n";
            print_1Dvector(new_loc);

            int tag;
            if (parallel_transfer) tag = 2;
            else tag = 3;
            //std::cout << "rank: " << rank << " parallel_transfer: " << parallel_transfer << "\n";
            //std::cout << "rank: " << rank << " tag: " << tag << "\n";
            
            if (j_old < 1) {/*communicate with proc to -x direction*/
                if (k_old < 1) {
                    //std::cout << "rank: " << rank << " ghost with proc to (-x,-y) direction \n\n";
                    /*
                    std::cout << "rank: " << rank << " comm with proc: " << proc_neighbors(rank,4) << "\n";
                    std::cout << "rank: " << rank << " comm with proc: " << proc_neighbors(rank,5) << "\n";
                    std::cout << "rank: " << rank << " comm with proc: " << proc_neighbors(rank,6) << "\n";
                    */
                    /*communicate with proc to (-x,-y) direction*/
                    if ((proc_neighbors(rank,6) == proc_neighbors(rank,5)) && (proc_neighbors(rank,4) == proc_neighbors(rank,5))) {std::cout << "rank: " << rank << " pass (-x,-y)\n\n";}     
                    else if ((proc_neighbors(rank,6) == rank) && (proc_neighbors(rank,4) != rank)) {
                        //std::cout << "rank: " << rank << " sending to: " << proc_neighbors(rank,4) << "\n";
                        new_proc = proc_neighbors(rank,4);
                        for (int idx=0; idx<new_loc.size(); idx++) {loc_buffer1[idx] = new_loc[idx];}
                        loc_buffer1[4] = rank;
                        loc_buffer1[5] = i_old;
                        loc_buffer1[6] = j_old;
                        loc_buffer1[7] = k_old;
                        loc_buffer1[8] = l_old;
                        MPI_Isend(
                            loc_buffer1.data(),        //Address of the message we are sending.
                            bufferlen,                  //Number of elements handled by that address.
                            MPI_INT,            //MPI_TYPE of the message we are sending.
                            new_proc,           //Rank of receiving process
                            tag,                  //Message Tag
                            MPI_COMM_WORLD,      //MPI Communicator
                            &request1 ); 
                    }
                    else if ((proc_neighbors(rank,4) == rank) && (proc_neighbors(rank,6) != rank)) {
                        //std::cout << "rank: " << rank << " sending to: " << proc_neighbors(rank,6) << "\n";
                        new_proc = proc_neighbors(rank,6);
                        for (int idx=0; idx<new_loc.size(); idx++) {loc_buffer1[idx] = new_loc[idx];}
                        loc_buffer1[4] = rank;
                        loc_buffer1[5] = i_old;
                        loc_buffer1[6] = j_old;
                        loc_buffer1[7] = k_old;
                        loc_buffer1[8] = l_old;
                        MPI_Isend(
                            loc_buffer1.data(),        //Address of the message we are sending.
                            bufferlen,                  //Number of elements handled by that address.
                            MPI_INT,            //MPI_TYPE of the message we are sending.
                            new_proc,           //Rank of receiving process
                            tag,                  //Message Tag
                            MPI_COMM_WORLD,      //MPI Communicator
                            &request1 ); 
                    }
                    else if ((proc_neighbors(rank,6) != rank) && (proc_neighbors(rank,4) != rank)) {
                        //std::cout << "rank: " << rank << " sending to: " << proc_neighbors(rank,6) << "\n";
                        new_proc = proc_neighbors(rank,6);
                        for (int idx=0; idx<new_loc.size(); idx++) {loc_buffer1[idx] = new_loc[idx];}
                        loc_buffer1[4] = rank;
                        loc_buffer1[5] = i_old;
                        loc_buffer1[6] = j_old;
                        loc_buffer1[7] = k_old;
                        loc_buffer1[8] = l_old;
                        MPI_Isend(
                            loc_buffer1.data(),        //Address of the message we are sending.
                            bufferlen,                  //Number of elements handled by that address.
                            MPI_INT,            //MPI_TYPE of the message we are sending.
                            new_proc,           //Rank of receiving process
                            tag,                  //Message Tag
                            MPI_COMM_WORLD,      //MPI Communicator
                            &request1 ); 
                        //std::cout << "rank: " << rank << " sending to: " << proc_neighbors(rank,5) << "\n";
                        new_proc = proc_neighbors(rank,5);
                        for (int idx=0; idx<new_loc.size(); idx++) {loc_buffer2[idx] = new_loc[idx];}
                        loc_buffer2[4] = rank;
                        loc_buffer2[5] = i_old;
                        loc_buffer2[6] = j_old;
                        loc_buffer2[7] = k_old;
                        loc_buffer2[8] = l_old;
                        MPI_Isend(loc_buffer2.data(), bufferlen, MPI_INT, new_proc, tag, MPI_COMM_WORLD, &request2);
                        //std::cout << "rank: " << rank << " sending to: " << proc_neighbors(rank,4) << "\n";
                        new_proc = proc_neighbors(rank,4);
                        for (int idx=0; idx<new_loc.size(); idx++) {loc_buffer3[idx] = new_loc[idx];}
                        loc_buffer3[4] = rank;
                        loc_buffer3[5] = i_old;
                        loc_buffer3[6] = j_old;
                        loc_buffer3[7] = k_old;
                        loc_buffer3[8] = l_old;
                        MPI_Isend(loc_buffer3.data(), bufferlen, MPI_INT, new_proc, tag, MPI_COMM_WORLD, &request3);
                    }
                }
                else {
                    /*communicate with proc to (-x) direction*/
                    //std::cout << "rank: " << rank << " ghost with proc to (-x) direction \n\n";
                    //std::cout << "rank: " << rank << " comm with proc: " << proc_neighbors(rank,4) << "\n";
                    new_proc = proc_neighbors(rank,4); // proc_neighbors indices follow: 0:+x, 1:(+x,+y), 2:+y, 3:(-x,+y), 4:-x, 5:(-x,-y), 6:-y, 7:(+x,-y)

                    if (new_proc == rank) {/*std::cout << "rank: " << rank << " pass (-x)\n\n";*/}
                    else {
                        //std::cout << "rank: " << rank << " sending to: " << proc_neighbors(rank,4) << "\n";
                        for (int idx=0; idx<new_loc.size(); idx++) {loc_buffer1[idx] = new_loc[idx];}
                        loc_buffer1[4] = rank;
                        loc_buffer1[5] = i_old;
                        loc_buffer1[6] = j_old;
                        loc_buffer1[7] = k_old;
                        loc_buffer1[8] = l_old;
                        MPI_Isend(loc_buffer1.data(), bufferlen, MPI_INT, new_proc, tag, MPI_COMM_WORLD, &request1);
                    }
                }
            }
            else if (j_old > (sublattice_dim[0] - 2)) {/*communicate with proc to +x direction*/
                if (k_old > (sublattice_dim[1] - 2)) {
                    
                    //std::cout << "rank: " << rank << " ghost with proc to (+x,+y) direction \n\n";
                    /*
                    std::cout << "rank: " << rank << " comm with proc: " << proc_neighbors(rank,0) << "\n";
                    std::cout << "rank: " << rank << " comm with proc: " << proc_neighbors(rank,1) << "\n";
                    std::cout << "rank: " << rank << " comm with proc: " << proc_neighbors(rank,2) << "\n";
                    */
                    /*communicate with proc to (+x,+y) direction*/
                    if ((proc_neighbors(rank,2) == proc_neighbors(rank,1)) && (proc_neighbors(rank,2) == proc_neighbors(rank,0))) {std::cout << "rank: " << rank << " pass (+x,+y)\n\n";}
                    else if ((proc_neighbors(rank,0) == rank) && (proc_neighbors(rank,2) != rank)) {
                        //std::cout << "rank: " << rank << " sending to: " << proc_neighbors(rank,2) << "\n";
                        new_proc = proc_neighbors(rank,2);
                        for (int idx=0; idx<new_loc.size(); idx++) {loc_buffer1[idx] = new_loc[idx];}
                        loc_buffer1[4] = rank;
                        loc_buffer1[5] = i_old;
                        loc_buffer1[6] = j_old;
                        loc_buffer1[7] = k_old;
                        loc_buffer1[8] = l_old;
                        MPI_Isend(
                            loc_buffer1.data(),        //Address of the message we are sending.
                            bufferlen,                  //Number of elements handled by that address.
                            MPI_INT,            //MPI_TYPE of the message we are sending.
                            new_proc,           //Rank of receiving process
                            tag,                  //Message Tag
                            MPI_COMM_WORLD,      //MPI Communicator
                            &request1 ); 
                    }
                    else if ((proc_neighbors(rank,2) == rank) && (proc_neighbors(rank,0) != rank)) {
                        //std::cout << "rank: " << rank << " sending to: " << proc_neighbors(rank,0) << "\n";
                        new_proc = proc_neighbors(rank,0);
                        for (int idx=0; idx<new_loc.size(); idx++) {loc_buffer1[idx] = new_loc[idx];}
                        loc_buffer1[4] = rank;
                        loc_buffer1[5] = i_old;
                        loc_buffer1[6] = j_old;
                        loc_buffer1[7] = k_old;
                        loc_buffer1[8] = l_old;
                        MPI_Isend(
                            loc_buffer1.data(),        //Address of the message we are sending.
                            bufferlen,                  //Number of elements handled by that address.
                            MPI_INT,            //MPI_TYPE of the message we are sending.
                            new_proc,           //Rank of receiving process
                            tag,                  //Message Tag
                            MPI_COMM_WORLD,      //MPI Communicator
                            &request1 ); 
                    }
                    else if ((proc_neighbors(rank,2) != rank) &&  (proc_neighbors(rank,0) != rank)) {
                        //std::cout << "rank: " << rank << " sending to: " << proc_neighbors(rank,2) << "\n";
                        new_proc = proc_neighbors(rank,2);
                        for (int idx=0; idx<new_loc.size(); idx++) {loc_buffer1[idx] = new_loc[idx];}
                        loc_buffer1[4] = rank;
                        loc_buffer1[5] = i_old;
                        loc_buffer1[6] = j_old;
                        loc_buffer1[7] = k_old;
                        loc_buffer1[8] = l_old;
                        MPI_Isend(
                            loc_buffer1.data(),        //Address of the message we are sending.
                            bufferlen,                  //Number of elements handled by that address.
                            MPI_INT,            //MPI_TYPE of the message we are sending.
                            new_proc,           //Rank of receiving process
                            tag,                  //Message Tag
                            MPI_COMM_WORLD,      //MPI Communicator
                            &request1 ); 
                        new_proc = proc_neighbors(rank,1);
                        //std::cout << "rank: " << rank << " sending to: " << proc_neighbors(rank,1) << "\n";
                        for (int idx=0; idx< new_loc.size(); idx++) {loc_buffer2[idx] = new_loc[idx];}
                        loc_buffer2[4] = rank;
                        loc_buffer2[5] = i_old;
                        loc_buffer2[6] = j_old;
                        loc_buffer2[7] = k_old;
                        loc_buffer2[8] = l_old;
                        MPI_Isend(loc_buffer2.data(), bufferlen, MPI_INT, new_proc, tag, MPI_COMM_WORLD, &request2);
                        new_proc = proc_neighbors(rank,0);
                        //std::cout << "rank: " << rank << " sending to: " << proc_neighbors(rank,0) << "\n";
                        for (int idx=0; idx< new_loc.size(); idx++) {loc_buffer3[idx] = new_loc[idx];}
                        loc_buffer3[4] = rank;
                        loc_buffer3[5] = i_old;
                        loc_buffer3[6] = j_old;
                        loc_buffer3[7] = k_old;
                        loc_buffer3[8] = l_old;
                        MPI_Isend(loc_buffer3.data(), bufferlen, MPI_INT, new_proc, tag, MPI_COMM_WORLD, &request3);
                    }
                }
                else {
                    //std::cout << "rank: " << rank << " ghost with proc to (+x) direction \n\n";
                    //std::cout << "rank: " << rank << " comm with proc: " << proc_neighbors(rank,0) << "\n";
                    /*communicate with proc to (+x) direction*/
                    
                    new_proc = proc_neighbors(rank,0); // proc_neighbors indices follow: 0:+x, 1:(+x,+y), 2:+y, 3:(-x,+y), 4:-x, 5:(-x,-y), 6:-y, 7:(+x,-y) 
                    //std::cout << "rank: " << rank << " sending to: " << proc_neighbors(rank,0) << "\n";
                    if (new_proc == rank) {/*std::cout << "rank: " << rank << " pass (+x)\n\n";*/}
                    else {
                        for (int idx=0; idx<new_loc.size(); idx++) {loc_buffer1[idx] = new_loc[idx];}
                        loc_buffer1[4] = rank;
                        loc_buffer1[5] = i_old;
                        loc_buffer1[6] = j_old;
                        loc_buffer1[7] = k_old;
                        loc_buffer1[8] = l_old;
                        MPI_Isend(loc_buffer1.data(), bufferlen, MPI_INT, new_proc, tag, MPI_COMM_WORLD, &request1);
                    }
                }
            }
            else if ( (k_old < 1) ) {/*communicate with proc to -y direction*/ 
                //std::cout << "rank: " << rank << " ghost with proc to (-y) direction \n\n";
                //std::cout << "rank: " << rank << " comm with proc: " << proc_neighbors(rank,6) << "\n";
                new_proc = proc_neighbors(rank,6); // proc_neighbors indices follow: 0:+x, 1:(+x,+y), 2:+y, 3:(-x,+y), 4:-x, 5:(-x,-y), 6:-y, 7:(+x,-y) 

                if (new_proc == rank) {/*std::cout << "rank: " << rank << " pass (-y)\n\n";*/}
                else {
                    //std::cout << "rank: " << rank << " sending to: " << proc_neighbors(rank,6) << "\n";
                    for (int idx=0; idx<new_loc.size(); idx++) {loc_buffer1[idx] = new_loc[idx];}
                    loc_buffer1[4] = rank;
                    loc_buffer1[5] = i_old;
                    loc_buffer1[6] = j_old;
                    loc_buffer1[7] = k_old;
                    loc_buffer1[8] = l_old; 
                    MPI_Isend(loc_buffer1.data(), bufferlen, MPI_INT, new_proc, tag, MPI_COMM_WORLD, &request1);
                }
            }
            else if ( (k_old > (sublattice_dim[1] - 2)) ) {/*communicate with proc to +y direction*/
                //std::cout << "rank: " << rank << " ghost with proc to (+y) direction \n\n";
                //std::cout << "rank: " << rank << " comm with proc: " << proc_neighbors(rank,2) << "\n";
                new_proc = proc_neighbors(rank,2); // proc_neighbors indices follow: 0:+x, 1:(+x,+y), 2:+y, 3:(-x,+y), 4:-x, 5:(-x,-y), 6:-y, 7:(+x,-y) 

                if (new_proc == rank) {/*std::cout << "rank: " << rank << " pass (+y)\n\n";*/}
                else {
                    //std::cout << "rank: " << rank << " sending to: " << proc_neighbors(rank,2) << "\n";
                    for (int idx=0; idx<new_loc.size(); idx++) {loc_buffer1[idx] = new_loc[idx];}
                    loc_buffer1[4] = rank;
                    loc_buffer1[5] = i_old;
                    loc_buffer1[6] = j_old;
                    loc_buffer1[7] = k_old;
                    loc_buffer1[8] = l_old;
                    MPI_Isend(loc_buffer1.data(), bufferlen, MPI_INT, new_proc, tag, MPI_COMM_WORLD, &request1);
                }
            }
        }


        /**
        * @brief Handles the reception of ghost site information from neighboring processes.
        *
        * This function updates the neighbor process lists based on the new location of a ghost site 
        * received in the `new_loc_buffer`. It accounts for both parallel and non-parallel transfers 
        * and modifies the appropriate neighbors according to the ghost site's old and new positions.
        *
        * @param new_loc_buffer A vector containing the new location and other related data of the ghost site.
        *                       Expected format:
        *                       - new_loc_buffer[0]: New lattice coordinate (vertex/bc site)
        *                       - new_loc_buffer[1]: New x-coordinate
        *                       - new_loc_buffer[2]: New y-coordinate
        *                       - new_loc_buffer[3]: New z-coordinate
        *                       - new_loc_buffer[4]: Old processor ID
        *                       - new_loc_buffer[5]: Old lattice coordinate (vertex/bc site)
        *                       - new_loc_buffer[6]: Old x-coordinate
        *                       - new_loc_buffer[7]: Old y-coordinate
        *                       - new_loc_buffer[8]: Old z-coordinate
        * @param parallel_transfer A boolean indicating whether the transfer is parallel.
        */
        void ghost_site_recieve(const std::vector<int>& new_loc_buffer, bool parallel_transfer) {
            //std::cout << "rank: " << rank << " ghost_site_recieve()\n\n";
            //
            assert( (((new_loc_buffer)[5] == 0) || ((new_loc_buffer)[5] == 1)) && (((new_loc_buffer)[0] == 0) || ((new_loc_buffer)[0] == 1)) );
            assert( (((new_loc_buffer)[6] >= 0) && ((new_loc_buffer)[6] < sublattice_dim[0])) && (((new_loc_buffer)[1] >= -1) && ((new_loc_buffer)[1] < (sublattice_dim[0] + 1))) );
            assert( (((new_loc_buffer)[7] >= 0) && ((new_loc_buffer)[7] < sublattice_dim[1])) && (((new_loc_buffer)[2] >= -1) && ((new_loc_buffer)[2] < (sublattice_dim[1] + 1))) );
            assert( (((new_loc_buffer)[8] >= 0) && ((new_loc_buffer)[8] < sublattice_dim[2])) && (((new_loc_buffer)[3] >= -1) && ((new_loc_buffer)[3] < (sublattice_dim[2] + 1))) );
            //
            int old_proc = (new_loc_buffer)[4];
            int i_old = (new_loc_buffer)[5];
            int j_old = (new_loc_buffer)[6];
            int k_old = (new_loc_buffer)[7];
            int l_old = (new_loc_buffer)[8];

            /*
            make cases for which proc neighbors to update depending on the location of the
            old site
            */
            int i = (new_loc_buffer)[0];
            int j = (((new_loc_buffer)[1] % sublattice_dim[0] + sublattice_dim[0]) % sublattice_dim[0]);
            int k = (((new_loc_buffer)[2] % sublattice_dim[1] + sublattice_dim[1]) % sublattice_dim[1]);
            int l = (((new_loc_buffer)[3] % sublattice_dim[2] + sublattice_dim[2]) % sublattice_dim[2]);
	        /*
            std::cout << "rank: " << rank <<  " ghost_site_recieve()\n";
            std::cout << "rank: " << rank << " old_proc " << old_proc << " \n";
            std::cout << "rank: " << rank << " new_loc [ " << i << " " << j << " " << k << " " << l << " ] \n";
            std::cout << "rank: " << rank << " old_loc [ " << i_old << " " << j_old << " " << k_old << " " << l_old << " ] \n\n";
            */
            size_t ipos = i_old; size_t ineg = !i_old;

            std::vector<size_t> x_dims = proc_pos_x_neighbors.size_vec;
            std::vector<size_t> y_dims = proc_pos_y_neighbors.size_vec;
            // ISSUE WITH OLD_PROC VARIABLE //

            if (!parallel_transfer) {
                /* check for parallel_transfer and if it's in the same processor on each if statement */
                //std::cout << "no parallel transfer\n";
                // issue with offset between corner site updates and edge updates (e.g. (-x,-y) vs -x)
                if ((j == (sublattice_dim[0]-1)) && (k == (sublattice_dim[1]-1))) {
                    //std::cout << "receive new ghost with proc with (-x,-y) \n";
                    if (proc_neighbors(rank,2) == rank) {
                        proc_neg_x_neighbors(0, (size_t)i, (size_t)0, (size_t)l) = 1;
                        proc_neg_y_neighbors(0, (size_t)i, (size_t)0, (size_t)l) = 1;  
                    }    
                    proc_neg_x_neighbors(0, (size_t)i, (size_t)(x_dims[2]-1), (size_t)(l)) = 1;
                    proc_pos_y_neighbors(0, (size_t)i, (size_t)(y_dims[2]-1), (size_t)(l)) = 1;
                    /*
                    std::cout << "rank: "  << rank << " proc_neg_x_neighbors new loc \n";
                    std::cout << "rank: "  << rank << " [ 0 " << i << " " << "0 " << l << " ]" << "\n";  
                    std::cout << "rank: "  << rank << " proc_neg_y_neighbors new loc \n";
                    std::cout << "rank: "  << rank << " [ 0 " << i << " " << "0 " << l << " ]" << "\n"; 
                    */  
                }
                if ((j == 0) && (k == 0)) {
                    //std::cout << "receive new ghost with proc with (+x,+y) \n";
                    if (proc_neighbors(rank,6) == rank) {
                        proc_pos_x_neighbors(0, (size_t)i, (size_t)(x_dims[2]-1), (size_t)(l)) = 1;
                        proc_pos_y_neighbors(0, (size_t)i, (size_t)(y_dims[2]-1), (size_t)(l)) = 1;
                    } 
                    proc_pos_x_neighbors(0, (size_t)i, (size_t)0, (size_t)l) = 1;
                    proc_neg_y_neighbors(0, (size_t)i, (size_t)0, (size_t)l) = 1;   
                    /*
                    std::cout << "rank: "  << rank << " proc_pos_x_neighbors new loc \n";
                    std::cout << "rank: "  << rank << " [ 0 " << i << " " << (x_dims[2]-1) << " " << l << " ]" << "\n";  
                    std::cout << "rank: "  << rank << " proc_pos_y_neighbors new loc \n";
                    std::cout << "rank: "  << rank << " [ 0 " << i << " " << (y_dims[2]-1) << " " << l << " ]" << "\n";  
                    */
                }
                else {
                    if ((k >= 0) && (k <= (sublattice_dim[1]-1)) && (j < 1)) {
                        //std::cout << "receive new ghost with proc with (+x) \n";
                        proc_pos_x_neighbors(0, (size_t)i, (size_t)(k), (size_t)(l)) = 1;  
                        //std::cout << "rank: "  << rank << " proc_pos_x_neighbors new loc \n";
                        //std::cout << "rank: "  << rank << " [ 0 " << i << " " << k << " " << l << " ]" << "\n";   
                    }                
                    if ((j >= 0) && (j <= (sublattice_dim[0]-1)) && (k < 1)) {
                        //std::cout << "receive new ghost with proc with (+y) \n";
                        proc_pos_y_neighbors(0, (size_t)i, (size_t)(j), (size_t)(l)) = 1; 
                         
                        //std::cout << "rank: "  << rank << " proc_pos_y_neighbors new loc \n";
                        //std::cout << "rank: "  << rank << " [ 0 " << i << " " << j << " " << l << " ]" << "\n";  
                        
                    }
                    if ((k >= 0) && (k <= (sublattice_dim[1]-1)) && (j > (sublattice_dim[0]-2))) {
                        //std::cout << "receive new ghost with proc with (-x) \n";
                        proc_neg_x_neighbors(0, (size_t)i, (size_t)(k+1), (size_t)(l)) = 1; 
                         
                        //std::cout << "rank: "  << rank << " proc_neg_x_neighbors new loc \n";
                        //std::cout << "rank: "  << rank << " [ 0 " << i << " " << (k+1) << " " << l << " ]" << "\n";  
                         
                    }
                    if ((j >= 0) && (j <= (sublattice_dim[0]-1)) && (k > (sublattice_dim[1]-2))) {
                        //std::cout << "receive new ghost with proc with (-y) \n";
                        proc_neg_y_neighbors(0, (size_t)i, (size_t)(j+1), (size_t)(l)) = 1; 
                         
                        //std::cout << "rank: "  << rank << " proc_neg_y_neighbors new loc \n";
                        //std::cout << "rank: "  << rank << " [ 0 " << i << " " << (j+1) << " " << l << " ]" << "\n";  
                        
                    }
                }
                
            }
            else {
                /* check for parallel_transfer and if it's in the same processor on each if statement */
                //std::cout << "yes parallel transfer\n";
                if ((j == (sublattice_dim[0]-1)) && (k == (sublattice_dim[1]-1))) {
                    if (rank == proc_neighbors(rank,5)) {
                        //std::cout << "receive (with parallel) new ghost with proc with (-x,-y) \n";
                        if (proc_neighbors(rank,2) == rank) {
                            proc_neg_x_neighbors(0, (size_t)i, (size_t)0, (size_t)l) = 1;
                            proc_neg_y_neighbors(0, (size_t)i, (size_t)0, (size_t)l) = 1;  
                        }    
                        proc_neg_x_neighbors(0, (size_t)i, (size_t)(x_dims[2]-1), (size_t)(l)) = 1;
                        proc_pos_y_neighbors(0, (size_t)i, (size_t)(y_dims[2]-1), (size_t)(l)) = 1;
                        /*
                        std::cout << "rank: "  << rank << " proc_neg_x_neighbors new loc \n";
                        std::cout << "rank: "  << rank << " [ 0 " << i << " " << "0 " << l << " ]" << "\n";  
                        std::cout << "rank: "  << rank << " proc_neg_y_neighbors new loc \n";
                        std::cout << "rank: "  << rank << " [ 0 " << i << " " << "0 " << l << " ]" << "\n";  
                        */
                        } 
                }
                if ((j == 0) && (k == 0)) {
                    if (rank == proc_neighbors(rank,1)) {
                        //std::cout << "receive (with parallel) new ghost with proc with (+x,+y) \n";
                        if (proc_neighbors(rank,6) == rank) {
                            proc_pos_x_neighbors(0, (size_t)i, (size_t)(x_dims[2]-1), (size_t)(l)) = 1;
                            proc_pos_y_neighbors(0, (size_t)i, (size_t)(y_dims[2]-1), (size_t)(l)) = 1;
                        } 
                        proc_pos_x_neighbors(0, (size_t)i, (size_t)0, (size_t)l) = 1;
                        proc_neg_y_neighbors(0, (size_t)i, (size_t)0, (size_t)l) = 1;   
 
                        /*
                        std::cout << "rank: "  << rank << " proc_pos_x_neighbors new loc \n";
                        std::cout << "rank: "  << rank << " [ 0 " << i << " " << (x_dims[2]-1) << " " << l << " ]" << "\n";  
                        std::cout << "rank: "  << rank << " proc_pos_y_neighbors new loc \n";
                        std::cout << "rank: "  << rank << " [ 0 " << i << " " << (y_dims[2]-1) << " " << l << " ]" << "\n"; 
                        */
 
                    }
                }
                if ((k >= 0) && (k <= (sublattice_dim[1]-1)) && (j == (sublattice_dim[0]-1))) {
                    if (rank == proc_neighbors(rank,0)) { 
                        //std::cout << "receive (with parallel) new ghost with proc with (-x) \n";
                        proc_neg_x_neighbors(0, (size_t)i, (size_t)(k+1), (size_t)(l)) = 1; 
                        
                        //std::cout << "rank: "  << rank << " proc_neg_x_neighbors new loc \n";
                        //std::cout << "rank: "  << rank << " [ 0 " << i << " " << (k+1) << " " << l << " ]" << "\n"; 
                         
                    }
                }
                if ((j >= 0) && (j <= (sublattice_dim[0]-1)) && (k == (sublattice_dim[1]-1))) {
                    if (rank == proc_neighbors(rank,2)) { 
                        //std::cout << "receive (with parallel) new ghost with proc with (-y) \n";
                        proc_neg_y_neighbors(0, (size_t)i, (size_t)(j+1), (size_t)(l)) = 1; 
                         
                        //std::cout << "rank: "  << rank << " proc_neg_y_neighbors new loc \n";
                        //std::cout << "rank: "  << rank << " [ 0 " << i << " " << (j+1) << " " << l << " ]" << "\n";
                         
                    }
                }
                if ((k >= 0) && (k <= (sublattice_dim[1]-1)) && (j == 0)) {
                    if (rank == proc_neighbors(rank,4)) { 
                        //std::cout << "receive (with parallel) new ghost with proc with (+x) \n";
                        proc_pos_x_neighbors(0, (size_t)i, (size_t)(k), (size_t)(l)) = 1; 
                         
                        //std::cout << "rank: "  << rank << " proc_pos_x_neighbors new loc \n";
                        //std::cout << "rank: "  << rank << " [ 0 " << i << " " << k << " " << l << " ]" << "\n"; 
                         
                    }
                }
                if ((j >= 0) && (j <= (sublattice_dim[0]-1)) && (k == 0)) {
                    if (rank == proc_neighbors(rank,6)) { 
                        //std::cout << "receive (with parallel) new ghost with proc with (+y) \n";
                        proc_pos_y_neighbors(0, (size_t)i, (size_t)(j), (size_t)(l)) = 1; 
                         
                        //std::cout << "rank: "  << rank << " proc_pos_y_neighbors new loc \n";
                        //std::cout << "rank: "  << rank << " [ 0 " << i << " " << j << " " << l << " ]" << "\n";
                         
                        }
                }
            }
            // issues with removing corner vs edge sites because more than one ghost corner site maps to same lattice site
            /* removing site */
            if ((j_old == (sublattice_dim[0]-1)) && (k_old == (sublattice_dim[1]-1))) {
                //std::cout << "receive old ghost with proc with (-x,-y) \n";
                proc_neg_x_neighbors(0, (size_t)i_old, (size_t)0, (size_t)l_old) = 0;
                proc_neg_y_neighbors(0, (size_t)i_old, (size_t)0, (size_t)l_old) = 0; 
                /*
                std::cout << "rank: "  << rank << " proc_neg_x_neighbors new loc \n";
                std::cout << "rank: "  << rank << " [ 0 " << i_old << " " << "0 " << l_old << " ]" << "\n";  
                std::cout << "rank: "  << rank << " proc_neg_y_neighbors new loc \n";
                std::cout << "rank: "  << rank << " [ 0 " << i_old << " " << "0 " << l_old << " ]" << "\n";  
                */
            }
            if ((j_old == 0) && (k_old == 0)) {
                //std::cout << "receive old ghost with proc with (+x,+y) \n";
                proc_pos_x_neighbors(0, (size_t)i_old, (size_t)(x_dims[2]-1), (size_t)(l_old)) = 0;
                proc_pos_y_neighbors(0, (size_t)i_old, (size_t)(y_dims[2]-1), (size_t)(l_old)) = 0; 
                /*
                std::cout << "rank: "  << rank << " proc_pos_x_neighbors new loc \n";
                std::cout << "rank: "  << rank << " [ 0 " << i_old << " " << (x_dims[2]-1) << " " << l_old << " ]" << "\n";  
                std::cout << "rank: "  << rank << " proc_pos_y_neighbors new loc \n";
                std::cout << "rank: "  << rank << " [ 0 " << i_old << " " << (y_dims[2]-1) << " " << l_old << " ]" << "\n";  
                */
            }
            if ((k_old >= 0) && (k_old <= (sublattice_dim[1] - 1))) {
                if (j_old == 0 ) { 
                    //std::cout << "receive old ghost with proc with (+x) \n";
                    proc_pos_x_neighbors(0, (size_t)i_old, (size_t)(k_old), (size_t)(l_old)) = 0;  
                    //std::cout << "rank: "  << rank << " proc_pos_x_neighbors new loc \n";
                    //std::cout << "rank: "  << rank << " [ 0 " << i_old << " " << (k_old) << " " << l_old << " ]" << "\n";  
                    }
                if (j_old == (sublattice_dim[0] - 1)) {
                    //std::cout << "receive old ghost with proc with (-x) \n"; 
                    proc_neg_x_neighbors(0, (size_t)i_old, (size_t)(k_old+1), (size_t)(l_old)) = 0;  
                    //std::cout << "rank: "  << rank << " proc_neg_x_neighbors new loc \n";
                    //std::cout << "rank: "  << rank << " [ 0 " << i_old << " " << (k_old+1) << " " << l_old << " ]" << "\n";  
                    }
            }
            if ((j_old >= 0) && (j_old <= (sublattice_dim[0] - 1))) {
                if (k_old == 0) { 
                    //std::cout << "receive old ghost with proc with (+y) \n";
                    proc_pos_y_neighbors(0, (size_t)i_old, (size_t)(j_old), (size_t)(l_old)) = 0; 
                    //std::cout << "rank: "  << rank << " proc_pos_y_neighbors new loc \n";
                    //std::cout << "rank: "  << rank << " [ 0 " << i_old << " " << (j_old) << " " << l_old << " ]" << "\n";   
                }
                if (k_old == (sublattice_dim[1] - 1)) { 
                    //std::cout << "receive old ghost with proc with (-y) \n";
                    proc_neg_y_neighbors(0, (size_t)i_old, (size_t)(j_old+1), (size_t)(l_old)) = 0; 
                    //std::cout << "rank: "  << rank << " proc_neg_y_neighbors new loc \n";
                    //std::cout << "rank: "  << rank << " [ 0 " << i_old << " " << (j_old+1) << " " << l_old << " ]" << "\n";  
                }
            }
            
            
           // std::cout << "done with ghost site receive\n";
            
        }


        /**
        * @brief Update ghost sites in the local domain of the processor.
        * 
        * This function checks the current location of ghost sites and updates 
        * their references based on the new location provided. It distinguishes 
        * between parallel and non-parallel transfer scenarios, adjusting 
        * neighbor references accordingly. The function uses the global rank 
        * of the processor to determine which ghost sites need to be updated 
        * and applies the changes in their respective positions based on the 
        * defined neighbor relationships.
        * 
        * @param i_old Old lattice coordinate (vertex/bc site) of site in ghost zone.
        * @param j_old Old x-coordinate of the ghost site.
        * @param k_old Old y-coordinate of the ghost site.
        * @param l_old Old z-coordinate of the ghost site.
        * @param new_loc A vector containing the new location coordinates: 
        *                {lattice type, x, y, z}.
        * @param parallel_transfer Boolean flag indicating whether the update 
        *                         occurs in parallel.
        */
        void ghost_site_self_reference(int i_old,int j_old,int k_old,int l_old, const std::vector<int>& new_loc, bool parallel_transfer) {
            //std::cout << "rank: " << rank << " enter ghost_site_self_reference\n";
            
            if (!(parallel_transfer)) {
                if ((new_loc[1] < 1) && (rank == proc_neighbors(rank,0))) {
                    proc_pos_x_neighbors(0, new_loc[0], (size_t)new_loc[2], (size_t)new_loc[3]) = 1;    
                    /*
                    std::cout << "rank: "  << rank << " proc_pos_x_neighbors new loc \n";
                    std::cout << "rank: "  << rank << " [ 0 1 " << (new_loc[2]) << " " << (new_loc[3]) << " ]" << "\n";  
                    */
                }
                if ((new_loc[1] > sublattice_dim[0] - 2) && (rank == proc_neighbors(rank,4))) {
                    proc_neg_x_neighbors(0, new_loc[0], (size_t)(new_loc[2]+1), (size_t)new_loc[3]) = 1;
                    /*
                    std::cout << "rank: "  << rank << " proc_neg_x_neighbors new loc \n"; 
                    std::cout << "rank: "  << rank << " [ 0 0 " << (new_loc[2]+1) << " " << (new_loc[3]) << " ]" << "\n";  
                    */
                }
                if ((new_loc[2] < 1) && (rank == proc_neighbors(rank,2))) {
                    proc_pos_y_neighbors(0, new_loc[0], (size_t)new_loc[1], (size_t)new_loc[3]) = 1;
                    /*
                    std::cout << "rank: "  << rank << " proc_pos_y_neighbors new loc \n"; 
                    std::cout << "rank: "  << rank << " [ 0 1 " << (new_loc[1]) << " " << (new_loc[3]) << " ]" << "\n"; 
                    */
                }
                if ((new_loc[2] > (sublattice_dim[1] - 2)) && (rank == proc_neighbors(rank,6))) {
                    proc_neg_y_neighbors(0, new_loc[0], (size_t)(new_loc[1]+1), (size_t)new_loc[3]) = 1;
                    /*
                    std::cout << "rank: "  << rank << " proc_neg_y_neighbors new loc \n"; 
                    std::cout << "rank: "  << rank << " [ 0 0 " << (new_loc[1]+1) << " " << (new_loc[3]) << " ]" << "\n"; 
                    */
                }
            }
            else {
                if ((new_loc[1] < 1) && (rank != proc_neighbors(rank,0))) {
                    proc_pos_x_neighbors(0, new_loc[0], (size_t)new_loc[2], (size_t)new_loc[3]) = 1; 
                    
                    //std::cout << "rank: "  << rank << " proc_pos_x_neighbors new loc \n";
                    //std::cout << "rank: "  << rank << " [ 0 " << (new_loc[0]) << " " << (new_loc[2]) << " " << (new_loc[3]) << " ]" << "\n";  
                    
                    if ((new_loc[2] == 0) && (rank == proc_neighbors(rank,2))) {
                        proc_pos_x_neighbors(0, (size_t)new_loc[0], (size_t)sublattice_dim[1], (size_t)new_loc[3]) = 1; 
                    }
                }
                if ((new_loc[1] > sublattice_dim[0] - 2) && (rank != proc_neighbors(rank,4))) {
                    proc_neg_x_neighbors(0, new_loc[0], (size_t)(new_loc[2]+1), (size_t)new_loc[3]) = 1;
                    //std::cout << "rank: "  << rank << " proc_neg_x_neighbors new loc \n"; 
                    //std::cout << "rank: "  << rank << " [ 0 " << (new_loc[0]) << " " << (new_loc[2]+1) << " " << (new_loc[3]) << " ]" << "\n";  
                    
                    if ((new_loc[2] == (size_t)(sublattice_dim[1] - 1)) && (rank == proc_neighbors(rank,6))) {
                        proc_neg_x_neighbors(0,  (size_t)new_loc[0], 0, (size_t)new_loc[3]) = 1; 
                    }
                }
                if ((new_loc[2] < 1) && (rank != proc_neighbors(rank,2))) {
                    proc_pos_y_neighbors(0, new_loc[0], (size_t)new_loc[1], (size_t)new_loc[3]) = 1; 
                    //std::cout << "rank: "  << rank << " proc_pos_y_neighbors new loc \n"; 
                    //std::cout << "rank: "  << rank << " [ 0 " << (new_loc[0]) << " " << (new_loc[1]) << " " << (new_loc[3]) << " ]" << "\n";  

                    if ((new_loc[1] == 0) && (rank == proc_neighbors(rank,0))) {
                        proc_pos_y_neighbors(0,  (size_t)new_loc[0],  (size_t)sublattice_dim[0], (size_t)new_loc[3]) = 1; 
                    }
                }
                if ((new_loc[2] > (sublattice_dim[1] - 2)) && (rank != proc_neighbors(rank,6))) {
                    proc_neg_y_neighbors(0, new_loc[0], (size_t)(new_loc[1]+1), (size_t)new_loc[3]) = 1; 
                    //std::cout << "rank: "  << rank << " proc_neg_y_neighbors new loc \n"; 
                    //std::cout << "rank: "  << rank << " [ 0 " << (new_loc[0]) << " " << (new_loc[1]+1) << " " << (new_loc[3]) << " ]" << "\n";  

                    if ((new_loc[1] == (size_t)(sublattice_dim[0] -1)) && (rank == proc_neighbors(rank,4))) {
                        proc_neg_y_neighbors(0,  (size_t)new_loc[0],  (size_t)0, (size_t)new_loc[3]) = 1; 
                    }
                }
            }

            if ((j_old < 1) && (proc_neighbors(rank,0) == rank)) {
                proc_pos_x_neighbors(0, i_old, (size_t)(k_old), (size_t)(l_old)) = 0; 
                //std::cout << "rank: "  << rank << " proc_pos_x_neighbors old loc \n"; 
                //std::cout << "rank: "  << rank << " [ 0 " << i_old << " " << (k_old) << " " << (l_old) << " ]" << "\n"; 
            }
            if ((j_old > (sublattice_dim[0] - 2)) && (proc_neighbors(rank,4) == rank)) {
                proc_neg_x_neighbors(0, i_old, (size_t)(k_old+1), (size_t)(l_old)) = 0; 
                //std::cout << "rank: "  << rank << " proc_neg_x_neighbors old loc \n"; 
                //std::cout << "rank: "  << rank << " [ 0 " << i_old << " " << (k_old+1) << " " << (l_old) << " ]" << "\n";   
            }
            if ((k_old < 1) && (proc_neighbors(rank,2) == rank)) {
                proc_pos_y_neighbors(0, i_old, (size_t)(j_old), (size_t)(l_old)) = 0;          
                //std::cout << "rank: "  << rank << " proc_pos_y_neighbors old loc \n"; 
                //std::cout << "rank: "  << rank << " [ 0 " << i_old << " " << (j_old) << " " << (l_old) << " ]" << "\n";   
            }
            if ((k_old > (sublattice_dim[1] - 2)) && (proc_neighbors(rank,6) == rank)) {
                proc_neg_y_neighbors(0, i_old, (size_t)(j_old+1), (size_t)(l_old)) = 0; 
                //std::cout << "rank: "  << rank << " proc_neg_y_neighbors old loc \n";
                //std::cout << "rank: "  << rank << " [ 0 " << i_old << " " << (j_old+1) << " " << (l_old) << " ]" << "\n";  
            }

            //std::cout << "rank: " << rank << " exit ghost_site_self_reference\n";
        }


        /**
        * @brief Updates the positions of atoms on the lattice according to the selected move.
        * 
        * This function updates the positions of vacancies and atoms in the lattice based on the 
        * selected move for a given index. It handles both parallel and non-parallel transfer cases 
        * while ensuring that the boundaries of the domain are respected.
        * 
        * @param idx The index of to access information about move being propagated from moves_vacs,
        *               moves_lattice, and other Matrix data structures
        * @param move_ticks The number of move that is being propagated
        * 
        * @note This function checks the boundaries of the domain against the total dimensions. 
        *       If equal and the move is within ghost sites, a modified move is created to find 
        *       a new location in ghost sites.
        */
        void new_update_lattice(int idx, int move_ticks) {

            bool parallel_transfer = false;
            if (moves_lattice(idx,0) == 5) { 
                prev_move_type.push_back(0);
                prev_move_type_ticks.push_back(move_ticks);
                store_move_info(idx, parallel_transfer, -1, 5);
                }
            else {
                std::vector<int> new_loc(4);
                for (int i=0; i< (int)new_loc.size(); i++) { new_loc[i] = moves_coords(idx,i); /* new location of vacancy */ } 
                std::vector<int> old_loc(3); // new location of vacancy
                int vacs_idx = moves_vacs(idx,0); // index of vacancy in master vector
                std::vector<int> shift(3); // shift of move
                int lattice = moves_lattice(idx,0); // lattice type of move
            
                for (int i=0; i<3; i++) {
                    old_loc[i] = (((new_loc[i+1] - moves_shifts(idx,i)) % sublattice_dim[i]) + sublattice_dim[i]) % sublattice_dim[i];
                    shift[i] = moves_shifts(idx,i);
                }

                int i_old;
                int j_old = old_loc[0];
                int k_old = old_loc[1];
                int l_old = old_loc[2];

                if (moves_lattice(idx,0) == 3) { i_old = 1; }
                else if  (moves_lattice(idx,0) == 2) { i_old = 0; }
                else if  (moves_lattice(idx,0) == 1) { i_old = 1; }
                else if  (moves_lattice(idx,0) == 0) { i_old = 0; }
                
                parallel_transfer = parallel_processes_check(i_old,j_old,k_old,l_old,idx,new_loc);
                
               //std::cout << "lattice rank: " << rank << " update_lattice()\n";
               //std::cout << "lattice rank: " << rank << " moves_lattice[idx][0]: " << moves_lattice(idx,0) << "\n";
               //std::cout << "lattice rank: " << rank << " old_loc: " << "[ " << i_old << " " << j_old << " " << k_old << " " << l_old << " ]\n";
               //std::cout << "lattice rank: " << rank << " new_loc: " << "[ " << new_loc[0] << " " <<  new_loc[1] << " " <<  new_loc[2] << " " <<  new_loc[3] << " ]\n";
                
                if (parallel_transfer) {
                    // switching occupancy for old and new site in lattice array // 
                   //std::cout << "lattice rank: " << rank << " move_ticks: " << move_ticks <<  " moves_lattice yes parallel_transfer \n";

                    if (moves_lattice(idx,0) == 3) {
                        bc_sites((size_t)0, (size_t)old_loc[0], (size_t)old_loc[1], (size_t)old_loc[2]) = 1;
                        vacancies((size_t)1, (size_t)old_loc[0], (size_t)old_loc[1], (size_t)old_loc[2]) = 0;
                    }
                    else if (moves_lattice(idx,0) == 2) {
                        vertex_sites((size_t)0, (size_t)old_loc[0], (size_t)old_loc[1], (size_t)old_loc[2]) = 1;
                        vacancies((size_t)0, (size_t)old_loc[0], (size_t)old_loc[1], (size_t)old_loc[2]) = 0;
                    }
                    else if (moves_lattice(idx,0) == 1) {
                        bc_sites((size_t)0, (size_t)old_loc[0], (size_t)old_loc[1], (size_t)old_loc[2]) = 1;
                        vacancies((size_t)1, (size_t)old_loc[0], (size_t)old_loc[1], (size_t)old_loc[2]) = 0;
                    }
                    else if (moves_lattice(idx,0) == 0) {
                        vertex_sites((size_t)0, (size_t)old_loc[0], (size_t)old_loc[1], (size_t)old_loc[2]) = 1;
                        vacancies((size_t)0, (size_t)old_loc[0], (size_t)old_loc[1], (size_t)old_loc[2]) = 0;
                    }
                    
                    prev_move_type.push_back(1);
                    store_move_info(idx, parallel_transfer, vacs_idx, lattice, shift, new_loc, {i_old, j_old, k_old, l_old});

                    vacancies_pos.remove_row(vacs_idx, rank);
                    moves_vacs.remove_row(idx, rank);
                    moves_lattice.remove_row(idx, rank);
                    moves_shifts.remove_row(idx, rank);
                    moves_coords.remove_row(idx, rank);

                    num_of_vacs --;

                }
                else {
                   //std::cout << "moves_lattice no parallel_transfer \n";
                    // adding vacancy corresponding to stripping move 
                    if (moves_lattice(idx,0) == 4) {
                        /*---
                        BC EDGE MOVES
                        ---*/
                        int new_site = bc_sites((size_t)0, (size_t)new_loc[1], (size_t)new_loc[2], (size_t)new_loc[3]);

                        // removing atom from lattice ###
                        bc_sites((size_t)0, (size_t)new_loc[1], (size_t)new_loc[2], (size_t)new_loc[3]) = (new_site ^ 1);

                        // adding vacancy to lattice ###
                        vacancies((size_t)1, (size_t)new_loc[1], (size_t)new_loc[2], (size_t)new_loc[3]) = new_site;
                        
                        // adjusting the size of data structures to account for new vacancy 
                        num_of_vacs ++;
                        vacancies_pos.reshape(num_of_vacs, 4, rank);
                        for (int i=0; i<4; i++) { vacancies_pos((num_of_vacs-1),i) = new_loc[i]; }
                        
                    }

                    // moving vacancy from bc site to bc site 
                    else if (moves_lattice(idx,0) == 3) {
                        /*---
                        BC EDGE MOVES
                        ---*/
                        int new_site = bc_sites((size_t)0, (size_t)new_loc[1], (size_t)new_loc[2], (size_t)new_loc[3]);
                        int old_site = bc_sites((size_t)0, (size_t)old_loc[0], (size_t)old_loc[1], (size_t)old_loc[2]);
                                    
                        // switching occupancy for old and new site in lattice array ###
                        bc_sites((size_t)0, (size_t)new_loc[1], (size_t)new_loc[2], (size_t)new_loc[3]) = old_site;
                        bc_sites((size_t)0, (size_t)old_loc[0], (size_t)old_loc[1], (size_t)old_loc[2]) = new_site;

                        // switching occupancy for old and new site in vacancy and mobileion arrays ###
                        vacancies((size_t)1, (size_t)new_loc[1], (size_t)new_loc[2], (size_t)new_loc[3]) = (old_site ^ 1);
                        vacancies((size_t)1, (size_t)old_loc[0], (size_t)old_loc[1], (size_t)old_loc[2]) = (new_site ^ 1);
                        
                        vacancies_pos(vacs_idx,0) = 1;    

                    }
                                
                    // moving vacancy from vertex site to vertex site 
                    else if (moves_lattice(idx,0) == 2) {
                        /*---
                        VERTEX EDGE MOVES
                        ---*/
                        int new_site = vertex_sites((size_t)0, (size_t)new_loc[1], (size_t)new_loc[2], (size_t)new_loc[3]);
                        int old_site = vertex_sites((size_t)0, (size_t)old_loc[0], (size_t)old_loc[1], (size_t)old_loc[2]);
                        
                        // switching occupancy for old and new site in lattice array ###
                        vertex_sites((size_t)0, (size_t)new_loc[1], (size_t)new_loc[2], (size_t)new_loc[3]) = old_site;
                        vertex_sites((size_t)0, (size_t)old_loc[0], (size_t)old_loc[1], (size_t)old_loc[2]) = new_site;
                        
                        // switching occupancy for old and new site in vacancy and mobileion arrays ###
                        vacancies((size_t)0, (size_t)new_loc[1], (size_t)new_loc[2], (size_t)new_loc[3]) = (old_site ^ 1);
                        vacancies((size_t)0, (size_t)old_loc[0], (size_t)old_loc[1], (size_t)old_loc[2]) = (new_site ^ 1);
                        
                        vacancies_pos(vacs_idx,0) = 0;

                    }
                                
                    // moving vacancy from bc site to vertex site 
                    else if (moves_lattice(idx,0) == 1) {
                        /*---
                        BC MOVES
                        ---*/
                        int new_site = vertex_sites((size_t)0, (size_t)new_loc[1], (size_t)new_loc[2], (size_t)new_loc[3]);
                        int old_site = bc_sites((size_t)0, (size_t)old_loc[0], (size_t)old_loc[1], (size_t)old_loc[2]);

                        // switching occupancy for old and new site in lattice array ###
                        vertex_sites((size_t)0, (size_t)new_loc[1], (size_t)new_loc[2], (size_t)new_loc[3]) = old_site;
                        bc_sites((size_t)0, (size_t)old_loc[0], (size_t)old_loc[1], (size_t)old_loc[2]) = new_site;

                        // switching occupancy for old and new site in vacancy and mobileion arrays ###
                        vacancies((size_t)0, (size_t)new_loc[1], (size_t)new_loc[2], (size_t)new_loc[3]) = (old_site ^ 1);
                        vacancies((size_t)1, (size_t)old_loc[0], (size_t)old_loc[1], (size_t)old_loc[2]) = (new_site ^ 1);

                        vacancies_pos(vacs_idx,0) = 0; 

                    }
                                                
                    // moving vacancy from vertex site to bc site 
                    else if (moves_lattice(idx,0) == 0) {
                        /*---
                        VERTEX MOVES
                        ---*/
                        // switching occupancy for old and new site in vacancy and mobileion arrays ###
                        int old_site = vertex_sites((size_t)0, (size_t)old_loc[0], (size_t)old_loc[1], (size_t)old_loc[2]);
                        int new_site = bc_sites((size_t)0, (size_t)new_loc[1], (size_t)new_loc[2], (size_t)new_loc[3]);

                        // switching occupancy for old and new site in lattice array ###
                        vertex_sites((size_t)0, (size_t)old_loc[0], (size_t)old_loc[1], (size_t)old_loc[2]) = new_site;
                        bc_sites((size_t)0, (size_t)new_loc[1], (size_t)new_loc[2], (size_t)new_loc[3]) = old_site;

                        // switching occupancy for old and new site in vacancy and mobileion arrays ###
                        vacancies((size_t)0, (size_t)old_loc[0], (size_t)old_loc[1], (size_t)old_loc[2]) = (new_site ^ 1);
                        vacancies((size_t)1, (size_t)new_loc[1], (size_t)new_loc[2], (size_t)new_loc[3]) = (old_site ^ 1); 

                        vacancies_pos(vacs_idx,0) = 1;

                    }
                    
                    // updating vector of positions of all vacancies
                    vacancies_pos(vacs_idx, 1) = new_loc[1];
                    vacancies_pos(vacs_idx, 2) = new_loc[2];
                    vacancies_pos(vacs_idx, 3) = new_loc[3]; 

                    prev_move_type.push_back(0);                
                    store_move_info(idx, parallel_transfer, vacs_idx, lattice, shift, new_loc, {i_old, j_old, k_old, l_old});
                } 
            
                prev_move_type_ticks.push_back(move_ticks);
                
                if ( ((new_loc[1] < 1) || (new_loc[1] > (sublattice_dim[0] - 2))) || ((new_loc[2] < 1) || (new_loc[2] > (sublattice_dim[1] - 2))) ||
                ((old_loc[0] < 1) || (old_loc[0] > (sublattice_dim[0] - 2))) || ((old_loc[1] < 1) || (old_loc[1] > (sublattice_dim[1] - 2))) ) 
                { /*checking to see if neighbor ghost sites need to be updated */
                    ghost_site_send(i_old,j_old,k_old,l_old,idx,new_loc,parallel_transfer);
                    ghost_site_self_reference(i_old,j_old,k_old,l_old,new_loc,parallel_transfer);
                }
            }
        }


        /**
        * @brief Calculates the time elapsed for a move.
        * 
        * This function generates a random double between 0 and 1, 
        * then calculates the elapsed time based on the cumulative rate.
        * 
        * @return double The calculated time elapsed for a move.
        */
        double new_random_times() {
            // creating a random double between 0 and 1
            unsigned int random = mt_obj();
            double random_double = ((random / (1.+ UINT32_MAX)) +  (1 / (1.+ UINT32_MAX)));
            
            //calculating the time elapsed
            int last_idx = (int)rate_cumsum.size() - 1;
            double time = ((-1/ rate_cumsum[last_idx]) * log(random_double));
            
            return time;
        }


        /**
        * @brief Communicates the total rate in each processor domain at each timestep.
        * 
        * This function gathers the maximum rate across all processes, 
        * broadcasts it to all processes, and creates a null move corresponding 
        * to the difference between the total rate (Rtot_i) and the maximum rate (Rmax).
        * 
        */
        void communicate_rates() {
            double max_i_rate = 0;
            int max_rate_idx;
            double max_rate = 0;
            std::vector<double> max_rates(num_procs);

            MPI_Request request;
            int end_idx;
            if (rate_cumsum.size() != 0) { max_i_rate = rate_cumsum[rate_cumsum.size()-1]; }
            else { max_i_rate = 0; }

            MPI_Gather(&max_i_rate, 1, MPI_DOUBLE, max_rates.data(), 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
                            
            if (rank == 0) {
                max_rate_idx = find_max_element(&max_rates);
                //std::cout << "rank: " << rank << " max_rate_idx " << max_rate_idx << "\n"; 
                //std::cout << "rank: " << rank << " max_rates \n";
                //print_1Dvector(max_rates);
                max_rate = max_rates[max_rate_idx];
            }
            
            MPI_Bcast(&max_rate, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
            
            if (rate_cumsum.size() != 0 ) {
                //std::cout << "rank: " << rank << " rate_cumsum[(rate_cumsum.size() - 1)] " << rate_cumsum[(rate_cumsum.size() - 1)]  << "\n";
                if (rate_cumsum[(rate_cumsum.size() - 1)] == max_rate) {//std::cout << "rank: " << rank << " has the max rate! \n";
                }
                else {
                    fflush(stdout);
                    end_idx = (int)moves_lattice.rows() + 1;
                    //std::cout << "rank: " << rank << " doesn't have the max rate - end_idx: " << end_idx << "\n";
                    moves_lattice.reshape(end_idx, 1, rank);
                    moves_lattice((end_idx-1),0) = 5;
                    rate_cumsum.push_back(max_rate);
                    moves_shifts.reshape(end_idx, 3, rank);
                    for (int i=0; i<3; i++) moves_shifts((end_idx-1), i) = 0;
                    moves_vacs.reshape(end_idx, 1, rank);
                    moves_vacs((end_idx-1),0) = -1;
                    moves_coords.reshape(end_idx, 4, rank);
                    for (int i=0; i<4; i++) moves_coords((end_idx-1) ,i) = -1;
                }

            }
            else {
                moves_lattice.reshape(1,1, rank);
                moves_lattice(0,0) = 5;
                rate_cumsum.push_back(max_rate);
                moves_shifts.reshape(1,3, rank);
                for (int i=0; i<3; i++) moves_shifts(0, i) = 0;
                moves_vacs.reshape(1,1, rank);
                moves_vacs(0,0) = -1;
                moves_coords.reshape(1,4, rank);
                for (int i=0; i<4; i++) moves_coords(0, i) = -1;
            }
        }


        /**
        * @brief Selects a random move from a vector of moves.
        * 
        * This function selects a random index in the vector of moves, 
        * with selection probability proportional to the rate constant 
        * corresponding to each move. It communicates the current rates 
        * and generates a random number to determine the selected index.
        * 
        * @return int The index of the selected move.
        */
        int get_idx() {
            //
            assert(rate_cumsum.size() == moves_coords.rows());
            //
            // creating a random double between 0 and 1
            communicate_rates();
            unsigned int random = mt_obj();
            double random_double = ((random / (1.+ UINT32_MAX)) +  (1 / (1.+ UINT32_MAX)));
            int min_idx;

            // accessing random element into cumulative sum array, 
            // probability of access proportional to the value at that point
            int last_idx = (int)rate_cumsum.size() - 1;
            if (last_idx == 0) { min_idx = 0; }
            else {
                double rand_pos = rate_cumsum[last_idx] * random_double;
                // CHECK FUNCTIONALITY OF searchsorted_recursive // 
                //min_idx = searchsorted_recursive(&rate_cumsum, rand_pos, 0, last_idx);
                min_idx = idx_to_insert(rate_cumsum, rand_pos);
            }

            return min_idx;
        }


        /**
        * @brief Communicates a boundary conflict to other processes.
        * 
        * This MPI routine informs all other processes to roll back two steps 
        * and restart the simulation run. It uses a non-blocking send to notify 
        * other ranks of the rollback action.
        */
        void comm_boundary_conflict() {
            int tag = 6;
            MPI_Request request;

            for (int new_proc = 0; new_proc < num_procs; new_proc++) { 
                if (new_proc != rank) 
                    MPI_Isend(NULL, 0, MPI_CHAR, new_proc, tag, MPI_COMM_WORLD, &request); 
            }
        }


        /**
        * @brief Truncates the vector of previous parallel moves.
        * 
        * This function removes old parallel moves from par_prev_moves, retaining 
        * only those moves from the previous two iterative steps. It updates 
        * par_prev associated vectors to ensure they reflect the current state of moves.
        * 
        * @param move_ticks The current number move in simulation
        */
        void remove_old_par_moves(int move_ticks) {
            //std::cout << "rank: " << rank << " entering remove_old_par_moves\n";

            std::vector<std::vector<int>>::iterator ptr1 = par_prev_moves.begin();
            std::vector<int>::iterator ptr2 = par_prev_idx.begin();
            std::vector<int> remove_idxs;
            /*
            std::cout << "rank: " << rank << " pre prev_move_type: \n";
            print_1Dvector(prev_move_type);
            std::cout << "rank: " << rank << " pre prev_move_type_ticks: \n";
            print_1Dvector(prev_move_type_ticks);
            */

            for (int i = 0; i < (int)par_move_ticks.size(); i++) {
                if (move_ticks > (par_move_ticks.at(i) + 1)) {
                    par_prev_moves.erase(ptr1);
                    remove_idxs.push_back(i);
                    par_prev_idx.erase(ptr2);
                }
            }
            
            std::vector<int>::iterator ptr3 = par_move_ticks.begin();
            for (int i = 0; i < (int)remove_idxs.size(); i++) {
                par_move_ticks.erase((ptr3 + (remove_idxs.at(i) - i)));
            }
            remove_idxs.clear();

            for (int i = 0; i < (int)prev_move_type_ticks.size(); i++) {
                if (move_ticks > (prev_move_type_ticks.at(i) + 1)) {
                    remove_idxs.push_back(i);
                }
            }

            std::vector<int>::iterator ptr4 = prev_move_type_ticks.begin();
            std::vector<int>::iterator ptr5 = prev_move_type.begin();

            for (int i = 0; i < (int)remove_idxs.size(); i++) {
                prev_move_type_ticks.erase((ptr4 + (remove_idxs.at(i) - i)));
                prev_move_type.erase((ptr5 + (remove_idxs.at(i) - i)));
            }

            //std::cout << "rank: " << rank << " post prev_move_type: \n";
            //print_1Dvector(prev_move_type);
            //std::cout << "rank: " << rank << " post prev_move_type_ticks: \n";
            //print_1Dvector(prev_move_type_ticks);
        }


        /** 
        * @brief Rolls back a move originating in the native process when an error is encountered in the simulation. 
        * 
        * This function handles two scenarios: one involving parallel transfer and the other without. 
        * It updates the relevant prev_moves vectors based on the previous state of the system. 
        * 
        * @param parallel_transfer A boolean indicating whether the move involved parallel transfer.
        */
        void reverse_move(bool parallel_transfer) {
            std::vector<int> new_vac;
            std::vector<int> old_vac;
            std::vector<int> move;
            int lattice;
            int idx;
            int i1; int i2; int i3; int i4;
            int vacs_idx;

            /*
            std::cout << "rank: " << rank << " prev_newlocs: \n";
            print_2Dvector(prev_newlocs);
            std::cout << "rank: " << rank << " prev_oldlocs \n";
            print_2Dvector(prev_oldlocs);
            std::cout << "rank: " << rank << " prev_moves: \n";
            print_2Dvector(prev_moves);
            std::cout << "rank: " << rank << " prev_lattice: \n";
            print_1Dvector(prev_lattice);
            */
                    
            new_vac = prev_newlocs.at((prev_newlocs.size() - 1));
            old_vac = prev_oldlocs.at((prev_oldlocs.size() - 1));
            move = prev_moves.at((prev_moves.size() - 1));
            lattice = prev_lattice.at((prev_lattice.size() - 1));
            vacs_idx  = prev_idxs.at((prev_idxs.size() - 1));

            if (prev_lattice.at((prev_lattice.size()-1)) == 5) {}
            else {
                if (parallel_transfer) {
                    //std::cout << "rank: " << rank << " parallel_transfer reverse \n"; 

                    i1 = old_vac[0];
                    i2 = old_vac[1];
                    i3 = old_vac[2];
                    i4 = old_vac[3];
                    /*
                    std::cout << "rank: " << rank << " prev_site (set 1): " << "[ " << old_vac[0] << " " << old_vac[1] << " " << old_vac[2] << " " << old_vac[3] << " ]\n"; 
                    std::cout << "rank: " << rank << " prev site status: " << vacancies(old_vac[0],old_vac[1],old_vac[2],old_vac[3]) << " \n"; 
                    std::cout << "rank: " << rank << " new_site (set 0): " << "[ " << new_vac[0] << " " << new_vac[1] << " " << new_vac[2] << " " << new_vac[3] << " ]\n";
                    std::cout << "rank: " << rank << " new site status: " << vacancies(new_vac[0],new_vac[1],new_vac[2],new_vac[3]) << " \n"; 
                    */
                    vacancies_pos.add_row(vacs_idx, rank);
                    vacancies_pos(vacs_idx,1) = (size_t)old_vac[1];
                    vacancies_pos(vacs_idx,2) = (size_t)old_vac[2];
                    vacancies_pos(vacs_idx,3) = (size_t)old_vac[3];

                    // moving vacancy from bc site to bc site 
                    if (lattice == 3) {
                        /*---
                        BC EDGE MOVES
                        ---*/
                                    
                        // switching occupancy for old and new site in lattice array ###
                        bc_sites((size_t)0, (size_t)old_vac[1], (size_t)old_vac[2], (size_t)old_vac[3]) = 0;

                        // switching occupancy for old and new site in vacancy and mobileion arrays ###
                        vacancies((size_t)1, (size_t)old_vac[1], (size_t)old_vac[2], (size_t)old_vac[3]) = 1;

                        vacancies_pos(vacs_idx,0) = 1;
                    }
                                    
                    // moving vacancy from vertex site to vertex site 
                    else if (lattice == 2) {
                        /*---
                        VERTEX EDGE MOVES
                        ---*/
                        
                        // switching occupancy for old and new site in lattice array ###
                        vertex_sites((size_t)0, (size_t)old_vac[1], (size_t)old_vac[2], (size_t)old_vac[3]) = 0;
                        
                        // switching occupancy for old and new site in vacancy and mobileion arrays ###
                        vacancies((size_t)0, (size_t)old_vac[1], (size_t)old_vac[2], (size_t)old_vac[3]) = 1;

                        vacancies_pos(vacs_idx,0) = 0;
                    }
                                    
                    // moving vacancy from bc site to vertex site 
                    else if (lattice == 1) {
                        /*---
                        BC MOVES
                        ---*/

                        // switching occupancy for old and new site in lattice array ###
                        bc_sites((size_t)0, (size_t)old_vac[1], (size_t)old_vac[2], (size_t)old_vac[3]) = 0;

                        // switching occupancy for old and new site in vacancy and mobileion arrays ###
                        vacancies((size_t)1, (size_t)old_vac[1], (size_t)old_vac[2], (size_t)old_vac[3]) = 1;

                        vacancies_pos(vacs_idx,0) = 1;
                    }
                                                    
                    // moving vacancy from vertex site to bc site 
                    else if (lattice == 0) {
                        /*---
                        VERTEX MOVES
                        ---*/

                        // switching occupancy for old and new site in lattice array ###
                        vertex_sites((size_t)0, (size_t)old_vac[1], (size_t)old_vac[2], (size_t)old_vac[3]) = 0;

                        // switching occupancy for old and new site in vacancy and mobileion arrays ###
                        vacancies((size_t)0, (size_t)old_vac[1], (size_t)old_vac[2], (size_t)old_vac[3]) = 1;

                        vacancies_pos(vacs_idx,0) = 0;
                    }
                    //std::cout << "rank: " << rank << " vacancies_pos(vacs_idx) [ " << vacancies_pos(vacs_idx,0) << " " << vacancies_pos(vacs_idx,1) << " " << vacancies_pos(vacs_idx,2) << " " << vacancies_pos(vacs_idx,3) << " ]\n"; 

                }
                else {
                    //std::cout << "rank: " << rank << " vacancies_pos(vacs_idx) [ " << vacancies_pos(vacs_idx,0) << " " << vacancies_pos(vacs_idx,1) << " " << vacancies_pos(vacs_idx,2) << " " << vacancies_pos(vacs_idx,3) << " ]\n"; 
                    //std::cout << "rank: " << rank << " one proc reverse \n"; 

                    i1 = old_vac[0];
                    i2 = old_vac[1];
                    i3 = old_vac[2];
                    i4 = old_vac[3];
                    
                    //std::cout << "rank: " << rank << " new_site (set 0): " << "[ " << new_vac[0] << " " << new_vac[1] << " " << new_vac[2] << " " << new_vac[3] << " ]\n";
                    //std::cout << "rank: " << rank << " new site status: " << vacancies(new_vac[0],new_vac[1],new_vac[2],new_vac[3]) << " \n"; 
                    
                    // moving vacancy from bc site to bc site 
                    if (lattice == 3) {
                        /*---
                        BC EDGE MOVES
                        ---*/
                        
                        //std::cout << "rank: " << rank << " prev_site (set 1): " << "[ 1" << " " << old_vac[1] << " " << old_vac[2] << " " << old_vac[3] << " ]\n"; 
                        //std::cout << "rank: " << rank << " prev site status: " << vacancies(1,old_vac[1],i3,old_vac[3]) << " \n"; 
                                        
                        // switching occupancy for old and new site in lattice array ###
                        bc_sites((size_t)0, (size_t)new_vac[1], (size_t)new_vac[2], (size_t)new_vac[3]) = 1;
                        bc_sites((size_t)0, (size_t)old_vac[1], (size_t)old_vac[2], (size_t)old_vac[3]) = 0;

                        // switching occupancy for old and new site in vacancy and mobileion arrays ###
                        vacancies((size_t)new_vac[0], (size

                    }
                                
                    // moving vacancy from vertex site to vertex site 
                    else if (lattice == 2) {
                        /*---
                        VERTEX EDGE MOVES
                        ---*/
                        //std::cout << "rank: " << rank << " prev_site (set 1): " << "[ 0" << " " << old_vac[1] << " " << old_vac[2] << " " << old_vac[3] << " ]\n"; 
                        //std::cout << "rank: " << rank << " prev site status: " << vacancies(0,old_vac[1],old_vac[2],old_vac[3]) << " \n"; 
                        
                        // switching occupancy for old and new site in lattice array ###
                        vertex_sites((size_t)0, (size_t)new_vac[1], (size_t)new_vac[2], (size_t)new_vac[3]) = 1;
                        vertex_sites((size_t)0, (size_t)old_vac[1], (size_t)old_vac[2], (size_t)old_vac[3]) = 0;
                        
                        // switching occupancy for old and new site in vacancy and mobileion arrays ###
                        vacancies((size_t)new_vac[0], (size_t)new_vac[1], (size_t)new_vac[2], (size_t)new_vac[3]) = 0;
                        vacancies((size_t)0, (size_t)old_vac[1], (size_t)old_vac[2], (size_t)old_vac[3]) = 1;
                        
                        vacancies_pos(vacs_idx,0) = 0;
                    }
                                
                    // moving vacancy from bc site to vertex site 
                    else if (lattice == 1) {
                        /*---
                        BC MOVES
                        ---*/
                        //std::cout << "rank: " << rank << " prev_site (set 1): " << "[ 1" << " " << old_vac[1] << " " << old_vac[2] << " " << old_vac[3] << " ]\n";
                        //std::cout << "rank: " << rank << " prev site status: " << vacancies(1,old_vac[1],old_vac[2],old_vac[3]) << " \n";

                        // switching occupancy for old and new site in lattice array ###
                        vertex_sites((size_t)0, (size_t)new_vac[1], (size_t)new_vac[2], (size_t)new_vac[3]) = 1;
                        bc_sites((size_t)0, (size_t)old_vac[1], (size_t)old_vac[2], (size_t)old_vac[3]) = 0;

                        // switching occupancy for old and new site in vacancy and mobileion arrays ###
                        vacancies((size_t)new_vac[0], (size_t)new_vac[1], (size_t)new_vac[2], (size_t)new_vac[3]) = 0;
                        vacancies((size_t)1, (size_t)old_vac[1], (size_t)old_vac[2], (size_t)old_vac[3]) = 1;

                        vacancies_pos(vacs_idx,0) = 1; 
                    }
                                                
                    // moving vacancy from vertex site to bc site 
                    else if (lattice == 0) {
                        /*---
                        VERTEX MOVES
                        ---*/
                        //std::cout << "rank: " << rank << " prev_site (set 1): " << "[ 0" << " " << old_vac[1] << " " << old_vac[2] << " " << old_vac[3] << " ]\n"; 
                        //std::cout << "rank: " << rank << " prev site status: " << vacancies(0,old_vac[1],old_vac[2],old_vac[3]) << " \n"; 

                        // switching occupancy for old and new site in lattice array ###
                        bc_sites((size_t)0, (size_t)new_vac[1], (size_t)new_vac[2], (size_t)new_vac[3]) = 1;
                        vertex_sites((size_t)0, (size_t)old_vac[1], (size_t)old_vac[2], (size_t)old_vac[3]) = 0;

                        // switching occupancy for old and new site in vacancy and mobileion arrays ###
                        vacancies((size_t)new_vac[0], (size_t)new_vac[1], (size_t)new_vac[2], (size_t)new_vac[3]) = 0;
                        vacancies((size_t)0, (size_t)old_vac[1], (size_t)old_vac[2], (size_t)old_vac[3]) = 1;

                        vacancies_pos(vacs_idx,0) = 0;
                    }

                    vacancies_pos(vacs_idx,1) = (size_t)old_vac[1];
                    vacancies_pos(vacs_idx,2) = (size_t)old_vac[2];
                    vacancies_pos(vacs_idx,3) = (size_t)old_vac[3];
                    
                    
                }
            }

            prev_moves.pop_back();
            prev_newlocs.pop_back();
            prev_oldlocs.pop_back();
            prev_lattice.pop_back();
            prev_idxs.pop_back(); 
                    
        }        


        /**
        * @brief Rolls back moves involving parallel communication when an error 
        * encountered in simulation.
        *
        * This function reverses the most recent move by updating the occupancy of
        * lattice sites and adjusting the status of vacancies. It handles restoring 
        * the original occupancy in the native processor - the removal of vacancy
        * from adjacent processor domain is handled by reverse_move() in adjacent 
        * rank process 
        *
        * It also manages the rollback of associated data structures used for 
        * tracking previous moves.
        */
        void reverse_move_parallel() {
            int vacs_idx;
            int i1; int i2; int i3; int i4;

            vacs_idx = par_prev_idx.at(par_prev_idx.size()-1),0;
            
            //std::cout << "rank: " << rank << " parallel_prev_moves: \n";
            //print_2Dvector(par_prev_moves);

            i1 = (par_prev_moves.at((par_prev_moves.size()-1)).at(0));
            i2 = (par_prev_moves.at((par_prev_moves.size()-1)).at(1));
            i3 = (par_prev_moves.at((par_prev_moves.size()-1)).at(2));
            i4 = (par_prev_moves.at((par_prev_moves.size()-1)).at(3));
            
            //std::cout << "rank: " << rank << " parallel site (set 0): " << "[ " << i1 << " " << i2 << " " << i3 << " " << i4 << " ]\n";
            //std::cout << "rank: " << rank << " parallel old site status: " << vacancies(i1,i2,i3,i4) << " \n"; 
            //std::cout << "rank: " << rank << " parallel vacancies_pos(vacs_idx): [ " << vacancies_pos(vacs_idx,0) << " " << vacancies_pos(vacs_idx,1) << " " << vacancies_pos(vacs_idx,2) << " " << vacancies_pos(vacs_idx,3) << " ]\n"; 

            if (i1 == 1) {
                // moving vacancy from bc site to bc site 
                bc_sites((size_t)0, (size_t)i2, (size_t)i3, (size_t)i4) = 1;
            }
            else if (i1 == 0) {
                // moving vacancy from vertex site to vertex site 
                vertex_sites((size_t)0, (size_t)i2, (size_t)i3, (size_t)i4) = 1;
            }

            vacancies(i1,i2,i3,i4) = 0;

            //std::cout << "rank: " << rank << " parallel remove vacancies_pos.rows(): " << vacancies_pos.rows() << "\n";
            //vacancies_pos.print();
            
            vacancies_pos.remove_row(vacs_idx, rank);

            num_of_vacs--;
            //std::cout << "rank: " << rank << " par_prev_moves: \n";
            //print_2Dvector(par_prev_moves);

            // emptying out rollback data structures //    
            par_move_ticks.pop_back();
            par_prev_moves.pop_back();
            par_prev_idx.pop_back();

            //std::cout << "rank: " << rank << " post parallel remove vacancies_pos.rows(): " << vacancies_pos.rows() << "\n";
            //vacancies_pos.print();

            //std::cout << "rank: " << rank << " EOL\n";
        }


        /**
        * @brief Decrements the simulation clock based on previously stored time values.
        *
        * This function rolls back the simulation clock by subtracting the recorded
        * time intervals from the previous moves. It also clears the list of previous 
        * time values after adjusting the clock.
        */
        void deincrement_time() {
            print_1Dvector(prev_times);
            for (int i=0; i < (int)prev_times.size(); i++) {
                t = t - prev_times[i];
            }

            prev_times.clear();
        }


        /**
        * @brief Stores information about a move originating in the native process.
        *
        * This function manages the storage of move information based on whether the
        * move involved interprocessor communication. It updates the locations of 
        * vacancies and tracks the shifts made during the move. Two different storage
        * paths are utilized depending on whether the move is a parallel transfer or 
        * a standard move.
        *
        * @param idx The index of the move being stored.
        * @param parallel_transfer Indicates if the move involves parallel transfer.
        * @param vac_idx The index of the vacancy being moved.
        * @param lattice Optional parameter for lattice type, (default is -1).
        * @param shift A vector representing the shift of the move (default is {0,0,0}).
        * @param new_loc A vector representing the new location of the vacancy (default is {0,0,0,0}).
        * @param old_loc A vector representing the old location of the vacancy (default is {0,0,0,0}).
        */
        void store_move_info(int idx,  bool parallel_transfer, int vac_idx, int lattice = -1, std::vector<int> shift = {0,0,0}, const std::vector<int>& new_loc = {0,0,0,0}, const std::vector<int>& old_loc = {0,0,0,0}) {
            
            // POTENTIALLY ISSUE WITH STORING MOVE LOCATIONS FROM UPDATE_LATTICE() //

            std::vector<int> old_vac(4); 
            std::vector<int> new_vac(4); 

            // print statement to check for errors
            if (vac_idx >= (int)vacancies_pos.rows()) {
                //std::cout << "rank: " << rank << " out of bounds access in vacancies_pos: vacancies_pos.rows(): " 
                //<< vacancies_pos.rows() << " vs moves_vacs(idx,0): " << moves_vacs(idx,0) << " for lattice type: " << moves_lattice(idx,0) << "\n";
            }

            // keeping track of new location if parallel transfer
            if (parallel_transfer) { 
                for (int i=0; i<(int)new_vac.size(); i++)  {
                    new_vac[i] = new_loc[i];
                    old_vac[i] = old_loc[i];
                }
            }
            // entering filler values if null move
            else if (vac_idx == -1) {
                for (int i=0; i<(int)shift.size(); i++)  {
                    shift[i] = -1;
                } 
                for (int i=0; i<(int)new_vac.size(); i++)  {
                    new_vac[i] = new_loc[i];
                    old_vac[i] = old_loc[i];
                }
            }
            // standard move
            else {
                for (int i=0; i<(int)new_vac.size(); i++)  {
                    new_vac[i] = new_loc[i];
                    old_vac[i] = old_loc[i];
                }
            }

            //std::cout << "rank: " << rank << " store vac: [ " << new_vac[0] << " " << new_vac[1] << " " << new_vac[2] << " " << new_vac[3] << " ]" << "\n";
            //std::cout << "rank: " << rank << " store move: [ " << shift[0] << " " << shift[1] << " " << shift[2] << " ]" << "\n";

            if ((int)prev_newlocs.size() >= 2) {

                prev_moves[0] = prev_moves[1];
                prev_moves[1] = shift;

                prev_newlocs[0] = prev_newlocs[1];
                prev_newlocs[1] = new_vac;

                prev_oldlocs[0] = prev_oldlocs[1];
                prev_oldlocs[1] = old_vac;

                prev_idxs[0] = prev_idxs[1];
                prev_idxs[1] = vac_idx;

                prev_lattice[0] = prev_lattice[1];
                prev_lattice[1] = lattice;
            }
            else {
                prev_moves.push_back(shift);
                prev_newlocs.push_back(new_vac);
                prev_oldlocs.push_back(old_vac);
                prev_idxs.push_back(vac_idx);
                prev_lattice.push_back(lattice);
            }

            //std::cout << "rank: " << rank << " post prev_moves: \n";
            //print_2Dvector(prev_moves);
            //std::cout << "rank: " << rank << " post prev_newlocs: \n";
            //print_2Dvector(prev_newlocs);
            //std::cout << "rank: " << rank << " post prev_oldlocs: \n";
            //print_2Dvector(prev_oldlocs);
            //std::cout << "rank: " << rank << " post prev_lattice: \n";
            //print_1Dvector(prev_lattice);

            //std::cout << "rank: " << rank << " prev_move_type: \n";
            //print_1Dvector(prev_move_type);
        }

        /**
        * @brief Stores information about a move involving parallel transfer from another process.
        *
        * This function captures details of a move that was transferred in parallel, 
        * updating the necessary vacancy and move tracking structures.
        *
        * @param parallel_buffer A vector containing the details of the vacancy involved in the move.
        * @param move_ticks The number of ticks associated with the move.
        * @param vac_idx The index of the vacancy that was moved.
        */
        void store_parallel_info(const std::vector<int>& parallel_buffer, int move_ticks, int vac_idx) {
            std::vector<int> vac(4); 
            vac[0] = parallel_buffer[0];
            vac[1] = parallel_buffer[1];
            vac[2] = parallel_buffer[2];
            vac[3] = parallel_buffer[3];
            
            par_prev_moves.push_back(vac);
            par_move_ticks.push_back(move_ticks);
            par_prev_idx.push_back(vac_idx);
        }

        /**
        * @brief Stores the time elapsed during the previous iteration of the KMC process.
        *
        * This function updates the time tracking structures by storing the time increment 
        * from the last Kinetic Monte Carlo (KMC) process iteration. If there are 
        * already two previous time entries, the oldest one is overwritten.
        *
        * @param timestep The time increment to be stored.
        */
        void store_time_incr(double timestep) {
            if ((int)prev_times.size() >= 2) {
                prev_times[1] = prev_times[0];
                prev_times[0] = timestep;
            }          
        }


        /**
        * @brief Helper function to reverse the previous two moves.
        *
        * This function checks the type of each previous move (in-lattice or parallel transfer)
        * and calls the appropriate reverse function to undo the moves. After reversing,
        * it clears the recorded previous move types and ticks.
        */
        void reverse_moves_wrapper() {
            std::cout << "rank: " << rank << " reverse_moves_wrapper: \n";
            std::cout << "rank: " << rank << " prev_move_type: \n";
            print_1Dvector(prev_move_type);
            std::cout << "rank: " << rank << " pre prev_move_type_ticks: \n";
            print_1Dvector(prev_move_type_ticks);
            std::cout << "rank: " << rank << " prev_newlocs: \n";
            print_2Dvector(prev_newlocs);

            for (int i=(prev_move_type.size()-1); i>=0; i--) { 
                std::cout << "rank: " << rank << " prev_move_type[i]: " << prev_move_type[i] << "\n";
                if (prev_move_type[i] == 0) reverse_move(false); // reversing in-lattice move with no parallel transfer
                if (prev_move_type[i] == 1) reverse_move(true);  // reversing move with sending parallel transfer
                if (prev_move_type[i] == 2) reverse_move_parallel(); // reversing move with receiving parallel transfer
            }

            prev_move_type.clear();
            prev_move_type_ticks.clear();
        }

        /**
        * @brief Receives a move from another process using non-blocking communication.
        *
        * This function updates the lattice with the received move details while checking for 
        * potential conflicts with existing vacancies. It extracts the new location and previous 
        * location data from the input buffer, performs necessary checks, and stores relevant 
        * information about the move.
        *
        * @param new_loc_buffer A vector containing the new location and metadata of the move.
        * @param move_ticks The number of ticks associated with the received move.
        * @return true if an interboundary conflict occurred, false otherwise.
        */
        bool recieve_move_parallel(const std::vector<int>& new_loc_buffer, int move_ticks) {
            //
            assert( (((new_loc_buffer)[5] == 0) || ((new_loc_buffer)[5] == 1)) && (((new_loc_buffer)[0] == 0) || ((new_loc_buffer)[0] == 1)) );
            assert( (((new_loc_buffer)[6] >= 0) && ((new_loc_buffer)[6] < sublattice_dim[0])) && (((new_loc_buffer)[1] >= -1) && ((new_loc_buffer)[1] < (sublattice_dim[0] + 1))) );
            assert( (((new_loc_buffer)[7] >= 0) && ((new_loc_buffer)[7] < sublattice_dim[1])) && (((new_loc_buffer)[2] >= -1) && ((new_loc_buffer)[2] < (sublattice_dim[1] + 1))) );
            assert( (((new_loc_buffer)[8] >= 0) && ((new_loc_buffer)[8] < sublattice_dim[2])) && (((new_loc_buffer)[3] >= -1) && ((new_loc_buffer)[3] < (sublattice_dim[2] + 1))) );
            //
            
            //std::cout << "rank: " << rank << " move_ticks: " << move_ticks << " parallel process new_loc\n";

            //std::cout << "rank: " << rank << " recieve_move_parallel()\n";

            int old_proc = (new_loc_buffer[4]);
            int i_old = (new_loc_buffer[5]);
            int j_old = (new_loc_buffer[6]);
            int k_old = (new_loc_buffer[7]);
            int l_old = (new_loc_buffer[8]);

            /*
            make cases for which proc neighbors to update depending on the location of the
            old site
            */
            int i = (new_loc_buffer[0]);
            int j = (new_loc_buffer[1]);
            int k = (new_loc_buffer[2]);
            int l = (new_loc_buffer[3]);

            bool interboundary_conflict = false;
            
            std::vector<size_t> x_dims = proc_pos_x_neighbors.size_vec;
            std::vector<size_t> y_dims = proc_pos_y_neighbors.size_vec;

            //std::cout << "rank: " << rank << " old  : [ " << i_old << " " << j_old << " " << k_old << " " << l_old << " ]\n";
            //std::cout << "rank: " << rank << " new loc: [ " << i << " " << j << " " << k << " " << l << " ]\n";
            
            
            if (vacancies((size_t)i, (size_t)j, (size_t)k, (size_t)l) == 1) {
                std::cout << "rank: " << rank << " ERROR: Vacancy overwriting already existing vacancy in interprocessor communication \n";
                interboundary_conflict = true;
                comm_boundary_conflict();
            }
            else {
                /* store previous two moves */

                size_t rows = vacancies_pos.rows();
                size_t cols = vacancies_pos.cols();
                vacancies_pos.reshape(rows+1, cols, rank);
                //std::cout << "rank: " << rank << " parallel recv site status: " << vacancies(i,j,k,l) << " \n"; 
                //std::cout << "rank: " << rank << " parallel recv vacancies_pos(vacs_idx): [ " << vacancies_pos(rows,0) << " " << vacancies_pos(rows,1) << " " << vacancies_pos(rows,2) << " " << vacancies_pos(rows,3) << " ]\n"; 

                
                store_parallel_info(new_loc_buffer, move_ticks, rows);
                
                vacancies((size_t)i, (size_t)j, (size_t)k, (size_t)l) = 1;
                if (i == 0) vertex_sites((size_t)0, (size_t)j, (size_t)k, (size_t)l) = 0;
                else if (i == 1) bc_sites((size_t)0, (size_t)j, (size_t)k, (size_t)l) = 0;
                //std::cout << "rank: " << rank << " rows: " << rows << "\n";
                //std::cout << "rank: " << rank << " vacancies_pos.rows(): " << vacancies_pos.rows() << "\n";
                //std::cout << "rank: " << rank << " vacancies_pos.cols(): " << vacancies_pos.cols() << "\n";

                vacancies_pos(rows,0) = i;
                vacancies_pos(rows,1) = j;
                vacancies_pos(rows,2) = k;
                vacancies_pos(rows,3) = l;
                // std::cout << "i: " << i << "  n j: " << j << " k: " << k << " l: " << l <<  "\n";

                num_of_vacs ++;
                prev_move_type.push_back(2);
                prev_move_type_ticks.push_back(move_ticks);

                //std::cout << "rank: " << rank << " post parallel recv site status: " << vacancies(i,j,k,l) << " \n"; 
                //std::cout << "rank: " << rank << " post parallel recv vacancies_pos(vacs_idx): [ " << vacancies_pos(rows,0) << " " << vacancies_pos(rows,1) << " " << vacancies_pos(rows,2) << " " << vacancies_pos(rows,3) << " ]\n"; 

            }

            return interboundary_conflict;
        }



        /**
        * @brief Receive communication (ghost site update, lattice site update, conflict) using blocking communication.
        *
        * This function handles the reception of messages from other processes,
        * including updates to ghost sites, lattice sites, and conflicts.
        *
        * @param move_ticks The current iteration for moves in the simulation.
        */
        void receive_parallel_comm_helper(int move_ticks) {

            int test_flag1 = 0;
            MPI_Status status1;
            std::vector<int> new_loc_buffer1(9);
            std::vector<int> new_loc_buffer2(9);
            std::vector<int> new_loc_buffer3(9);
            

            std::vector<bool> stop_ghost_Array(num_procs, 0);
            std::vector<bool> stop_par_Array(num_procs, 0);
            std::vector<bool> stop_rev_array(num_procs, 0);
            std::vector<bool> stop_conflict_array(num_procs, 0);
            stop_ghost_Array[rank] = 1;
            stop_par_Array[rank] = 1;
            stop_rev_array[rank] = 1;
            stop_conflict_array[rank] = 1;
            bool stop_ghost = 0;
            bool stop_par = 0;
            bool stop_rev = 0;
            bool need_reverse = 0;
            bool stop_conflict = 0;
            bool interboundary_conflict = false;

            while ((!stop_par) || (!stop_ghost)) {
                //std::cout << "rank: " << rank  << " stop_ghost_Array " << "\n";
                //print_1Dvector(stop_ghost_Array);

                //std::cout << "rank: " << rank  << " stop_par_Array " << "\n";
                //print_1Dvector(stop_par_Array);

                stop_par = 1;
                stop_ghost = 1;
                stop_rev = 1;

                status1.MPI_TAG = 0;
                
                MPI_Probe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status1);

                // test for both message 1 and message 2 simulatneously in if statement
                if (status1.MPI_TAG == 1) {
                    MPI_Recv(new_loc_buffer1.data(), 9, MPI_INT, status1.MPI_SOURCE, 1, MPI_COMM_WORLD, &status1);
                    interboundary_conflict = recieve_move_parallel(new_loc_buffer1, move_ticks);
                }
                
                if (status1.MPI_TAG == 2) {
                    MPI_Recv(new_loc_buffer2.data(), 9, MPI_INT, status1.MPI_SOURCE, 2, MPI_COMM_WORLD, &status1);
                    ghost_site_recieve(new_loc_buffer2, true);
                }
                
                if (status1.MPI_TAG == 3) {
                    MPI_Recv(new_loc_buffer3.data(), 9, MPI_INT, status1.MPI_SOURCE, 3, MPI_COMM_WORLD, &status1);
                    ghost_site_recieve(new_loc_buffer3, false);
                }
                
                if (status1.MPI_TAG == par_done_tag) {
                    MPI_Recv(NULL, 0, MPI_CHAR, status1.MPI_SOURCE, par_done_tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                    stop_par_Array[status1.MPI_SOURCE] = 1;
                }
                
                for (int j = 0; j < (int)stop_par_Array.size(); j++) {
                    stop_par &= stop_par_Array[j];
                }
                
                if (status1.MPI_TAG == ghost_done_tag) {
                    MPI_Recv(NULL, 0, MPI_CHAR, status1.MPI_SOURCE, ghost_done_tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                    stop_ghost_Array[status1.MPI_SOURCE] = 1;
                }
                
                for (int j = 0; j < (int)stop_ghost_Array.size(); j++) {
                    stop_ghost &= stop_ghost_Array[j];
                }
            }

            MPI_Request request1;
            for (int new_proc=0; new_proc<num_procs; new_proc++) { if (new_proc != rank) MPI_Isend( NULL, 0, MPI_CHAR, new_proc, conflict_done_flag, MPI_COMM_WORLD, &request1); }
                
            //std::cout << "rank: " << rank  << " exited 1 \n";

            remove_old_par_moves(move_ticks);

            MPI_Barrier(MPI_COMM_WORLD);
            
            //std::cout << "rank: " << rank  << " pre conflict \n";
            while ((!stop_conflict)) {
                stop_conflict = 1;
                status1.MPI_TAG = 0;
                MPI_Probe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status1);
                
                if (status1.MPI_TAG == 6) {
                    MPI_Recv(NULL, 0, MPI_CHAR, status1.MPI_SOURCE, 6, MPI_COMM_WORLD, &status1);
                    need_reverse = 1;
                    //std::cout << "rank: " << rank  << " stop_conflict_array, status1.MPI_SOURCE: " << status1.MPI_SOURCE << "\n";
                    //print_1Dvector(stop_conflict_array);

                }
                if (status1.MPI_TAG == conflict_done_flag) {
                    MPI_Recv(NULL, 0, MPI_CHAR, status1.MPI_SOURCE, conflict_done_flag, MPI_COMM_WORLD, &status1);
                    stop_conflict_array[status1.MPI_SOURCE] = 1;
                    //std::cout << "rank: " << rank  << " stop_conflict_array, status1.MPI_SOURCE: " << status1.MPI_SOURCE << "\n";
                    //print_1Dvector(stop_conflict_array);

                }
                for (int j = 0; j < (int)stop_conflict_array.size(); j++) {
                    stop_conflict &= stop_conflict_array[j];
                }
            }
            //std::cout << "rank: " << rank  << " past conflict \n";

            MPI_Barrier(MPI_COMM_WORLD);
            
            if ((need_reverse) || (interboundary_conflict)) {
                deincrement_time();
                reverse_moves_wrapper();
            }
        }

        /**
        * @brief Wrapper method containing initialization of all variables in the system, timer, and calls to update state.
        *
        * This method initializes necessary variables and manages the iterative
        * process of the Kinetic Monte Carlo simulation. It handles time limits,
        * updates the state of the system, and logs results to specified data 
        * structures and folders.
        *
        * @param time_lim The time limit for the simulation.
        * @param start The starting time point of the simulation.
        * @param folder The folder path for output data.
        * @param iteration The current iteration of the simulation.
        * @param rates_i An index representing the rate being processed.
        * @return A structure containing the results of the KMC simulation.
        */
        lattice_return_struct new_kmc_iterator(double time_lim, std::chrono::system_clock::time_point start, std::string folder, int iteration, int rates_i) {
            fprintf(stdout, "%s", "beginning kmc iterations \n\n"); 

            // INITIALIZING VARIABLES PRIOR TO BEGINNING FIRST KMC STEP //

            std::chrono::system_clock::time_point end; // current (real) clock time
            std::chrono::duration<double> elapsed_seconds; // elapesed (simulated) time in simulation

            std::vector<double> timesteps; // time elapsed at each step
            std::vector<int> move_counts(6); // each type of move propogated by simulation 
            std::vector<double> time_count(6); // time elapsed by each type of move
            int rand_idx; // index of move selected
            double timestep;
            std::vector<double> all_times; // vector containing trajectory of time elapsed by each type of move
            parallel_get_actions(); // updating list of moves in system

            int move_ticks = 0;
            double old_time;

            int prev_num_vacs = 0;
            int curr_num_vacs = 0;
            bool get_rank = true;

            std::ostringstream ss;


            Matrix<int> only_vacancies = vacancies.nonzero();
            std::cout << "rank: " << rank << " move_ticks: " << move_ticks << " vacancies_pos.rows(): " << vacancies_pos.rows() << " vs only_vacancies.rows(): " << only_vacancies.rows() << "\n";
            
            MPI_Barrier(MPI_COMM_WORLD);
            if (rank == 0) {
                std::cout << "rank: " << rank << " vacancies_pos.rows: " << vacancies_pos.rows() << "\n";
                std::cout << "rank: " << rank << " num of vacancies: " << only_vacancies.rows() << "\n";
                std::cout << "rank: " << rank << " vacancies_pos: \n";
                vacancies_pos.print();
                std::cout << "rank: " << rank << " only_vacancies: \n";
                only_vacancies.print();
                if (vacancies_pos.rows() != only_vacancies.rows()) { 
                    std::cout << "rank: " << rank << " vacancies_pos: \n";
                    vacancies_pos.print();
                    std::cout << "rank: " << rank << " only_vacancies: \n";
                    only_vacancies.print();
                    //std::cout << "rank: " << rank << " vacancies_pos.rows(): " << vacancies_pos.rows() << " mismatch with only_vacancies.rows(): " << only_vacancies.rows() << "\n"; 
                    //Matrix<int> difference = comparison(vacancies_pos, only_vacancies);
                    //std::cout << "rank: " << rank << " difference: \n";
                    //difference.print();
                    exit(0);
                }
            }
            MPI_Barrier(MPI_COMM_WORLD);
            if (rank == 1) {
                std::cout << "rank: " << rank << " vacancies_pos.rows: " << vacancies_pos.rows() << "\n";
                std::cout << "rank: " << rank << " num of vacancies: " << only_vacancies.rows() << "\n";
                std::cout << "rank: " << rank << " vacancies_pos: \n";
                vacancies_pos.print();
                std::cout << "rank: " << rank << " only_vacancies: \n";
                only_vacancies.print();
                if (vacancies_pos.rows() != only_vacancies.rows()) { 
                    std::cout << "rank: " << rank << " vacancies_pos: \n";
                    vacancies_pos.print();
                    std::cout << "rank: " << rank << " only_vacancies: \n";
                    only_vacancies.print();
                    //std::cout << "rank: " << rank << " vacancies_pos.rows(): " << vacancies_pos.rows() << " mismatch with only_vacancies.rows(): " << only_vacancies.rows() << "\n"; 
                    //Matrix<int> difference = comparison(vacancies_pos, only_vacancies);
                    //std::cout << "rank: " << rank << " difference: \n";
                    //difference.print();
                    //exit(0);
                }
            }
            MPI_Barrier(MPI_COMM_WORLD);
            if (rank == 2) {
                std::cout << "rank: " << rank << " vacancies_pos.rows: " << vacancies_pos.rows() << "\n";
                std::cout << "rank: " << rank << " num of vacancies: " << only_vacancies.rows() << "\n";
                std::cout << "rank: " << rank << " vacancies_pos: \n";
                vacancies_pos.print();
                std::cout << "rank: " << rank << " only_vacancies: \n";
                only_vacancies.print();
                if (vacancies_pos.rows() != only_vacancies.rows()) { 
                    std::cout << "rank: " << rank << " vacancies_pos: \n";
                    vacancies_pos.print();
                    std::cout << "rank: " << rank << " only_vacancies: \n";
                    only_vacancies.print();
                    //std::cout << "rank: " << rank << " vacancies_pos.rows(): " << vacancies_pos.rows() << " mismatch with only_vacancies.rows(): " << only_vacancies.rows() << "\n"; 
                    //Matrix<int> difference = comparison(vacancies_pos, only_vacancies);
                    //std::cout << "rank: " << rank << " difference: \n";
                    //difference.print();
                    //exit(0);
                }
            }
            MPI_Barrier(MPI_COMM_WORLD);
            if (rank == 3) {
                std::cout << "rank: " << rank << " vacancies_pos.rows: " << vacancies_pos.rows() << "\n";
                std::cout << "rank: " << rank << " num of vacancies: " << only_vacancies.rows() << "\n";
                std::cout << "rank: " << rank << " vacancies_pos: \n";
                vacancies_pos.print();
                std::cout << "rank: " << rank << " only_vacancies: \n";
                only_vacancies.print();
                if (vacancies_pos.rows() != only_vacancies.rows()) { 
                    std::cout << "rank: " << rank << " vacancies_pos: \n";
                    vacancies_pos.print();
                    std::cout << "rank: " << rank << " only_vacancies: \n";
                    only_vacancies.print();
                    //std::cout << "rank: " << rank << " vacancies_pos.rows(): " << vacancies_pos.rows() << " mismatch with only_vacancies.rows(): " << only_vacancies.rows() << "\n"; 
                    //Matrix<int> difference = comparison(vacancies_pos, only_vacancies);
                    //std::cout << "rank: " << rank << " difference: \n";
                    //difference.print();
                    //exit(0);
                }
            }
            
            MPI_Barrier(MPI_COMM_WORLD);

            while (move_ticks < 250000) {
                std::cout << "rank: " << rank  << " move_ticks: " << move_ticks  << "\n";


                if ( rate_cumsum.size() != 0) {
                    end = std::chrono::system_clock::now(); 
                    rand_idx = get_idx();
                    new_update_lattice(rand_idx, move_ticks);
                    move_counts[moves_lattice[rand_idx][0]] ++;
                    timestep = new_random_times();
                    time_count[moves_lattice[rand_idx][0]] += timestep;
                    store_time_incr(timestep);
                }
                else {
                    std::cout << "rank: " << rank << " only communicate_rates " << "\n";
                    communicate_rates(); 
                    timestep = new_random_times();
                }
                MPI_Barrier(MPI_COMM_WORLD);

                MPI_Request request1;
                for (int new_proc=0; new_proc<num_procs; new_proc++) { if (new_proc != rank) MPI_Isend( NULL, 0, MPI_CHAR, new_proc, par_done_tag, MPI_COMM_WORLD, &request1); }
                MPI_Request request2;
                for (int new_proc=0; new_proc<num_procs; new_proc++) { if (new_proc != rank) MPI_Isend( NULL, 0, MPI_CHAR, new_proc, ghost_done_tag, MPI_COMM_WORLD, &request2); }
                
                
                std::cout << "rank: " << rank << " t: " << t << "\n";
                
                MPI_Barrier(MPI_COMM_WORLD);
                
                
                receive_parallel_comm_helper(move_ticks);

                MPI_Barrier(MPI_COMM_WORLD);
               
                if (move_ticks % 1000 == 0) {
                    Matrix<int> only_vacancies = vacancies.nonzero(); // configuration of vacancies at current timestep
                    write_output_parallel(only_vacancies, total_dims, folder, proc_dims[0], proc_dims[1], num_procs, rank, iteration, rates_i, move_ticks, t, get_rank);
                }
            
                only_vacancies = vacancies.nonzero();
                //std::cout << "rank: " << rank << " move_ticks: " << move_ticks << " vacancies_pos.rows(): " << vacancies_pos.rows() << " vs only_vacancies.rows(): " << only_vacancies.rows() << "\n";
                //std::cout << "rank: " << rank << " test to see running\n";
                
                
                MPI_Barrier(MPI_COMM_WORLD);


                parallel_get_actions();
                
                fflush(stdout);
                MPI_Barrier(MPI_COMM_WORLD);

                //std::cout << "rank: " << rank << " t: " << t << "\n";
               
                
                curr_num_vacs = sum_vacs_allprocs(only_vacancies, total_dims, proc_dims[0], proc_dims[1], num_procs, rank);

                if (rank == 0) std::cout << "Total number of vacancies in simulation -- curr_num_vacs: " << curr_num_vacs << "\n"; 
                if ((rank == 0) && (prev_num_vacs != curr_num_vacs)) {
                    std::cout << "Change in total number of vacancies in simulation -- curr_num_vacs: " << curr_num_vacs << " prev_num_vacs: " << prev_num_vacs << " move_ticks: " << move_ticks << "\n"; 
                    if (move_ticks != 0) {
                        MPI_Barrier(MPI_COMM_WORLD);
                        exit(0);
                    }
                }
                prev_num_vacs = curr_num_vacs;

                elapsed_seconds = end-start;
                t += timestep;
                move_ticks ++;
                old_time = t;                

                // terminating simulation after real-time limit reached
                if (elapsed_seconds.count() >= 172800) {
                    t = time_lim + 1;
                    std::cout << "elapsed_seconds: " << elapsed_seconds.count() << "\n";
                    std::cout << "exiting\n";
                    std::cout << "t: " << t << "\n";
                    print_1Dvector(move_counts);
                    break;
                }
                MPI_Barrier(MPI_COMM_WORLD);
            }

            all_times.push_back(t);

            std::cout << "rank: " << rank << " move_ticks: " << move_ticks << "\n";

            lattice_return_struct output_vals(move_counts, time_count, all_times);

            MPI_Barrier(MPI_COMM_WORLD);
            return output_vals;
        }
};

/*---------------------------------------------------------------------------*/

/**
 * @brief Creates entries in the rate catalog corresponding to all configurations of "len" sites with m vacancies.
 *
 * This function generates all binary strings of length "len" with exactly m zeros,
 * representing configurations of vacancies. The output is a vector of integers,
 * where each integer corresponds to a specific binary configuration.
 *
 * @param len The total length of the binary string (total number of sites).
 * @param m The number of zeros in the binary string (vacancies).
 * @param size The number of atomic types in the system used to create the encoding.
 * @return A vector of integers representing all binary configurations of length "len"
 *         with m zeros.
 */
std::vector<int> bin_m_zeros(int len, int m, int size) {
    int num_combos = NCR(len, m);
    std::vector<int> vec(num_combos);
    int smallest_bin = 0;
    unsigned int t=0;
    unsigned int v=0;
    unsigned int w=0;

    for (int i=0; i<(len-m); i++) {
        smallest_bin += exp_int((int)(size), i);
    }
    v = smallest_bin;

    for (int j=0; j<(num_combos-1); j++) {
        t = (v | (v - 1)) + 1;
        w = t | ((((t & -t) / (v & -v)) >> 1) - 1);
        vec[j] = v;
        v = w;
    }

    vec[((int)vec.size() - 1)] = v;

    return vec;
}

/**
 * @brief Converts an inputted base-m string to an integer.
 *
 * This function takes a string representation of a number in base-m,
 * where each character represents a digit in that base, and converts it
 * to an integer.
 *
 * @param config The input string representing the number in base-m.
 * @param size The base (m) of the input string.
 * @return The integer value of the input string in base-m.
 */
int base_m_to_int(std::string config, int size) {

    int result = 0;
    std::vector<char> toks = split_by_char(config);
    char tok;


    for (int i=0; i<(int)toks.size(); i++) {
        tok = toks[i];
        result += (int)((tok) - '0') * exp_int(size, i);
    }
        
    return result;
}

/**
 * @brief Reads a rate catalog file and creates a catalog of rates, allowed moves,
 *        and configurations corresponding to regions.
 *
 * This function processes an input file that contains information about atom types,
 * configurations, and associated energies. It constructs and returns a structured
 * catalog (Matrix) that includes allowed configurations and rates for migration.
 *
 * @param catalogfile The path to the rate catalog file.
 * @param atype_list A list of atom types specified for use in the catalog.
 * @return A `ratecatalog_struct` containing the catalogs of configurations,
 *         energies, and regions.
 *
 * @note The expected file format is:
 *       - Line 0: Atom types (e.g., "1:sodium, 2:LLZO, etc...")
 *       - Line 1: Number of configurations (e.g., "#configs: m")
 *       - Subsequent lines: Configurations and corresponding energies.
 */
ratecatalog_struct updated_create_ratecatalog(std::string catalogfile, std::vector<int> atype_list) {
    
    //FILE FORMAT:
    //line 0: 1:sodium, 2:LLZO, etc... 
    //line 1: #configs: m
    //line 1: 1020420001:" migrationE #config of atom type 1,0,2,0,4 for 
    //                                "northmost" atom, clockwise
    //...
    //line n: 
    std::fstream cat_file;

    cat_file.open(catalogfile);
    std::vector<std::string> lines;
    std::string line;
    std::string output;

    if (cat_file.is_open()) {
        while ( getline (cat_file,line) )
        {
            lines.push_back(line);
        }
        cat_file.close();
    }

    // getting atom types from line 0 ###
    std::string typeline = lines[0];

    std::vector<std::string> types = tokenizer(typeline, " ");

    types = slice_1Dvec_str_inp(types, 1, (int)types.size());


    std::string atype;
    std::vector<std::string> key_and_val;
    std::map<int, std::string> rate_typedict;

    for (int i=0; i<(int)types.size(); i++) {
        atype = types[i];
        std::cout << "type " << atype << "\n";
        key_and_val = tokenizer(atype, ":");
        rate_typedict[std::stoi(key_and_val[0])] = key_and_val[1];
    }

    // checking that the atom types specified in input file is subset of types specified in rate catalog ###
    if  ((atype_list).size() < rate_typedict.size()) {
        output = "Atom type not found in rate catalog";
        printf("%s", output.c_str());
    }
    
    lines = slice_1Dvec_str_inp(lines, 1, (int)lines.size());

    // making arrays for configs and energies ###
    // dft_energies = np.zeros((2, num_of_configs))
    std::vector< std::vector< std::vector<double> > > all_energies;
    std::vector< std::vector<int> > all_configs;

    // getting atom config energies into catalog ###
    int catalog = 0;
    std::vector<int> idxs;
    int lines_read = 0;
    int j;
    int k;
    int encoded_config;
    int num_of_configs;
    std::vector<int> configs;
    std::vector< std::vector<double> > dft_energies;
    std::vector<std::string> toks;
    std::string tok;
    std::vector< std::string> atom_strs;
    std::string atom;
    std::string configline;
    std::vector<std::string> atom_info;
    std::string atom_type;
    int atom_idx;
    std::vector<std::string> region_num_info;
    int region_num=0;
    std::vector<Region> regions;
    std::vector<int> unsorted_idxs;


    for (int i=0; i<(int)lines.size(); i++) {
        configline = lines[i];
        //std::cout << "configline: " << configline << "\n";

        if (configline.find("stop catalog") != std::string::npos) {
            
            auto comparator = [configs](int idx1, int idx2) {
                return configs[idx1] < configs[idx2];
            };

            std::sort(unsorted_idxs.begin(), unsorted_idxs.end(), comparator);
            configs = reorder_inp(configs, unsorted_idxs);

            dft_energies[0] = (reorder_inp(dft_energies[0], unsorted_idxs));
            dft_energies[1] = (reorder_inp(dft_energies[1], unsorted_idxs));

            all_configs.push_back(configs);
            all_energies.push_back(dft_energies);

            catalog ++;

            for (int i1=0; i1<(int)dft_energies.size(); i1++) {
                for (int i2=0; i2<(int)dft_energies[0].size(); i2++) {  
                    dft_energies[i1][i2] = 0;
                }
            }
        }
        
        else if (configline == " ") {}

        else if (configline == "\n") {}
                    
        else if (configline.find("starting catalogs") != std::string::npos) {
        }

        else if (configline.find("/") != std::string::npos) {}

        else if (configline.find("start catalog") != std::string::npos) {
            j = 0;
        }

        else if (configline.find("count") != std::string::npos) {
            // getting num of configs###
            toks = tokenizer(configline, ": ");
            num_of_configs = std::stoi(toks[1]);        
            configs.resize(num_of_configs,0);
            idxs.resize(num_of_configs,0);
            dft_energies.resize(2);
            dft_energies[0].resize(num_of_configs,0);
            dft_energies[1].resize(num_of_configs,0);
            unsorted_idxs.resize(num_of_configs,0);

            for (int idx=0; idx<(int)unsorted_idxs.size(); idx++) {unsorted_idxs[idx] = idx;}
        }

        else if (configline.find("vac") != std::string::npos) {}

        else if (configline.find("config:") != std::string::npos) {
            toks = tokenizer(configline, " ");
            toks = slice_1Dvec_str_inp(toks, 1, (int)toks.size());
            k = j;

            for (int l=0; l<(int)toks.size(); l++) {
                tok = toks[l];
                encoded_config = base_m_to_int(tok, (int)atype_list.size()); 
                configs[k] = encoded_config;
                k ++;
            }
        }
            
        else if (configline.find("L:") != std::string::npos) {
            toks = tokenizer(configline, " ");
            toks = slice_1Dvec_str_inp(toks, 1, (int)toks.size());
            k = j;

            for (int l=0; l<(int)toks.size(); l++) {
                tok = toks[l];
                dft_energies[0][k] = std::stof(tok);
                k ++;
            }
        }

        else if (configline.find("R:") != std::string::npos) {
            toks = tokenizer(configline, " ");
            toks = slice_1Dvec_str_inp(toks, 1, (int)toks.size());

            k = j;
            for (int l=0; l<(int)toks.size(); l++) {
                tok = toks[l];
                dft_energies[1][k] = std::stof(tok);
                k ++;
            }

            j = k;
        }

        else if (configline.find("finishing catalogs") != std::string::npos) {break;}
        lines_read ++;
    }

    lines = slice_1Dvec_str_inp(lines, lines_read, (int)lines.size());

    // generating and assigning energies for various migration directions in all regions ###

    std::vector< std::vector<double> > atom_e;
    std::vector< std::vector< std::vector<double> > > region_e;
    std::vector< std::vector< std::vector< std::vector<double> > > > regions_catalog;
    lines_read = 0;
    std::vector<double> catalog_out;


    for (int i=0; i<(int)lines.size(); i++) {
        
        configline = lines[i];
        //std::cout << "configline: " << configline << "\n";

        if (configline.find("stop regions") != std::string::npos) {
            lines.pop_back();
            break;
        }

        else if ((str_isalpha(configline)) || (configline == "\n")) {}

        
        else if (configline.find("starting regions") != std::string::npos) {}


        else if (configline.find("region #") != std::string::npos) {
            region_num_info = tokenizer(configline, " ");
            region_num = std::stoi(region_num_info[2]);
            std::cout << "region_num: " << region_num << "\n";
            if (region_num > 1) {regions_catalog.push_back(region_e);}
            region_e = vect_create_3D_float(0,0,0);
        }

        else if (configline.find("atom") != std::string::npos) {
            atom_strs = tokenizer(configline, " ");
            atom = atom_strs[1];
            atom_info = tokenizer(atom, ":");
            atom_type = atom_info[0];
            atom_idx = std::stoi(atom_info[1]);

            if (atom_idx > 1) {region_e.push_back(atom_e);}

            atom_e = vect_create_2D_float(0,0,0);
        }

        else if ((configline.find("ver:") != std::string::npos ) || (configline.find("bc:") != std::string::npos) || (configline.find("ver_edge:") != std::string::npos) || (configline.find("bc_edge:") != std::string::npos)) {
            toks = tokenizer(configline, ": ");
            catalog_out = create_vec_1D_float((int)toks.size());

            toks = slice_1Dvec_str_inp(toks, 1, (int)toks.size());

            for (j=0; j<(int)toks.size(); j++) {catalog_out[j] = std::stof(toks[j]);}

            catalog_out = slice_1Dvec_float(catalog_out, 1, (int)catalog_out.size());
            atom_e.push_back(catalog_out);
        }

        else if (configline.find("finishing regions") != std::string::npos ) {
            region_e.push_back(atom_e);
            regions[(region_num-1)].energies = region_e;
            regions_catalog.push_back(region_e);
        }

        lines_read ++;
    }

    
    ratecatalog_struct output_vals(all_configs, all_energies, regions_catalog, region_num);
    

    return output_vals;
}

/**
 * @brief Creates a region object based on the provided information.
 *
 * This function constructs a `Region` object using data extracted from the 
 * input vector. It supports different region types, such as "GB" (Grain Boundary)
 * and "BLOCK" (rectangular prism). The parameters for the region are parsed
 * from the input string vector.
 *
 * @param info A vector of strings containing region information:
 *             - info[0]: Region ID (e.g., "id:1")
 *             - info[1]: Region type (e.g., "GB" or "BLOCK")
 *             - Subsequent elements contain parameters relevant to the region type.
 * @return A pointer to the newly created `Region` object.
 *
 * @note The expected format for `info` elements varies based on the region type.
 *       - For "GB": 
 *         - info[2-4]: First set of parameters (3 integers)
 *         - info[5-7]: Second set of parameters (3 integers)
 *       - For "BLOCK":
 *         - info[2-3]: X-dimension parameters (2 integers)
 *         - info[4-5]: Y-dimension parameters (2 integers)
 *         - info[6-7]: Z-dimension parameters (2 integers)
 */
Region* add_region(std::vector<std::string> info) {
    std::cout << "adding region \n";
    int id = std::stoi(tokenizer(info[0], ":")[0]); // region id number
    std::vector< std::vector<int> > params = vect_create_2D(2,3);
    std::string reg_type = info[1]; // region type

    if (info[1] == "GB") {
        // case of grain boundary region
        params[0][0] = std::stoi(tokenizer(info[2], ":")[1]);
        params[0][1] = std::stoi(tokenizer(info[3], ":")[1]);
        params[0][2] = std::stoi(tokenizer(info[4], ":")[1]);
        
        params[1][0] = std::stoi(tokenizer(info[5], ":")[1]);
        params[1][1] = std::stoi(tokenizer(info[6], ":")[1]);
        params[1][2] = std::stoi(tokenizer(info[7], ":")[1]);
    }
        
    if (info[1] == "BLOCK") {
        // case of region defined as rectangular prism (block)
        params[0][0] = std::stoi(tokenizer(info[2], ":")[1]);
        params[1][0] = std::stoi(tokenizer(info[3], ":")[1]);
        
        params[0][1] = std::stoi(tokenizer(info[4], ":")[1]);
        params[1][1] = std::stoi(tokenizer(info[5], ":")[1]);

        params[0][2] = std::stoi(tokenizer(info[6], ":")[1]);
        params[1][2] = std::stoi(tokenizer(info[7], ":")[1]);
    }
        
    Region* new_region = new Region(id, reg_type, params);
    return new_region;
}

/**
 * @brief Populates the FourDArr `sites` with values corresponding to custom regions from an input file.
 *
 * This function reads a custom regions file and updates the specified 
 * `FourDArr` structure with region IDs based on the coordinates specified 
 * in the file. The function processes each region from the provided list 
 * and updates the 4D array for each lattice position.
 *
 * @param sites A pointer to a `FourDArr` object to be populated with region IDs.
 * @param regions A vector of pointers to `Region` objects representing the regions to be processed.
 * @param custom_reg_idx An index specifying which custom region to draw; currently not used in implementation.
 * @param dim A vector of integers representing the dimensions of the simulation cell.
 * @param infile_name The name of the input file containing custom region data.
 *
 * @note The input file is expected to contain lines of data formatted for 
 *       each region's lattice position, which are parsed and used to update 
 *       the `FourDArr`. The function checks that the coordinates do not exceed 
 *       the bounds of the simulation cell dimensions.
 *
 * @warning If coordinates exceed simulation cell bounds, no updates will be made 
 *          for those entries, but no exception will be thrown. 
 */
void custom_draw_regions(FourDArr* sites, std::vector<Region*> regions, int custom_reg_idx, std::vector<int> dim, std::string infile_name) {
    std::cout << "drawing regions \n";
    Region* region;

    std::fstream in_file;
    std::vector<std::string> lines;
    std::string line;
    std::string output;
    int read_idx = 0;
    std::string lattice_pos; int x; int y; int z;

    for (int i=0; i<(int)regions.size(); i++) {  
        region = regions[i];
        std::cout << "opening custom region file \n";
        in_file.open(infile_name);

        if (in_file.is_open()) {
            std::cout << "custom region file open\n";
            while ( getline (in_file,line) )
            {
                
                reg_line_struct tuple_out = parse_reg_line(line);
                lattice_pos = tuple_out.get_latice_type(); 
                x = tuple_out.get_x(); 
                y = tuple_out.get_y();  
                z = tuple_out.get_z(); 
                
                if ( (x >= dim[0]) || (y >= dim[1]) || (z >= dim[2]) ) {
                    //printf("ERROR: region site exceed simulation cell bounds");
                    //throw std::exception();
                }
                else {
                    //std::cout << region->id << "\n";
                    //std::cout << x << y << z << "\n";
                    (*sites)(0,x,y,z) = region->id;
                    (*sites)(1,x,y,z) = region->id;
                }
            }
            in_file.close();
        }
    }
}

/**
 * @brief Populates the `FourDArr` `sites` with values corresponding to defined regions.
 *
 * This function iterates through a vector of `Region` objects and updates the 
 * provided `FourDArr` structure based on the type of each region (either 
 * "GB" for grain boundaries or "BLOCK" for rectangular prisms). For grain 
 * boundary regions, the function computes the appropriate coordinates based on 
 * the slopes and shifts defined for the region. For block regions, the function 
 * uses the lower and upper bounds to determine the coordinates to populate.
 *
 * @param sites A pointer to a `FourDArr` object where the region IDs will be assigned.
 * @param regions A vector of pointers to `Region` objects that define the regions to draw.
 * @param dim A vector of integers representing the dimensions of the simulation cell.
 *
 * @note The function computes the coordinates for each region based on its 
 *       specific parameters and assigns the corresponding region ID to those 
 *       coordinates in the `FourDArr`.
 *
 * @warning Ensure that the coordinates calculated do not exceed the bounds of 
 *          the `FourDArr`, as this implementation does not include checks for 
 *          out-of-bounds access.
 */
void draw_regions(FourDArr* sites, std::vector<Region*> regions, std::vector<int> dim ) {
    Region* region;
    std::vector<int> lo(3);
    std::vector<int> hi(3);
    int x_ceil;
    int x_floor;
    std::vector<int> start; 
    std::vector<int> end;
    std::vector< std::vector<int> > coords;
    std::vector<int> values;
    std::cout << "drawing regions \n";

    // looping over region objects
    for (int i=0; i<(int)regions.size(); i++) {
        
        region = regions[i];

        if (region->type == "GB") {
            // case of grain boundary region
            for (int y=0; y<(int)(dim[1]); y++) {
                //std::cout << "y: " << y << "\n";
                x_ceil = (int)(ceil((double)(region->slopes[0]/region->slopes[1]) * ( y - (double)(region->shifts[1]) ) + (double)region->shifts[0] )) % dim[0];
                x_floor = (int)(floor((double)(region->slopes[0]/region->slopes[1]) * ( y - (double)(region->shifts[1]) ) + (double)region->shifts[0] )) % dim[0];

                start = {0, x_ceil, y, 0}; end = {0, x_ceil, y, dim[1]};
                coords = FourD_idxs(start, end);
                sites->assign_idxs(coords, create_vec_1D((int)coords.size(), region->id)); 

                start = {0, x_floor, y, 0}; end = {0, x_floor, y, dim[1]};
                coords = FourD_idxs(start, end);
                sites->assign_idxs(coords, create_vec_1D((int)coords.size(), region->id));

                start = {1, x_ceil, y, 0}; end = {1, x_ceil, y, dim[1]};
                coords = FourD_idxs(start, end);
                sites->assign_idxs(coords, create_vec_1D((int)coords.size(), region->id));

                start = {1, x_floor, y, 0}; end = {1, x_floor, y, dim[1]};
                coords = FourD_idxs(start, end);
                sites->assign_idxs(coords, create_vec_1D((int)coords.size(), region->id));
            } 
        }

        else if (region->type == "BLOCK") {
            // case of region defined as rectangular prism (block)
            lo[0] = floor(region->lowerbound[0]);
            lo[1] = floor(region->lowerbound[1]);
            lo[2] = floor(region->lowerbound[2]);

            hi[0] = ceil(region->upperbound[0]);
            hi[1] = ceil(region->upperbound[1]);
            hi[2] = ceil(region->upperbound[2]);

            start = {0, lo[0], lo[1], lo[2]}; end = {0, hi[0], hi[1], hi[2]};
            coords = FourD_idxs(start, end);
            values = create_vec_1D((int)coords.size(),region->id);
            sites->assign_idxs(coords, values);
            start = {1, lo[0], lo[1], lo[2]}; end = {1, hi[0], hi[1], hi[2]};
            coords = FourD_idxs(start, end);
            sites->assign_idxs(coords, values);
        }
    }
}

/**
 * @brief Initializes region objects from input data and populates a `FourDArr` with region sites.
 *
 * This wrapper function reads lines from an input file to create `Region` objects
 * and initializes their corresponding site entries in a `FourDArr`. The function
 * handles both regions of predefined shape (block or slab) and custom regions based 
 * defined by a coordinate file depending whether a custom region input file is provided. 
 * If the number of regions specified exceeds those defined in the input, additional 
 * block regions are generated with default parameters.
 *
 * @param lines A vector of strings containing lines read from the input file.
 * @param dims A vector of integers representing the dimensions of the simulation space.
 * @param num_regions The total number of regions to initialize.
 * @param region_infile The name of the input file for custom regions (can be empty).
 *
 * @return add_reg_struct A structure containing information about the initialized regions,
 *         the number of lines read, and the populated `FourDArr` which indicates which sites
 *         are specially-defined regions.
 *
 * @note The function initializes Region objects based upon the 'num_regions' information 
 *       and region definitions in the regions input file. It also ensures that any 
 *       additional regions needed are created with default block parameters.
 *
 * @warning Make sure that the input lines are properly formatted to avoid runtime errors 
 *          when parsing region definitions. Also, ensure that the `region_infile` 
 *          provided is valid if custom regions are required.
 */
add_reg_struct
init_regions(std::vector<std::string> lines, 
std::vector<int> dims, int num_regions, std::string region_infile) {

    int read_idx = 0;
    std::vector<Region*> regions;
    std::vector<std::string> region_info;
    std::string curr_line = lines[read_idx];
    std::cout << "curr_line: " << lines[read_idx] << "\n";
    std::cout << "num_regions: " << num_regions << "\n";
    FourDArr* temp_region_sites = new FourDArr(2, dims[0], dims[1], dims[2]);

    // initiliazing all entries as 0s
    for (int i=0; i<2; i++) {
        for (int j=0; j<dims[0]; j++) {
            for (int k=0; k<dims[1]; k++) {
                for (int l=0; l<dims[2]; l++) {
                    (*temp_region_sites)(i,j,k,l) = 0;
                }
            }
        }
    }

    // reading input file and initializing regions
    while (curr_line.find("regions end") == std::string::npos) {
        region_info = tokenizer(curr_line, " ");
        Region* region = add_region(region_info);
        regions.push_back(region);
        read_idx ++;
        curr_line = lines[read_idx];
    }

    std::string arbitrary_reigon;
    if (num_regions > (int)regions.size()) {
        std::cout << "adding regions\n";

        for (int i=(int)regions.size(); i <= num_regions; i++) {
            
            arbitrary_reigon += std::to_string(i + (int)regions.size() + 1);
            arbitrary_reigon += ": BLOCK xmin:0 xmax:";
            arbitrary_reigon += std::to_string(dims[0]);
            arbitrary_reigon += " ymin:0 ymax:";
            arbitrary_reigon += std::to_string(dims[1]);
            arbitrary_reigon += " zmin:0 zmax:";
            arbitrary_reigon += std::to_string(dims[2]);
            arbitrary_reigon += "\n";
            region_info = {arbitrary_reigon};

            Region* region = add_region(region_info);
            regions.push_back(region);
        }
    } 

    // entering region id in sites correspoding to pre-defined regions
    std::cout << "regions.size(): " << regions.size() << "\n";
    std::cout << "dims: [" << dims[0] << " " << dims[1] << " " << dims[2] << "]\n";
    std::cout << "region_infile: " << region_infile << "\n";

    if (region_infile.empty()) {
        std::cout << "draw_region: \n";
        draw_regions(temp_region_sites, regions, dims);
    }
    else {
        std::cout << "custom_draw_region: \n";
        custom_draw_regions(temp_region_sites, regions, num_regions, dims, region_infile);
    }

    add_reg_struct returnval(read_idx, regions, temp_region_sites);

    return returnval;
}

/**
 * @brief Wrapper function for populating the lattice.
 *
 * This function reads input files to populate FourDArr and FourDBoolArr 
 * data structures for vacancy and atomic positions corresponding to the 
 * initial configuration of the simulation. It also defines regions within
 * the simulation cell and builds a rate catalog for bulk and pre-defined 
 * regions.
 *
 * @param infile_name The name of the input file to read.
 * @param catalogfile_name The name of the rate catalog file.
 * @param region_infile The name of the region input file.
 * @param vertex_rate The rate for vertex interactions.
 * @param edge_rate The rate for edge interactions.
 * @param reg_rates A vector of vectors containing region rates.
 * @param total_dims A vector containing the total dimensions of the lattice.
 * @param chunk_bounds A vector of vectors defining the bounds of the chunks.
 * @param rank The rank of the current process in a parallel computation.
 * @param procs A vector containing the number of processes in each dimension.
 *
 * @return A pointer to the populated Lattice object.
 */
Lattice* populate_lattice(std::string infile_name, std::string catalogfile_name, std::string region_infile, double vertex_rate, double edge_rate, 
std::vector<std::vector<double>> reg_rates, std::vector<int> total_dims, std::vector<std::vector<int>> chunk_bounds, int rank, std::vector<int> procs) {

    std::fstream in_file;
    in_file.open(infile_name);


    std::vector<std::string> lines;
    std::string line;
    std::string output;
    int read_idx = 0;

    int nprocs = (int)(procs[0] * procs[1]);
    int min_x = chunk_bounds[0][0];
    int max_x = chunk_bounds[0][1];
    int min_y = chunk_bounds[1][0];
    int max_y = chunk_bounds[1][1];
    int min_z = chunk_bounds[2][0];
    int max_z = chunk_bounds[2][1];

    std::cout << "nprocs: " << nprocs << "\n";
    std::cout << "chunk_bounds: \n";
    print_2Dvector(chunk_bounds);

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
    std::vector<int> dims_int(3);

    for (int i=0; i<(int)dims_int.size(); i++) {
        dims_int[i] = std::stoi(dims_str[i+1]); 
    }

    //parsing through tokens except line label (this is why i=1)    

    // reading in geo type ###
    std::vector<bool> periodic = {true, true, true};
    
    // reading in types of atoms mapped onto lattice ###
    std::string atypes = lines[read_idx];
    read_idx ++;
    std::vector<std::string> atypes_str = tokenizer(atypes, " ");
    std::string a;            
    
    //parsing through tokens except line label (this is why i=1)
    std::vector<std::string> a_type_values;
    std::vector<int> a_type_keys;

    for (int i=0; i<(int)(atypes_str.size()-1); i++) {
        a = atypes_str[i+1];
        std::vector< std::string > tokens = tokenizer(a,":");
        int key = std::stoi(tokens[0]); std::string value = tokens[1];
        //std::cout << "key " << key <<"\n";
        //std::cout << "value " << value <<"\n";
        a_type_values.push_back(value);
        a_type_keys.push_back(key);
    } 

    // gettiing num of regions
    int num_regions = 0;
    std::vector<std::string> num_regions_info;

    if (lines[read_idx].find("num_regions") == std::string::npos) {
        std::cout << "ERROR: no number of regions specified" << "\n";
        exit(0);
    }
    else { 
        num_regions_info = tokenizer(lines[read_idx], " ");
        num_regions = std::stoi(num_regions_info[1]); 
    }

    read_idx ++;

    // reading in information about regions // 
    // reading in region dimensions and constants //
    std::string substring = "regions begin";
    std::vector<Region*> temp_regions;
    FourDArr* temp_region_sites;
    int incriment;

    std::cout << "start initializing regions\n";

    if (lines[read_idx].find(substring) != std::string::npos) {
        read_idx ++;
        add_reg_struct regions_tuple = init_regions(slice_1Dvec_str(lines, read_idx, (int)lines.size()), dims_int, num_regions, region_infile); 
        std::cout << "region!\n";
        incriment = regions_tuple.get_idx(); temp_regions = regions_tuple.get_regions(); temp_region_sites = regions_tuple.get_region_sites();
    }
    else {
        printf("ERROR: regions section mis-formatted in geometry file (check for extra newlines)");
        throw std::exception();
    }

    read_idx = read_idx + incriment;
    read_idx ++;

    // reading in atoms, along with their type and coordinate //
    std::tuple<std::string, double, double, double, int> tuple_out;
    std::string lattice_pos;
    double x_raw; double y_raw; double z_raw; int x_unmod; int y_unmod; int z_unmod;
    int x; int y; int z; int x_idx; int y_idx; int z_idx;

    int atomtype;
    int vacancies_count = 0;

    dims_int[0] = chunk_bounds[0][1] - chunk_bounds[0][0];
    dims_int[1] = chunk_bounds[1][1] - chunk_bounds[1][0];
    dims_int[2] = chunk_bounds[2][1] - chunk_bounds[2][0];

    FourDBoolArr* temp_vacancies = new FourDBoolArr(2, (size_t)dims_int[0], (size_t)dims_int[1], (size_t)dims_int[2]);
    FourDBoolArr* temp_vertex_sites = new FourDBoolArr(1, (size_t)dims_int[0], (size_t)dims_int[1], (size_t)dims_int[2]);
    FourDBoolArr* temp_bc_sites = new FourDBoolArr(1, (size_t)dims_int[0], (size_t)dims_int[1], (size_t)dims_int[2]);

    FourDBoolArr* temp_proc_neg_x_neighbors = new FourDBoolArr(1, 2, (size_t)(dims_int[1]+1), (size_t)(dims_int[2]));
    FourDBoolArr* temp_proc_neg_y_neighbors = new FourDBoolArr(1, 2, (size_t)(dims_int[0]+1), (size_t)(dims_int[2]));
    FourDBoolArr* temp_proc_pos_x_neighbors = new FourDBoolArr(1, 2, (size_t)(dims_int[1]+1), (size_t)(dims_int[2]));
    FourDBoolArr* temp_proc_pos_y_neighbors = new FourDBoolArr(1, 2, (size_t)(dims_int[0]+1), (size_t)(dims_int[2]));

    std::vector<size_t> vacs_size_tuple = (*temp_vacancies).size_vec;


    for (size_t i=0; i<2; i++) {
        for (size_t j=0; j<(size_t)dims_int[0]; j++) {
            for (size_t k=0; k<(size_t)dims_int[1]; k++) {
                for (size_t l=0; l<(size_t)dims_int[2]; l++) {
                    if (i == 0) {
                        (*temp_vertex_sites)(0,j,k,l) = 1;
                    }
                    else if (i == 1) {
                        (*temp_bc_sites)(0,j,k,l) = 1;
                    }
                    (*temp_vacancies)(i,j,k,l) = 0;
                }
            }
        }
    }

    int xhi_edge; int xlo_edge;
    int yhi_edge; int ylo_edge; 
    
    int xbound_lo; int xbound_hi;
    int ybound_lo; int ybound_hi;

    xlo_edge = (((chunk_bounds[0][0]-1) % total_dims[0] + total_dims[0]) % total_dims[0]);
    xhi_edge = (((chunk_bounds[0][1]) % total_dims[0] + total_dims[0]) % total_dims[0]);
    ylo_edge = (((chunk_bounds[1][0]-1) % total_dims[1] + total_dims[1]) % total_dims[1]);
    yhi_edge = (((chunk_bounds[1][1]) % total_dims[1] + total_dims[1]) % total_dims[1]);


    std::vector<size_t> vertex_size_tup = (*temp_vertex_sites).size_vec;
    std::vector<size_t> bc_size_tup = (*temp_bc_sites).size_vec;
    std::vector<size_t> vacancies_size_tup = (*temp_vacancies).size_vec;

    
    std::vector< std::vector<size_t> > coords; 
    std::vector< std::vector<size_t> > out_coords;


    std::vector<size_t> x_dims = (*temp_proc_neg_x_neighbors).size_vec;
    std::vector<size_t> y_dims = (*temp_proc_neg_y_neighbors).size_vec;
    std::vector<size_t> x_dims_pos = (*temp_proc_pos_x_neighbors).size_vec;
    std::vector<size_t> y_dims_pos = (*temp_proc_pos_y_neighbors).size_vec;
    std::cout << "x_dims: ";
    print_1Dvector(x_dims);
    std::cout << "y_dims: ";
    print_1Dvector(y_dims);
    std::cout << "x_dims_pos: ";
    print_1Dvector(x_dims_pos);
    std::cout << "y_dims_pos: ";
    print_1Dvector(y_dims_pos);


    /* making procs neighbors */
    int total_procs = procs[0] * procs[1];
    x_idx = 0; y_idx = 0;
    int curr_x=0; int curr_y=0;
    Matrix<int> all_procs((size_t)procs[0], (size_t)procs[1]);
    Matrix<int> temp_proc_neighbors((size_t)total_procs, (size_t)8);
    std::vector<size_t> new_idxs;

    for (int rank_i=0; rank_i<total_procs; rank_i++) {
        x_idx = rank_i % procs[0];
        y_idx = floor(rank_i / procs[0]);

        all_procs[x_idx][y_idx] = rank_i;
        //if (rank == rank_i) { curr_x = x_idx; curr_y = y_idx; }
    }
    
    std::cout << "all_procs: \n";
    all_procs.print();

    for (int rank_i=0; rank_i<total_procs; rank_i++) {
        curr_x = rank_i % procs[0];
        curr_y = floor(rank_i / procs[0]);
        //std::cout << "curr_x: " << curr_x << " curr_y: " << curr_y << "\n";

        new_idxs = mod_with_bounds(curr_x + 1, curr_y, procs[0], procs[1]);
        temp_proc_neighbors[rank_i][0] = all_procs[new_idxs[0]][new_idxs[1]];
        new_idxs = mod_with_bounds(curr_x + 1, curr_y + 1, procs[0], procs[1]);
        temp_proc_neighbors[rank_i][1] = all_procs[new_idxs[0]][new_idxs[1]];
        new_idxs = mod_with_bounds(curr_x, curr_y + 1, procs[0], procs[1]);
        temp_proc_neighbors[rank_i][2] = all_procs[new_idxs[0]][new_idxs[1]];
        new_idxs = mod_with_bounds(curr_x - 1, curr_y + 1, procs[0], procs[1]);
        temp_proc_neighbors[rank_i][3] = all_procs[new_idxs[0]][new_idxs[1]];
        new_idxs = mod_with_bounds(curr_x - 1, curr_y, procs[0], procs[1]);
        temp_proc_neighbors[rank_i][4] = all_procs[new_idxs[0]][new_idxs[1]];
        new_idxs = mod_with_bounds(curr_x - 1, curr_y - 1, procs[0], procs[1]);
        temp_proc_neighbors[rank_i][5] = all_procs[new_idxs[0]][new_idxs[1]];
        new_idxs = mod_with_bounds(curr_x, curr_y - 1, procs[0], procs[1]);
        temp_proc_neighbors[rank_i][6] = all_procs[new_idxs[0]][new_idxs[1]];
        new_idxs = mod_with_bounds(curr_x + 1, curr_y - 1, procs[0], procs[1]);
        temp_proc_neighbors[rank_i][7] = all_procs[new_idxs[0]][new_idxs[1]];   
    }



    for (int i=read_idx; i<(int)lines.size(); i++) {
        std::cout << "rank: " << rank << " lines[i]: " << lines[i] << "\n";
        line_struct tuple_out = parse_line(lines[i]);
        // adjust the x_idx and y_idx to be for negative arrays vs positive arrays (+1 vs +0 on idx) rather than bc vs v
        lattice_pos = tuple_out.get_latice_pos(); x_raw = tuple_out.get_x(); y_raw = tuple_out.get_y(); 
        z_raw = tuple_out.get_z(); atomtype = tuple_out.get_atype();
        atomtype = (int)(atomtype);
        // PERHAPS SEPARATE IF STATEMENTS HERE AS WELL
        x_unmod = floor(x_raw);
        y_unmod = floor(y_raw);
        z_unmod = floor(z_raw);

        if (((procs[0] != 1) || (procs[1] != 1)) && (atomtype == 0)) {
            

            if (lattice_pos == "v") {
                x = mod_with_bounds(x_raw, dims_int[0]);
                y = mod_with_bounds(y_raw, dims_int[1]);
                z = mod_with_bounds(z_raw, dims_int[2]);
                z_idx = z;
            }
            else {
                x = mod_with_bounds((x_raw - 0.5), dims_int[0]);
                y = mod_with_bounds((y_raw - 0.5), dims_int[1]);
                z = mod_with_bounds((z_raw - 0.5), dims_int[2]);
                //CHANGE WHEN TO USE + 1 IN IDX (ONLY FOR NEG NEIGHBOR ARRAYS)
                z_idx = z;
            }

            //std::cout << "lattice_pos: " << lattice_pos << " x: " << x << " y: " << y << " z: " << z << "\n";
            //std::cout << "xlo_edge: " << xlo_edge << " xhi_edge: " << xhi_edge << " ylo_edge: " << ylo_edge << " yhi_edge: " << yhi_edge << "\n";

            //ADD IF STATEMENT TO CHECK IF CHUNK_BOUNDS[i] == TOTAL_DIMS[i]
            
            /* first layer ghost sites for (111) and (100) moves */
            if ((y_unmod >= (chunk_bounds[1][0] - 1)) && (y_unmod < (chunk_bounds[1][1]) )) {
                if ((x_unmod == xlo_edge) && (lattice_pos == "bc")) {
                    
                    std::cout << "neg x neigh 2 rank: " << rank  << " y_idx: " << y_idx << " z_idx: " << z_idx << "\n";
                    std::cout << "neg x neigh 2 rank: " << rank << " x_unmod: " << x_unmod << " y_unmod: " << y_unmod << " z_unmod: " << z_unmod << "\n";
                    y_idx = mod_with_bounds(((y_raw - 0.5) - chunk_bounds[1][0] + 1), total_dims[1]);
                    (*temp_proc_neg_x_neighbors)(0,1,y_idx,z_idx) = 1;
                    
                    //filling in other corner to preserve corner periodic boundary conditions
                    if ((y_idx == x_dims[2]-1) && (temp_proc_neighbors(rank,2) == rank)) { 
                        std::cout << "neg x neigh 2 rank corner: " << rank << " x_dims[2]-1: " << (x_dims[2]-1) << "\n";
                        (*temp_proc_neg_x_neighbors)(0,1,0,z_idx) = 1; }

                    if ((y_unmod == ylo_edge)) {
                        std::cout << "neg y neigh 2 rank: " << rank << " x_idx: " << x_idx << " z_idx: " << z_idx << "\n";
                        std::cout << "neg y neigh 2 rank: " << rank << " x_unmod: " << x_unmod << " y_unmod: " << y_unmod << " z_unmod: " << z_unmod << "\n";
                        x_idx = mod_with_bounds(((x_raw - 0.5) - chunk_bounds[0][0] + 1), total_dims[0]);
                        (*temp_proc_neg_y_neighbors)(0,1,x_idx,z_idx) = 1;
                    }   
                }
            }
            if ((x_unmod >= (chunk_bounds[0][0] - 1)) && (x_unmod < (chunk_bounds[0][1]) )) {
                if ((y_unmod == ylo_edge)) { 
                    if (lattice_pos == "bc") {
                        std::cout << "neg y neigh 2 rank: " << rank << " x_idx: " << x_idx << " z_idx: " << z_idx << "\n";
                        std::cout << "neg y neigh 2 rank: " << rank << " x_unmod: " << x_unmod << " y_unmod: " << y_unmod << " z_unmod: " << z_unmod << "\n";
                        x_idx = mod_with_bounds(((x_raw - 0.5) - chunk_bounds[0][0] + 1), total_dims[0]);
                        (*temp_proc_neg_y_neighbors)(0,1,x_idx,z_idx) = 1;

                        //filling in other corner to preserve corner periodic boundary conditions
                        if ((x_idx == y_dims[2]-1) && (temp_proc_neighbors(rank,0) == rank)) { 
                            std::cout << "pos x neigh 2 rank corner: " << rank << " y_dims[2]-1: " << (y_dims[2]-1) << "\n";
                            (*temp_proc_neg_y_neighbors)(0,1,0,z_idx) = 1; }
                    }
                }
            }

            /* second layer ghost sites for (100) moves */
            if ((y_unmod >= (chunk_bounds[1][0])) && (y_unmod < (chunk_bounds[1][1] + 1) )) {
                if ((x_unmod == xhi_edge) && (lattice_pos == "bc")) {
                    std::cout << "pos x neigh 100 rank: " << rank << " y_idx: " << y_idx << " z_idx: " << z_idx << "\n";
                    std::cout << "pos x neigh 100 rank: " << rank << " x_unmod: " << x_unmod << " y_unmod: " << y_unmod << " z_unmod: " << z_unmod << "\n";
                    y_idx = mod_with_bounds((y_raw - chunk_bounds[1][0]), total_dims[1]);
                    (*temp_proc_pos_x_neighbors)(0,1,y_idx,z_idx) = 1;

                    //filling in other corner to preserve corner periodic boundary conditions
                    if ((y_idx == x_dims_pos[2]-1) && (temp_proc_neighbors(rank,2) == rank)) { 
                        std::cout << "pos x neigh 100 rank corner: " << rank << " x_dims[2]-1: " << (x_dims[2]-1) << "\n";
                        (*temp_proc_pos_x_neighbors)(0,1,0,z_idx) = 1; }
                }
            }
            if ((x_unmod >= (chunk_bounds[0][0])) && (x_unmod < (chunk_bounds[0][1] + 1) )) {
                if ((y_unmod == yhi_edge) && (lattice_pos == "bc")) {
                    
                    std::cout << "pos y neigh 100 rank: " << rank << " x_idx: " << x_idx << " z_idx: " << z_idx << "\n";
                    std::cout << "pos y neigh 100 rank: " << rank << " x_unmod: " << x_unmod << " y_unmod: " << y_unmod << " z_unmod: " << z_unmod << "\n";
                    x_idx = mod_with_bounds((x_raw - chunk_bounds[0][0]), total_dims[0]);
                    (*temp_proc_pos_y_neighbors)(0,1,x_idx,z_idx) = 1;

                    //filling in other corner to preserve corner periodic boundary conditions
                    if ((x_idx == 0) && (temp_proc_neighbors(rank,0) == rank)) { 
                        std::cout << "pos y neigh 100 rank corner: " << rank << " y_dims[2]-1: " << (y_dims[2]-1) << "\n";
                        (*temp_proc_pos_y_neighbors)(0,1,y_dims_pos[2]-1,z_idx) = 1; }
                }
            }

            /* first layer ghost sites for (111) and (100) moves */
            if ((y_unmod >= (chunk_bounds[1][0])) && (y_unmod < (chunk_bounds[1][1] + 1) )) {
                if ((x_unmod == xhi_edge) && (lattice_pos == "v")) {
                    
                    y_idx = mod_with_bounds((y_raw - chunk_bounds[1][0]), total_dims[1]);
                    //std::cout << "pos x neigh 2 rank: " << rank << " x: " << x << " y: " << y << " z: " << z << "\n";
                    std::cout << "pos x neigh 1 rank: " << rank << " y_idx: " << y_idx << " z_idx: " << z_idx << "\n";
                    std::cout << "pos x neigh 1 rank: " << rank << " x_unmod: " << x_unmod << " y_unmod: " << y_unmod << " z_unmod: " << z_unmod << "\n";
                    (*temp_proc_pos_x_neighbors)(0,0,y_idx,z_idx) = 1;

                    //filling in other corner to preserve corner periodic boundary conditions
                    if ((y_idx == 0) && (temp_proc_neighbors(rank,2) == rank)) { 
                        std::cout << "pos x neigh 1 rank corner: " << rank << " x_dims[2]-1: " << (x_dims[2]-1) << "\n";
                        (*temp_proc_pos_x_neighbors)(0,0,(size_t)(x_dims_pos[2]-1),z_idx) = 1; 
                        }
                    
                    if ((y_unmod == yhi_edge)) {
                        std::cout << "pos y neigh 1 rank: " << rank <<  " x_idx: " << x_idx << " z_idx: " << z_idx << "\n";
                        std::cout << "pos y neigh 1 rank: " << rank << " x_unmod: " << x_unmod << " y_unmod: " << y_unmod << " z_unmod: " << z_unmod << "\n";
                        x_idx = mod_with_bounds((x_raw - chunk_bounds[0][0]), total_dims[0]);
                        (*temp_proc_pos_y_neighbors)(0,0,x_idx,z_idx) = 1;
                    }
                }
            }
            if ((x_unmod >= (chunk_bounds[0][0])) && (x_unmod < (chunk_bounds[0][1] + 1) )) {
                if ((y_unmod == yhi_edge)) {
                    
                    std::cout << "pos y neigh 1 rank: " << rank <<  " x_idx: " << x_idx << " z_idx: " << z_idx << "\n";
                    std::cout << "pos y neigh 1 rank: " << rank << " x_unmod: " << x_unmod << " y_unmod: " << y_unmod << " z_unmod: " << z_unmod << "\n";
                    x_idx = mod_with_bounds((x_raw - chunk_bounds[0][0]), total_dims[0]);    
                    (*temp_proc_pos_y_neighbors)(0,0,x_idx,z_idx) = 1; 

                    //filling in other corner to preserve corner periodic boundary conditions
                    if ((x_idx == 0) && (temp_proc_neighbors(rank,0) == rank)) { 
                        std::cout << "pos y neigh 1 rank corner: " << rank << " y_dims[2]-1: " << (y_dims[2]-1) << "\n";
                        (*temp_proc_pos_y_neighbors)(0,0,(size_t)(y_dims_pos[2]-1),z_idx) = 1; 
                        }
                }
            }

            /* second layer ghost sites for (100) moves */
            if ((y_unmod >= (chunk_bounds[1][0] - 1)) && (y_unmod < (chunk_bounds[1][1]) )) {
                if ((x_unmod == xlo_edge) && (lattice_pos == "v")) {

                    std::cout << "neg x neigh 100 rank: " << rank <<  " y_idx: " << y_idx << " z_idx: " << z_idx << "\n";
                    std::cout << "neg x neigh 100 rank: " << rank << " x_unmod: " << x_unmod << " y_unmod: " << y_unmod << " z_unmod: " << z_unmod << "\n";
                    y_idx = mod_with_bounds(((y_raw - 0.5) - chunk_bounds[1][0] + 1), total_dims[1]);
                    (*temp_proc_neg_x_neighbors)(0,0,y_idx,z_idx) = 1;

                    //filling in other corner to preserve corner periodic boundary conditions
                    //if (y_idx == x_dims[2]-1) { (*temp_proc_neg_x_neighbors)(0,0,0,z_idx) = 1; }
                }
            }
            if ((x_unmod >= (chunk_bounds[0][0] - 1)) && (x_unmod < (chunk_bounds[0][1]) )) {
                if ((y_unmod == ylo_edge) && (lattice_pos == "v")) {
                    
                    std::cout << "neg y neigh 100 rank: " << rank <<  " x_idx: " << x_idx << " z_idx: " << z_idx << "\n";
                    std::cout << "neg y neigh 100 rank: " << rank << " x_unmod: " << x_unmod << " y_unmod: " << y_unmod << " z_unmod: " << z_unmod << "\n";
                    x_idx = mod_with_bounds(((x_raw - 0.5) - chunk_bounds[0][0] + 1), total_dims[0]);
                    (*temp_proc_neg_y_neighbors)(0,0,x_idx,z_idx) = 1;

                    //filling in other corner to preserve corner periodic boundary conditions
                    //if (x_idx == y_dims[2]-1) { (*temp_proc_pos_y_neighbors)(0,0,0,z_idx) = 1; }
                }
            }
        }
        
        if ((x_unmod >= chunk_bounds[0][0]) && (x_unmod < chunk_bounds[0][1]) && 
        (y_unmod >= chunk_bounds[1][0]) && (y_unmod < chunk_bounds[1][1]) &&
        (z_unmod >= chunk_bounds[2][0]) && (z_unmod < chunk_bounds[2][1]) ) {
            if (lattice_pos == "v") {
                x = mod_with_bounds(x_unmod, dims_int[0]);
                y = mod_with_bounds(y_unmod, dims_int[1]);
                z = mod_with_bounds(z_unmod, dims_int[2]);

                if (atomtype == 0) {
                    (*temp_vertex_sites)(0,(size_t)x,(size_t)y,(size_t)z) = 0;
                    (*temp_vacancies)(0,(size_t)x,(size_t)y,(size_t)z) = 1;
                    //std::cout << "rank: " << rank << " x_unmod: " << x_unmod << " y_unmod: " << y_unmod << " z_unmod: " << z_unmod << "\n";
                    //std::cout << "rank: " << rank << " x: " << x << " y: " << y << " z: " << z << "\n";
                    if ( is_in(coords, {0,(size_t)x,(size_t)y,(size_t)z})) {
                        out_coords.push_back({0,(size_t)x,(size_t)y,(size_t)z});
                    }
                    else { coords.push_back({0,(size_t)x,(size_t)y,(size_t)z}); }
                    vacancies_count ++;
                }
                else if (atomtype == 1) {
                    (*temp_vertex_sites)(0,x,y,z) = 1;
                }
                else {
                    printf("Unrecognized atom type");
                    throw std::exception();
                }
            }
            else {
                x = mod_with_bounds((x_unmod), dims_int[0]);
                y = mod_with_bounds((y_unmod), dims_int[1]);
                z = mod_with_bounds((z_unmod), dims_int[2]);

                if (atomtype == 0) {
                    (*temp_bc_sites)(0,(size_t)x,(size_t)y,(size_t)z) = 0;
                    (*temp_vacancies)(1,(size_t)x,(size_t)y,(size_t)z) = 1;
                    //std::cout << "rank: " << rank << " x_unmod: " << x_unmod << " y_unmod: " << y_unmod << " z_unmod: " << z_unmod << "\n";
                    //std::cout << "rank: " << rank << " x: " << x << " y: " << y << " z: " << z << "\n";
                    if (is_in(coords, {1,(size_t)x,(size_t)y,(size_t)z})) {
                        out_coords.push_back({1,(size_t)x,(size_t)y,(size_t)z});
                    }
                    else { coords.push_back({1,(size_t)x,(size_t)y,(size_t)z}); }
                    vacancies_count ++;
                }
                else if (atomtype == 1) {
                    (*temp_bc_sites)(0,x,y,z) = 1;
                }
                else {
                    printf("Unrecognized atom type");
                    throw std::exception();
                }
            }
        }
    }

    if (catalogfile_name != "None") {
        ratecatalog_struct cat_tuple = updated_create_ratecatalog(catalogfile_name, a_type_keys);
        std::cout << "cat_tuple out \n";
        std::vector< std::vector< std::vector<double> > > reg_energies = cat_tuple.get_energies();
        int region_num = cat_tuple.get_region_num();
    }
    else {
        std::cout << "ERROR: No rate catalog file provided";
        exit(0);
    }

    std::cout << "read_catalog \n";

    //creating & sorting rate catalog
    std::vector< std::vector< std::vector<double> > > new_energies;
    std::vector< std::vector< std::vector< std::vector<double> > > > reg_energies;
    std::vector< std::vector<int> > new_configs;
    std::vector<int> diag_configs;
    std::vector<int> lateral_configs;
    std::vector< std::vector<double> > diag_E;
    std::vector< std::vector<double> > lateral_E;
    std::vector<int> temp_vec;
    std::vector< std::vector<double> > temp_vec_2D;
    std::vector<int> nums_vacs {4,5,6,7,8};
    std::vector<double> vacs_Es {0.02, 0.02, 0.02, 0.02, 0.02};
    int reg_num = 0;

    Matrix<double>* temp_111_catalog = new Matrix<double>((size_t)2, (size_t)exp_int(2, 8));
    Matrix<double>* temp_100_catalog = new Matrix<double>((size_t)2, (size_t)exp_int(2, 14));

    Matrix<int> configs_111((size_t)1, (size_t)exp_int(2, 8));
    Matrix<int> configs_100((size_t)1, (size_t)exp_int(2, 14)); 
    
    for (int i=0; i<9; i++) {
        temp_vec = bin_m_zeros(8, i, (int)a_type_values.size());
        diag_configs.insert(diag_configs.end(), temp_vec.begin(), temp_vec.end() );
        if (i==0) {
            temp_vec_2D = vect_create_2D_float(2, NCR(8,i), vertex_rate);
            diag_E = temp_vec_2D;
        }
        else {
            temp_vec_2D = vect_create_2D_float(2, NCR(8,i), vertex_rate);
            diag_E[0].insert(diag_E[0].end(), temp_vec_2D[0].begin(), temp_vec_2D[0].end());
            diag_E[1].insert(diag_E[1].end(), temp_vec_2D[1].begin(), temp_vec_2D[1].end());
        }
        /*
        else if (i>=3) {
            temp_vec_2D = vect_create_2D_float(2, NCR(8,i), 0.116);
            diag_E[0].insert(diag_E[0].end(), temp_vec_2D[0].begin(), temp_vec_2D[0].end());
            diag_E[1].insert(diag_E[1].end(), temp_vec_2D[1].begin(), temp_vec_2D[1].end());
        }
        */
    }

    new_configs.push_back(diag_configs);
    new_energies.push_back(diag_E);

    for (int i=0; i<15; i++) {
        temp_vec = bin_m_zeros(14, i, (int)a_type_values.size());
        lateral_configs.insert(lateral_configs.end(), temp_vec.begin(), temp_vec.end());
        if (i==0) {
            temp_vec_2D = vect_create_2D_float(2, NCR(14,i), edge_rate);
            lateral_E = temp_vec_2D;
        }
        else {
            temp_vec_2D = vect_create_2D_float(2, NCR(14,i), edge_rate);
            lateral_E[0].insert(lateral_E[0].end(), temp_vec_2D[0].begin(), temp_vec_2D[0].end());
            lateral_E[1].insert(lateral_E[1].end(), temp_vec_2D[1].begin(), temp_vec_2D[1].end());
        }
    }

    new_configs.push_back(lateral_configs);
    new_energies.push_back(lateral_E);

    std::vector<int> unsorted_idxs_8bit = arange(0,(int)new_configs[0].size(),1);
    auto comparator_8bit = [new_configs](int idx1, int idx2) {
                return new_configs[0][idx1] < new_configs[0][idx2];
        };

    std::sort(unsorted_idxs_8bit.begin(), unsorted_idxs_8bit.end(), comparator_8bit);
    new_configs[0] = reorder_inp(new_configs[0], unsorted_idxs_8bit);
    new_energies[0] = reorder_inp(new_energies[0], unsorted_idxs_8bit);

    std::vector<int> unsorted_idxs_14bit = arange(0,(int)new_configs[1].size(),1);
    auto comparator_14bit = [new_configs](int idx1, int idx2) {
                return new_configs[1][idx1] < new_configs[1][idx2];
        };

    std::sort(unsorted_idxs_14bit.begin(), unsorted_idxs_14bit.end(), comparator_14bit);
    new_configs[1] = reorder_inp(new_configs[1], unsorted_idxs_14bit);
    new_energies[1] = reorder_inp(new_energies[1], unsorted_idxs_14bit);

    for (int i=0; i<(int)new_configs[0].size(); i++) {
       configs_111[0][i] = new_configs[0][i]; 
    }
    for (int i=0; i<(int)new_configs[1].size(); i++) {
       configs_100[0][i] = new_configs[1][i]; 
    }
    for (int i=0; i<(int)new_energies[0].size(); i++) {
        for (int j=0; j<(int)new_energies[0][i].size(); j++) {
            (*temp_111_catalog)(i,j) = new_energies[0][i][j];
        }
    }
    for (int i=0; i<(int)new_energies[1].size(); i++) {
        for (int j=0; j<(int)new_energies[1][i].size(); j++) {
            (*temp_100_catalog)(i,j) = new_energies[1][i][j];
        }
    }

    // intialzing lattice, basis vectors, vacancies, mobile ions, and fixed //
    // atoms based upon dimensions //
    
    int num_x_neigh = dims_int[0] + 1;
    int num_y_neigh = dims_int[1] + 1;

    Lattice* new_lattice = new Lattice(dims_int[0], dims_int[1], dims_int[2], vacancies_count, num_regions, nprocs, num_x_neigh, num_y_neigh, total_dims[0], total_dims[1], total_dims[2]);

    for (size_t i=0; i<2; i++) {
        for (size_t j=0; j<(size_t)dims_int[0]; j++) {
            for (size_t k=0; k<(size_t)dims_int[1]; k++) {
                for (size_t l=0; l<(size_t)dims_int[2]; l++) {
                    if (i == 0) {
                        new_lattice->vertex_sites(0,j,k,l) = (*temp_vertex_sites)(0,j,k,l);
                    }
                    else if (i == 1) {
                        new_lattice->bc_sites(0,j,k,l) = (*temp_bc_sites)(0,j,k,l);
                    }

                    new_lattice->region_sites(i,j,k,l) = (*temp_region_sites)(i,j,k,l);
                    new_lattice->vacancies(i,j,k,l) = (*temp_vacancies)(i,j,k,l);
                }
            }
        }
    }

    size_t xlen; size_t ylen; size_t zlen;
    
    x_dims = new_lattice->proc_pos_x_neighbors.size_vec;
    y_dims = new_lattice->proc_pos_y_neighbors.size_vec;
    ylen = x_dims[2];
    zlen = x_dims[3];
    for (size_t i=0; i<2; i++) {
        for (size_t j=0; j<ylen; j++) {
            for (size_t k=0; k<zlen; k++) {
                new_lattice->proc_neg_x_neighbors(0,i,j,k) = (*temp_proc_neg_x_neighbors)(0,i,j,k);
                new_lattice->proc_pos_x_neighbors(0,i,j,k) = (*temp_proc_pos_x_neighbors)(0,i,j,k);
            }
        }
    }


    xlen = y_dims[2];
    zlen = y_dims[3];
    for (size_t i=0; i<2; i++) {
        for (size_t j=0; j<xlen; j++) {
            for (size_t k=0; k<zlen; k++) {
                new_lattice->proc_neg_y_neighbors(0,i,j,k) = (*temp_proc_neg_y_neighbors)(0,i,j,k);
                new_lattice->proc_pos_y_neighbors(0,i,j,k) = (*temp_proc_pos_y_neighbors)(0,i,j,k);
            }
        }
    }
    
    new_lattice->rank = rank;
    new_lattice->configs_111 = configs_111;
    new_lattice->configs_100 = configs_100;

    size_t rows = new_lattice->ratecatalog_111.rows();
    size_t cols = new_lattice->ratecatalog_111.cols();
    for (size_t i=0; i<rows; i++) {
        for (size_t j=0; j<cols; j++) {
            new_lattice->ratecatalog_111(i,j) = (*temp_111_catalog)(i,j);
        }
    }
    
    rows = new_lattice->ratecatalog_100.rows();
    cols = new_lattice->ratecatalog_100.cols();
    for (size_t i=0; i<rows; i++) {
        for (size_t j=0; j<cols; j++) {
            new_lattice->ratecatalog_100(i,j) = (*temp_100_catalog)(i,j);
        }
    }
    
    rows = new_lattice->regionrates_100_L.rows();
    cols = new_lattice->regionrates_100_L.cols();
    for (size_t i=0; i<rows; i++) {
        for (size_t j=0; j<cols; j++) {
            new_lattice->regionrates_100_L(i,j) = 1e-10;
            new_lattice->regionrates_100_R(i,j) = 1e-10 ;
        }
    }

    rows = new_lattice->regionrates_111_L.rows();
    cols = new_lattice->regionrates_111_L.cols();
    for (size_t i=0; i<rows; i++) {
        std::cout << i << "\n";
        std::cout << reg_rates[i][0] << "\n";
        std::cout << reg_rates[i][1] << "\n";
        for (size_t j=0; j<cols; j++) {
            new_lattice->regionrates_111_L(i,j) = reg_rates[i][0];
            new_lattice->regionrates_111_R(i,j) = reg_rates[i][1];
        }
    }

    new_lattice->rates = create_vec_1D_float(vacancies_count);

    for (int i=0; i<(int)a_type_values.size(); i++) {new_lattice->a_types[a_type_keys[i]] = a_type_values[i];}

    new_lattice->region_energies = reg_energies;    
    new_lattice->chunk_bounds = chunk_bounds;
    new_lattice->proc_dims = procs;

    Matrix<int> nonzero_vacs = new_lattice->vacancies.nonzero();

    for (int i=0; i<(int)nonzero_vacs.rows(); i++) {
        for (int j=0; j<(int)nonzero_vacs.cols(); j++) {
            new_lattice->vacancies_pos(i,j) = nonzero_vacs[i][j];
        }
    }

    std::cout << "rank: " << rank << " nonzero_vacs: \n";
    nonzero_vacs.print();
    
    new_lattice->proc_neighbors = temp_proc_neighbors;
   
    std::cout << "proc_neighbors: \n";
    new_lattice->proc_neighbors.print();

    MPI_Barrier(MPI_COMM_WORLD);
    if (rank == 0) {
        std::cout << "rank: " << rank << " new_lattice->proc_neg_x_neighbors: \n";
        new_lattice->proc_neg_x_neighbors.print();

        std::cout << "rank: " << rank << " new_lattice->proc_pos_x_neighbors: \n";
        new_lattice->proc_pos_x_neighbors.print();

        std::cout << "rank: " << rank << " new_lattice->proc_neg_y_neighbors: \n";
        new_lattice->proc_neg_y_neighbors.print();

        std::cout << "rank: " << rank << " new_lattice->proc_pos_y_neighbors: \n";
        new_lattice->proc_pos_y_neighbors.print();
    }

    MPI_Barrier(MPI_COMM_WORLD);
    if (rank == 1) {
        std::cout << "rank: " << rank << " new_lattice->proc_neg_x_neighbors: \n";
        new_lattice->proc_neg_x_neighbors.print();

        std::cout << "rank: " << rank << " new_lattice->proc_pos_x_neighbors: \n";
        new_lattice->proc_pos_x_neighbors.print();

        std::cout << "rank: " << rank << " new_lattice->proc_neg_y_neighbors: \n";
        new_lattice->proc_neg_y_neighbors.print();

        std::cout << "rank: " << rank << " new_lattice->proc_pos_y_neighbors: \n";
        new_lattice->proc_pos_y_neighbors.print();
    }

    MPI_Barrier(MPI_COMM_WORLD);
    if (rank == 2) {
        std::cout << "rank: " << rank << " new_lattice->proc_neg_x_neighbors: \n";
        new_lattice->proc_neg_x_neighbors.print();

        std::cout << "rank: " << rank << " new_lattice->proc_pos_x_neighbors: \n";
        new_lattice->proc_pos_x_neighbors.print();

        std::cout << "rank: " << rank << " new_lattice->proc_neg_y_neighbors: \n";
        new_lattice->proc_neg_y_neighbors.print();

        std::cout << "rank: " << rank << " new_lattice->proc_pos_y_neighbors: \n";
        new_lattice->proc_pos_y_neighbors.print();
    }

    MPI_Barrier(MPI_COMM_WORLD);
    if (rank == 3) {
        std::cout << "rank: " << rank << " new_lattice->proc_neg_x_neighbors: \n";
        new_lattice->proc_neg_x_neighbors.print();

        std::cout << "rank: " << rank << " new_lattice->proc_pos_x_neighbors: \n";
        new_lattice->proc_pos_x_neighbors.print();

        std::cout << "rank: " << rank << " new_lattice->proc_neg_y_neighbors: \n";
        new_lattice->proc_neg_y_neighbors.print();

        std::cout << "rank: " << rank << " new_lattice->proc_pos_y_neighbors: \n";
        new_lattice->proc_pos_y_neighbors.print();
    }
    
    new_lattice->void_threshold = 2;
    new_lattice->void_barrier = 2761667.608;


    delete temp_111_catalog;
    delete temp_100_catalog;
    delete temp_region_sites;

    delete temp_vacancies;
    delete temp_vertex_sites;
    delete temp_bc_sites;

    delete temp_proc_neg_x_neighbors;
    delete temp_proc_neg_y_neighbors;
    delete temp_proc_pos_x_neighbors;
    delete temp_proc_pos_y_neighbors;


    new_lattice->rate_cumsum.resize(14*nonzero_vacs.rows());

    return new_lattice;
}
/*---------------------------------------------------------------------------*/
