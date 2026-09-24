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
#include <set>

#include "hpp_files_allrates/math_func.hpp"
#include "hpp_files_allrates/vec_func.hpp"

// See definition below (after class Lattice) for full documentation.
// value: 1 to add this call's contribution to the target cell (reference count += 1), 0 to
// remove it (reference count -= 1, clamped at 0). See the comment on the Lattice ghost-array
// members for why these are counts, not booleans -- a single (w,x_idx,y_idx,z) cell can be the
// legitimate geometric target of more than one distinct vacancy near a process-grid corner.
void set_ghost_position_impl(int w, int x_unmod, int y_unmod, int z_unmod, int value,
                              const std::vector<int>& total_dims,
                              const std::vector<std::vector<int>>& chunk_bounds,
                              FourDArr& neg_x_arr, FourDArr& pos_x_arr,
                              FourDArr& neg_y_arr, FourDArr& pos_y_arr,
                              int trace_rank = -1, int trace_tick = -1, int trace_remote_proc = -2,
                              const std::string& trace_tag = "");

/*------------------------------------------------------------------------------------*/
 /*! \brief A class for storing a simulation cell of atoms and propogating moves around the
 lattice */

class Lattice {

    std::map<int, double> rate_typedict;

    public:
        std::vector<Region*> regions;
        std::vector<int> sublattice_dim;
        std::vector<int> total_dims;
        std::map<int, std::string> a_types;
        FourDArr vertex_sites;
        FourDBoolArr vacancies;
        FourDArr bc_sites;
        FourDArr region_sites;
        int num_procs;
        double t;
        // These 4 arrays used to be plain booleans ("is there a ghost vacancy here"). Near a
        // process-grid corner (see set_ghost_position_impl), the same (w,x_idx,y_idx,z) cell can
        // be the legitimate geometric target of two DIFFERENT, physically-distinct vacancies at
        // once (e.g. one rank's own stationary vacancy plus another rank's vacancy transiently
        // crossing through/being reversed near the same corner). A boolean can't tell "cleared
        // because my contribution left" from "still occupied by something else's contribution",
        // so one source's clear could wipe out another's still-valid mark. These are now reference
        // counts (see set_ghost_position_impl): incremented per contributing vacancy, decremented
        // (clamped at 0) when one leaves; "ghost present" is count > 0, which every existing
        // boolean-style read site (if (proc_pos_x_neighbors(...))) still gets for free from int's
        // implicit bool conversion.
        FourDArr proc_neg_x_neighbors;
        FourDArr proc_neg_y_neighbors;
        FourDArr proc_pos_x_neighbors;
        FourDArr proc_pos_y_neighbors;
        Matrix<int> proc_neighbors;
        std::vector<int> proc_dims;
        std::vector<double> probs;
        std::vector<double> rates;
        std::vector<double> rate_cumsum;

        Matrix<int> diag_directions;
        Matrix<int> edge_directions;
        Matrix<int> moves_coords; 
        Matrix<int> moves_coords_unmod; 
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
        std::vector< std::vector<int> > prev_newlocs_unmod;
        std::vector< std::vector<int> > prev_oldlocs;
        std::vector<int> prev_lattice;
        std::vector<int> prev_idxs;
        std::vector< std::vector<int> > par_prev_oldlocs;
        std::vector< std::vector<int> > par_prev_newlocs;
        std::vector<int> par_move_ticks;
        std::vector<int> par_prev_idx;
        std::vector<int> prev_move_type;
        std::vector<int> prev_move_type_ticks;
        int ghost_done_tag;
        int par_done_tag;
        int conflict_done_tag;
        int conflict_ghost_done_tag;
        int rank;
        // Diagnostic only: best-effort current move_ticks for GHOST_WRITE source-tracking logs
        // (see set_ghost_position_impl's trace_tag param). Updated at the top of the handful of
        // functions that originate live ghost-array writes; not meant to be exact between updates.
        int ghost_trace_tick = -1;
        int void_threshold;
        double void_barrier;
        double terrace_barrier_100;
        double terrace_barrier_111;
        double void_gb_diss_barrier;
        double bulk_migration_111;
        double bulk_migration_100;
        double void_E;
        double interface_E;
        double temperature;
        double interface_barrier;
        double voidsurface_E_below_bulk;
        double total_cost;
        double system_energy;
        int last_idx_chosen;
        int last_newNN;
        int last_currNN;
        std::vector<int> last_oldloc;
        std::vector<int> last_newloc;
        int solo_vacs;
        int adaptive_gb_id;
            

        Lattice(int xdim, int ydim, int zdim, int num_vacancies, int num_regions, int number_procs, int xtot, int ytot, int ztot, std::vector<Region*> regs_in, int rank_in):
            regions(regs_in),
            vertex_sites((size_t)1, (size_t)xdim, (size_t)ydim, (size_t)zdim),
            vacancies((size_t)2, (size_t)xdim, (size_t)ydim, (size_t)zdim),
            bc_sites((size_t)1, (size_t)xdim, (size_t)ydim, (size_t)zdim),
            region_sites((size_t)2, (size_t)xdim, (size_t)ydim, (size_t)zdim),
            num_procs((size_t)number_procs),

            proc_neg_y_neighbors(2, (size_t)(xdim+4), (size_t)(2), (size_t)(zdim)), 
            proc_neg_x_neighbors(2, (size_t)(2), (size_t)(ydim+4), (size_t)(zdim)), 
            proc_pos_y_neighbors(2, (size_t)(xdim+4), (size_t)(2), (size_t)(zdim)), 
            proc_pos_x_neighbors(2, (size_t)(2), (size_t)(ydim+4), (size_t)(zdim)), 
            diag_directions((size_t)8, (size_t)3),
            edge_directions((size_t)8, (size_t)3),

            proc_neighbors((size_t)number_procs, (size_t)8),
            moves_coords((size_t)(14 * num_vacancies), 4),
            moves_coords_unmod((size_t)(14 * num_vacancies), 4),
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
            //mt_obj((unsigned int)(std::chrono::high_resolution_clock::now().time_since_epoch().count() + rank_in)),
            mt_obj((unsigned int)725863834569),
            x_rand(0, (size_t)xdim),
            y_rand(0, (size_t)ydim)

            {
                // FourDArr's backing store is malloc'd, not zero-initialized -- unlike the old
                // FourDBoolArr these replaced, whose std::vector-backed storage was. Must zero
                // explicitly so initial ghost reference counts start at 0, not garbage.
                proc_neg_x_neighbors.zero();
                proc_neg_y_neighbors.zero();
                proc_pos_x_neighbors.zero();
                proc_pos_y_neighbors.zero();
                // Diagnostic only: directly verifies .zero() actually leaves every cell at 0
                // immediately after construction, before populate_lattice's seed loop runs --
                // ruling in/out uninitialized malloc'd memory as the source of unexplained nonzero
                // ghost cells that no traced write (seed or live) ever accounts for.
                std::cout << "rank: " << rank_in << " POST_ZERO_CHECK neg_x nonzero: " << proc_neg_x_neighbors.nonzero_elems().rows()
                          << " neg_y nonzero: " << proc_neg_y_neighbors.nonzero_elems().rows()
                          << " pos_x nonzero: " << proc_pos_x_neighbors.nonzero_elems().rows()
                          << " pos_y nonzero: " << proc_pos_y_neighbors.nonzero_elems().rows() << "\n";
                sublattice_dim = {xdim,ydim,zdim};
                total_dims = {xtot,ytot,ztot};
                t = 0;
                num_of_moves = 0;
                num_of_vacs = num_vacancies;
                ghost_done_tag = 10;
                par_done_tag = 11;
                conflict_done_tag = 12;
                conflict_ghost_done_tag = 13;
                rank = rank_in;
                system_energy = 0;
            }
        /**
        * @brief Wrapper function to assign rates to region-specific rate catalogs.
        * 
        * @param temp_regions A vector of pointers to Region objects representing different regions.
        * @param misc_rates A constant reference to a vector of double values representing miscellaneous rates.
        */
        void assign_region_rates_wrapper(std::vector<Region*>& temp_regions, std::vector<double>& misc_rates) {
            std::cout << "temp_regions.size(): " <<temp_regions.size() << "\n";
            std::cout << "misc_rates:  ";
            print_1Dvector( misc_rates);

            for (int i=0; i<(int)temp_regions.size(); i++) {
                std::cout << "we in the looooooop i: " << i << "\n";
                assign_region_rates(temp_regions[i], misc_rates, i);
            }
        }

        /**
        * @brief Assigns rates to region-specific rate catalogs.
        * 
        * @param region A constant pointer to a Region object whose rates will be assigned.
        * @param misc_rates A constant reference to a vector of double values representing miscellaneous rates.
        * @param i An integer index representing the position of the region in the rate catalog.
        */
        void assign_region_rates(const Region* region, const std::vector<double>& misc_rates, int i) {
            std::cout << "assign_region_rates() \n";
            
            int cols = regionrates_100_L.cols();
            for (size_t j=0; j<cols; j++) {
                regionrates_100_L(i,j) = misc_rates[1];
                regionrates_100_R(i,j) = misc_rates[1];
            }

            cols = regionrates_111_L.cols();
            for (size_t j=0; j<cols; j++) {
                regionrates_111_L(i,j) = region->rates[0];
                regionrates_111_R(i,j) = region->rates[1];
            }         
        }

        /**
        * @brief Checks if an adjacent site is unoccupied for a move.
        *
        * This function determines the number of nearest-neighbor vacancies of a given site
        * in a specified lattice configuration.
        *
        * @param i First coordinate of the site.
        * @param j Second coordinate of the site.
        * @param k Third coordinate of the site.
        * @param l Fourth coordinate of the site.
        * @param direc_sign Directional sign indicator.
        * @param s Direction index.
        * @param lattice Type of lattice structure.
        * @return int Number of nearest-neighbor vacancies.
        */
        int get_NNcountofNN(int i, int j, int k, int l, int direc_sign, int s, int lattice) {
            int i1; int i2; int i3; int i4; int direc_sign_NN;
            int i1_NN; int i2_NN; int i3_NN; int i4_NN;
            int i1_NN_unmod; int i2_NN_unmod; int i3_NN_unmod; int i4_NN_unmod;

            std::vector<size_t> x_dims = proc_pos_x_neighbors.size_vec;
            std::vector<size_t> y_dims = proc_pos_y_neighbors.size_vec;
            
            if ((lattice == 2) || (lattice == 3)) {
                if (lattice == 2) { direc_sign_NN = 1; }
                else if (lattice == 3) { direc_sign_NN = 1; }
                i1 = i;
                i2 = (((j + edge_directions[s][0]) % sublattice_dim[0] + sublattice_dim[0]) % sublattice_dim[0]);
                i3 = (((k + edge_directions[s][1]) % sublattice_dim[1] + sublattice_dim[1]) % sublattice_dim[1]);
                i4 = (((l + edge_directions[s][2]) % sublattice_dim[2] + sublattice_dim[2]) % sublattice_dim[2]); 
            }
            else if ((lattice == 0) || (lattice == 1)) {
                if (lattice == 0) { i1 = 1; direc_sign_NN =  1; }
                else if (lattice == 1) { i1 = 0; direc_sign_NN = -1; }
                i2 = (((j + direc_sign * diag_directions[s][0]) % sublattice_dim[0] + sublattice_dim[0]) % sublattice_dim[0]);
                i3 = (((k + direc_sign * diag_directions[s][1]) % sublattice_dim[1] + sublattice_dim[1]) % sublattice_dim[1]);
                i4 = (((l + direc_sign * diag_directions[s][2]) % sublattice_dim[2] + sublattice_dim[2]) % sublattice_dim[2]);
            }

            int NN_count = 0;
            
            for (int s2=0; s2 < (int)diag_directions.rows(); s2++) {
                i1_NN_unmod = !i1;
                i2_NN_unmod = (i2 + direc_sign_NN * diag_directions[s2][0]);
                i3_NN_unmod = (i3 + direc_sign_NN * diag_directions[s2][1]);
                i4_NN_unmod = (i4 + direc_sign_NN * diag_directions[s2][2]);

                
                i1_NN = !i1;
                i2_NN = (((i2 + direc_sign_NN * diag_directions[s2][0]) % sublattice_dim[0] + sublattice_dim[0]) % sublattice_dim[0]);
                i3_NN = (((i3 + direc_sign_NN * diag_directions[s2][1]) % sublattice_dim[1] + sublattice_dim[1]) % sublattice_dim[1]);
                i4_NN = (((i4 + direc_sign_NN * diag_directions[s2][2]) % sublattice_dim[2] + sublattice_dim[2]) % sublattice_dim[2]);

                if (((lattice == 0) || (lattice == 1)) && (i4 == 0) && (i1 == 0) && (diag_directions[s2][2] == 1)) {/* checking for leftmost non-periodic boundary along z-axis*/}
                else if (((lattice == 0) || (lattice == 1)) && (i4 == (int)(sublattice_dim[2]-1)) && (i1 == 1) && (diag_directions[s2][2] == 1)) {/* checking for rightmost non-periodic boundary along z-axis*/}
                else if (((lattice == 2) || (lattice == 3)) && (i4 == 0) && (edge_directions[s][2] == -1)) {}  
                else if (((lattice == 2) || (lattice == 3)) && (i4 == (int)(sublattice_dim[2]-1)) && (edge_directions[s][2] == 1)) {}
                else if ((i1_NN == i) && (i2_NN == j) && (i3_NN == k) && (i4_NN == l)) {}

                else if ((i == 0) && (j == 0) && (diag_directions[s][0] == 1)) {/*check with proc to -x direction*/
                    
                    if ((k == 0) && (diag_directions[s][1] == 1)) {
                        /*check neighbor -x,-y array */
                        if ((proc_neighbors(rank,5) == rank) && (!check_move_free(i1,i2,i3,i4,-1,s,0))) {NN_count ++;}
                        else if ((proc_neighbors(rank,5) != rank) && (proc_neg_x_neighbors(0, 1, 0, (size_t)(i4_NN)))) {NN_count ++;}
                    }
                    else {
                        /*check neighbor -x array */
                        if ((proc_neighbors(rank,4) == rank) && (!check_move_free(i1,i2,i3,i4,-1,s,0))) {NN_count ++;}
                        else if ((proc_neighbors(rank,4) != rank) && (proc_neg_x_neighbors(0, 1, (size_t)(i3_NN+1), (size_t)(i4_NN)))) {NN_count ++;}
                    }
                }
                
                else if ((i == 1) && (j == (sublattice_dim[0] - 1)) && (diag_directions[s][0] == 1)) {/*check with proc to +x direction*/ 
                    
                    if ((k == ((sublattice_dim[1] - 1))) && (diag_directions[s][1] == 1)) {
                        /*check neighbor (+x,+y) array */
                        if (((proc_neighbors(rank,1) == rank)) && (!check_move_free(i1,i2,i3,i4,1,s,1))) {NN_count ++;}
                        else if (((proc_neighbors(rank,1) != rank)) && (proc_pos_x_neighbors(0, 0, (size_t)(x_dims[2]-1), (size_t)(i4_NN)))) {NN_count ++;}
                    }
                    else {
                        /*check neighbor (+x) array */
                        if (((proc_neighbors(rank,0) == rank)) && (!check_move_free(i1,i2,i3,i4,1,s,1))) {NN_count ++;}
                        else if (((proc_neighbors(rank,0) != rank)) && (proc_pos_x_neighbors(0, 0, (size_t)(i3_NN), (size_t)(i4_NN)))) {NN_count ++;}
                    }
                }

                else if ((i == 0) && (k == 0) && (diag_directions[s][1] == 1)) {/*check with proc to -y direction*/
                    
                    if ((j == 0) && (diag_directions[s][0] == 1)) {
                        /*check neighbor (-x,-y) array */
                        if ((proc_neighbors(rank,5) == rank) && (!check_move_free(i1,i2,i3,i4,-1,s,0))) {NN_count ++;}
                        else if ((proc_neighbors(rank,5) != rank) && (proc_neg_y_neighbors(0, 1, 0, (size_t)(i4_NN)))) {NN_count ++;}
                    }
                    else {
                        if (((proc_neighbors(rank,6) == rank)) && (!check_move_free(i1,i2,i3,i4,-1,s,0))) {NN_count ++;}
                        else if (((proc_neighbors(rank,6) != rank)) && (proc_neg_y_neighbors(0, 1, (size_t)(i3_NN+1), (size_t)(i4_NN)))) {NN_count ++;}
                    }
                }

                else if ((i == 1) && (k == (sublattice_dim[1] - 1)) && (diag_directions[s][1] == 1)) {/*check with proc to +y direction*/

                    if ((j == ((sublattice_dim[0] - 1))) && (diag_directions[s][0] == 1)) {
                        /*check neighbor +x,+y array */
                        if ((proc_neighbors(rank,1) == rank) && (!check_move_free(i1,i2,i3,i4,1,s,1))) {NN_count ++;}
                        else if (((proc_neighbors(rank,1) != rank)) && (proc_pos_y_neighbors(0, 0, (size_t)(y_dims[2]-1), (size_t)(i4_NN)))) {NN_count ++;}
                    }
                    else {
                        /*check neighbor +y array */
                        if ((proc_neighbors(rank,2) == rank) && (!check_move_free(i1,i2,i3,i4,1,s,1))) {NN_count ++;}
                        else if ((proc_neighbors(rank,2) != rank) && (proc_pos_y_neighbors(0, 0, (size_t)(i2_NN), (size_t)(i4_NN)))) {NN_count ++;}
                    }
                }

                else if (vacancies(i1_NN,i2_NN,i3_NN,i4_NN)) {NN_count++;} // std::cout << "incriment \n";}
            }

            return NN_count;
        }

        /*
        subroutine for adding move information to matrices stored as attributes of Lattice struc
        */
        int add_move(int i, int j, int k, int l, int curr_move_num, int direc_sign, int s, int idx, int lattice, int curr_NN, int new_NN) {
            double rate = 0; 
            std::vector<int> old_loc = {i,j,k,l};

            if ((lattice == 2) || (lattice == 3)) {
                moves_coords(curr_move_num,0) = i;
                moves_coords(curr_move_num,1) = (((j + edge_directions[s][0]) % sublattice_dim[0] + sublattice_dim[0]) % sublattice_dim[0]);
                moves_coords(curr_move_num,2) = (((k + edge_directions[s][1]) % sublattice_dim[1] + sublattice_dim[1]) % sublattice_dim[1]);
                moves_coords(curr_move_num,3) = (((l + edge_directions[s][2]) % sublattice_dim[2] + sublattice_dim[2]) % sublattice_dim[2]); 
                moves_coords_unmod(curr_move_num,0) = i;
                moves_coords_unmod(curr_move_num,1) = (j + edge_directions[s][0]);
                moves_coords_unmod(curr_move_num,2) = (k + edge_directions[s][1]);
                moves_coords_unmod(curr_move_num,3) = (l + edge_directions[s][2]);                           
                moves_shifts(curr_move_num,0) = edge_directions[s][0];
                moves_shifts(curr_move_num,1) = edge_directions[s][1];
                moves_shifts(curr_move_num,2) = edge_directions[s][2];
                moves_lattice(curr_move_num,0) = lattice;
                moves_vacs(curr_move_num,0) = idx;

                // getting rate corresponding to move
                if ((curr_NN >= void_threshold) && (new_NN >= void_threshold)) { rate = terrace_barrier_100; }  
                else {  rate = new_get_rateconstants(old_loc, moves_shifts[curr_move_num], moves_lattice(curr_move_num,0), curr_NN, new_NN); }
            }
            else if ((lattice == 0) || (lattice == 1)) {
                if (lattice == 0) { 
                    moves_coords(curr_move_num,0) = 1;
                    moves_coords_unmod(curr_move_num,0) = 1; 
                }
                else if (lattice == 1) { 
                    moves_coords(curr_move_num,0) = 0;
                    moves_coords_unmod(curr_move_num,0) = 0; 
                } 
                moves_coords(curr_move_num,1) = (((j + direc_sign * diag_directions[s][0]) % sublattice_dim[0] + sublattice_dim[0]) % sublattice_dim[0]);
                moves_coords(curr_move_num,2) = (((k + direc_sign * diag_directions[s][1]) % sublattice_dim[1] + sublattice_dim[1]) % sublattice_dim[1]);
                moves_coords(curr_move_num,3) = (((l + direc_sign * diag_directions[s][2]) % sublattice_dim[2] + sublattice_dim[2]) % sublattice_dim[2]); 
                moves_coords_unmod(curr_move_num,1) = (j + direc_sign * diag_directions[s][0]);
                moves_coords_unmod(curr_move_num,2) = (k + direc_sign * diag_directions[s][1]);
                moves_coords_unmod(curr_move_num,3) = (l + direc_sign * diag_directions[s][2]);                                   
                moves_shifts(curr_move_num,0) = direc_sign * diag_directions[s][0];
                moves_shifts(curr_move_num,1) = direc_sign * diag_directions[s][1];
                moves_shifts(curr_move_num,2) = direc_sign * diag_directions[s][2]; 
                moves_lattice(curr_move_num,0) = lattice;
                moves_vacs(curr_move_num,0) = idx;

                // getting rate corresponding to move
                rate = get_rateconstants_Elandscape_interface_GB(old_loc, moves_shifts[curr_move_num], lattice, curr_NN, new_NN);
                
            }
            


            //if (rank == 0) std::cout << "rate: " << rate << "\n";
            if (rate == -1) {curr_move_num --;} //case for when move not available in rate catalog
            else {
                if (curr_move_num == 0) {rate_cumsum[curr_move_num] = rate;}
                else {rate_cumsum[curr_move_num] = rate + rate_cumsum[curr_move_num-1];}
            }

            curr_move_num ++;

            return curr_move_num;
        }

        /*
        routine to check if a adjacent site is unoccupied for move
        */
        bool check_move_free(int i, int j, int k, int l, int direc_sign, int s, int lattice) {
            int i1; int i2; int i3; int i4;
            
            if ((lattice == 2) || (lattice == 3)) {
                i1 = i;
                i2 = (((j + edge_directions[s][0]) % sublattice_dim[0] + sublattice_dim[0]) % sublattice_dim[0]);
                i3 = (((k + edge_directions[s][1]) % sublattice_dim[1] + sublattice_dim[1]) % sublattice_dim[1]);
                i4 = (((l + edge_directions[s][2]) % sublattice_dim[2] + sublattice_dim[2]) % sublattice_dim[2]); 
            }
            else if ((lattice == 0) || (lattice == 1)) {
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
            int new_i=0; int new_j=0; int new_k=0; int new_l=0;
            int unmod_new_j=0; int unmod_new_k=0; int unmod_new_l=0;
            int NN_count = 0;

            std::vector<size_t> x_dims = proc_pos_x_neighbors.size_vec;
            std::vector<size_t> y_dims = proc_pos_y_neighbors.size_vec;

            for (int s=0; s < (int)diag_directions.rows(); s++) {

                if (i == 0) {
                    new_i = 1;
                    unmod_new_j = (((j - diag_directions[s][0]) % sublattice_dim[0] + sublattice_dim[0]) % sublattice_dim[0]);
                    unmod_new_k = (((k - diag_directions[s][1]) % sublattice_dim[1] + sublattice_dim[1]) % sublattice_dim[1]);
                    unmod_new_l = (((l - diag_directions[s][2]) % sublattice_dim[2] + sublattice_dim[2]) % sublattice_dim[2]);
                    new_j = (((j - diag_directions[s][0]) % sublattice_dim[0] + sublattice_dim[0]) % sublattice_dim[0]);
                    new_k = (((k - diag_directions[s][1]) % sublattice_dim[1] + sublattice_dim[1]) % sublattice_dim[1]);
                    new_l = (((l - diag_directions[s][2]) % sublattice_dim[2] + sublattice_dim[2]) % sublattice_dim[2]);
                }
                else if (i == 1) {
                    new_i = 0;
                    unmod_new_j = (((j + diag_directions[s][0]) % sublattice_dim[0] + sublattice_dim[0]) % sublattice_dim[0]);
                    unmod_new_k = (((k + diag_directions[s][1]) % sublattice_dim[1] + sublattice_dim[1]) % sublattice_dim[1]);
                    unmod_new_l = (((l + diag_directions[s][2]) % sublattice_dim[2] + sublattice_dim[2]) % sublattice_dim[2]);
                    new_j = (((j + diag_directions[s][0]) % sublattice_dim[0] + sublattice_dim[0]) % sublattice_dim[0]);
                    new_k = (((k + diag_directions[s][1]) % sublattice_dim[1] + sublattice_dim[1]) % sublattice_dim[1]);
                    new_l = (((l + diag_directions[s][2]) % sublattice_dim[2] + sublattice_dim[2]) % sublattice_dim[2]);
                }

                if ((i == 0) && (l == 0) && (diag_directions[s][2] == 1)) {/*check with proc to -z direction*/}
                    
                else if ((l == (sublattice_dim[2]-1)) && (diag_directions[s][2] == 1)) {/*check with proc to +z direction*/}

                else if ((unmod_new_j < 0)) {/*check with proc to -x direction*/   
                    /*check neighbor -x array */
                    //if ((proc_neighbors(rank,4) == rank) && (!check_move_free(i,j,k,l,-1,s,0))) {NN_count ++;}
                    if ((proc_neighbors(rank,4) != rank) && (proc_neg_x_neighbors((size_t)new_i, (size_t)(unmod_new_j%2), (size_t)(unmod_new_k+2), (size_t)(new_l)))) {NN_count ++;}
                
                }                
                else if ((unmod_new_j > (sublattice_dim[0] - 1))) {/*check with proc to +x direction*/ 
                
                    /*check neighbor (+x) array */
                    //if (((proc_neighbors(rank,0) == rank)) && (!check_move_free(i,j,k,l,1,s,1))) {NN_count ++;}
                    if (((proc_neighbors(rank,0) != rank)) && (proc_pos_x_neighbors((size_t)new_i, (size_t)(unmod_new_j%2), (size_t)(unmod_new_k+2), (size_t)(new_l)))) {NN_count ++;}
                
                }

                else if ((unmod_new_k < 0)) { /*check with proc to -y direction*/
                    /*check with proc to -y direction*/
                    //if (((proc_neighbors(rank,6) == rank)) && (!check_move_free(i,j,k,l,-1,s,0))) {NN_count ++;}
                    if (((proc_neighbors(rank,6) != rank)) && (proc_neg_y_neighbors((size_t)new_i, (size_t)(unmod_new_j+2), (size_t)(unmod_new_k%2), (size_t)new_l))) {NN_count ++;}
                    
                }

                else if ((unmod_new_k > (sublattice_dim[1] - 1))) {/*check with proc to +y direction*/
                    /*check neighbor +y array */
                    //if ((proc_neighbors(rank,2) == rank) && (!check_move_free(i,j,k,l,1,s,1))) {NN_count ++;}
                    if ((proc_neighbors(rank,2) != rank) && (proc_pos_y_neighbors((size_t)new_i, (size_t)(unmod_new_j+2), (size_t)(unmod_new_k%2), (size_t)(new_l)))) {NN_count ++;}
                }

                else {
                    if ((i == 0) && (vacancies(1, (((j - diag_directions[s][0]) % sublattice_dim[0] + sublattice_dim[0]) % sublattice_dim[0]), (((k - diag_directions[s][1]) % sublattice_dim[1] + sublattice_dim[1]) % sublattice_dim[1]), (((l - diag_directions[s][2]) % sublattice_dim[2] + sublattice_dim[2]) % sublattice_dim[2])) == 1)) {
                        // checking that vertex site -> bc site move has new site occupied by atom
                        NN_count ++;
                    }
                    else if ((i == 1) && (vacancies(0, (((j + diag_directions[s][0]) % sublattice_dim[0] + sublattice_dim[0]) % sublattice_dim[0]), (((k + diag_directions[s][1]) % sublattice_dim[1] + sublattice_dim[1]) % sublattice_dim[1]), (((l + diag_directions[s][2]) % sublattice_dim[2] + sublattice_dim[2]) % sublattice_dim[2])) == 1)) {
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
            
            int curr_move_num = 0; // total number of moves at this current timestep  
            int vacs_on_interface = 0; // vacancies at last z-index of lattice (used to calculate rate for stripping)
            int num_interface_sites = sublattice_dim[0] * sublattice_dim[1]; // total number of sites at last z-index 
            std::fill(rate_cumsum.begin(), rate_cumsum.end(), 0); // zeroing out entries to rates cumulative sum array
            solo_vacs = 0;

            rate_cumsum.resize((int)moves_coords.rows());
            int i=0; int j=0; int k=0; int l=0; 
            int new_i=0; int new_j=0; int new_k=0; int new_l=0;
            int unmod_new_j=0; int unmod_new_k=0; int unmod_new_l=0;

            std::vector<size_t> x_dims = proc_pos_x_neighbors.size_vec;
            std::vector<size_t> y_dims = proc_pos_y_neighbors.size_vec;

            int bulk_vacs_count = 0;
            int interface_vacs_count = 0;
            int NN_vac=0;
            int NN_vac_test=0;
            int NN_newsite=0;
            system_energy = 0;
            
                        
            // looping over all vacancies in system
            for (int idx=0; idx < (int)vacancies_pos.rows(); idx++) {
                // position in lattice of vacancy 
                i = vacancies_pos(idx,0);
                j = vacancies_pos(idx,1);
                k = vacancies_pos(idx,2);
                l = vacancies_pos(idx,3);

                if ((curr_move_num  >= ((int)moves_shifts.rows() - 20))) {
                    // resizing data structures to accommodate all moves 

                    int newsize = 2 * ((int)moves_shifts.rows() + (int)vacancies_pos.rows() * 14);
                    rate_cumsum.resize(newsize);
                    moves_coords.reshape(newsize, 4, rank);
                    moves_coords_unmod.reshape(newsize, 4, rank);
                    moves_shifts.reshape(newsize, 3, rank);
                    moves_lattice.reshape(newsize, 1, rank);
                    moves_vacs.reshape(newsize, 1, rank);
                }

                int reg_id = region_sites(i, j, k, l);
                NN_vac = get_NN_count(vacancies_pos[idx], i); 
                //NN_vac = get_NN_count_2NNshell(vacancies_pos[idx], i); 
                
                double E_initial = 0;                               
                int curr_NN_SE = 0;
                
                // determining if in region or solid electrolyte region
                if (reg_id != 0) { }
                else if (l == (sublattice_dim[2]-1)) { curr_NN_SE = 1; }
                
                // getting initial site energy
                if (reg_id != 0) { E_initial = regions[(reg_id-1)]->e_below_bulk; }
                else if (curr_NN_SE != 0) { 
                    if ((NN_vac >= void_threshold)) {  E_initial = void_E; }
                    else { E_initial = interface_E; } 
                }              
                else {                                                                                                                                                                                                                                                                                                                                     
                    if ((NN_vac >= void_threshold)) { E_initial = void_E; }
                    else { E_initial = 0;}
                }

                system_energy += E_initial;
                
                if (NN_vac < void_threshold) {
                    if (l > (sublattice_dim[2] - 2)) { 

                        if (rank == 0) {
                            std::cout << "\n rank: " << rank << " NN_vac: " << NN_vac << " interface vac i: " << i << " j: " << j << " k: " << k << " l: " << l << "\n";
                            NN_vac_test = get_NN_count(vacancies_pos[idx], i, true); 
                        }

                        if (rank == 0) std::cout << "interface found: [ "  << i << " " << j << " " << k << " " << l << " ] NN_vac: " << NN_vac <<  "\n";
                        interface_vacs_count ++; 
                    }
                    else {
                        if (rank == 0) {
                            std::cout << "\n rank: " << rank << " NN_vac: " << NN_vac << " bulk vac i: " << i << " j: " << j << " k: " << k << " l: " << l << "\n";
                            NN_vac_test = get_NN_count(vacancies_pos[idx], i, true); 
                        }

                        if (rank == 0) std::cout << "bulk found: [ "  << i << " " << j << " " << k << " " << l << " ] \n";
                        bulk_vacs_count ++; 
                    }
                }

                // finding all moves along the {111} family of vectors
                for (int s=0; s < (int)diag_directions.rows(); s++) {
                    
                    if (i == 0) {
                        unmod_new_j = (j - diag_directions[s][0]);
                        unmod_new_k = (k - diag_directions[s][1]);
                        unmod_new_l = (l - diag_directions[s][2]);
                        new_i = 1;
                        new_j = (((unmod_new_j) % sublattice_dim[0] + sublattice_dim[0]) % sublattice_dim[0]);
                        new_k = (((unmod_new_k) % sublattice_dim[1] + sublattice_dim[1]) % sublattice_dim[1]);
                        new_l = (((unmod_new_l) % sublattice_dim[2] + sublattice_dim[2]) % sublattice_dim[2]);
                    }

                    else if (i == 1) {
                        new_i = 0;
                        unmod_new_j = (j + diag_directions[s][0]);
                        unmod_new_k = (k + diag_directions[s][1]);
                        unmod_new_l = (l + diag_directions[s][2]);
                        new_j = (((unmod_new_j) % sublattice_dim[0] + sublattice_dim[0]) % sublattice_dim[0]);
                        new_k = (((unmod_new_k) % sublattice_dim[1] + sublattice_dim[1]) % sublattice_dim[1]);
                        new_l = (((unmod_new_l) % sublattice_dim[2] + sublattice_dim[2]) % sublattice_dim[2]);
                    }

                    std::vector<int> newsite_coord = {new_i, new_j, new_k, new_l};
                    std::vector<int> oldsite_coord = {i, j, k, l};
                    NN_newsite = get_NN_count(newsite_coord, new_i, oldsite_coord, true);

                    if ((l == 0) && (i == 0) && (diag_directions[s][2] == 1)) {/* checking for leftmost non-periodic boundary along z-axis*/}
                    
                    else if ((l == (int)(sublattice_dim[2]-1)) && (i == 1) && (diag_directions[s][2] == 1)) {/* checking for rightmost non-periodic boundary along z-axis*/}
                    
                    else if (unmod_new_j < 0) {/*communicate with proc to -x direction*/
                        
                        if (unmod_new_k < 0) {
                            /*check neighbor -x,-y array */
                            if ((proc_neighbors(rank,5) == rank) && check_move_free(i,j,k,l,-1,s,0)) {
                                curr_move_num = add_move(i,j,k,l,curr_move_num,-1,s,idx,0,NN_vac, NN_newsite);
                            }
                            else if ((proc_neighbors(rank,5) != rank) && (!proc_neg_x_neighbors(new_i, mod_with_bounds(new_j + chunk_bounds[0][0], 2), (size_t)(unmod_new_k+2), (size_t)(new_l)))) {
                                curr_move_num = add_move(i,j,k,l,curr_move_num,-1,s,idx,0,NN_vac, NN_newsite);
                            }
                        }
                        else {
                            /*check neighbor -x array */
                            if ((proc_neighbors(rank,4) == rank) && check_move_free(i,j,k,l,-1,s,0)) {
                                curr_move_num = add_move(i,j,k,l,curr_move_num,-1,s,idx,0,NN_vac, NN_newsite);
                            }
                            else if ((proc_neighbors(rank,4) != rank) && (!proc_neg_x_neighbors(new_i, mod_with_bounds(new_j + chunk_bounds[0][0], 2), (size_t)(new_k+2), (size_t)(new_l)))) {
                                curr_move_num = add_move(i,j,k,l,curr_move_num,-1,s,idx,0,NN_vac, NN_newsite);
                            }
                        }
                    }
                    
                    else if (j > (sublattice_dim[0] - 1)) {/*communicate with proc to +x direction*/ 
                        
                        if (k > ((sublattice_dim[1] - 1))) {
                            /*check neighbor +x,+y array */
                            if (((proc_neighbors(rank,1) == rank)) && check_move_free(i,j,k,l,1,s,1)) {
                                curr_move_num = add_move(i,j,k,l,curr_move_num,1,s,idx,1,NN_vac, NN_newsite);
                            }
                            else if (((proc_neighbors(rank,1) != rank)) && (!proc_pos_x_neighbors(new_i, mod_with_bounds(new_j + chunk_bounds[0][0], 2), (size_t)(unmod_new_k+2), (size_t)(new_l)))) {
                                curr_move_num = add_move(i,j,k,l,curr_move_num,1,s,idx,1,NN_vac, NN_newsite);
                            }
                        }
                        else {
                            /*check neighbor +x array */
                            if (((proc_neighbors(rank,0) == rank)) && check_move_free(i,j,k,l,1,s,1)) {
                                curr_move_num = add_move(i,j,k,l,curr_move_num,1,s,idx,1,NN_vac, NN_newsite);
                            }
                            else if (((proc_neighbors(rank,0) != rank)) && (!proc_pos_x_neighbors(new_i, mod_with_bounds(new_j + chunk_bounds[0][0], 2), (size_t)(new_k+2), (size_t)(new_l)))) {
                                curr_move_num = add_move(i,j,k,l,curr_move_num,1,s,idx,1,NN_vac, NN_newsite);
                            }
                        }
                    }

                    else if ((i == 0) && (k == 0) && (diag_directions[s][1] == 1)) {/*communicate with proc to -y direction*/
                        
                        if (((proc_neighbors(rank,6) == rank)) && check_move_free(i,j,k,l,-1,s,0)) {
                            curr_move_num = add_move(i,j,k,l,curr_move_num,-1,s,idx,0,NN_vac, NN_newsite);
                        }
                        else if (((proc_neighbors(rank,6) != rank)) && (!proc_neg_y_neighbors(new_i, (size_t)(new_j+2), mod_with_bounds(new_k + chunk_bounds[1][0], 2), (size_t)(new_l)))) {
                            curr_move_num = add_move(i,j,k,l,curr_move_num,-1,s,idx,0,NN_vac, NN_newsite);
                        }
                    }

                    else if ((i == 1) && (k == (sublattice_dim[1] - 1)) && (diag_directions[s][1] == 1)) {/*communicate with proc to +y direction*/
                        
                        if ((proc_neighbors(rank,2) == rank) && check_move_free(i,j,k,l,1,s,1)) {
                            curr_move_num = add_move(i,j,k,l,curr_move_num,1,s,idx,1,NN_vac, NN_newsite);
                        }
                        else if ((proc_neighbors(rank,2) != rank) && (!proc_pos_y_neighbors(new_i, (size_t)(new_j+2), mod_with_bounds(new_k + chunk_bounds[1][0], 2), (size_t)(new_l)))) {
                            curr_move_num = add_move(i,j,k,l,curr_move_num,1,s,idx,1,NN_vac, NN_newsite);
                        }
                    }

                    else {
                        if ((i == 0) && (vacancies(1, (((j - diag_directions[s][0]) % sublattice_dim[0] + sublattice_dim[0]) % sublattice_dim[0]), (((k - diag_directions[s][1]) % sublattice_dim[1] + sublattice_dim[1]) % sublattice_dim[1]), (((l - diag_directions[s][2]) % sublattice_dim[2] + sublattice_dim[2]) % sublattice_dim[2])) == 0)) {
                            // checking that vertex site -> bc site move has new site occupied by atom
                            curr_move_num = add_move(1,j,k,l,curr_move_num,-1,s,idx,0,NN_vac, NN_newsite);
                        }
                        else if ((i == 1) && (vacancies(0, (((j + diag_directions[s][0]) % sublattice_dim[0] + sublattice_dim[0]) % sublattice_dim[0]), (((k + diag_directions[s][1]) % sublattice_dim[1] + sublattice_dim[1]) % sublattice_dim[1]), (((l + diag_directions[s][2]) % sublattice_dim[2] + sublattice_dim[2]) % sublattice_dim[2])) == 0)) {
                            // checking that bc site -> vertex site move has new site occupied by atom
                            curr_move_num = add_move(0,j,k,l,curr_move_num,1,s,idx,1,NN_vac, NN_newsite);
                        }
                    }
                }

                // finding all moves along the {100} family of vectors
                for (int s=0; s < (int)edge_directions.rows(); s++) {

                    new_j = (((j + edge_directions[s][0]) % sublattice_dim[0] + sublattice_dim[0]) % sublattice_dim[0]);
                    new_k = (((k + edge_directions[s][1]) % sublattice_dim[1] + sublattice_dim[1]) % sublattice_dim[1]);
                    new_l = (((l + edge_directions[s][2]) % sublattice_dim[2] + sublattice_dim[2]) % sublattice_dim[2]);

                    unmod_new_j = (j + edge_directions[s][0]);
                    unmod_new_k = (k + edge_directions[s][1]);
                    unmod_new_l = (l + edge_directions[s][2]);

                    std::vector<int> edge_newsite_coord = {i, new_j, new_k, new_l};
                    std::vector<int> edge_oldsite_coord = {i, j, k, l};
                    NN_newsite = get_NN_count(edge_newsite_coord, i, edge_oldsite_coord, true);

                    if ((l == 0) && (edge_directions[s][2] == -1)) {}

                    else if ((l == (int)(sublattice_dim[2]-1)) && (edge_directions[s][2] == 1)) {}

                    else if (unmod_new_j < 0) {/*communicate with proc to -x direction*/
                        if ((proc_neighbors(rank,4) == rank) && check_move_free(i,j,k,l,1,s,(i+2))) {
                            curr_move_num = add_move(i,j,k,l,curr_move_num,1,s,idx,(i+2),NN_vac, NN_newsite);
                        }       
                        else if ((proc_neighbors(rank,4) != rank) && !(proc_neg_x_neighbors(i, mod_with_bounds(new_j + chunk_bounds[0][0], 2), (size_t)(new_k+2),(size_t)(new_l)))) {
                            curr_move_num = add_move(i,j,k,l,curr_move_num,1,s,idx,(i+2),NN_vac, NN_newsite);
                        }                    
                    }
                    else if (unmod_new_j > 0)  {/*communicate with proc to +x direction*/
                        if ((proc_neighbors(rank,0) == rank) && check_move_free(i,j,k,l,1,s,(i+2))) {
                            curr_move_num = add_move(i,j,k,l,curr_move_num,1,s,idx,(i+2),NN_vac, NN_newsite);
                        }
                        else if ((proc_neighbors(rank,0) != rank) && !(proc_pos_x_neighbors(i, mod_with_bounds(new_j + chunk_bounds[0][0], 2), (size_t)(new_k+2),(size_t)(new_l)))) {
                            curr_move_num = add_move(i,j,k,l,curr_move_num,1,s,idx,(i+2),NN_vac, NN_newsite);
                        }                    
                    }
                    else if (unmod_new_k < 0) {/*communicate with proc to -y direction*/
                        if ((proc_neighbors(rank,6) == rank) && check_move_free(i,j,k,l,1,s,(i+2))) {
                            curr_move_num = add_move(i,j,k,l,curr_move_num,1,s,idx,(i+2),NN_vac, NN_newsite);
                        }       
                        if ((proc_neighbors(rank,6) != rank) && !(proc_neg_y_neighbors(i, (size_t)(new_j+2), mod_with_bounds(new_k + chunk_bounds[1][0], 2), (size_t)(new_l)))) {
                            curr_move_num = add_move(i,j,k,l,curr_move_num,1,s,idx,(i+2),NN_vac, NN_newsite);
                        }                    
                    }
                   else if (unmod_new_k > 0) {/*communicate with proc to +y direction*/
                        if ((proc_neighbors(rank,2) == rank) && check_move_free(i,j,k,l,1,s,(i+2))) {
                            curr_move_num = add_move(i,j,k,l,curr_move_num,1,s,idx,(i+2),NN_vac, NN_newsite);
                        }     
                        if ((proc_neighbors(rank,2) != rank) && !(proc_pos_y_neighbors(i, (size_t)(new_j+2), mod_with_bounds(new_k + chunk_bounds[1][0], 2), (size_t)(new_l)))) {
                            curr_move_num = add_move(i,j,k,l,curr_move_num,1,s,idx,(i+2),NN_vac, NN_newsite);
                        }                    
                    }
                    else if (vacancies(i, (((j + edge_directions[s][0]) % sublattice_dim[0] + sublattice_dim[0]) % sublattice_dim[0]), 
                        (((k + edge_directions[s][1]) % sublattice_dim[1] + sublattice_dim[1]) % sublattice_dim[1]), 
                        (((l + edge_directions[s][2]) % sublattice_dim[2] + sublattice_dim[2]) % sublattice_dim[2])) == 0) {
                        // checking that vertex site -> vertex site or bc site -> bc site move has new site occupied by atom
                        curr_move_num = add_move(i,j,k,l,curr_move_num,1,s,idx,(i+2),NN_vac, NN_newsite);
                    }
                }
            }

            // UPDATING SIZE OF DATA STRUCTURES CONTIANING COORDINATES AND RATES OF MOVES
            num_of_moves = curr_move_num;
            rate_cumsum.resize(num_of_moves);
            moves_vacs.reshape(num_of_moves, 1, rank);
            moves_coords.reshape(num_of_moves, 4, rank);
            moves_coords_unmod.reshape(num_of_moves, 4, rank);
            moves_shifts.reshape(num_of_moves, 3, rank);
            moves_lattice.reshape(num_of_moves, 1, rank);

            //if (rank==0){
            std::cout << "rank: " << rank << " bulk_vacs_count: " << bulk_vacs_count << "\n";
            std::cout << "rank: " << rank << "  interface_vacs_count: " << interface_vacs_count << "\n";
            
            MPI_Barrier(MPI_COMM_WORLD);
            //exit(0);

        }


        bool crossing_boundary(int unmod_new_i, int unmod_new_j, int unmod_new_k, int unmod_new_l, int s, int rank) {
            
            if (unmod_new_j < 0) {/*communicate with proc to -x direction*/
                
                if (unmod_new_k < 0) {
                    /*check neighbor -x,-y array */
                    if (proc_neighbors(rank,5) != rank) { return true; }
                }
                else {
                    /*check neighbor -x array */
                    if (proc_neighbors(rank,4) != rank) { return true; }
                }
            }
            
            else if (unmod_new_j > (sublattice_dim[0] - 1)) {/*communicate with proc to +x direction*/ 
                
                if (unmod_new_k > (sublattice_dim[1] - 1)) {
                    /*check neighbor +x,+y array */
                    if (proc_neighbors(rank,1) != rank) { return true; }
                }
                else {
                    /*check neighbor +x array */
                    if (proc_neighbors(rank,0) != rank) { return true; }
                }
            }

            else if (unmod_new_k < 0) {/*communicate with proc to -y direction*/
                
                if (unmod_new_j < 0) {
                    /*check neighbor -x,-y array */
                    if (proc_neighbors(rank,5) != rank) { return true; }
                }
                else {
                    if (proc_neighbors(rank,6) != rank) { return true; }
                }
            }

            else if (unmod_new_k > (sublattice_dim[1] - 1)) {/*communicate with proc to +y direction*/
                
                if (unmod_new_j > (sublattice_dim[0] - 1)) {
                    /*check neighbor +x,+y array */
                    if (proc_neighbors(rank,1) != rank) { return true; }
                }
                else {
                    /*check neighbor +y array */
                    if (proc_neighbors(rank,2) != rank) { return true; }
                }
            }

            return false;
        }


        bool check_for_vacancy_boundary(int unmod_new_i, int unmod_new_j, int unmod_new_k, int unmod_new_l, int new_i, int new_j, int new_k, int new_l, int s, int rank) {

            std::vector<size_t> x_dims = proc_pos_x_neighbors.size_vec;
            std::vector<size_t> y_dims = proc_pos_y_neighbors.size_vec;

            if (unmod_new_j < 0) {/*communicate with proc to -x direction*/
                
                if (unmod_new_k < 0) {
                    /*check neighbor -x,-y array */
                    if (proc_neg_x_neighbors(new_i, mod_with_bounds(new_j + chunk_bounds[0][0], 2), (size_t)(unmod_new_k+2), (size_t)(new_l))) {
                        return true;
                    }
                }
                else if (unmod_new_k > (sublattice_dim[1]-1)) {
                    if (proc_neg_x_neighbors(new_i, mod_with_bounds(new_j + chunk_bounds[0][0], 2), (size_t)(unmod_new_k+2), (size_t)(new_l)))  {
                        return true;
                    }
                }
                else {
                    /*check neighbor -x array */
                    if (proc_neg_x_neighbors(new_i, mod_with_bounds(new_j + chunk_bounds[0][0], 2), (size_t)(unmod_new_k+2), (size_t)(new_l))) {
                        return true;
                    }
                }
            }
            
            else if ((unmod_new_j > (sublattice_dim[0] - 1))) {/*check with proc to +x direction*/ 
                /*check neighbor (+x) array */
                
                if (unmod_new_k < 0) {
                    if (proc_pos_x_neighbors(new_i, mod_with_bounds(new_j + chunk_bounds[0][0], 2), (size_t)(unmod_new_k+2), (size_t)(new_l))) {
                        return true;
                    }
                }
                else if (unmod_new_k > (sublattice_dim[1]-1)) {
                    if (proc_pos_x_neighbors(new_i, mod_with_bounds(new_j + chunk_bounds[0][0], 2), (size_t)(unmod_new_k+2), (size_t)(new_l))) {
                        return true;
                    }
                }
                else {
                    if (proc_pos_x_neighbors(new_i, mod_with_bounds(new_j + chunk_bounds[0][0], 2), (size_t)(unmod_new_k+2), (size_t)(new_l))) {
                        return true;
                    }
                }

            }

            else if ((unmod_new_k < 0)) { /*check with proc to -y direction*/
                /*check with proc to -y direction*/
                if (proc_neg_y_neighbors(new_i, (size_t)(unmod_new_j+2), mod_with_bounds(new_k + chunk_bounds[1][0], 2), (size_t)new_l)) {
                    return true;
                }
                
            }

            else if ((unmod_new_k > (sublattice_dim[1] - 1))) {/*check with proc to +y direction*/
                /*check neighbor +y array */
                if (proc_pos_y_neighbors(new_i, (size_t)(unmod_new_j+2), mod_with_bounds(new_k + chunk_bounds[1][0], 2), (size_t)(new_l))) {
                    return true;
                }
            }

            return false;
        }   


        double get_E_of_NN_void_in_reg(std::vector<int>& init_vec, std::vector<int>& dest_vec, int lattice, bool in_initial_state, bool debug=false) {       
            double NN_count = 0; double total_E = 0;
            int i = init_vec[0]; int j = init_vec[1]; int k = init_vec[2]; int l = init_vec[3];
            int dest_i1 = dest_vec[0]; int dest_i2 = dest_vec[1]; int dest_i3 = dest_vec[2]; int dest_i4 = dest_vec[3];

            int i1; int i2; int i3; int i4; int direc_sign_NN;
            int i1_unmod; int i2_unmod; int i3_unmod; int i4_unmod;
            int i1_NN; int i2_NN; int i3_NN; int i4_NN;
            int i1_NN_unmod; int i2_NN_unmod; int i3_NN_unmod; int i4_NN_unmod;
            if (debug) std::cout << "lattice: " << lattice << " i: " << i << " j: " << j << " k: " << k << " l: " << l << "\n";   
            std::set< std::vector<int> > used_vecs;
            int reg_id;

            for (int s1=0; s1 < (int)diag_directions.rows(); s1++) {
                NN_count = 0;
                
                // getting coordinates of NN of initial site 
                if (i == 0) { i1 = 1; direc_sign_NN = -1; }
                else if (i == 1) { i1 = 0; direc_sign_NN = 1; }
                i2 = (((j + direc_sign_NN * diag_directions[s1][0]) % sublattice_dim[0] + sublattice_dim[0]) % sublattice_dim[0]);
                i3 = (((k + direc_sign_NN * diag_directions[s1][1]) % sublattice_dim[1] + sublattice_dim[1]) % sublattice_dim[1]);
                i4 = (((l + direc_sign_NN * diag_directions[s1][2]) % sublattice_dim[2] + sublattice_dim[2]) % sublattice_dim[2]);
                
                i1_unmod = i1;
                i2_unmod = (j + direc_sign_NN * diag_directions[s1][0]);
                i3_unmod = (k + direc_sign_NN * diag_directions[s1][1]);
                i4_unmod = (l + direc_sign_NN * diag_directions[s1][2]);
                //if (rank == 0) std::cout << "rank: " << rank << " i1: " << i1 << " i2: " << i2 << " i3: " << i3 << " i4: " << i4 << "\n";
                //if (rank == 0) std::cout << "rank: " << rank << " i1_unmod: " << i1_unmod << " i2_unmod: " << i2_unmod << " i3_unmod: " << i3_unmod << " i4_unmod: " << i4_unmod << "\n";
                if ((i1 == i) && (i2 == j) && (i3 == k) && (i4 == l)) { }
                else if ((i1 == dest_i1) && (i2 == dest_i2) && (i3 == dest_i3) && (i4 == dest_i4)) { }
                if (i4_unmod < 0 ) { /* checking for leftmost non-periodic boundary along z-axis*/}
                else if (i4_unmod > (int)(sublattice_dim[2]-1)) {  /* checking for rightmost non-periodic boundary along z-axis*/}                                  
                else { 
                    if (crossing_boundary(i1_unmod,i2_unmod,i3_unmod,i4_unmod,s1,rank)) {
                        if (used_vecs.count({i1_unmod,i2_unmod,i3_unmod,i4_unmod})) {}
                        else {
                            used_vecs.insert({i1_unmod,i2_unmod,i3_unmod,i4_unmod});
                            if (check_for_vacancy_boundary(i1,i2_unmod,i3_unmod,i4_unmod,i1,i2,i3,i4,s1,rank)) {
                                
                                if (debug) std::cout << "i1: " << i1 << " i2: " << i2 << " i3: " << i3 << " i4: " << i4 << "\n";  
                                for (int s2=0; s2 < (int)diag_directions.rows(); s2++) {

                                    // getting coordinates of NN of NN
                                    if (i1 == 0) { direc_sign_NN = -1; }
                                    else if (i1 == 1) { direc_sign_NN = 1; }
                                    i1_NN = !i1_unmod;
                                    i2_NN = (((i2_unmod + direc_sign_NN * diag_directions[s2][0]) % sublattice_dim[0] + sublattice_dim[0]) % sublattice_dim[0]);
                                    i3_NN = (((i3_unmod + direc_sign_NN * diag_directions[s2][1]) % sublattice_dim[1] + sublattice_dim[1]) % sublattice_dim[1]);
                                    i4_NN = (((i4_unmod + direc_sign_NN * diag_directions[s2][2]) % sublattice_dim[2] + sublattice_dim[2]) % sublattice_dim[2]);

                                    i1_NN_unmod = !i1_unmod;
                                    i2_NN_unmod = (i2_unmod + direc_sign_NN * diag_directions[s2][0]);
                                    i3_NN_unmod = (i3_unmod + direc_sign_NN * diag_directions[s2][1]);
                                    i4_NN_unmod = (i4_unmod + direc_sign_NN * diag_directions[s2][2]);
                                    /*
                                        if (rank == 0) std::cout << "rank: " << rank << " i1_NN: " << i1_NN << " i2_NN: " << i2_NN << " i3_NN: " << i3_NN << " i4_NN: " << i4_NN << "\n";
                                        if (rank == 0) std::cout << "rank: " << rank << " i1_NN_unmod: " << i1_NN_unmod << " i2_NN_unmod: " << i2_NN_unmod << " i3_NN_unmod: " << i3_NN_unmod << " i4_NN_unmod: " << i4_NN_unmod << "\n";
                                    */
                                    if (i4_NN_unmod < 0) { /* checking for leftmost non-periodic boundary along z-axis*/}
                                    else if (i4_NN_unmod > (sublattice_dim[2]-1)) {  /* checking for rightmost non-periodic boundary along z-axis*/}
                                    else if (crossing_boundary(i1_NN_unmod,i2_NN_unmod,i3_NN_unmod,i4_NN_unmod,s2,rank)) {
                                        if ((i1_NN == i) && (i2_NN == j) && (i3_NN == k) && (i4_NN == l)) {
                                            if (in_initial_state) {
                                                //std::cout<< "in final state accept\n";
                                                NN_count ++;}
                                        }
                                        else if ((i1_NN == dest_i1) && (i2_NN == dest_i2) && (i3_NN == dest_i3) && (i4_NN == dest_i4)) {
                                            if (!in_initial_state) {
                                                //std::cout<< "in initial state accept\n";
                                                NN_count ++;}
                                        }
                                        else if (check_for_vacancy_boundary(i1_NN_unmod,i2_NN_unmod,i3_NN_unmod,i4_NN_unmod,i1_NN,i2_NN,i3_NN,i4_NN,s2,rank)) { NN_count ++; }
                                        else { /* no vacancy */ }
                                    }                                
                                    else { 
                                        if ((i1_NN == i) && (i2_NN == j) && (i3_NN == k) && (i4_NN == l)) {
                                            if (in_initial_state) {
                                                //std::cout<< "in final state accept\n";
                                                NN_count ++;}
                                        }
                                        else if ((i1_NN == dest_i1) && (i2_NN == dest_i2) && (i3_NN == dest_i3) && (i4_NN == dest_i4)) {
                                            if (!in_initial_state) {
                                                //std::cout<< "in initial state accept\n";
                                                NN_count ++;}
                                        }
                                        else if (vacancies(i1_NN,i2_NN,i3_NN,i4_NN)) {
                                            //std::cout << "accepted \n";    
                                            NN_count++;
                                        }
                                    }
                                }
                                //std::cout << "NN_count: " << NN_count << "\n";
                                reg_id = region_sites(i1, i2, i3, i4);

                                if (NN_count >= void_threshold) {
                                    if (i4 == (sublattice_dim[2]-1)) { total_E += void_E; }
                                    else if (reg_id != 0) { }
                                    else { total_E += void_E; }
                                }
                                else if (i4 == (sublattice_dim[2]-1)) { total_E += interface_E; }
                                
                            }
                        }
                    }  
                    else if (vacancies(i1_unmod,i2_unmod,i3_unmod,i4_unmod)) {
                        
                        if (used_vecs.count({i1,i2,i3,i4})) {}
                        else {
                            used_vecs.insert({i1,i2,i3,i4});
                            //std::cout << "loop 1 used_vecs: ";
                            //print_set(used_vecs);
                            if (debug) std::cout << "i1: " << i1 << " i2: " << i2 << " i3: " << i3 << " i4: " << i4 << "\n";  
                            
                            for (int s2=0; s2 < (int)diag_directions.rows(); s2++) {
                                // getting coordinates of NN of NN
                                if (i1 == 0) { direc_sign_NN = -1; }
                                else if (i1 == 1) { direc_sign_NN = 1; }
                                i1_NN = !i1_unmod;
                                i2_NN = (((i2_unmod + direc_sign_NN * diag_directions[s2][0]) % sublattice_dim[0] + sublattice_dim[0]) % sublattice_dim[0]);
                                i3_NN = (((i3_unmod + direc_sign_NN * diag_directions[s2][1]) % sublattice_dim[1] + sublattice_dim[1]) % sublattice_dim[1]);
                                i4_NN = (((i4_unmod + direc_sign_NN * diag_directions[s2][2]) % sublattice_dim[2] + sublattice_dim[2]) % sublattice_dim[2]);

                                i1_NN_unmod = !i1_unmod;
                                i2_NN_unmod = (i2_unmod + direc_sign_NN * diag_directions[s2][0]);
                                i3_NN_unmod = (i3_unmod + direc_sign_NN * diag_directions[s2][1]);
                                i4_NN_unmod = (i4_unmod + direc_sign_NN * diag_directions[s2][2]);
                                /*
                                    if (rank == 0) std::cout << "rank: " << rank << " i1_NN: " << i1_NN << " i2_NN: " << i2_NN << " i3_NN: " << i3_NN << " i4_NN: " << i4_NN << "\n";
                                    if (rank == 0) std::cout << "rank: " << rank << " i1_NN_unmod: " << i1_NN_unmod << " i2_NN_unmod: " << i2_NN_unmod << " i3_NN_unmod: " << i3_NN_unmod << " i4_NN_unmod: " << i4_NN_unmod << "\n";
                                */
                                //std::cout << "loc_NN: [ " << i1_NN << " " << i2_NN << " " << i3_NN << " " << i4_NN << " ]\n";    
                                if ((i4 == 0) && (i1 == 0) && (diag_directions[s2][2] == 1)) { /* checking for leftmost non-periodic boundary along z-axis*/}
                                else if ((i4 == (int)(sublattice_dim[2]-1)) && (i1 == 1) && (diag_directions[s2][2] == 1)) {  /* checking for rightmost non-periodic boundary along z-axis*/}
                                else if (crossing_boundary(i1_NN_unmod,i2_NN_unmod,i3_NN_unmod,i4_NN_unmod,s2,rank)) {
                                    if ((i1_NN == i) && (i2_NN == j) && (i3_NN == k) && (i4_NN == l)) {
                                        if (in_initial_state) {
                                            //std::cout<< "in final state accept\n";
                                            NN_count ++;}
                                    }
                                    else if ((i1_NN == dest_i1) && (i2_NN == dest_i2) && (i3_NN == dest_i3) && (i4_NN == dest_i4)) {
                                        if (!in_initial_state) {
                                            //std::cout<< "in initial state accept\n";
                                            NN_count ++;}
                                    }                                    
                                    else if (check_for_vacancy_boundary(i1_NN_unmod,i2_NN_unmod,i3_NN_unmod,i4_NN_unmod,i1_NN,i2_NN,i3_NN,i4_NN,s2,rank)) { NN_count ++; }
                                    else { /* no vacancy */ }
                                }                                
                                else { 
                                    if ((i1_NN == i) && (i2_NN == j) && (i3_NN == k) && (i4_NN == l)) {
                                        if (in_initial_state) {
                                            //std::cout<< "in final state accept\n";
                                            NN_count ++;}
                                    }
                                    else if ((i1_NN == dest_i1) && (i2_NN == dest_i2) && (i3_NN == dest_i3) && (i4_NN == dest_i4)) {
                                        if (!in_initial_state) {
                                            //std::cout<< "in initial state accept\n";
                                            NN_count ++;}
                                    }
                                    else if (vacancies(i1_NN_unmod,i2_NN_unmod,i3_NN_unmod,i4_NN_unmod)) {
                                        //std::cout << "accepted \n";    
                                        NN_count++;
                                    }
                                }
                            }
                                
                            
                            //std::cout << "NN_count: " << NN_count << "\n";
                            reg_id = region_sites(i1, i2, i3, i4);

                            if (NN_count >= void_threshold) {
                                if (i4 == (sublattice_dim[2]-1)) { total_E += void_E; }
                                else if (reg_id != 0) { }
                                else { total_E += void_E; }
                            }
                            else if (i4 == (sublattice_dim[2]-1)) { total_E += interface_E; }

                        }                 
                    }
                }
            }

            for (int s1=0; s1 < (int)diag_directions.rows(); s1++) {
                NN_count = 0;

                // getting coordinates of NN of initial site 
                if (dest_i1 == 0) { i1 = 1; direc_sign_NN = -1; }
                else if (dest_i1 == 1) { i1 = 0; direc_sign_NN = 1; }
                i2 = (((dest_i2 + direc_sign_NN * diag_directions[s1][0]) % sublattice_dim[0] + sublattice_dim[0]) % sublattice_dim[0]);
                i3 = (((dest_i3 + direc_sign_NN * diag_directions[s1][1]) % sublattice_dim[1] + sublattice_dim[1]) % sublattice_dim[1]);
                i4 = (((dest_i4 + direc_sign_NN * diag_directions[s1][2]) % sublattice_dim[2] + sublattice_dim[2]) % sublattice_dim[2]);
                
                i1_unmod = dest_i1;
                i2_unmod = (dest_i2 + direc_sign_NN * diag_directions[s1][0]);
                i3_unmod = (dest_i3 + direc_sign_NN * diag_directions[s1][1]);
                i4_unmod = (dest_i4 + direc_sign_NN * diag_directions[s1][2]);
                /*
                    if (rank == 0) std::cout << "rank: " << rank << " second loop i1: " << i1 << " i2: " << i2 << " i3: " << i3 << " i4: " << i4 << "\n";
                    if (rank == 0) std::cout << "rank: " << rank << " second loop i1_unmod: " << i1_unmod << " i2_unmod: " << i2_unmod << " i3_unmod: " << i3_unmod << " i4_unmod: " << i4_unmod << "\n";
                    if (rank == 0) std::cout << "rank: " << rank << " ((i1 == i) && (i2 == j) && (i3 == k) && (i4 == l)): " << ((i1 == i) && (i2 == j) && (i3 == k) && (i4 == l)) << "\n";
                    if (rank == 0) std::cout << "rank: " << rank << " ((i1 == dest_i1) && (i2 == dest_i2) && (i3 == dest_i3) && (i4 == dest_i4)): " << ((i1 == dest_i1) && (i2 == dest_i2) && (i3 == dest_i3) && (i4 == dest_i4)) << "\n";
                    if (rank == 0) std::cout << "rank: " << rank << " (crossing_boundary(i1_unmod,i2_unmod,i3_unmod,i4_unmod,s1,rank)): " << crossing_boundary(i1_unmod,i2_unmod,i3_unmod,i4_unmod,s1,rank) << "\n";
                */
                if ((i1 == i) && (i2 == j) && (i3 == k) && (i4 == l)) { }
                else if ((i1 == dest_i1) && (i2 == dest_i2) && (i3 == dest_i3) && (i4 == dest_i4)) { }
                if (i4_unmod < 0) { /* checking for leftmost non-periodic boundary along z-axis*/}
                else if (i4_unmod > (int)(sublattice_dim[2]-1)) {  /* checking for rightmost non-periodic boundary along z-axis*/}                                  
                else { 
                    if (crossing_boundary(i1_unmod,i2_unmod,i3_unmod,i4_unmod,s1,rank)) {
                        //if (rank == 0) std::cout << "rank: " << rank << " crossing boundary outer\n";
                        if (used_vecs.count({i1_unmod,i2_unmod,i3_unmod,i4_unmod})) {}
                        else {
                            used_vecs.insert({i1_unmod,i2_unmod,i3_unmod,i4_unmod});
                            if (check_for_vacancy_boundary(i1_unmod,i2_unmod,i3_unmod,i4_unmod,i1,i2,i3,i4,s1,rank)) {
                                //std::cout << "loop 1 used_vecs: ";
                                //print_set(used_vecs);
                                if (debug) std::cout << "i1: " << i1 << " i2: " << i2 << " i3: " << i3 << " i4: " << i4 << "\n";  
                                for (int s2=0; s2 < (int)diag_directions.rows(); s2++) {

                                    // getting coordinates of NN of NN
                                    if (i1 == 0) { direc_sign_NN = -1; }
                                    else if (i1 == 1) { direc_sign_NN = 1; }
                                    i1_NN = !i1_unmod;
                                    i2_NN = (((i2_unmod + direc_sign_NN * diag_directions[s2][0]) % sublattice_dim[0] + sublattice_dim[0]) % sublattice_dim[0]);
                                    i3_NN = (((i3_unmod + direc_sign_NN * diag_directions[s2][1]) % sublattice_dim[1] + sublattice_dim[1]) % sublattice_dim[1]);
                                    i4_NN = (((i4_unmod + direc_sign_NN * diag_directions[s2][2]) % sublattice_dim[2] + sublattice_dim[2]) % sublattice_dim[2]);

                                    i1_NN_unmod = !i1_unmod;
                                    i2_NN_unmod = (i2_unmod + direc_sign_NN * diag_directions[s2][0]);
                                    i3_NN_unmod = (i3_unmod + direc_sign_NN * diag_directions[s2][1]);
                                    i4_NN_unmod = (i4_unmod + direc_sign_NN * diag_directions[s2][2]);
                                    /*
                                        if (rank == 0) std::cout << "rank: " << rank << " second loop i1_NN: " << i1_NN << " i2_NN: " << i2_NN << " i3_NN: " << i3_NN << " i4_NN: " << i4_NN << "\n";
                                        if (rank == 0) std::cout << "rank: " << rank << " second loop i1_NN_unmod: " << i1_NN_unmod << " i2_NN_unmod: " << i2_NN_unmod << " i3_NN_unmod: " << i3_NN_unmod << " i4_NN_unmod: " << i4_NN_unmod << "\n";
                                    */
                                    //std::cout << "loc_NN: [ " << i1_NN << " " << i2_NN << " " << i3_NN << " " << i4_NN << " ]\n";    
                                    if (i4_NN_unmod < 0) { /* checking for leftmost non-periodic boundary along z-axis*/}
                                    else if (i4_NN_unmod > (int)(sublattice_dim[2]-1)) {  /* checking for rightmost non-periodic boundary along z-axis*/}
                                    else if (crossing_boundary(i1_NN_unmod,i2_NN_unmod,i3_NN_unmod,i4_NN_unmod,s2,rank)) {
                                        // if (rank == 0) std::cout << "rank: " << rank << " second loop crossing boundary \n";
                                        if ((i1_NN == i) && (i2_NN == j) && (i3_NN == k) && (i4_NN == l)) {
                                            if (in_initial_state) {
                                                //std::cout<< "in final state accept\n";
                                                NN_count ++;}
                                        }
                                        else if ((i1_NN == dest_i1) && (i2_NN == dest_i2) && (i3_NN == dest_i3) && (i4_NN == dest_i4)) {
                                            if (!in_initial_state) {
                                                //std::cout<< "in initial state accept\n";
                                                NN_count ++;}
                                        }                                    
                                        else if (check_for_vacancy_boundary(i1_NN_unmod,i2_NN_unmod,i3_NN_unmod,i4_NN_unmod,i1_NN,i2_NN,i3_NN,i4_NN,s2,rank)) { NN_count ++; }
                                        else { /* no vacancy */ }
                                    }                                
                                    else { 
                                        // if (rank == 0) std::cout << "rank: " << rank << " second loop else statement \n";
                                        if ((i1_NN == i) && (i2_NN == j) && (i3_NN == k) && (i4_NN == l)) {
                                            if (in_initial_state) {
                                                //std::cout<< "in final state accept\n";
                                                NN_count ++;}
                                        }
                                        else if ((i1_NN == dest_i1) && (i2_NN == dest_i2) && (i3_NN == dest_i3) && (i4_NN == dest_i4)) {
                                            if (!in_initial_state) {
                                                //std::cout<< "in initial state accept\n";
                                                NN_count ++;}
                                        }
                                        else if (vacancies(i1_NN_unmod,i2_NN_unmod,i3_NN_unmod,i4_NN_unmod)) {
                                            //std::cout << "accepted \n";    
                                            NN_count++;
                                        }
                                    }
                                }
                                //std::cout << "NN_count: " << NN_count << "\n";
                                reg_id = region_sites(i1, i2, i3, i4);

                                if (NN_count >= void_threshold) {
                                    if (i4 == (sublattice_dim[2]-1)) { total_E += void_E; }
                                    else if (reg_id != 0) { }
                                    else { total_E += void_E; }
                                }
                                else if (i4 == (sublattice_dim[2]-1)) { total_E += interface_E; }                                
                            }
                        }
                    }  
                    else if (vacancies(i1_unmod,i2_unmod,i3_unmod,i4_unmod)) {     
                        // if (rank == 0) std::cout << "rank: " << rank << " else outer\n";                   
                        if (used_vecs.count({i1_unmod,i2_unmod,i3_unmod,i4_unmod})) {}
                        else {
                            used_vecs.insert({i1_unmod,i2_unmod,i3_unmod,i4_unmod});
                            //std::cout << "loop 1 used_vecs: ";
                            //print_set(used_vecs);
                            if (debug) std::cout << "i1: " << i1 << " i2: " << i2 << " i3: " << i3 << " i4: " << i4 << "\n";  
                            
                            for (int s2=0; s2 < (int)diag_directions.rows(); s2++) {

                                // getting coordinates of NN of NN
                                if (i1 == 0) { direc_sign_NN = -1; }
                                else if (i1 == 1) { direc_sign_NN = 1; }
                                i1_NN = !i1_unmod;
                                i2_NN = (((i2_unmod + direc_sign_NN * diag_directions[s2][0]) % sublattice_dim[0] + sublattice_dim[0]) % sublattice_dim[0]);
                                i3_NN = (((i3_unmod + direc_sign_NN * diag_directions[s2][1]) % sublattice_dim[1] + sublattice_dim[1]) % sublattice_dim[1]);
                                i4_NN = (((i4_unmod + direc_sign_NN * diag_directions[s2][2]) % sublattice_dim[2] + sublattice_dim[2]) % sublattice_dim[2]);

                                i1_NN_unmod = !i1_unmod;
                                i2_NN_unmod = (i2_unmod + direc_sign_NN * diag_directions[s2][0]);
                                i3_NN_unmod = (i3_unmod + direc_sign_NN * diag_directions[s2][1]);
                                i4_NN_unmod = (i4_unmod + direc_sign_NN * diag_directions[s2][2]);
                                /*
                                    if (rank == 0) std::cout << "rank: " << rank << " second loop i1_NN 2 : " << i1_NN << " i2_NN: " << i2_NN << " i3_NN: " << i3_NN << " i4_NN: " << i4_NN << "\n";
                                    if (rank == 0) std::cout << "rank: " << rank << " second loop i1_NN_unmod 2: " << i1_NN_unmod << " i2_NN_unmod: " << i2_NN_unmod << " i3_NN_unmod: " << i3_NN_unmod << " i4_NN_unmod: " << i4_NN_unmod << "\n";
                                */
                                //std::cout << "loc_NN: [ " << i1_NN << " " << i2_NN << " " << i3_NN << " " << i4_NN << " ]\n";    
                                if (i4_NN_unmod < 0) { /* checking for leftmost non-periodic boundary along z-axis*/}
                                else if (i4_NN_unmod > (int)(sublattice_dim[2]-1)) {  /* checking for rightmost non-periodic boundary along z-axis*/}
                                else if (crossing_boundary(i1_NN_unmod,i2_NN_unmod,i3_NN_unmod,i4_NN_unmod,s2,rank)) {
                                    // if (rank == 0) std::cout << "rank: " << rank << " second loop crossing boundary 2\n";
                                    if ((i1_NN == i) && (i2_NN == j) && (i3_NN == k) && (i4_NN == l)) {
                                        if (in_initial_state) {
                                            //std::cout<< "in final state accept\n";
                                            NN_count ++;}
                                    }
                                    else if ((i1_NN == dest_i1) && (i2_NN == dest_i2) && (i3_NN == dest_i3) && (i4_NN == dest_i4)) {
                                        if (!in_initial_state) {
                                            //std::cout<< "in initial state accept\n";
                                            NN_count ++;}
                                    }                                    
                                    else if (check_for_vacancy_boundary(i1_NN_unmod,i2_NN_unmod,i3_NN_unmod,i4_NN_unmod,i1_NN,i2_NN,i3_NN,i4_NN,s2,rank)) { NN_count ++; }
                                    else { /* no vacancy */ }
                                }                                
                                else { 
                                    // if (rank == 0) std::cout << "rank: " << rank << " second loop else 2: \n";
                                    if ((i1_NN == i) && (i2_NN == j) && (i3_NN == k) && (i4_NN == l)) {
                                        if (in_initial_state) {
                                            //std::cout<< "in final state accept\n";
                                            NN_count ++;}
                                    }
                                    else if ((i1_NN == dest_i1) && (i2_NN == dest_i2) && (i3_NN == dest_i3) && (i4_NN == dest_i4)) {
                                        if (!in_initial_state) {
                                            //std::cout<< "in initial state accept\n";
                                            NN_count ++;}
                                    }
                                    else if (vacancies(i1_NN,i2_NN,i3_NN,i4_NN)) {
                                        //std::cout << "accepted \n";    
                                        NN_count++;
                                    }
                                }
                            }
                                
                            
                            //std::cout << "NN_count: " << NN_count << "\n";
                            reg_id = region_sites(i1, i2, i3, i4);

                            if (NN_count >= void_threshold) {
                                if (i4 == (sublattice_dim[2]-1)) { total_E += void_E; }
                                else if (reg_id != 0) { }
                                else { total_E += void_E; }
                            }
                            else if (i4 == (sublattice_dim[2]-1)) { total_E += interface_E; }

                        }                 
                    }
                }
            }
            
            return total_E;
        }


        /**
        * @brief Checks energy of neighboring site with and without occupancy in current site.
        *
        * This function determines the number of nearest-neighbor vacancies of a given site
        * in a specified lattice configuration.
        *
        * @param i First coordinate of the site.
        * @param j Second coordinate of the site.
        * @param k Third coordinate of the site.
        * @param l Fourth coordinate of the site.
        * @param direc_sign Directional sign indicator.
        * @param s Direction index.
        * @param lattice Type of lattice structure.
        * @return int Number of nearest-neighbor vacancies.
        */
        double get_E_of_NN_wrapper(int i, int j, int k, int l, int* shift, int lattice) {
            
            int dest_i1; int dest_i2; int dest_i3; int dest_i4;
            
            if ((lattice == 2) || (lattice == 3)) {
                dest_i1 = i;
                dest_i2 = (((j + shift[0]) % sublattice_dim[0] + sublattice_dim[0]) % sublattice_dim[0]);
                dest_i3 = (((k + shift[1]) % sublattice_dim[1] + sublattice_dim[1]) % sublattice_dim[1]);
                dest_i4 = (((l + shift[2]) % sublattice_dim[2] + sublattice_dim[2]) % sublattice_dim[2]); 
            }
            else if ((lattice == 0) || (lattice == 1)) {
                if (lattice == 0) { dest_i1 = 1; }
                else if (lattice == 1) { dest_i1 = 0; }
                dest_i2 = (((j + shift[0]) % sublattice_dim[0] + sublattice_dim[0]) % sublattice_dim[0]);
                dest_i3 = (((k + shift[1]) % sublattice_dim[1] + sublattice_dim[1]) % sublattice_dim[1]);
                dest_i4 = (((l + shift[2]) % sublattice_dim[2] + sublattice_dim[2]) % sublattice_dim[2]);
            }

            std::vector<int> init_coords{i,j,k,l};
            std::vector<int> final_coords = {dest_i1,dest_i2,dest_i3,dest_i4};
            
            
            double init_E = get_E_of_NN_void_in_reg(init_coords, final_coords, lattice, true);
            double final_E = get_E_of_NN_void_in_reg(init_coords, final_coords, lattice, false);
            double energy_difference = final_E - init_E;
            //double energy_difference = delta_E_site_i + delta_E_site_f;
            //std::cout << "energy_difference: " << energy_difference << "\n";

            return energy_difference;
        }

        double get_rateconstants_Elandscape_interface_GB(std::vector<int> coord, int* shift, int lattice, int curr_NN, int new_NN) {  
            //std::cout << "entering get_rateconstants()\n";
            double rate = -1;
            double E_initial = 0;
            double E_final = 0;

            int i_new; 
            if ((lattice == 0) || (lattice == 3)) {i_new = 1;}
            else {i_new = 0;}
            int j_new = (((coord[1] + shift[0]) % sublattice_dim[0] + sublattice_dim[0]) % sublattice_dim[0]);
            int k_new = (((coord[2] + shift[1]) % sublattice_dim[1] + sublattice_dim[1]) % sublattice_dim[1]);
            int l_new = (((coord[3] + shift[2]) % sublattice_dim[2] + sublattice_dim[2]) % sublattice_dim[2]);
            
            int reg_id = region_sites(coord[0], coord[1], coord[2], coord[3]);
            int new_reg_id = region_sites(i_new, j_new, k_new, l_new);
            
            int curr_NN_SE = 0;
            int new_NN_SE = 0;

            // determining if in region or solid electrolyte region
            if (reg_id != 0) { }
            else if (coord[3] == (sublattice_dim[2]-1)) { curr_NN_SE = 1; }
            if (new_reg_id != 0) { }
            else if (l_new == (sublattice_dim[2]-1)) { new_NN_SE = 1; }

            
            // getting initial site energy
            if (reg_id != 0) { E_initial = regions[(reg_id-1)]->e_below_bulk; }
            else if (curr_NN_SE != 0) { 
                if ((curr_NN >= void_threshold)) { E_initial = void_E; }
                else { E_initial = interface_E; } }              
            else {                                                                                                                                                                                                                                                                                                                                     
                if ((curr_NN >= void_threshold)) { E_initial = void_E; }
                else { E_initial = 0;}
            }

            // getting final site energy
            if (new_reg_id != 0) { E_final += regions[(new_reg_id-1)]->e_below_bulk; }
            else if (new_NN_SE != 0) { 
                if ((new_NN >= void_threshold)) { E_final = void_E; }
                else { E_final = interface_E; }
             }   
            else {                                                                                                                                                                                                                                                                                                                                     
                if ((new_NN >= void_threshold)) { E_final = void_E; }
                else { E_final = 0;}
            }


            // getting change in energy of nearest-neighbor sites
            double delta_E_NN = get_E_of_NN_wrapper(coord[0], coord[1], coord[2], coord[3], shift, lattice);

            double barrier = 0;
            if ((new_NN >= void_threshold) && (curr_NN >= void_threshold)) { 
                if ((lattice == 0) || (lattice == 1)) { barrier = terrace_barrier_111; }
                else if ((lattice == 2) || (lattice == 3)) { barrier = terrace_barrier_100; }
            }
            else if ((new_NN_SE) || (curr_NN_SE)) { barrier = interface_barrier; }
            else { 
                if ((lattice == 0) || (lattice == 1)) { barrier = bulk_migration_111; }
                else if ((lattice == 2) || (lattice == 3)) { barrier = bulk_migration_100; }                
            }
            //barrier = bulk_migration_111;
            
            // energy difference 
            double delta_endpoints = E_final - E_initial;
            double neighbor_deltaE = delta_E_NN + delta_endpoints;
                        
            /*if (rank == 0) {
                std::cout << "rank: " << rank << " i: " << coord[0] << " j: " << coord[1] << " k: " << coord[2] << " l: " << coord[3] << "\n";
                std::cout << "rank: " << rank << " i_new: " << i_new << " j_new: " << j_new << " k_new: " << k_new << " l_new: " << l_new << "\n";
                std::cout << "rank: " << rank << " neighbor_deltaE: " << neighbor_deltaE <<  " barrier: " << barrier << "\n";
                std::cout << "rank: " << rank << " curr_NN: " << curr_NN <<  " new_NN: " << new_NN << "\n";
            }*/


            if (neighbor_deltaE >= 0) { rate = 5e12 * std::exp( -(neighbor_deltaE + barrier) * (1 / (8.6173e-5 * temperature)));  }
            else { rate = 5e12 * std::exp( -(barrier) * (1 / (8.6173e-5 * temperature))); }


            //system_energy += E_initial;

            return rate;
        }

        double delta_E_init_to_final(int* coord, int* shift, int lattice, int curr_NN, int new_NN, bool update_lattice_check=false) {
            //std::cout << "entering get_rateconstants()\n";
            
            double E_initial = 0;
            double E_final = 0;
            //std::cout << "curr_NN: " << curr_NN << "new_NN: " << new_NN << "\n";

            int i_new; 

            if ((lattice == 0) || (lattice == 3)) {i_new = 1;}
            else {i_new = 0;}

            int j_new = (((coord[1] + shift[0]) % sublattice_dim[0] + sublattice_dim[0]) % sublattice_dim[0]);
            int k_new = (((coord[2] + shift[1]) % sublattice_dim[1] + sublattice_dim[1]) % sublattice_dim[1]);
            int l_new = (((coord[3] + shift[2]) % sublattice_dim[2] + sublattice_dim[2]) % sublattice_dim[2]);
            
            int reg_id = region_sites(coord[0], coord[1], coord[2], coord[3]);
            int new_reg_id = region_sites(i_new, j_new, k_new, l_new);
            
            int curr_NN_SE = 0;
            int new_NN_SE = 0;

            // determining if in region or solid electrolyte region
            if (reg_id != 0) { }
            else if (coord[3] == (sublattice_dim[2]-1)) { curr_NN_SE = 1; }
            if (new_reg_id != 0) { }
            else if (l_new == (sublattice_dim[2]-1)) { new_NN_SE = 1; }
            
            // getting initial site energy
            if (reg_id != 0) { E_initial = regions[(reg_id-1)]->e_below_bulk; }
            else if (curr_NN_SE != 0) { 
                if ((curr_NN >= void_threshold)) { E_initial = void_E; }
                else { E_initial = interface_E; } }              
            else {
                if ((curr_NN >= void_threshold)) { E_initial = void_E; }
                else { E_initial = 0;}
            }

            // getting final site energy
            if (new_reg_id != 0) { E_final += regions[(new_reg_id-1)]->e_below_bulk; }
            else if (new_NN_SE != 0) { 
                if ((new_NN >= void_threshold)) { E_final = void_E; }
                else { E_final = interface_E; }
             }   
            else {            
                if ((new_NN >= void_threshold)) { E_final = void_E; }
                else { E_final = 0;}
            }
                    
            // getting change in energy of nearest-neighbor sites
            double delta_E_NN = get_E_of_NN_wrapper(coord[0], coord[1], coord[2], coord[3], shift, lattice);            

            // energy difference 
            double delta_endpoints = E_final - E_initial;
            double neighbor_deltaE = delta_E_NN + delta_endpoints;
            
            if ((update_lattice_check)) {
                std::cout << "rank: " << rank << " delta_endpoints: " << delta_endpoints << "\n";
                std::cout << "rank: " << rank << " delta_E_NN: " << delta_E_NN << "\n";
                std::cout << "rank: " << rank << " proc_E_cost: " << neighbor_deltaE << "\n";
            }
            return neighbor_deltaE;
        }


        /**
        * @brief Method for obtaining the number of vacancies in the nearest neighbor shell of a vacancy.
        *
        * @param vac Pointer to an array representing the vacancy coordinates.
        * @param lattice The lattice type.
        * @return The number of vacancies in the nearest neighbor shell.
        */
        template <typename T>
        int get_NN_count(T vac, int lattice_idx, T exclude_loc_vac, bool exclude_loc_bool) {
            
            int new_i=0; int new_j=0; int new_k=0; int new_l=0;
            int unmod_new_j=0; int unmod_new_k=0; int unmod_new_l=0;
            int NN_count = 0;

            std::vector<size_t> x_dims = proc_pos_x_neighbors.size_vec;
            std::vector<size_t> y_dims = proc_pos_y_neighbors.size_vec;

            // moving vacancy from bc site to vertex site
            for (int i=0; i<(int)diag_directions.rows(); i++) {
                  
                if (lattice_idx == 1) {
                    new_i = 0;
                    new_j = (((vac[1] + diag_directions[i][0]) % sublattice_dim[0] + sublattice_dim[0]) % sublattice_dim[0]);
                    new_k = (((vac[2] + diag_directions[i][1]) % sublattice_dim[1] + sublattice_dim[1]) % sublattice_dim[1]);
                    new_l = (((vac[3] + diag_directions[i][2]) % sublattice_dim[2] + sublattice_dim[2]) % sublattice_dim[2]);
                    
                    unmod_new_j = (vac[1] + diag_directions[i][0]);
                    unmod_new_k = (vac[2] + diag_directions[i][1]);
                    unmod_new_l = (vac[3] + diag_directions[i][2]);                    
                   
                }

                else if (lattice_idx == 0) {
                
                    new_i = 1;
                    new_j = (((vac[1] - diag_directions[i][0]) % sublattice_dim[0] + sublattice_dim[0]) % sublattice_dim[0]);
                    new_k = (((vac[2] - diag_directions[i][1]) % sublattice_dim[1] + sublattice_dim[1]) % sublattice_dim[1]);
                    new_l = (((vac[3] - diag_directions[i][2]) % sublattice_dim[2] + sublattice_dim[2]) % sublattice_dim[2]);

                    unmod_new_j = (vac[1] - diag_directions[i][0]);
                    unmod_new_k = (vac[2] - diag_directions[i][1]);
                    unmod_new_l = (vac[3] - diag_directions[i][2]);                    
                }

                if ((exclude_loc_bool) && (new_i == exclude_loc_vac[0]) 
                    && (new_j == exclude_loc_vac[1]) 
                    && (new_k == exclude_loc_vac[2]) 
                    && (new_l == exclude_loc_vac[3])) {

                }
            
                else if ((vac[3] == (int)(sublattice_dim[2]-1)) && (diag_directions[i][2] == 1)) {/* checking for rightmost non-periodic boundary along z-axis*/}
                
                else if ((vac[3] == 0) && (diag_directions[i][2] == 1)) {/* checking for leftmost non-periodic boundary along z-axis*/}
                    
                else if ((unmod_new_j < 0)) {/*check with proc to -x direction*/   
                    if (unmod_new_k < 0) {
                        if ((proc_neighbors(rank,5) != rank) && (proc_neg_x_neighbors(new_i, mod_with_bounds(new_j + chunk_bounds[0][0], 2), (size_t)(unmod_new_k+2), (size_t)(new_l)))) {NN_count ++;}
                    }
                    else if (unmod_new_k > (sublattice_dim[1]-1)) {
                        if ((proc_neighbors(rank,3) != rank) && (proc_neg_x_neighbors(new_i, mod_with_bounds(new_j + chunk_bounds[0][0], 2), (size_t)(unmod_new_k+2), (size_t)(new_l)))) {NN_count ++;}
                    }
                    else {
                        /*check neighbor -x array */
                        if ((proc_neighbors(rank,4) != rank) && (proc_neg_x_neighbors(new_i, mod_with_bounds(new_j + chunk_bounds[0][0], 2), (size_t)(new_k+2), (size_t)(new_l)))) {NN_count ++;}
                    }
                }     

                else if ((unmod_new_j > (sublattice_dim[0] - 1))) {/*check with proc to +x direction*/ 
                    /*check neighbor (+x) array */
                    
                    if (unmod_new_k < 0) {
                        if ((proc_neighbors(rank,7) != rank) && (proc_pos_x_neighbors(new_i, mod_with_bounds(new_j + chunk_bounds[0][0], 2), (size_t)(unmod_new_k+2), (size_t)(new_l)))) {NN_count ++;}
                    }
                    else if (unmod_new_k > (sublattice_dim[1]-1)) {
                        if ((proc_neighbors(rank,1) != rank) && (proc_pos_x_neighbors(new_i, mod_with_bounds(new_j + chunk_bounds[0][0], 2), (size_t)(unmod_new_k+2), (size_t)(new_l)))) {NN_count ++;}
                    }
                    else if (((proc_neighbors(rank,0) != rank)) && (proc_pos_x_neighbors(new_i, mod_with_bounds(new_j + chunk_bounds[0][0], 2), (size_t)(new_k+2), (size_t)(new_l)))) {NN_count ++;}

                }

                else if ((unmod_new_k < 0)) { /*check with proc to -y direction*/
                    /*check with proc to -y direction*/
                    if (((proc_neighbors(rank,6) != rank)) && (proc_neg_y_neighbors(new_i, (size_t)(new_j+2), mod_with_bounds(new_k + chunk_bounds[1][0], 2), (size_t)new_l))) {NN_count ++;}
                    
                }

                else if ((unmod_new_k > (sublattice_dim[1] - 1))) {/*check with proc to +y direction*/
                    /*check neighbor +y array */
                    if ((proc_neighbors(rank,2) != rank) && (proc_pos_y_neighbors(new_i, (size_t)(new_j+2), mod_with_bounds(new_k + chunk_bounds[1][0], 2), (size_t)(new_l)))) {NN_count ++;}
                }
                
                else { 
                    NN_count += vacancies(new_i, new_j, new_k, new_l);
                }   

            }

            return NN_count;
        }

        template <typename T>
        int get_NN_count(T vac, int lattice_idx, bool debug=false) {
                            
            int new_i=0; int new_j=0; int new_k=0; int new_l=0;
            int unmod_new_j=0; int unmod_new_k=0; int unmod_new_l=0;
            int NN_count = 0;

            std::vector<size_t> x_dims = proc_pos_x_neighbors.size_vec;
            std::vector<size_t> y_dims = proc_pos_y_neighbors.size_vec;

            if (debug) std::cout << " rank: " << rank << " sublattice_dim: [ " << sublattice_dim[0] << " " << sublattice_dim[1] << " " << sublattice_dim[2] << " " << sublattice_dim[3] << "]\n";
                
            // moving vacancy from bc site to vertex site
            for (int i=0; i<(int)diag_directions.rows(); i++) {
                
                if (lattice_idx == 1) {
                    new_i = 0;
                    new_j = (((vac[1] + diag_directions[i][0]) % sublattice_dim[0] + sublattice_dim[0]) % sublattice_dim[0]);
                    new_k = (((vac[2] + diag_directions[i][1]) % sublattice_dim[1] + sublattice_dim[1]) % sublattice_dim[1]);
                    new_l = (((vac[3] + diag_directions[i][2]) % sublattice_dim[2] + sublattice_dim[2]) % sublattice_dim[2]);
                    
                    unmod_new_j = (vac[1] + diag_directions[i][0]);
                    unmod_new_k = (vac[2] + diag_directions[i][1]);
                    unmod_new_l = (vac[3] + diag_directions[i][2]);                    
                   
                }

                else if (lattice_idx == 0) {
                
                    new_i = 1;
                    new_j = (((vac[1] - diag_directions[i][0]) % sublattice_dim[0] + sublattice_dim[0]) % sublattice_dim[0]);
                    new_k = (((vac[2] - diag_directions[i][1]) % sublattice_dim[1] + sublattice_dim[1]) % sublattice_dim[1]);
                    new_l = (((vac[3] - diag_directions[i][2]) % sublattice_dim[2] + sublattice_dim[2]) % sublattice_dim[2]);

                    unmod_new_j = (vac[1] - diag_directions[i][0]);
                    unmod_new_k = (vac[2] - diag_directions[i][1]);
                    unmod_new_l = (vac[3] - diag_directions[i][2]);                    
                }

                if (debug) std::cout << " rank: " << rank << " new_loc: [ " << new_i << " " << new_j << " " << new_k << " " << new_l << "]\n";
                if (debug) std::cout << " rank: " << rank << " new_loc_unmod: [ " << new_i << " " << unmod_new_j << " " << unmod_new_k << " " << unmod_new_l << "]\n";
                            
                if (unmod_new_l > (int)(sublattice_dim[2]-1)) {/* checking for rightmost non-periodic boundary along z-axis*/}
                
                else if (unmod_new_l < 0) {/* checking for leftmost non-periodic boundary along z-axis*/}
                    
                else if ((unmod_new_j < 0)) {/*check with proc to -x direction*/   
                    if (unmod_new_k < 0) {
                        if ((proc_neighbors(rank,5) != rank) && (proc_neg_x_neighbors(new_i, mod_with_bounds(new_j + chunk_bounds[0][0], 2), (size_t)(unmod_new_k+2), (size_t)(new_l)))) {NN_count ++;}
                    }
                    else if (unmod_new_k > (sublattice_dim[1]-1)) {
                        if ((proc_neighbors(rank,3) != rank) && (proc_neg_x_neighbors(new_i, mod_with_bounds(new_j + chunk_bounds[0][0], 2), (size_t)(unmod_new_k+2), (size_t)(new_l)))) {NN_count ++;}
                    }
                    else {
                        /*check neighbor -x array */
                        if ((proc_neighbors(rank,4) != rank) && (proc_neg_x_neighbors(new_i, mod_with_bounds(new_j + chunk_bounds[0][0], 2), (size_t)(new_k+2), (size_t)(new_l)))) {NN_count ++;}
                    }
                }     

                else if ((unmod_new_j > (sublattice_dim[0] - 1))) {/*check with proc to +x direction*/ 
                    /*check neighbor (+x) array */
                    
                    if (unmod_new_k < 0) {
                        if ((proc_neighbors(rank,7) != rank) && (proc_pos_x_neighbors(new_i, mod_with_bounds(new_j + chunk_bounds[0][0], 2), (size_t)(unmod_new_k+2), (size_t)(new_l)))) {NN_count ++;}
                    }
                    else if (unmod_new_k > (sublattice_dim[1]-1)) {
                        if ((proc_neighbors(rank,1) != rank) && (proc_pos_x_neighbors(new_i, mod_with_bounds(new_j + chunk_bounds[0][0], 2), (size_t)(unmod_new_k+2), (size_t)(new_l)))) {NN_count ++;}
                    }
                    else if (((proc_neighbors(rank,0) != rank)) && (proc_pos_x_neighbors(new_i, mod_with_bounds(new_j + chunk_bounds[0][0], 2), (size_t)(new_k+2), (size_t)(new_l)))) {NN_count ++;}

                }

                else if ((unmod_new_k < 0)) { /*check with proc to -y direction*/
                    /*check with proc to -y direction*/
                    if (((proc_neighbors(rank,6) != rank)) && (proc_neg_y_neighbors(new_i, (size_t)(new_j+2), mod_with_bounds(new_k + chunk_bounds[1][0], 2), (size_t)new_l))) {NN_count ++;}
                    
                }

                else if ((unmod_new_k > (sublattice_dim[1] - 1))) {/*check with proc to +y direction*/
                    /*check neighbor +y array */
                    if ((proc_neighbors(rank,2) != rank) && (proc_pos_y_neighbors(new_i, (size_t)(new_j+2), mod_with_bounds(new_k + chunk_bounds[1][0], 2), (size_t)(new_l)))) {NN_count ++;}
                }
                
                else { 
                    if (debug) std::cout << " rank: " << rank << " in same proc \n";
                    NN_count += vacancies(new_i, new_j, new_k, new_l);
                }    

            }
                     
            return NN_count;
        }

        template <typename T>
        int get_NN_count_2NNshell(T vac, int lattice_idx, T exclude_loc_vac, bool exclude_loc_bool) {
            int count = 0;
            
            int new_i1; int new_i2; int new_i3; int new_i4;
            // moving vacancy from bc site to vertex site
            if (lattice_idx == 1) {
                for (int i=0; i<(int)diag_directions.rows(); i++) {
                    new_i1 = !(vac[0]);
                    new_i2 = (((vac[1] + diag_directions[i][0]) % sublattice_dim[0] + sublattice_dim[0]) % sublattice_dim[0]);
                    new_i3 = (((vac[2] + diag_directions[i][1]) % sublattice_dim[1] + sublattice_dim[1]) % sublattice_dim[1]);
                    new_i4 = (((vac[3] + diag_directions[i][2]) % sublattice_dim[2] + sublattice_dim[2]) % sublattice_dim[2]);
                    
                    //if ((vac[3] == 0) && (diag_directions[i][2] == 1)) {/* checking for leftmost non-periodic boundary along z-axis*/}
                    if ((vac[3] == (int)(sublattice_dim[2]-1)) && (diag_directions[i][2] == 1)) {/* checking for rightmost non-periodic boundary along z-axis*/}
                    else { 
                        if ((exclude_loc_bool) && (new_i1 == exclude_loc_vac[0]) 
                            && (new_i2 == exclude_loc_vac[1]) 
                            && (new_i3 == exclude_loc_vac[2]) 
                            && (new_i4 == exclude_loc_vac[3])) {}
                        else { count += vacancies(0, new_i2, new_i3, new_i4); }
                    }                        
                }
                for (int i=0; i<(int)edge_directions.rows(); i++) {
                    new_i1 = vac[0];
                    new_i2 = (((vac[1] + edge_directions[i][0]) % sublattice_dim[0] + sublattice_dim[0]) % sublattice_dim[0]);
                    new_i3 = (((vac[2] + edge_directions[i][1]) % sublattice_dim[1] + sublattice_dim[1]) % sublattice_dim[1]);
                    new_i4 = (((vac[3] + edge_directions[i][2]) % sublattice_dim[2] + sublattice_dim[2]) % sublattice_dim[2]);
                    
                    if ((vac[3] == 0) && (edge_directions[i][2] == 1)) {/* checking for leftmost non-periodic boundary along z-axis*/}
                    if ((vac[3] == (int)(sublattice_dim[2]-1)) && (edge_directions[i][2] == 1)) {/* checking for rightmost non-periodic boundary along z-axis*/}
                    else { 
                        
                        if ((exclude_loc_bool) && (new_i1 == exclude_loc_vac[0]) 
                            && (new_i2 == exclude_loc_vac[1]) 
                            && (new_i3 == exclude_loc_vac[2]) 
                            && (new_i4 == exclude_loc_vac[3])) {}
                        else { count += vacancies(1, new_i2, new_i3, new_i4); }
                    }
                }
            }

            // moving vacancy from vertex site to bc site     
            else if (lattice_idx == 0) {
                for (int i=0; i<(int)diag_directions.rows(); i++) {
                    new_i1 = !(vac[0]);
                    new_i2 = (((vac[1] - diag_directions[i][0]) % sublattice_dim[0] + sublattice_dim[0]) % sublattice_dim[0]);
                    new_i3 = (((vac[2] - diag_directions[i][1]) % sublattice_dim[1] + sublattice_dim[1]) % sublattice_dim[1]);
                    new_i4 = (((vac[3] - diag_directions[i][2]) % sublattice_dim[2] + sublattice_dim[2]) % sublattice_dim[2]);
                                        
                    if ((vac[3] == 0) && (diag_directions[i][2] == 1)) {/* checking for leftmost non-periodic boundary along z-axis*/}
                    // else if ((vac[3] == (int)(lattice_dim[2]-1)) && (diag_directions[i][2] == 1)) {/* checking for rightmost non-periodic boundary along z-axis*/}
                    else { 
                        if ((exclude_loc_bool) && (new_i1 == exclude_loc_vac[0]) 
                            && (new_i2 == exclude_loc_vac[1]) 
                            && (new_i3 == exclude_loc_vac[2]) 
                            && (new_i4 == exclude_loc_vac[3])) {}
                        else { count += vacancies(new_i1, new_i2, new_i3, new_i4); }
                    }
                
                }
                for (int i=0; i<(int)edge_directions.rows(); i++) {
                    new_i1 = vac[0];
                    new_i2 = (((vac[1] + edge_directions[i][0]) % sublattice_dim[0] + sublattice_dim[0]) % sublattice_dim[0]);
                    new_i3 = (((vac[2] + edge_directions[i][1]) % sublattice_dim[1] + sublattice_dim[1]) % sublattice_dim[1]);
                    new_i4 = (((vac[3] + edge_directions[i][2]) % sublattice_dim[2] + sublattice_dim[2]) % sublattice_dim[2]);
                    
                    if ((vac[3] == 0) && (edge_directions[i][2] == 1)) {/* checking for leftmost non-periodic boundary along z-axis*/}
                    else if ((vac[3] == (int)(sublattice_dim[2]-1)) && (edge_directions[i][2] == 1)) {/* checking for rightmost non-periodic boundary along z-axis*/}
                    else { 
                        if ((exclude_loc_bool) && (new_i1 == exclude_loc_vac[0]) 
                            && (new_i2 == exclude_loc_vac[1]) 
                            && (new_i3 == exclude_loc_vac[2]) 
                            && (new_i4 == exclude_loc_vac[3])) {}
                       else { count += vacancies(0, new_i2, new_i3, new_i4); }
                    }
                }
            }
            
            return count;
        }

        template <typename T>
        int get_NN_count_2NNshell(T vac, int lattice_idx) {

            int count = 0;
            
            // moving vacancy from bc site to vertex site
            if (lattice_idx == 1) {
                for (int i=0; i<(int)diag_directions.rows(); i++) {
                    
                    //if ((vac[3] == 0) && (diag_directions[i][2] == 1)) {/* checking for leftmost non-periodic boundary along z-axis*/}
                    if ((vac[3] == (int)(sublattice_dim[2]-1)) && (diag_directions[i][2] == 1)) {/* checking for rightmost non-periodic boundary along z-axis*/}
                    else {
                        count += vacancies(0, (((vac[1] + diag_directions[i][0]) % sublattice_dim[0] + sublattice_dim[0]) % sublattice_dim[0]), 
                                            (((vac[2] + diag_directions[i][1]) % sublattice_dim[1] + sublattice_dim[1]) % sublattice_dim[1]), 
                                            (((vac[3] + diag_directions[i][2]) % sublattice_dim[2] + sublattice_dim[2]) % sublattice_dim[2])); 
                    }                        
                }
                for (int i=0; i<(int)edge_directions.rows(); i++) {
                    
                    if ((vac[3] == 0) && (edge_directions[i][2] == 1)) {/* checking for leftmost non-periodic boundary along z-axis*/}
                    if ((vac[3] == (int)(sublattice_dim[2]-1)) && (edge_directions[i][2] == 1)) {/* checking for rightmost non-periodic boundary along z-axis*/}
                    else {     
                        count += vacancies(1, (((vac[1] + edge_directions[i][0]) % sublattice_dim[0] + sublattice_dim[0]) % sublattice_dim[0]), 
                                            (((vac[2] + edge_directions[i][1]) % sublattice_dim[1] + sublattice_dim[1]) % sublattice_dim[1]), 
                                            (((vac[3] + edge_directions[i][2])% sublattice_dim[2] + sublattice_dim[2]) % sublattice_dim[2])); 
                    }
                }
                
            }

            // moving vacancy from vertex site to bc site     
            else if (lattice_idx == 0) {
                for (int i=0; i<(int)diag_directions.rows(); i++) {
                                        
                    if ((vac[3] == 0) && (diag_directions[i][2] == 1)) {/* checking for leftmost non-periodic boundary along z-axis*/}
                    // else if ((vac[3] == (int)(lattice_dim[2]-1)) && (diag_directions[i][2] == 1)) {/* checking for rightmost non-periodic boundary along z-axis*/}
                    else { 
                        count += vacancies(1, (((vac[1] - diag_directions[i][0]) % sublattice_dim[0] + sublattice_dim[0]) % sublattice_dim[0]), 
                                            (((vac[2] - diag_directions[i][1]) % sublattice_dim[1] + sublattice_dim[1]) % sublattice_dim[1]), 
                                            (((vac[3] - diag_directions[i][2])% sublattice_dim[2] + sublattice_dim[2]) % sublattice_dim[2]));        
                    }
                
                }
                for (int i=0; i<(int)edge_directions.rows(); i++) {
                    
                    if ((vac[3] == 0) && (edge_directions[i][2] == 1)) {/* checking for leftmost non-periodic boundary along z-axis*/}
                    else if ((vac[3] == (int)(sublattice_dim[2]-1)) && (edge_directions[i][2] == 1)) {/* checking for rightmost non-periodic boundary along z-axis*/}
                    else { 
                        count += vacancies(0, (((vac[1] + edge_directions[i][0]) % sublattice_dim[0] + sublattice_dim[0]) % sublattice_dim[0]), 
                                                (((vac[2] + edge_directions[i][1]) % sublattice_dim[1] + sublattice_dim[1]) % sublattice_dim[1]), 
                                                (((vac[3] + edge_directions[i][2])% sublattice_dim[2] + sublattice_dim[2]) % sublattice_dim[2])); 
                    
                    }
                }
            }

            //std::cout << "count: " << count << "\n";
            return count;
        }


        /**
        * @brief Method for obtaining the number of vacancies in the nearest neighbor shells of a vacancies.
        *
        * @param vac Type T data structure (2-dimensional) containing all vacancy coordinates.
        * @return The number of vacancies in the nearest neighbor shell for each vacancy.
        */
        template <typename T>
        std::vector<int> get_all_NN(T vacs, int shells) {
            int size = vacs.size();
            std::vector<int> vac_NN(size);
            int vac_count;

            if (shells == 1) {
                for ( int i=0; i<(int)vacs.size(); i++ ) {
                    vac_count = get_NN_count(vacs[i], vacs[i][0]); 
                    vac_NN[i] = vac_count;
                }
            }
            else if (shells == 2) {
                for ( int i=0; i<(int)vacs.size(); i++ ) {
                    vac_count = get_NN_count_2NNshell(vacs[i], vacs[i][0]); 
                    vac_NN[i] = vac_count;
                }
            }
            else { 
                std::cout << "ERROR: wrong number of shells in get_all_NN()\n";
                exit(0); 
            }

            return vac_NN;
        }


        /**
        * @brief Method for obtaining the number of adaptive GB sites in the nearest neighbor shell of a vacancy.
        *
        * @param vac Pointer to an array representing the vacancy coordinates.
        * @param lattice The lattice type.
        * @return The number of vacancies in the nearest neighbor shell.
        */
        int get_adaptivesites_NN_count(std::vector<int> vac, int lattice) {
            int count = 0;
            int reg_id = 0;
            
            
            // moving vacancy from bc site to vertex site
            if (lattice == 1) {
                for (int i=0; i<(int)diag_directions.rows(); i++) {
                    if ((vac[3] == 0) && (diag_directions[i][2] == 1)) {/* checking for leftmost non-periodic boundary along z-axis*/}
                    else if ((vac[3] == (int)(sublattice_dim[2]-1)) && (diag_directions[i][2] == 1)) {/* checking for rightmost non-periodic boundary along z-axis*/}
                    else { 
                        reg_id = region_sites(0, 
                            (((vac[1] + diag_directions[i][0]) % sublattice_dim[0] + sublattice_dim[0]) % sublattice_dim[0]), 
                            (((vac[2] + diag_directions[i][1]) % sublattice_dim[1] + sublattice_dim[1]) % sublattice_dim[1]), 
                            (((vac[3] + diag_directions[i][2]) % sublattice_dim[2] + sublattice_dim[2]) % sublattice_dim[2]));
                            
                        if (reg_id == adaptive_gb_id) { count ++;}
                        else if ((reg_id != 0) && (regions[(reg_id-1)]->is_gb)) { count ++;}
                    }
                    
                }
            }

            // moving vacancy from vertex site to bc site     
            else if (lattice == 0) {
                for (int i=0; i<(int)diag_directions.rows(); i++) {
                    if ((vac[3] == 0) && (diag_directions[i][2] == 1)) {/* checking for leftmost non-periodic boundary along z-axis*/}
                    else if ((vac[3] == (int)(sublattice_dim[2]-1)) && (diag_directions[i][2] == 1)) {/* checking for rightmost non-periodic boundary along z-axis*/}
                    else { 
                        reg_id = region_sites(1, 
                            (((vac[1] - diag_directions[i][0]) % sublattice_dim[0] + sublattice_dim[0]) % sublattice_dim[0]), 
                            (((vac[2] - diag_directions[i][1]) % sublattice_dim[1] + sublattice_dim[1]) % sublattice_dim[1]), 
                            (((vac[3] - diag_directions[i][2]) % sublattice_dim[2] + sublattice_dim[2]) % sublattice_dim[2])); 
                            
                        if (reg_id == adaptive_gb_id) { count ++;}
                        else if ((reg_id != 0) && (regions[(reg_id-1)]->is_gb)) { count ++;}
                    }
                }
            }
            std::cout << "adaptive_count: " << count << "\n";
            
            return count;
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
        double new_get_rateconstants(std::vector<int>& coord, int* shift, int lattice, int curr_NN, int new_NN) {  
            
            double rate = -1; 
            int LR_idx;  // index corresponding to direction of movement in lattice (left/right)
            int idx = 0;

            // return -1 if move not allowed
            if (idx == -1) {
                return -1;
            }
 
            int reg_id = region_sites(coord[0], coord[1], coord[2], coord[3]);
            // DETERMINING DIRECTION OF MOVE //


            // moving vacancy from vertex site to bc site
            if (reg_id != 0) { 
                
                if (regions[(reg_id-1)]->bias == "X") {
                    if (lattice == 0) {
                        // moving vacancy from vertex site to bc site
                        if (shift[0] == 0) {LR_idx = 0;}
                        else {LR_idx = 1;}
                    }
                        
                    else if (lattice == 1) {
                        // moving vacancy from bc site to vertex site
                        if (shift[0] == 1) {LR_idx = 0;}
                        else {LR_idx = 1;}
                    }

                    else if ((lattice == 2) || (lattice == 3)) {
                        // moving vacancy from  vertex site to vertex site OR bc site to bc site
                        if (shift[0] == 1) {LR_idx = 0;}
                        else if (shift[0] == -1) {LR_idx = 1;}
                        else LR_idx = -1;
                    }
                }
                else if (regions[(reg_id-1)]->bias == "Y") {
                    if (lattice == 0) {
                        // moving vacancy from vertex site to bc site
                        if (shift[1] == 0) {LR_idx = 0;}
                        else {LR_idx = 1;}
                    }
                        
                    else if (lattice == 1) {
                        // moving vacancy from bc site to vertex site
                        if (shift[1] == 1) {LR_idx = 0;}
                        else {LR_idx = 1;}
                    }

                    else if ((lattice == 2) || (lattice == 3)) {
                        // moving vacancy from  vertex site to vertex site OR bc site to bc site
                        if (shift[1] == 1) {LR_idx = 0;}
                        else if (shift[1] == -1) {LR_idx = 1;}
                        else LR_idx = -1;
                    }
                }
                else if (regions[(reg_id-1)]->bias == "Z") {
                    if (lattice == 0) {
                        // moving vacancy from vertex site to bc site
                        if (shift[2] == 0) {LR_idx = 0;}
                        else {LR_idx = 1;}
                    }
                        
                    else if (lattice == 1) {
                        // moving vacancy from bc site to vertex site
                        if (shift[2] == 1) {LR_idx = 0;}
                        else {LR_idx = 1;}
                    }

                    else if ((lattice == 2) || (lattice == 3)) {
                        // moving vacancy from  vertex site to vertex site OR bc site to bc site
                        if (shift[2] == 1) {LR_idx = 0;}
                        else if (shift[2] == -1) {LR_idx = 1;}
                        else LR_idx = -1;
                    }
                }
                else {
                    std::cout << "ERROR: invalid directional bias" << "\n";
                    exit(0);
                }

                if ((lattice == regions[(reg_id-1)]->interface_i) && (regions[(reg_id-1)]->interface) && 
                    (shift[regions[(reg_id-1)]->interface_dim] != 0)) {
                    
                    rate = regions[(reg_id-1)]->interface_100_rate;
                }
                else if (regions[(reg_id-1)]->random) {
                    if ( regions[(reg_id-1)]->get_rate(coord[0], coord[1], coord[2], coord[3], LR_idx) == -1) { rate = -1; }
                    
                    else if ((curr_NN >= void_threshold)) { //else if ((curr_NN >= void_threshold) && (new_NN < void_threshold)) { 
                        rate = void_gb_diss_barrier; 
                    }  
                    else {
                        if (LR_idx == 1) rate = regions[(reg_id-1)]->get_rate(coord[0], coord[1], coord[2], coord[3], LR_idx);
                        if (LR_idx == 0) rate = regions[(reg_id-1)]->get_rate(coord[0], coord[1], coord[2], coord[3], LR_idx);
                    }
                }
                else {
                    //if ((curr_NN >= void_threshold) && (new_NN < void_threshold)) { 
                    if ((curr_NN >= void_threshold)) {
                        rate = void_gb_diss_barrier; 
                    }  
                    else if ((lattice == 0) || (lattice == 1)) {
                        if (LR_idx == 1) {rate = regionrates_111_L[(reg_id-1)][idx];}
                        else if (LR_idx == 0) {rate = regionrates_111_R[(reg_id-1)][idx];}                        
                    }
                    
                    else if ((lattice == 2) || (lattice == 3)) {
                        if (LR_idx == 1) {rate = regionrates_100_L[(reg_id-1)][idx];}
                        else if (LR_idx == 0) {rate = regionrates_100_R[(reg_id-1)][idx];}
                        else if (LR_idx == -1) {rate = ratecatalog_100[0][idx];}
                    }
                    else {
                        std::string str_output = "Error: invalid lattice type in search_catalog()";
                        printf("%s", str_output.c_str());
                        throw std::exception();
                    }
                }
            }            
            else {
                if (lattice == 0) {
                    // moving vacancy from vertex site to bc site
                    if (shift[2] == 0) {LR_idx = 0;}
                    else {LR_idx = 1;}
                }
                    
                else if (lattice == 1) {
                    // moving vacancy from bc site to vertex site
                    if (shift[2] == 1) {LR_idx = 0;}
                    else {LR_idx = 1;}
                }

                else if ((lattice == 2) || (lattice == 3)) {
                    // moving vacancy from  vertex site to vertex site OR bc site to bc site
                    if (shift[2] == 1) {LR_idx = 0;}
                    else {LR_idx = 1;}
                }

                // if ((curr_NN >= void_threshold) && (new_NN < void_threshold)) { rate = void_barrier; }  
                if ((curr_NN >= void_threshold)) { rate = void_barrier; }  
                else {
                    // in case of no pre-defined region, use bulk rate constants
                    if ((lattice == 1) || (lattice == 0)) {rate = ratecatalog_111[LR_idx][idx];}
                    else if ((lattice == 2) || (lattice == 3)) {rate = ratecatalog_100[0][idx];}
                }
            }
            
            /*
            if ( ((rank == 5) && ((coord[1] == 0) || (coord[2] == 0)))
                || ((rank == 6) && ((coord[1] == (sublattice_dim[0]-1)) || (coord[2] == 0)))
                || ((rank == 9) && ((coord[1] == 0) || (coord[2] == (sublattice_dim[1]-1))))
                || ((rank == 10) && ((coord[1] == (sublattice_dim[0]-1)) || (coord[2] == (sublattice_dim[1]-1))))) {
                
                // std::cout << "rank: " << rank << " [ " << coord[0] << " " << coord[1] << " " << coord[2] << " " << coord[3] << " ]\n";
                // std::cout << "rank: " << rank << " rate: " << rate << " curr_NN: " << curr_NN << " new_NN: " << new_NN << "\n\n";
                // std::cout << " sublattice_dim[0]: " << sublattice_dim[0] << " sublattice_dim[1]: " << sublattice_dim[1] << "\n"; 
            }
            */

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
            //
            assert((vac[0] == 0) || (vac[0] == 1));
            assert( ((vac[1] >= 0) && (vac[1] < sublattice_dim[0])) );
            assert( ((vac[2] >= 0) && (vac[2] < sublattice_dim[1])) );
            assert( ((vac[3] >= 0) && (vac[3] < sublattice_dim[2])) );
            //
            
            int sum = 0;
            int m = (int)a_types.size(); // number of atom types in system

            // moving vacancy from bc site to bc site 
            if (lattice == 3) {
                for (int i=0; i<(int)diag_directions.rows(); i++) {
                    sum += exp_int(m,i) * vertex_sites(0, (((vac[0] - diag_directions[i][0]) % sublattice_dim[0] + sublattice_dim[0]) % sublattice_dim[0]), (((vac[1] - diag_directions[i][1]) % sublattice_dim[1] + sublattice_dim[1]) % sublattice_dim[1]), (((vac[2] - diag_directions[i][2]) % sublattice_dim[2] + sublattice_dim[2]) % sublattice_dim[2]));
                }
                for (int i=0; i<(int)edge_directions.rows(); i++) {
                    sum += exp_int(m,i) *  bc_sites(0, (((vac[0] + edge_directions[i][0]) % sublattice_dim[0] + sublattice_dim[0]) % sublattice_dim[0]), (((vac[1] + edge_directions[i][1]) % sublattice_dim[1] + sublattice_dim[1]) % sublattice_dim[1]), (((vac[2] + edge_directions[i][2]) % sublattice_dim[2] + sublattice_dim[2]) % sublattice_dim[2]));
                }
            }

            // moving vacancy from vertex site to vertex site 
            else if (lattice == 2) {
                for (int i=0; i<(int)diag_directions.rows(); i++) {
                    sum +=  exp_int(m,i) * bc_sites(0, (((vac[0] + diag_directions[i][0])  % sublattice_dim[0] + sublattice_dim[0]) % sublattice_dim[0]), (((vac[1] + diag_directions[i][1]) % sublattice_dim[1] + sublattice_dim[1]) % sublattice_dim[1]), (((vac[2] + diag_directions[i][2]) % sublattice_dim[2] + sublattice_dim[2]) % sublattice_dim[2]));
                }
                for (int i=0; i<(int)edge_directions.rows(); i++) {
                    sum += exp_int(m,i) *  vertex_sites(0, (((vac[0] + edge_directions[i][0]) % sublattice_dim[0] + sublattice_dim[0]) % sublattice_dim[0]), (((vac[1] + edge_directions[i][1]) % sublattice_dim[1] + sublattice_dim[1]) % sublattice_dim[1]), (((vac[2] + edge_directions[i][2]) % sublattice_dim[2] + sublattice_dim[2]) % sublattice_dim[2]));
                }
            }

            // moving vacancy from bc site to vertex site
            else if (lattice == 1) {
                for (int i=0; i<(int)diag_directions.rows(); i++) {
                    sum += exp_int(m,i) * vertex_sites(0, (((vac[0] - diag_directions[i][0]) % sublattice_dim[0] + sublattice_dim[0]) % sublattice_dim[0]), (((vac[1] - diag_directions[i][1]) % sublattice_dim[1] + sublattice_dim[1]) % sublattice_dim[1]), (((vac[2] - diag_directions[i][2]) % sublattice_dim[2] + sublattice_dim[2]) % sublattice_dim[2]));
                }
            }

            // moving vacancy from vertex site to bc site     
            else if (lattice == 0) {
                for (int i=0; i<(int)diag_directions.rows(); i++) {
                    sum += exp_int(m,i) * bc_sites(0, (((vac[0] + diag_directions[i][0]) % sublattice_dim[0] + sublattice_dim[0]) % sublattice_dim[0]), (((vac[1] + diag_directions[i][1]) % sublattice_dim[1] + sublattice_dim[1]) % sublattice_dim[1]), (((vac[2] + diag_directions[i][2])% sublattice_dim[2] + sublattice_dim[2]) % sublattice_dim[2]));
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
        * @param new_loc_unmod The same new coordinates before wraparound, i.e. they may fall outside
        *                [0, sublattice_dim) by one cell. Used for the ghost-array writes, which need to
        *                know which side of the boundary was crossed.
        *
        * @return The rank of the process the vacancy was transferred to if the move crosses the boundary
        *         of the processor domain, otherwise -1.
        *
        * @note The function uses MPI non-blocking communication (MPI_Isend) to send data to neighboring
        *       processors. It handles multiple ghost regions and different data structures for various
        *       lattice configurations.
        */

        int parallel_processes_check(int i_old, int j_old, int k_old, int l_old, int move_idx, const std::vector<int>& new_loc, const std::vector<int>& new_loc_unmod) {
            //
            assert( ((i_old == 0) || (i_old == 1)) && ((new_loc[0] == 0) || (new_loc[0] == 1)) );
            assert( ((j_old >= 0) && (j_old < sublattice_dim[0])) && ((new_loc[1] >= 0) && (new_loc[1] < sublattice_dim[0])) );
            assert( ((k_old >= 0) && (k_old < sublattice_dim[1])) && ((new_loc[2] >= 0) && (new_loc[2] < sublattice_dim[1])) );
            assert( ((l_old >= 0) && (l_old < sublattice_dim[2])) && ((new_loc[3] >= 0) && (new_loc[3] < sublattice_dim[2])) );
            //
            //std::cout << "rank: " << rank  << " parallel i_old: " << i_old << " j_old: " << j_old << " k_old: " << k_old << " l_old: " << l_old << "\n";
            //std::cout << "rank: " << rank   << " parallel new_loc[0]: " << new_loc[0] << " new_loc[1]: " << new_loc[1] << " new_loc[2]: " << new_loc[2] << " new_loc[3]: " << new_loc[3] << "\n";
            
            int new_proc;
            int bufferlen = 9;
            
            std::vector<int> loc_buffer1(bufferlen); // new location of vacancy
            std::vector<int> loc_buffer2(bufferlen); // new location of vacancy
            std::vector<int> loc_buffer3(bufferlen); // new location of vacancy

            MPI_Request request1;
            MPI_Request request2;
            MPI_Request request3;


            int i_new = new_loc[0];
            int j_new = new_loc[1];
            int k_new = new_loc[2];
            int l_new = new_loc[3];

            int i_new_unmod = new_loc_unmod[0];
            int j_new_unmod = new_loc_unmod[1];
            int k_new_unmod = new_loc_unmod[2];
            int l_new_unmod = new_loc_unmod[3];
            int i = i_new; int l = l_new;

            std::vector<size_t> x_dims = proc_pos_x_neighbors.size_vec;
            std::vector<size_t> y_dims = proc_pos_y_neighbors.size_vec;

            if (j_old == 0) {
                if (k_old == 0) {
                    
                    if ((j_new == (sublattice_dim[0] - 1) ) && (k_new == (sublattice_dim[1] - 1) )) {  //proc 5 (-x,-y)
                        
                        new_proc = proc_neighbors(rank,5);
                        for (int idx=0; idx<new_loc.size(); idx++) {loc_buffer1[idx] = new_loc[idx];}
                        loc_buffer1[4] = rank;
                        loc_buffer1[5] = i_old;
                        loc_buffer1[6] = j_old;
                        loc_buffer1[7] = k_old;
                        loc_buffer1[8] = l_old;
                        MPI_Isend(loc_buffer1.data(), bufferlen, MPI_INT, new_proc, 1, MPI_COMM_WORLD, &request1);
                        MPI_Wait(&request1, MPI_STATUS_IGNORE);
                    
                        proc_neg_y_neighbors((size_t)i, (size_t)(j_new_unmod+2), mod_with_bounds(k_new_unmod + chunk_bounds[1][0], 2), l) = 1;
                        proc_neg_x_neighbors((size_t)i, mod_with_bounds(j_new_unmod + chunk_bounds[0][0], 2), (size_t)(k_new_unmod+2), l) = 1;

                        return new_proc;
                    }
                    else if (k_new == (sublattice_dim[1] - 1)) { //proc 6 (-y)
                        new_proc = proc_neighbors(rank,6); // proc_neighbors indices follow: 0:+x, 1:(+x,+y), 2:+y, 3:(-x,+y), 4:-x, 5:(-x,-y), 6:-y, 7:(+x,-y) 

                        for (int idx=0; idx<new_loc.size(); idx++) {loc_buffer1[idx] = new_loc[idx];}
                        loc_buffer1[4] = rank;
                        loc_buffer1[5] = i_old;
                        loc_buffer1[6] = j_old;
                        loc_buffer1[7] = k_old;
                        loc_buffer1[8] = l_old; 
                        
                        MPI_Isend(loc_buffer1.data(), bufferlen, MPI_INT, new_proc, 1, MPI_COMM_WORLD, &request1);

                        MPI_Wait(&request1, MPI_STATUS_IGNORE);
                            
                        proc_neg_y_neighbors((size_t)i, (size_t)(j_new_unmod+2), mod_with_bounds(k_new_unmod + chunk_bounds[1][0], 2), l) = 1;

                        return new_proc;
                    }
                    else if (j_new == (sublattice_dim[0] - 1)) { //proc 4 (-x)
                        new_proc = proc_neighbors(rank,4); // proc_neighbors indices follow: 0:+x, 1:(+x,+y), 2:+y, 3:(-x,+y), 4:-x, 5:(-x,-y), 6:-y, 7:(+x,-y)
                        
                        for (int idx=0; idx<new_loc.size(); idx++) {loc_buffer1[idx] = new_loc[idx];}
                        loc_buffer1[4] = rank;
                        loc_buffer1[5] = i_old;
                        loc_buffer1[6] = j_old;
                        loc_buffer1[7] = k_old;
                        loc_buffer1[8] = l_old; 
                        MPI_Isend(loc_buffer1.data(), bufferlen, MPI_INT, new_proc, 1, MPI_COMM_WORLD, &request1);                             
                        MPI_Wait(&request1, MPI_STATUS_IGNORE);
                    
                        proc_neg_x_neighbors((size_t)i, mod_with_bounds(j_new_unmod + chunk_bounds[0][0], 2), (size_t)(k_new_unmod+2), l) = 1;

                        return new_proc;                        
                    }
                }
                else if ((k_old > 0) && (k_old < (sublattice_dim[1] - 1))) {
                    
                    if (j_new == (sublattice_dim[0] - 1)) { //proc 4 (-x)
                        new_proc = proc_neighbors(rank,4); // proc_neighbors indices follow: 0:+x, 1:(+x,+y), 2:+y, 3:(-x,+y), 4:-x, 5:(-x,-y), 6:-y, 7:(+x,-y)
                        
                        for (int idx=0; idx<new_loc.size(); idx++) {loc_buffer1[idx] = new_loc[idx];}
                        loc_buffer1[4] = rank;
                        loc_buffer1[5] = i_old;
                        loc_buffer1[6] = j_old;
                        loc_buffer1[7] = k_old;
                        loc_buffer1[8] = l_old; 
                        MPI_Isend(loc_buffer1.data(), bufferlen, MPI_INT, new_proc, 1, MPI_COMM_WORLD, &request1);
                        MPI_Wait(&request1, MPI_STATUS_IGNORE);

                        // Investigating a 10-tick ghost mismatch traced to this exact write: it's a
                        // raw assignment (=1), not the increment every other write path uses for
                        // these reference-counted arrays -- if the cell already holds a nonzero
                        // count from an unrelated contribution, this silently discards it instead
                        // of adding to it. Read-only: prints, does not change, the pre-write value.
                        {
                            size_t dbg_xi = mod_with_bounds(j_new_unmod + chunk_bounds[0][0], 2);
                            size_t dbg_yi = (size_t)(k_new_unmod+2);
                            if ((i==1) && (dbg_xi==0) && (dbg_yi==7) && (l==119)) {
                                std::cout << "PPC_RAWSET rank: " << rank << " arr: neg_x w: " << i
                                          << " x_idx: " << dbg_xi << " y_idx: " << dbg_yi << " z: " << l
                                          << " PRE_WRITE_VALUE: " << proc_neg_x_neighbors((size_t)i, dbg_xi, dbg_yi, l)
                                          << " (about to be hard-set to 1 here)\n";
                            }
                        }
                        proc_neg_x_neighbors((size_t)i, mod_with_bounds(j_new_unmod + chunk_bounds[0][0], 2), (size_t)(k_new_unmod+2), l) = 1;

                        return new_proc;
                    }
                }
                else if (k_old == (sublattice_dim[1] - 1)) {
                    
                    if (j_new == (sublattice_dim[0] - 1)) {  //proc 4 (-x)
                        new_proc = proc_neighbors(rank,4); // proc_neighbors indices follow: 0:+x, 1:(+x,+y), 2:+y, 3:(-x,+y), 4:-x, 5:(-x,-y), 6:-y, 7:(+x,-y)
                        
                        for (int idx=0; idx<new_loc.size(); idx++) {loc_buffer1[idx] = new_loc[idx];}
                        loc_buffer1[4] = rank;
                        loc_buffer1[5] = i_old;
                        loc_buffer1[6] = j_old;
                        loc_buffer1[7] = k_old;
                        loc_buffer1[8] = l_old; 
                        MPI_Isend(loc_buffer1.data(), bufferlen, MPI_INT, new_proc, 1, MPI_COMM_WORLD, &request1);                             
                        MPI_Wait(&request1, MPI_STATUS_IGNORE);
                    
                        proc_neg_x_neighbors((size_t)i, mod_with_bounds(j_new_unmod + chunk_bounds[0][0], 2), (size_t)(k_new_unmod+2), l) = 1;

                        return new_proc;
                    }
                    else if (k_new == 0) { // proc 2 (+y)
                        new_proc = proc_neighbors(rank,2); // proc_neighbors indices follow: 0:+x, 1:(+x,+y), 2:+y, 3:(-x,+y), 4:-x, 5:(-x,-y), 6:-y, 7:(+x,-y) 

                        for (int idx=0; idx<new_loc.size(); idx++) {loc_buffer1[idx] = new_loc[idx];}
                        loc_buffer1[4] = rank;
                        loc_buffer1[5] = i_old;
                        loc_buffer1[6] = j_old;
                        loc_buffer1[7] = k_old;
                        loc_buffer1[8] = l_old;
                        
                        MPI_Isend(loc_buffer1.data(), bufferlen, MPI_INT, new_proc, 1, MPI_COMM_WORLD, &request1);

                        MPI_Wait(&request1, MPI_STATUS_IGNORE);
                        
                        proc_pos_y_neighbors((size_t)i, (size_t)(j_new_unmod+2), mod_with_bounds(k_new_unmod + chunk_bounds[1][0], 2), l) = 1;

                        return new_proc;
                    }
                }
            }
            else if (j_old == sublattice_dim[0] - 1) {
                if (k_old == 0) {
                    if (k_new == (sublattice_dim[1] - 1)) { //proc 6 (-y)
                        new_proc = proc_neighbors(rank,6); // proc_neighbors indices follow: 0:+x, 1:(+x,+y), 2:+y, 3:(-x,+y), 4:-x, 5:(-x,-y), 6:-y, 7:(+x,-y) 

                        for (int idx=0; idx<new_loc.size(); idx++) {loc_buffer1[idx] = new_loc[idx];}
                        loc_buffer1[4] = rank;
                        loc_buffer1[5] = i_old;
                        loc_buffer1[6] = j_old;
                        loc_buffer1[7] = k_old;
                        loc_buffer1[8] = l_old; 
                        
                        MPI_Isend(loc_buffer1.data(), bufferlen, MPI_INT, new_proc, 1, MPI_COMM_WORLD, &request1);

                        MPI_Wait(&request1, MPI_STATUS_IGNORE);
                            
                        proc_neg_y_neighbors((size_t)i, (size_t)(j_new_unmod+2), mod_with_bounds(k_new_unmod + chunk_bounds[1][0], 2), l) = 1;

                        return new_proc;
                    }
                    else if (j_new == 0) { //proc 0 (+x)
                        new_proc = proc_neighbors(rank,0); // proc_neighbors indices follow: 0:+x, 1:(+x,+y), 2:+y, 3:(-x,+y), 4:-x, 5:(-x,-y), 6:-y, 7:(+x,-y) 
                        
                        for (int idx=0; idx<new_loc.size(); idx++) {loc_buffer1[idx] = new_loc[idx];}
                        loc_buffer1[4] = rank;
                        loc_buffer1[5] = i_old;
                        loc_buffer1[6] = j_old;
                        loc_buffer1[7] = k_old;
                        loc_buffer1[8] = l_old;
                        MPI_Isend(loc_buffer1.data(), bufferlen, MPI_INT, new_proc, 1, MPI_COMM_WORLD, &request1);

                        MPI_Wait(&request1, MPI_STATUS_IGNORE);
                        
                        proc_pos_x_neighbors((size_t)i, mod_with_bounds(j_new_unmod + chunk_bounds[0][0], 2), (size_t)(k_new_unmod+2), l) = 1;

                        return new_proc;
                    }
                }
                else if ((k_old > 0) && (k_old < (sublattice_dim[1] - 1))) {
                    if (j_new == 0) {//proc 0 (+x)
                        new_proc = proc_neighbors(rank,0); // proc_neighbors indices follow: 0:+x, 1:(+x,+y), 2:+y, 3:(-x,+y), 4:-x, 5:(-x,-y), 6:-y, 7:(+x,-y) 
                        
                        for (int idx=0; idx<new_loc.size(); idx++) {loc_buffer1[idx] = new_loc[idx];}
                        loc_buffer1[4] = rank;
                        loc_buffer1[5] = i_old;
                        loc_buffer1[6] = j_old;
                        loc_buffer1[7] = k_old;
                        loc_buffer1[8] = l_old;
                        MPI_Isend(loc_buffer1.data(), bufferlen, MPI_INT, new_proc, 1, MPI_COMM_WORLD, &request1);

                        MPI_Wait(&request1, MPI_STATUS_IGNORE);
                        
                        proc_pos_x_neighbors((size_t)i, mod_with_bounds(j_new_unmod + chunk_bounds[0][0], 2), (size_t)(k_new_unmod+2), l) = 1;

                        return new_proc;
                    }
                }
                else if (k_old == (sublattice_dim[1] - 1)) {
                    
                    if ((j_new == 0) && (k_new == 0)) {  //proc 1 (+x,+y)
                        
                        new_proc = proc_neighbors(rank,1);
                        for (int idx=0; idx<new_loc.size(); idx++) {loc_buffer1[idx] = new_loc[idx];}
                        loc_buffer1[4] = rank;
                        loc_buffer1[5] = i_old;
                        loc_buffer1[6] = j_old;
                        loc_buffer1[7] = k_old;
                        loc_buffer1[8] = l_old;
                        MPI_Isend(loc_buffer1.data(), bufferlen, MPI_INT, new_proc, 1, MPI_COMM_WORLD, &request1);

                        MPI_Wait(&request1, MPI_STATUS_IGNORE);                            
                    
                        proc_pos_x_neighbors((size_t)i, mod_with_bounds(j_new_unmod + chunk_bounds[0][0], 2), (size_t)(k_new_unmod+2), l) = 1;
                        proc_pos_y_neighbors((size_t)i, (size_t)(j_new_unmod+2), mod_with_bounds(k_new_unmod + chunk_bounds[1][0], 2), l) = 1;
                        
                        return new_proc;
                    }
                    if (j_new == 0) {  //proc 0 (+x)
                        new_proc = proc_neighbors(rank,0); // proc_neighbors indices follow: 0:+x, 1:(+x,+y), 2:+y, 3:(-x,+y), 4:-x, 5:(-x,-y), 6:-y, 7:(+x,-y) 
                        
                        for (int idx=0; idx<new_loc.size(); idx++) {loc_buffer1[idx] = new_loc[idx];}
                        loc_buffer1[4] = rank;
                        loc_buffer1[5] = i_old;
                        loc_buffer1[6] = j_old;
                        loc_buffer1[7] = k_old;
                        loc_buffer1[8] = l_old;
                        MPI_Isend(loc_buffer1.data(), bufferlen, MPI_INT, new_proc, 1, MPI_COMM_WORLD, &request1);

                        MPI_Wait(&request1, MPI_STATUS_IGNORE);
                        
                        proc_pos_x_neighbors((size_t)i, mod_with_bounds(j_new_unmod + chunk_bounds[0][0], 2), (size_t)(k_new_unmod+2), l) = 1;

                        return new_proc;
                    }
                    else if (k_new == 0) { // proc 2 (+y)
                        new_proc = proc_neighbors(rank,2); // proc_neighbors indices follow: 0:+x, 1:(+x,+y), 2:+y, 3:(-x,+y), 4:-x, 5:(-x,-y), 6:-y, 7:(+x,-y) 

                        for (int idx=0; idx<new_loc.size(); idx++) {loc_buffer1[idx] = new_loc[idx];}
                        loc_buffer1[4] = rank;
                        loc_buffer1[5] = i_old;
                        loc_buffer1[6] = j_old;
                        loc_buffer1[7] = k_old;
                        loc_buffer1[8] = l_old;
                        
                        MPI_Isend(loc_buffer1.data(), bufferlen, MPI_INT, new_proc, 1, MPI_COMM_WORLD, &request1);

                        MPI_Wait(&request1, MPI_STATUS_IGNORE);
                        
                        proc_pos_y_neighbors((size_t)i, (size_t)(j_new_unmod+2), mod_with_bounds(k_new_unmod + chunk_bounds[1][0], 2), l) = 1;

                        return new_proc;
                    }
                }
            }
            else if ((k_old == 0) && (k_new == (sublattice_dim[1]-1))) { //proc 6 (-y)
                new_proc = proc_neighbors(rank,6); // proc_neighbors indices follow: 0:+x, 1:(+x,+y), 2:+y, 3:(-x,+y), 4:-x, 5:(-x,-y), 6:-y, 7:(+x,-y) 

                for (int idx=0; idx<new_loc.size(); idx++) {loc_buffer1[idx] = new_loc[idx];}
                loc_buffer1[4] = rank;
                loc_buffer1[5] = i_old;
                loc_buffer1[6] = j_old;
                loc_buffer1[7] = k_old;
                loc_buffer1[8] = l_old; 
                
                MPI_Isend(loc_buffer1.data(), bufferlen, MPI_INT, new_proc, 1, MPI_COMM_WORLD, &request1);

                MPI_Wait(&request1, MPI_STATUS_IGNORE);
                    
                proc_neg_y_neighbors((size_t)i, (size_t)(j_new_unmod+2), mod_with_bounds(k_new_unmod + chunk_bounds[1][0], 2), l) = 1;
                
                return new_proc;
            }
            else if ((k_old == (sublattice_dim[1] - 1)) && (k_new == 0)) { // proc 2 (+y)
                new_proc = proc_neighbors(rank,2); // proc_neighbors indices follow: 0:+x, 1:(+x,+y), 2:+y, 3:(-x,+y), 4:-x, 5:(-x,-y), 6:-y, 7:(+x,-y) 

                for (int idx=0; idx<new_loc.size(); idx++) {loc_buffer1[idx] = new_loc[idx];}
                loc_buffer1[4] = rank;
                loc_buffer1[5] = i_old;
                loc_buffer1[6] = j_old;
                loc_buffer1[7] = k_old;
                loc_buffer1[8] = l_old;
                
                MPI_Isend(loc_buffer1.data(), bufferlen, MPI_INT, new_proc, 1, MPI_COMM_WORLD, &request1);

                MPI_Wait(&request1, MPI_STATUS_IGNORE);
                
                proc_pos_y_neighbors((size_t)i, (size_t)(j_new_unmod+2), mod_with_bounds(k_new_unmod + chunk_bounds[1][0], 2), l) = 1;

                return new_proc;
            }

            //std::cout << "rank: " << rank  << "se vuelve nada \n";

            return -1;
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
        * @param new_loc Vector containing the new locations of the ghost sites.
        * @param parallel_transfer Boolean flag indicating if the transfer is parallel.
        */
        /**
        * @brief Direction-dispatch core of ghost_site_send, parameterized on which position
        * (old or new) drives the boundary-proximity branch selection.
        *
        * ghost_site_send used to always branch on (j_old, k_old) alone, on the assumption that
        * a move's old and new positions are always close enough that only the old position's
        * proximity to a boundary needs checking. That assumption breaks for longer-range move
        * types (this lattice has 14 distinct move types per vacancy, not just nearest-neighbor
        * hops): a move can have an old position that's comfortably interior while its new
        * position lands right in/near the ghost margin -- in that case branching on (j_old,k_old)
        * alone never notifies any neighbor about the new position at all. ghost_site_send now
        * calls this dispatch once for (j_old,k_old) and once for the new position's local
        * coordinates, so either position's boundary proximity is enough to trigger the right
        * neighbor notification(s). The transmitted buffer content (old position + new_loc_unmod)
        * is unchanged either way -- only which branch(es) fire depends on (j_pos, k_pos).
        */
        // already_sent: tracks which target ranks ghost_site_send has already messaged for THIS
        // move (shared across both of its dispatch calls -- see ghost_site_send). ghost_site_send
        // dispatches once keyed off the old position and once off the new position so that either
        // endpoint's boundary proximity is enough to trigger notification; when old and new are
        // both near the same edge/corner, both calls resolve to the same target rank(s), and
        // without this guard that rank received (and applied) the identical update twice. That
        // was harmless back when the live ghost arrays were plain booleans (setting/clearing twice
        // is idempotent) but corrupts the reference counts they are now (see the ghost-array
        // members' comment) by double-incrementing/decrementing. send_ghost_msg below is a thin
        // wrapper around the existing per-direction send blocks that consults/updates this set
        // instead of sending unconditionally; which branch fires for which (j_pos,k_pos), and
        // which proc_neighbors slot each targets, is unchanged.
        void ghost_site_send_dispatch(int j_pos, int k_pos, int i_old, int j_old, int k_old, int l_old,
                                       const std::vector<int>& new_loc_unmod, int tag,
                                       std::set<int>& already_sent) {
            int bufferlen = 10; // 9 location/old-proc fields + 1 explicit direction index (loc_buffer[9])

            std::vector<int> loc_buffer1(bufferlen); // new location of vacancy
            std::vector<int> loc_buffer2(bufferlen); // new location of vacancy
            std::vector<int> loc_buffer3(bufferlen); // new location of vacancy

            MPI_Request request1;
            MPI_Request request2;
            MPI_Request request3;

            auto send_ghost_msg = [&](int target_proc, int side_idx, std::vector<int>& loc_buffer, MPI_Request& req) {
                if ((target_proc == rank) || already_sent.count(target_proc)) { return; }
                already_sent.insert(target_proc);
                for (int idx=0; idx<(int)new_loc_unmod.size(); idx++) {loc_buffer[idx] = new_loc_unmod[idx];}
                loc_buffer[4] = rank;
                loc_buffer[5] = i_old;
                loc_buffer[6] = j_old;
                loc_buffer[7] = k_old;
                loc_buffer[8] = l_old;
                loc_buffer[9] = side_idx; // side_idx: direction (from sender's frame) this message was sent in
                MPI_Isend(loc_buffer.data(), bufferlen, MPI_INT, target_proc, tag, MPI_COMM_WORLD, &req);
                MPI_Wait(&req, MPI_STATUS_IGNORE);
            };

            // proc_neighbors indices follow: 0:+x, 1:(+x,+y), 2:+y, 3:(-x,+y), 4:-x, 5:(-x,-y), 6:-y, 7:(+x,-y)
            if (j_pos < 2) {/*communicate with proc to -x direction*/
                if (k_pos < 2) {
                    /*communicate with proc to (-x,-y) direction*/
                    if ((proc_neighbors(rank,6) == proc_neighbors(rank,5)) && (proc_neighbors(rank,4) == proc_neighbors(rank,5))) {}
                    else if ((proc_neighbors(rank,6) == rank) && (proc_neighbors(rank,4) != rank)) {
                        send_ghost_msg(proc_neighbors(rank,4), 4, loc_buffer1, request1);
                    }
                    else if ((proc_neighbors(rank,4) == rank) && (proc_neighbors(rank,6) != rank)) {
                        send_ghost_msg(proc_neighbors(rank,6), 6, loc_buffer1, request1);
                    }
                    else if ((proc_neighbors(rank,6) != rank) && (proc_neighbors(rank,4) != rank)) {
                        send_ghost_msg(proc_neighbors(rank,6), 6, loc_buffer1, request1);
                        send_ghost_msg(proc_neighbors(rank,5), 5, loc_buffer2, request2);
                        send_ghost_msg(proc_neighbors(rank,4), 4, loc_buffer3, request3);
                    }
                }
                else if (k_pos > (sublattice_dim[1] - 3)) {
                    /*communicate with proc to (-x,+y) direction -- was previously missing entirely:
                    this whole branch fell into the "pure -x" else below, which only ever notifies
                    proc_neighbors(rank,4) and never the +y or diagonal (-x,+y) neighbors, silently
                    dropping their notification for any move near this corner.*/
                    if ((proc_neighbors(rank,4) == proc_neighbors(rank,3)) && (proc_neighbors(rank,2) == proc_neighbors(rank,3))) {}
                    else if ((proc_neighbors(rank,2) == rank) && (proc_neighbors(rank,4) != rank)) {
                        send_ghost_msg(proc_neighbors(rank,4), 4, loc_buffer1, request1);
                    }
                    else if ((proc_neighbors(rank,4) == rank) && (proc_neighbors(rank,2) != rank)) {
                        send_ghost_msg(proc_neighbors(rank,2), 2, loc_buffer1, request1);
                    }
                    else if ((proc_neighbors(rank,2) != rank) && (proc_neighbors(rank,4) != rank)) {
                        send_ghost_msg(proc_neighbors(rank,2), 2, loc_buffer1, request1);
                        send_ghost_msg(proc_neighbors(rank,3), 3, loc_buffer2, request2);
                        send_ghost_msg(proc_neighbors(rank,4), 4, loc_buffer3, request3);
                    }
                }
                else {
                    /*communicate with proc to (-x) direction*/
                    send_ghost_msg(proc_neighbors(rank,4), 4, loc_buffer1, request1);
                }
            }
            else if (j_pos > (sublattice_dim[0] - 3)) {/*communicate with proc to +x direction*/
                if (k_pos > (sublattice_dim[1] - 3)) {
                    /*communicate with proc to (+x,+y) direction*/
                    if ((proc_neighbors(rank,2) == proc_neighbors(rank,1)) && (proc_neighbors(rank,2) == proc_neighbors(rank,0))) {}
                    else if ((proc_neighbors(rank,0) == rank) && (proc_neighbors(rank,2) != rank)) {
                        send_ghost_msg(proc_neighbors(rank,2), 2, loc_buffer1, request1);
                    }
                    else if ((proc_neighbors(rank,2) == rank) && (proc_neighbors(rank,0) != rank)) {
                        send_ghost_msg(proc_neighbors(rank,0), 0, loc_buffer1, request1);
                    }
                    else if ((proc_neighbors(rank,2) != rank) &&  (proc_neighbors(rank,0) != rank)) {
                        send_ghost_msg(proc_neighbors(rank,2), 2, loc_buffer1, request1);
                        send_ghost_msg(proc_neighbors(rank,1), 1, loc_buffer2, request2);
                        send_ghost_msg(proc_neighbors(rank,0), 0, loc_buffer3, request3);
                    }
                }
                else if (k_pos < 2) {
                    /*communicate with proc to (+x,-y) direction -- was previously missing entirely:
                    this is our flagship case (a move at j_pos near +x and k_pos near -y simultaneously
                    fell into the "pure +x" else below, which never notifies proc_neighbors(rank,6),
                    the true -y neighbor -- exactly the missing-clear bug traced on rank2's tick-342
                    move at local (53,0,121).*/
                    if ((proc_neighbors(rank,0) == proc_neighbors(rank,7)) && (proc_neighbors(rank,6) == proc_neighbors(rank,7))) {}
                    else if ((proc_neighbors(rank,6) == rank) && (proc_neighbors(rank,0) != rank)) {
                        send_ghost_msg(proc_neighbors(rank,0), 0, loc_buffer1, request1);
                    }
                    else if ((proc_neighbors(rank,0) == rank) && (proc_neighbors(rank,6) != rank)) {
                        send_ghost_msg(proc_neighbors(rank,6), 6, loc_buffer1, request1);
                    }
                    else if ((proc_neighbors(rank,0) != rank) && (proc_neighbors(rank,6) != rank)) {
                        send_ghost_msg(proc_neighbors(rank,0), 0, loc_buffer1, request1);
                        send_ghost_msg(proc_neighbors(rank,7), 7, loc_buffer2, request2);
                        send_ghost_msg(proc_neighbors(rank,6), 6, loc_buffer3, request3);
                    }
                }
                else {
                    /*communicate with proc to (+x) direction*/
                    send_ghost_msg(proc_neighbors(rank,0), 0, loc_buffer1, request1);
                }
            }
            else if ( (k_pos < 2) ) {/*communicate with proc to -y direction*/
                send_ghost_msg(proc_neighbors(rank,6), 6, loc_buffer1, request1);
            }
            else if ( (k_pos > (sublattice_dim[1] - 3)) ) {/*communicate with proc to +y direction*/
                send_ghost_msg(proc_neighbors(rank,2), 2, loc_buffer1, request1);
            }
        }

        void ghost_site_send(int i_old, int j_old, int k_old, int l_old, const std::vector<int>& new_loc_unmod, const std::vector<int>& shift, int parallel_transfer, bool reverse=false) {
            int tag;
            if (reverse) {
                if (parallel_transfer != -1) tag = 4; // with parallel transfer
                else tag = 5; // no parallel transfer
            }
            else {
                if (parallel_transfer != -1) tag = 2;
                else tag = 3;
            }

            std::cout << "rank: " << rank << " GSS_TRACE tag: " << tag << " reverse: " << reverse << " parallel_transfer: " << parallel_transfer
                      << " old(i,j,k,l): (" << i_old << "," << j_old << "," << k_old << "," << l_old << ")"
                      << " new_loc_unmod[0..3]: (" << new_loc_unmod[0] << "," << new_loc_unmod[1] << "," << new_loc_unmod[2] << "," << new_loc_unmod[3] << ")\n";

            // Shared across both dispatch calls below so a target rank that both the old and new
            // position resolve to only gets notified once -- see ghost_site_send_dispatch's
            // already_sent comment for why a duplicate here isn't just wasted bandwidth: it now
            // corrupts the receiving rank's ghost reference counts.
            std::set<int> already_sent;
            ghost_site_send_dispatch(j_old, k_old, i_old, j_old, k_old, l_old, new_loc_unmod, tag, already_sent);

            // Also dispatch on the new position's own local coordinates, so a move whose new
            // position lands near/in the ghost margin gets its neighbor(s) notified even when the
            // old position was comfortably interior (see ghost_site_send_dispatch's docstring).
            if ((new_loc_unmod[1] != j_old) || (new_loc_unmod[2] != k_old)) {
                ghost_site_send_dispatch(new_loc_unmod[1], new_loc_unmod[2], i_old, j_old, k_old, l_old, new_loc_unmod, tag, already_sent);
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
        // is_reverse: diagnostic only -- distinguishes tag 2/3 (forward moves) from tag 4/5
        // (reverse_move()/reverse_move_parallel() undoing a move) for GHOST_WRITE trace logs.
        // NOTE: the parallel_transfer parameter here does NOT indicate reverse vs. forward --
        // it's just "did this message carry a parallel transfer" (true for tag 2 AND tag 4, false
        // for tag 3 AND tag 5); is_reverse is what actually distinguishes them.
        void ghost_site_recieve(const std::vector<int>& new_loc_buffer, bool parallel_transfer, bool is_reverse) {
            //
            // Bounds on new_loc_buffer[5..8] ("old position") mirror the +/-1 slack already
            // granted to new_loc_buffer[0..3] ("new position") rather than requiring strict
            // in-range values: reverse_move() legitimately puts the raw (pre-wrap) arrival
            // coordinate of the move being undone into this slot, so it can carry the same
            // one-cell overshoot a genuine forward-move new-position value can.
            assert( (((new_loc_buffer)[5] == 0) || ((new_loc_buffer)[5] == 1)) && (((new_loc_buffer)[0] == 0) || ((new_loc_buffer)[0] == 1)) );
            assert( (((new_loc_buffer)[6] >= -1) && ((new_loc_buffer)[6] < (sublattice_dim[0] + 1))) && (((new_loc_buffer)[1] >= -1) && ((new_loc_buffer)[1] < (sublattice_dim[0] + 1))) );
            assert( (((new_loc_buffer)[7] >= -1) && ((new_loc_buffer)[7] < (sublattice_dim[1] + 1))) && (((new_loc_buffer)[2] >= -1) && ((new_loc_buffer)[2] < (sublattice_dim[1] + 1))) );
            assert( (((new_loc_buffer)[8] >= -1) && ((new_loc_buffer)[8] < (sublattice_dim[2] + 1))) && (((new_loc_buffer)[3] >= -1) && ((new_loc_buffer)[3] < (sublattice_dim[2] + 1))) );
            //
            int old_proc = (new_loc_buffer)[4];
            int i_old = (new_loc_buffer)[5];
            int j_old = (new_loc_buffer)[6];
            int k_old = (new_loc_buffer)[7];
            int l_old = (new_loc_buffer)[8];

            if (rank == 0) {
                std::cout << "GSR rank: " << rank << " old_proc: " << old_proc << " parallel_transfer: " << parallel_transfer
                          << " old(i,j,k,l): (" << i_old << "," << j_old << "," << k_old << "," << l_old << ")"
                          << " new_loc_buffer[0..3]: (" << new_loc_buffer[0] << "," << new_loc_buffer[1] << "," << new_loc_buffer[2] << "," << new_loc_buffer[3] << ")\n";
            }

            /*
            make cases for which proc neighbors to update depending on the location of the
            old site
            */
            int i_new = (new_loc_buffer)[0];
            int j_new = (((new_loc_buffer)[1] % sublattice_dim[0] + sublattice_dim[0]) % sublattice_dim[0]);
            int k_new = (((new_loc_buffer)[2] % sublattice_dim[1] + sublattice_dim[1]) % sublattice_dim[1]);
            int l_new = (((new_loc_buffer)[3] % sublattice_dim[2] + sublattice_dim[2]) % sublattice_dim[2]);


            int i_new_unmod = (new_loc_buffer)[0];
            int j_new_unmod = (new_loc_buffer)[1];
            int k_new_unmod = (new_loc_buffer)[2];
            int l_new_unmod = (new_loc_buffer)[3];

            int i = i_new;
            int j = j_new;
            int k = k_new;
            int l = l_new;

            size_t ipos = i_old; size_t ineg = !i_old;

            std::vector<size_t> x_dims = proc_pos_x_neighbors.size_vec;
            std::vector<size_t> y_dims = proc_pos_y_neighbors.size_vec;

            // Global (whole-simulation) chunk origin of old_proc, derived the same way the driver
            // derives every rank's chunk_bounds from its rank index -- lets us convert old_proc's
            // own local (unpadded) coordinates into global coordinates, then hand them to
            // set_ghost_position (the single canonical rule for where a global coordinate belongs
            // in THIS rank's ghost arrays, already used by populate_lattice/ghost_site_self_reference).
            // This replaces a hand-rolled direction table that indexed using old_proc's LOCAL
            // coordinate parity (j_old%2/k_old%2) -- which only agrees with the canonical GLOBAL-
            // coordinate parity when old_proc's chunk offset happens to be even, and silently
            // disagreed (missed/misplaced ghost bits) whenever it was odd.
            //
            // j_old/k_old/j_new_unmod/k_new_unmod are old_proc's own raw vacancy-position
            // coordinates (0 maps directly to old_proc's chunk_bounds[dim][0], same convention
            // ghost_site_self_reference uses via "+ chunk_bounds[dim][0]") -- NOT indices into the
            // +2-padded ghost-tracking arrays, so no "-2" belongs here; that padding only applies
            // once set_ghost_position itself maps a global coordinate into those arrays.
            int sender_x_idx_grid = old_proc % proc_dims[0];
            int sender_y_idx_grid = old_proc / proc_dims[0];
            int sender_x_start = (total_dims[0] / proc_dims[0]) * sender_x_idx_grid;
            int sender_y_start = (total_dims[1] / proc_dims[1]) * sender_y_idx_grid;

            int global_x_old = sender_x_start + j_old;
            int global_y_old = sender_y_start + k_old;
            int global_z_old = l_old;

            // The tag2/tag3 (parallel_transfer true/false) distinction only affects how many
            // neighboring ranks ghost_site_send notifies (corner cases get multiple messages) --
            // it doesn't change what THIS rank needs to do with the new position it was sent, so
            // both cases place it the same way via set_ghost_position. (Previously the "true"
            // branch used a separate, stale, 1-cell-margin-era hand-rolled block that treated
            // old_proc's local coordinates as if they were already this rank's own local frame.)
            int global_x_new = sender_x_start + j_new_unmod;
            int global_y_new = sender_y_start + k_new_unmod;
            int global_z_new = l_new_unmod;
            // A "new" position identical to the "old" one is this message's way of saying
            // "clear-only, no corresponding arrival" (see reverse_move_parallel()'s use of this):
            // that rank is deleting a vacancy it received via parallel transfer, and the position
            // it's undone back to belongs to a DIFFERENT, remote rank -- there's no single frame
            // (sender_x_start/sender_y_start are old_proc's, i.e. the rank sending THIS message)
            // that both i1..i4 (old_proc's own local coords) and that remote position could
            // correctly share, and that remote rank's own reverse_move() already restores its own
            // ghost tracking for it independently. Skipping the set here avoids computing a
            // meaningless "arrival" site from coordinates that were never in old_proc's frame.
            bool clear_only = (i == i_old) && (j_new_unmod == j_old) && (k_new_unmod == k_old) && (l_new_unmod == l_old);
            if (!clear_only) {
                set_ghost_position(i, global_x_new, global_y_new, global_z_new, 1,
                                    proc_neg_x_neighbors, proc_pos_x_neighbors,
                                    proc_neg_y_neighbors, proc_pos_y_neighbors,
                                    std::string(is_reverse ? "GSR_rev_set_new(" : "GSR_fwd_set_new(") + "parallel_transfer=" + std::to_string((int)parallel_transfer) + ")", old_proc);
            }

            // issues with removing corner vs edge sites because more than one ghost corner site maps to same lattice site
            /* removing site */
            set_ghost_position(i_old, global_x_old, global_y_old, global_z_old, 0,
                                proc_neg_x_neighbors, proc_pos_x_neighbors,
                                proc_neg_y_neighbors, proc_pos_y_neighbors,
                                std::string(is_reverse ? "GSR_rev_clear_old(" : "GSR_fwd_clear_old(") + "parallel_transfer=" + std::to_string((int)parallel_transfer) + ")", old_proc);
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
        // NOTE: new_loc must be the wrapped, in-range local sublattice coordinate (as produced by
        // moves_coords / the modulo-wrapped old_loc/i_old,j_old,k_old,l_old convention) -- NOT the
        // raw "_unmod" value, which can be -1 or sublattice_dim as a boundary-crossing sentinel.
        // set_ghost_position requires GLOBAL coordinates, so local coordinates are converted here
        // via chunk_bounds[dim][0] before being passed in.
        // source_tag: diagnostic-only label distinguishing which caller triggered this self-update
        // (forward move vs. the two reversal paths) -- forwarded into the GHOST_WRITE trace logs.
        void ghost_site_self_reference(int i_old,int j_old,int k_old,int l_old, const std::vector<int>& new_loc, int parallel_transfer, const std::string& source_tag = "GSSR") {

            if (rank == 0) {
                std::cout << "GSSR rank: " << rank << " parallel_transfer: " << parallel_transfer
                          << " old(i,j,k,l): (" << i_old << "," << j_old << "," << k_old << "," << l_old << ")"
                          << " new_loc[0..3]: (" << new_loc[0] << "," << new_loc[1] << "," << new_loc[2] << "," << new_loc[3] << ")\n";
            }

            // Clear this process's own ghost-tracking entry for the vacated (old) position.
            // Uses direct geometry (see set_ghost_position) rather than a rank/direction lookup,
            // so it can't be thrown off by a neighboring rank occupying multiple direction slots.
            set_ghost_position(i_old, j_old + chunk_bounds[0][0], k_old + chunk_bounds[1][0], l_old + chunk_bounds[2][0], 0,
                                proc_neg_x_neighbors, proc_pos_x_neighbors,
                                proc_neg_y_neighbors, proc_pos_y_neighbors,
                                source_tag + "_clear_old", -1);

            // Only mark the new position as a ghost site here if it's still on this process
            // (parallel_transfer == -1); a crossing move's new position belongs to a different
            // process, whose own ghost tracking is updated on the receiving side instead.
            if (parallel_transfer == -1) {
                set_ghost_position(new_loc[0], new_loc[1] + chunk_bounds[0][0], new_loc[2] + chunk_bounds[1][0], new_loc[3] + chunk_bounds[2][0], 1,
                                    proc_neg_x_neighbors, proc_pos_x_neighbors,
                                    proc_neg_y_neighbors, proc_pos_y_neighbors,
                                    source_tag + "_set_new", -1);
            }
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
            ghost_trace_tick = move_ticks; // diagnostic only, see GHOST_WRITE tracing
            int parallel_transfer = -1;

            if (moves_lattice(idx,0) == 5) { 
                //std::cout << "rank: " << rank <<  "moves_lattice(idx,0): " << moves_lattice(idx,0) << "\n";
                
                prev_move_type.push_back(0);
                prev_move_type_ticks.push_back(move_ticks);
                store_move_info(idx, parallel_transfer, -1, 5);
                }
            else {
                std::vector<int> new_loc(4);
                std::vector<int> new_loc_unmod(4);
                for (int i=0; i< (int)new_loc.size(); i++) { new_loc[i] = moves_coords(idx,i); /* new location of vacancy */ } 
                for (int i=0; i< (int)new_loc.size(); i++) { new_loc_unmod[i] = moves_coords_unmod(idx,i); /* new location of vacancy */ } 
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

                int new_i = 0;
                double energy_cost;

                if (*moves_lattice[idx] == 0) { new_i = 1; i_old = 0; }
                else if (*moves_lattice[idx] == 1)  { new_i = 0; i_old = 1; }
                else if (*moves_lattice[idx] == 2)  { new_i = 0; i_old = 0; }
                else if (*moves_lattice[idx] == 3)  { new_i = 1; i_old = 1; }
                else if (*moves_lattice[idx] == 4)  { new_i = new_loc[0]; i_old = new_loc[0]; }

                last_newloc = {new_loc[0], new_loc[1], new_loc[2], new_loc[3]};
                last_oldloc = {i_old, old_loc[0], old_loc[1], old_loc[2]};

                std::vector<int> temp_oldloc = {i_old, old_loc[0], old_loc[1], old_loc[2]};
                std::vector<int> temp_newloc = {new_i, new_loc[1], new_loc[2], new_loc[3]};
            
                int old_loc_arr[4]; 
                old_loc_arr[0] = i_old; old_loc_arr[1] = old_loc[0]; old_loc_arr[2] = old_loc[1]; old_loc_arr[3] = old_loc[2];
                
                std::vector<int> new_loc_arr = {new_i, new_loc[1], new_loc[2], new_loc[3]};        
                
                last_currNN = get_NN_count(temp_oldloc, i_old);
                last_newNN = get_NN_count(temp_newloc, new_i, temp_oldloc, true);

                energy_cost = delta_E_init_to_final(old_loc_arr, moves_shifts[idx], moves_lattice[idx][0], last_currNN, last_newNN, true);
                total_cost += energy_cost;

                parallel_transfer = parallel_processes_check(i_old,j_old,k_old,l_old,idx,new_loc,new_loc_unmod);
                
                std::cout << "rank: " << rank <<  " i_old: " << i_old << " j_old: " << old_loc[0] << " k_old: " << old_loc[1] << " l_old: " << old_loc[2] << "\n";
                std::cout << "rank: " << rank <<  " i_new: " << new_i << " j_new: " << new_loc[1] << " k_new: " << new_loc[2] << " l_new: " << new_loc[3] << "\n";
                std::cout << "rank: " << rank << " parallel_transfer: " << parallel_transfer << "\n";
                std::cout << "rank: " << rank << " last_curr_NN: " << last_currNN << "\n";
                std::cout << "rank: " << rank << " last_new_NN: " << last_newNN << "\n"; 
                
                if (rank==0) std::cout << "rank: " << rank << " energy_cost: " << energy_cost << "\n";
                
                if (parallel_transfer != -1) {
                    // switching occupancy for old and new site in lattice array // 

                    if (moves_lattice(idx,0) == 3) {
                        if (vacancies((size_t)1, (size_t)old_loc[0], (size_t)old_loc[1], (size_t)old_loc[2]) != 1) {
                            std::cout << "rank: " << rank << " move_ticks: " << move_ticks << "ERROR: par_proc send identifying wrong site with vacancies_pos or moves_coords or moves_vacs in parallel move etc \n";
                            std::cout << "i_old: 1" << " j_old: " << old_loc[0] << " k_old: " << old_loc[1] << " l_old: " << old_loc[2] << "\n";
                            std::cout << "vacancies_pos  i: " << vacancies_pos(vacs_idx, 0) << " j: " << vacancies_pos(vacs_idx, 1) << " k: " << vacancies_pos(vacs_idx, 2) << " l: " << vacancies_pos(vacs_idx, 3) << "\n";
                            
                            Matrix<int> only_vacancies = vacancies.nonzero(rank); // configuration of vacancies at current timestep
                            std::cout << "rank: " << rank << " vac_nonzero: \n";
                            only_vacancies.print();
                            std::cout << "rank: " << rank << " vacancies_pos: \n";
                            vacancies_pos.print();

                            Matrix<int> unequal_elems_mat1 = comparison(only_vacancies, vacancies_pos);
                            Matrix<int> unequal_elems_mat2 = comparison(vacancies_pos, only_vacancies);

                            std::cout << "rank: " << rank << " unequal_elems_mat1: \n";
                            unequal_elems_mat1.print();
                            std::cout << "rank: " << rank << " unequal_elems_mat2: \n";
                            unequal_elems_mat2.print();
                            
                            exit(0);
                        }
                        bc_sites((size_t)0, (size_t)old_loc[0], (size_t)old_loc[1], (size_t)old_loc[2]) = 1;
                        vacancies((size_t)1, (size_t)old_loc[0], (size_t)old_loc[1], (size_t)old_loc[2]) = 0;
                    }
                    else if (moves_lattice(idx,0) == 2) {
                        if (vacancies((size_t)0, (size_t)old_loc[0], (size_t)old_loc[1], (size_t)old_loc[2]) != 1) {
                            std::cout << "rank: " << rank << " move_ticks: " << move_ticks << "ERROR: par_proc identifying wrong site with vacancies_pos or moves_coords or moves_vacs parallel move etc \n";
                            std::cout << "i_old: 0" << " j_old: " << old_loc[0] << " k_old: " << old_loc[1] << " l_old: " << old_loc[2] << "\n";
                            std::cout << "vacancies_pos  i: " << vacancies_pos(vacs_idx, 0) << " j: " << vacancies_pos(vacs_idx, 1) << " k: " << vacancies_pos(vacs_idx, 2) << " l: " << vacancies_pos(vacs_idx, 3) << "\n";
                            
                            Matrix<int> only_vacancies = vacancies.nonzero(rank); // configuration of vacancies at current timestep
                            std::cout << "rank: " << rank << " vac_nonzero: \n";
                            only_vacancies.print();
                            std::cout << "rank: " << rank << " vacancies_pos: \n";
                            vacancies_pos.print();

                            Matrix<int> unequal_elems_mat1 = comparison(only_vacancies, vacancies_pos);
                            Matrix<int> unequal_elems_mat2 = comparison(vacancies_pos, only_vacancies);

                            std::cout << "rank: " << rank << " unequal_elems_mat1: \n";
                            unequal_elems_mat1.print();
                            std::cout << "rank: " << rank << " unequal_elems_mat2: \n";
                            unequal_elems_mat2.print();
                            
                            exit(0);
                        }
                        vertex_sites((size_t)0, (size_t)old_loc[0], (size_t)old_loc[1], (size_t)old_loc[2]) = 1;
                        vacancies((size_t)0, (size_t)old_loc[0], (size_t)old_loc[1], (size_t)old_loc[2]) = 0;
                    }
                    else if (moves_lattice(idx,0) == 1) {
                        if (vacancies((size_t)1, (size_t)old_loc[0], (size_t)old_loc[1], (size_t)old_loc[2]) != 1) {
                            std::cout << "rank: " << rank << " move_ticks: " << move_ticks << "ERROR: par_proc identifying wrong site with vacancies_pos or moves_coords or moves_vacs parallel move etc \n";
                            std::cout << "i_old: 1" << " j_old: " << old_loc[0] << " k_old: " << old_loc[1] << " l_old: " << old_loc[2] << "\n";
                            std::cout << "vacancies_pos  i: " << vacancies_pos(vacs_idx, 0) << " j: " << vacancies_pos(vacs_idx, 1) << " k: " << vacancies_pos(vacs_idx, 2) << " l: " << vacancies_pos(vacs_idx, 3) << "\n";
                            
                            Matrix<int> only_vacancies = vacancies.nonzero(rank); // configuration of vacancies at current timestep
                            std::cout << "rank: " << rank << " vac_nonzero: \n";
                            only_vacancies.print();
                            std::cout << "rank: " << rank << " vacancies_pos: \n";
                            vacancies_pos.print();

                            Matrix<int> unequal_elems_mat1 = comparison(only_vacancies, vacancies_pos);
                            Matrix<int> unequal_elems_mat2 = comparison(vacancies_pos, only_vacancies);

                            std::cout << "rank: " << rank << " unequal_elems_mat1: \n";
                            unequal_elems_mat1.print();
                            std::cout << "rank: " << rank << " unequal_elems_mat2: \n";
                            unequal_elems_mat2.print();

                            exit(0);
                        }
                        bc_sites((size_t)0, (size_t)old_loc[0], (size_t)old_loc[1], (size_t)old_loc[2]) = 1;
                        vacancies((size_t)1, (size_t)old_loc[0], (size_t)old_loc[1], (size_t)old_loc[2]) = 0;
                    }
                    else if (moves_lattice(idx,0) == 0) {
                        if (vacancies((size_t)0, (size_t)old_loc[0], (size_t)old_loc[1], (size_t)old_loc[2]) != 1) {
                            std::cout << "rank: " << rank << " move_ticks: " << move_ticks << "ERROR: par_proc identifying wrong site with vacancies_pos or moves_coords or moves_vacs parallel move etc \n";
                            std::cout << "i_old: 0" << " j_old: " << old_loc[0] << " k_old: " << old_loc[1] << " l_old: " << old_loc[2] << "\n";
                            std::cout << "vacancies_pos  i: " << vacancies_pos(vacs_idx, 0) << " j: " << vacancies_pos(vacs_idx, 1) << " k: " << vacancies_pos(vacs_idx, 2) << " l: " << vacancies_pos(vacs_idx, 3) << "\n";
                                
                            Matrix<int> only_vacancies = vacancies.nonzero(rank); // configuration of vacancies at current timestep
                            std::cout << "rank: " << rank << " vac_nonzero: \n";
                            only_vacancies.print();
                            std::cout << "rank: " << rank << " vacancies_pos: \n";
                            vacancies_pos.print();

                            Matrix<int> unequal_elems_mat1 = comparison(only_vacancies, vacancies_pos);
                            Matrix<int> unequal_elems_mat2 = comparison(vacancies_pos, only_vacancies);

                            std::cout << "rank: " << rank << " unequal_elems_mat1: \n";
                            unequal_elems_mat1.print();
                            std::cout << "rank: " << rank << " unequal_elems_mat2: \n";
                            unequal_elems_mat2.print();

                            exit(0);
                        }
                        vertex_sites((size_t)0, (size_t)old_loc[0], (size_t)old_loc[1], (size_t)old_loc[2]) = 1;
                        vacancies((size_t)0, (size_t)old_loc[0], (size_t)old_loc[1], (size_t)old_loc[2]) = 0;
                    }
                    
                    if ((rank == 2) && (move_ticks <= 25)) {
                        std::cout << "rank: " << rank << " PRE_STORE_MOVE_INFO_READBACK move_ticks: " << move_ticks
                                  << " proc_neg_y_neighbors(1,49,1,124): " << proc_neg_y_neighbors(1,49,1,124) << "\n";
                    }
                    prev_move_type.push_back(1);
                    store_move_info(idx, parallel_transfer, vacs_idx, lattice, shift, new_loc, {i_old, j_old, k_old, l_old}, new_loc_unmod);
                    if ((rank == 2) && (move_ticks <= 25)) {
                        std::cout << "rank: " << rank << " POST_STORE_MOVE_INFO_READBACK move_ticks: " << move_ticks
                                  << " proc_neg_y_neighbors(1,49,1,124): " << proc_neg_y_neighbors(1,49,1,124) << "\n";
                    }

                    vacancies_pos.remove_row(vacs_idx, rank);
                    moves_vacs.remove_row(idx, rank);
                    moves_lattice.remove_row(idx, rank);
                    moves_shifts.remove_row(idx, rank);
                    moves_coords.remove_row(idx, rank);
                    moves_coords_unmod.remove_row(idx, rank);
                    if ((rank == 2) && (move_ticks <= 25)) {
                        std::cout << "rank: " << rank << " POST_REMOVE_ROWS_READBACK move_ticks: " << move_ticks
                                  << " proc_neg_y_neighbors(1,49,1,124): " << proc_neg_y_neighbors(1,49,1,124) << "\n";
                    }

                    num_of_vacs --;

                }
                else {
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
                        
                        if ( (vacancies((size_t)1, (size_t)new_loc[1], (size_t)new_loc[2], (size_t)new_loc[3]) != 0 ) || 
                        (bc_sites((size_t)0, (size_t)old_loc[0], (size_t)old_loc[1], (size_t)old_loc[2]) != 0 ) ||
                        (bc_sites((size_t)0, (size_t)new_loc[1], (size_t)new_loc[2], (size_t)new_loc[3]) != 1 ) ) {
                            std::cout << "rank: " << rank << " move_ticks: " << move_ticks << " ERROR: identifying wrong site with vacancies_pos or moves_coords or moves_vacs non-parallel move etc \n";

                            std::cout << "moves_lattice(idx,0): " << moves_lattice(idx,0) << "\n";
                            std::cout << "vacancies((size_t)1, (size_t)old_loc[1], (size_t)old_loc[2], (size_t)old_loc[3]): " << vacancies((size_t)0, (size_t)old_loc[1], (size_t)old_loc[2], (size_t)old_loc[3]) << "\n";
                            std::cout << "vacancies((size_t)1, (size_t)new_loc[1], (size_t)new_loc[2], (size_t)new_loc[3]): " << vacancies((size_t)0, (size_t)new_loc[1], (size_t)new_loc[2], (size_t)new_loc[3]) << "\n";
                            std::cout << "bc_sites((size_t)0, (size_t)old_loc[0], (size_t)old_loc[1], (size_t)old_loc[2]): " << bc_sites((size_t)0, (size_t)old_loc[0], (size_t)old_loc[1], (size_t)old_loc[2]) << "\n";
                            std::cout << "bc_sites((size_t)0, (size_t)new_loc[1], (size_t)new_loc[2], (size_t)new_loc[3]): " << bc_sites((size_t)0, (size_t)new_loc[1], (size_t)new_loc[2], (size_t)new_loc[3]) << "\n";
                                                        
                            std::cout << "i_old: 1" << " j_old: " << old_loc[0] << " k_old: " << old_loc[1] << " l_old: " << old_loc[2] << "\n";
                            std::cout << "i_new: 1" << " j_new: " << new_loc[1] << " k_new: " << new_loc[2] << " l_new: " << new_loc[3] << "\n";
                            std::cout << "vacancies_pos  i: " << vacancies_pos(vacs_idx, 0) << " j: " << vacancies_pos(vacs_idx, 1) << " k: " << vacancies_pos(vacs_idx, 2) << " l: " << vacancies_pos(vacs_idx, 3) << "\n";
                            
                            Matrix<int> only_vacancies = vacancies.nonzero(rank); // configuration of vacancies at current timestep
                            std::cout << "rank: " << rank << " vac_nonzero: \n";
                            only_vacancies.print();
                            std::cout << "rank: " << rank << " vacancies_pos: \n";
                            vacancies_pos.print();
                            
                            Matrix<int> unequal_elems_mat1 = comparison(only_vacancies, vacancies_pos);
                            Matrix<int> unequal_elems_mat2 = comparison(vacancies_pos, only_vacancies);

                            std::cout << "rank: " << rank << " unequal_elems_mat1: \n";
                            unequal_elems_mat1.print();
                            std::cout << "rank: " << rank << " unequal_elems_mat2: \n";
                            unequal_elems_mat2.print();

                            exit(0);
                        }
                        
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

                        if ((vacancies((size_t)0, (size_t)new_loc[1], (size_t)new_loc[2], (size_t)new_loc[3]) != 0 ) || 
                        (vertex_sites((size_t)0, (size_t)old_loc[0], (size_t)old_loc[1], (size_t)old_loc[2]) != 0 ) ||
                        (vertex_sites((size_t)0, (size_t)new_loc[1], (size_t)new_loc[2], (size_t)new_loc[3]) != 1 )) {
                            std::cout << "rank: " << rank << " move_ticks: " << move_ticks << " ERROR: identifying wrong site with vacancies_pos or moves_coords or moves_vacs non-parallel move etc \n";
                            
                            std::cout << "moves_lattice(idx,0): " << moves_lattice(idx,0) << "\n";
                            std::cout << "vacancies((size_t)0, (size_t)old_loc[1], (size_t)old_loc[2], (size_t)old_loc[3]): " << vacancies((size_t)0, (size_t)old_loc[1], (size_t)old_loc[2], (size_t)old_loc[3]) << "\n";
                            std::cout << "vacancies((size_t)0, (size_t)new_loc[1], (size_t)new_loc[2], (size_t)new_loc[3]): " << vacancies((size_t)0, (size_t)new_loc[1], (size_t)new_loc[2], (size_t)new_loc[3]) << "\n";
                            std::cout << "vertex_sites((size_t)0, (size_t)old_loc[0], (size_t)old_loc[1], (size_t)old_loc[2]): " << vertex_sites((size_t)0, (size_t)old_loc[0], (size_t)old_loc[1], (size_t)old_loc[2]) << "\n";
                            std::cout << "vertex_sites((size_t)0, (size_t)new_loc[1], (size_t)new_loc[2], (size_t)new_loc[3]): " << vertex_sites((size_t)0, (size_t)new_loc[1], (size_t)new_loc[2], (size_t)new_loc[3]) << "\n";
                            
                            std::cout << "i_old: 0" << " j_old: " << old_loc[0] << " k_old: " << old_loc[1] << " l_old: " << old_loc[2] << "\n";
                            std::cout << "i_new: 0" << " j_new: " << new_loc[1] << " k_new: " << new_loc[2] << " l_new: " << new_loc[3] << "\n";
                            std::cout << "vacancies_pos  i: " << vacancies_pos(vacs_idx, 0) << " j: " << vacancies_pos(vacs_idx, 1) << " k: " << vacancies_pos(vacs_idx, 2) << " l: " << vacancies_pos(vacs_idx, 3) << "\n";
                            
                            Matrix<int> only_vacancies = vacancies.nonzero(rank); // configuration of vacancies at current timestep
                            std::cout << "rank: " << rank << " vac_nonzero: \n";
                            only_vacancies.print();
                            std::cout << "rank: " << rank << " vacancies_pos: \n";
                            vacancies_pos.print();

                            Matrix<int> unequal_elems_mat1 = comparison(only_vacancies, vacancies_pos);
                            Matrix<int> unequal_elems_mat2 = comparison(vacancies_pos, only_vacancies);

                            std::cout << "rank: " << rank << " unequal_elems_mat1: \n";
                            unequal_elems_mat1.print();
                            std::cout << "rank: " << rank << " unequal_elems_mat2: \n";
                            unequal_elems_mat2.print();

                            exit(0);
                        }

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
                        if ((vacancies((size_t)0, (size_t)new_loc[1], (size_t)new_loc[2], (size_t)new_loc[3]) != 0 ) || 
                        (bc_sites((size_t)0, (size_t)old_loc[0], (size_t)old_loc[1], (size_t)old_loc[2]) != 0 ) ||
                        (vertex_sites((size_t)0, (size_t)new_loc[1], (size_t)new_loc[2], (size_t)new_loc[3]) != 1 )) {
                            std::cout << "rank: " << rank << " move_ticks: " << move_ticks << " ERROR: identifying wrong site with vacancies_pos or moves_coords or moves_vacs non-parallel move etc \n";

                            std::cout << "moves_lattice(idx,0): " << moves_lattice(idx,0) << "\n";
                            std::cout << "vacancies((size_t)1, (size_t)old_loc[0], (size_t)old_loc[1], (size_t)old_loc[2]): " << vacancies((size_t)1, (size_t)old_loc[0], (size_t)old_loc[1], (size_t)old_loc[2]) << "\n";
                            std::cout << "vacancies((size_t)0, (size_t)new_loc[1], (size_t)new_loc[2], (size_t)new_loc[3]): " << vacancies((size_t)0, (size_t)new_loc[1], (size_t)new_loc[2], (size_t)new_loc[3]) << "\n";
                            std::cout << "bc_sites((size_t)0, (size_t)old_loc[0], (size_t)old_loc[1], (size_t)old_loc[2]): " << bc_sites((size_t)0, (size_t)old_loc[0], (size_t)old_loc[1], (size_t)old_loc[2]) << "\n";
                            std::cout << "vertex_sites((size_t)0, (size_t)new_loc[1], (size_t)new_loc[2], (size_t)new_loc[3]): " << vertex_sites((size_t)0, (size_t)new_loc[1], (size_t)new_loc[2], (size_t)new_loc[3]) << "\n";
                            
                            std::cout << "i_old: 1" << " j_old: " << old_loc[0] << " k_old: " << old_loc[1] << " l_old: " << old_loc[2] << "\n";
                            std::cout << "i_new: 0" << " j_new: " << new_loc[1] << " k_new: " << new_loc[2] << " l_new: " << new_loc[3] << "\n";
                            std::cout << "vacancies_pos  i: " << vacancies_pos(vacs_idx, 0) << " j: " << vacancies_pos(vacs_idx, 1) << " k: " << vacancies_pos(vacs_idx, 2) << " l: " << vacancies_pos(vacs_idx, 3) << "\n";

                            Matrix<int> only_vacancies = vacancies.nonzero(rank); // configuration of vacancies at current timestep
                            std::cout << "rank: " << rank << " vac_nonzero: \n";
                            only_vacancies.print();
                            std::cout << "rank: " << rank << " vacancies_pos: \n";
                            vacancies_pos.print();

                            Matrix<int> unequal_elems_mat1 = comparison(only_vacancies, vacancies_pos);
                            Matrix<int> unequal_elems_mat2 = comparison(vacancies_pos, only_vacancies);

                            std::cout << "rank: " << rank << " unequal_elems_mat1: \n";
                            unequal_elems_mat1.print();
                            std::cout << "rank: " << rank << " unequal_elems_mat2: \n";
                            unequal_elems_mat2.print();

                            exit(0);
                        }
                        
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
                        if ((vacancies((size_t)1, (size_t)new_loc[1], (size_t)new_loc[2], (size_t)new_loc[3]) != 0 ) || 
                        (vertex_sites((size_t)0, (size_t)old_loc[0], (size_t)old_loc[1], (size_t)old_loc[2]) != 0 ) ||
                        (bc_sites((size_t)0, (size_t)new_loc[1], (size_t)new_loc[2], (size_t)new_loc[3]) != 1 )) {
                            std::cout << "rank: " << rank << " move_ticks: " << move_ticks << " ERROR: identifying wrong site with vacancies_pos or moves_coords or moves_vacs non-parallel move etc \n";

                            std::cout << "moves_lattice(idx,0): " << moves_lattice(idx,0) << "\n";
                            std::cout << "vacancies((size_t)0, (size_t)old_loc[0], (size_t)old_loc[1], (size_t)old_loc[2]): " << vacancies((size_t)0, (size_t)old_loc[0], (size_t)old_loc[1], (size_t)old_loc[2]) << "\n";
                            std::cout << "vacancies((size_t)1, (size_t)new_loc[1], (size_t)new_loc[2], (size_t)new_loc[3]): " << vacancies((size_t)1, (size_t)new_loc[1], (size_t)new_loc[2], (size_t)new_loc[3]) << "\n";
                            std::cout << "vertex_sites((size_t)0, (size_t)old_loc[0], (size_t)old_loc[1], (size_t)old_loc[2]): " << vertex_sites((size_t)0, (size_t)old_loc[0], (size_t)old_loc[1], (size_t)old_loc[2]) << "\n";
                            std::cout << "bc_sites((size_t)0, (size_t)new_loc[1], (size_t)new_loc[2], (size_t)new_loc[3]): " << bc_sites((size_t)0, (size_t)new_loc[1], (size_t)new_loc[2], (size_t)new_loc[3]) << "\n";

                            std::cout << "i_old: 1" << " j_old: " << old_loc[0] << " k_old: " << old_loc[1] << " l_old: " << old_loc[2] << "\n";
                            std::cout << "i_old: 0" << " j_old: " << old_loc[0] << " k_old: " << old_loc[1] << " l_old: " << old_loc[2] << "\n";
                            std::cout << "i_new: 1" << " j_new: " << new_loc[1] << " k_new: " << new_loc[2] << " l_new: " << new_loc[3] << "\n";
                            std::cout << "vacancies_pos  i: " << vacancies_pos(vacs_idx, 0) << " j: " << vacancies_pos(vacs_idx, 1) << " k: " << vacancies_pos(vacs_idx, 2) << " l: " << vacancies_pos(vacs_idx, 3) << "\n";
                            
                            Matrix<int> only_vacancies = vacancies.nonzero(rank); // configuration of vacancies at current timestep
                            std::cout << "rank: " << rank << " vac_nonzero: \n";
                            only_vacancies.print();
                            std::cout << "rank: " << rank << " vacancies_pos: \n";
                            vacancies_pos.print();
                            
                            Matrix<int> unequal_elems_mat1 = comparison(only_vacancies, vacancies_pos);
                            Matrix<int> unequal_elems_mat2 = comparison(vacancies_pos, only_vacancies);

                            std::cout << "rank: " << rank << " unequal_elems_mat1: \n";
                            unequal_elems_mat1.print();
                            std::cout << "rank: " << rank << " unequal_elems_mat2: \n";
                            unequal_elems_mat2.print();

                            exit(0);
                        }

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
                    store_move_info(idx, parallel_transfer, vacs_idx, lattice, shift, new_loc, {i_old, j_old, k_old, l_old}, new_loc_unmod);
                
                } 
            
                prev_move_type_ticks.push_back(move_ticks);
                
                if ( ((new_loc[1] < 2) || (new_loc[1] > (sublattice_dim[0] - 3))) || ((new_loc[2] < 2) || (new_loc[2] > (sublattice_dim[1] - 3))) ||
                ((old_loc[0] < 2) || (old_loc[0] > (sublattice_dim[0] - 3))) || ((old_loc[1] < 2) || (old_loc[1] > (sublattice_dim[1] - 3))) )
                { /*checking to see if neighbor ghost sites need to be updated */
                    ghost_site_send(i_old,j_old,k_old,l_old,new_loc_unmod,shift,parallel_transfer);
                    ghost_site_self_reference(i_old,j_old,k_old,l_old,new_loc,parallel_transfer,"GSSR_forward");
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
            double time = 0;
            if (rank == 0) {
                unsigned int random = mt_obj();
                double random_double = ((random / (1.+ UINT32_MAX)) +  (1 / (1.+ UINT32_MAX)));
                
                //calculating the time elapsed
                int last_idx = (int)rate_cumsum.size() - 1;
                time = ((-1/ rate_cumsum[last_idx]) * log(random_double));
            }

            MPI_Bcast(&time, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
            
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
            
            MPI_Reduce(&max_i_rate, &max_rate, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
            MPI_Bcast(&max_rate, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
            
            if (rate_cumsum.size() != 0 ) {
                if (rate_cumsum[(rate_cumsum.size() - 1)] == max_rate) {
                }
                else {
                    fflush(stdout);
                    end_idx = (int)moves_lattice.rows() + 1;
                    moves_lattice.reshape(end_idx, 1, rank);
                    moves_lattice((end_idx-1),0) = 5;
                    rate_cumsum.push_back(max_rate);
                    moves_shifts.reshape(end_idx, 3, rank);
                    for (int i=0; i<3; i++) moves_shifts((end_idx-1), i) = 0;
                    moves_vacs.reshape(end_idx, 1, rank);
                    moves_vacs((end_idx-1),0) = -1;
                    moves_coords.reshape(end_idx, 4, rank);
                    for (int i=0; i<4; i++) moves_coords((end_idx-1) ,i) = -1;
                    moves_coords_unmod.reshape(end_idx, 4, rank);
                    for (int i=0; i<4; i++) moves_coords_unmod((end_idx-1) ,i) = -1;
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
                moves_coords_unmod.reshape(1,4, rank);
                for (int i=0; i<4; i++) moves_coords_unmod(0, i) = -1;
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
        
            for (int new_proc=0; new_proc<num_procs; new_proc++) { 
                if (new_proc != rank) {
                    MPI_Isend(NULL, 0, MPI_CHAR, new_proc, tag, MPI_COMM_WORLD, &request ); 
                    MPI_Wait(&request, MPI_STATUS_IGNORE);
                }    
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
            if (((rank == 0) || (rank == 1) || (rank == 2)) && (move_ticks <= 20)) {
                std::cout << "ROPM_PRE rank: " << rank << " move_ticks: " << move_ticks
                          << " prev_move_type.size(): " << prev_move_type.size()
                          << " par_prev_idx.size(): " << par_prev_idx.size() << "\n";
                for (int dbg_i = 0; dbg_i < (int)prev_move_type.size(); dbg_i++) {
                    std::cout << "ROPM_PRE rank: " << rank << " move_ticks: " << move_ticks
                              << " prev_move_type[" << dbg_i << "]: " << prev_move_type[dbg_i]
                              << " recorded_tick: " << prev_move_type_ticks[dbg_i] << "\n";
                }
                for (int dbg_i = 0; dbg_i < (int)par_prev_idx.size(); dbg_i++) {
                    std::cout << "ROPM_PRE rank: " << rank << " move_ticks: " << move_ticks
                              << " par_prev_idx[" << dbg_i << "]: " << par_prev_idx[dbg_i]
                              << " par_move_ticks: " << par_move_ticks[dbg_i]
                              << " par_prev_newlocs: (" << par_prev_newlocs[dbg_i][0] << "," << par_prev_newlocs[dbg_i][1] << "," << par_prev_newlocs[dbg_i][2] << "," << par_prev_newlocs[dbg_i][3] << ")\n";
                }
            }

            std::vector<std::vector<int>>::iterator ptrnew = par_prev_newlocs.begin();
            std::vector<std::vector<int>>::iterator ptrold = par_prev_oldlocs.begin();
            std::vector<int>::iterator ptr2 = par_prev_idx.begin();
            std::vector<int> remove_idxs;

            for (int i=0; i < (int)par_move_ticks.size(); i++) {
                if (move_ticks > (par_move_ticks.at(i) + 1)) {
                    par_prev_newlocs.erase(ptrnew);
                    par_prev_oldlocs.erase(ptrold);
                    remove_idxs.push_back(i);
                    par_prev_idx.erase(ptr2);
                }
            }
            
            std::vector<int>::iterator ptr3 = par_move_ticks.begin();
            for (int i=0; i < (int)remove_idxs.size(); i++) {
                par_move_ticks.erase((ptr3 + (remove_idxs.at(i) - i)));
            }

            remove_idxs.clear();
            for (int i=0; i < (int)prev_move_type_ticks.size(); i++) {
                if (move_ticks > (prev_move_type_ticks.at(i) + 1)) {
                    remove_idxs.push_back(i);
                }
            }

            std::vector<int>::iterator ptr4 = prev_move_type_ticks.begin();
            std::vector<int>::iterator ptr5 = prev_move_type.begin();

            for (int i=0; i < (int)remove_idxs.size(); i++) {
                prev_move_type_ticks.erase((ptr4 + (remove_idxs.at(i) - i)));
                prev_move_type.erase((ptr5 + (remove_idxs.at(i) - i)));
            }
        }

        /* 
        rolling back move originating in native process when error encountered 
        in simulation - one case for involving parallel transfer, one for no
        parallel transfer
        */
        void reverse_move(int parallel_transfer) {
            std::vector<int> new_vac;
            std::vector<int> new_vac_unmod;
            std::vector<int> old_vac;
            std::vector<int> move;
            int lattice;
            int idx;
            int i1; int i2; int i3; int i4;
            int vacs_idx;

            // Defensive: reverse_moves_wrapper's local_reversals_done guard is supposed to ensure
            // this ring buffer (see its comment) always has an entry left whenever this is called;
            // this only trips if that invariant is ever broken elsewhere (a prior instance was a
            // stray double pop_back() in reverse_move_parallel() draining a sibling ring buffer
            // out of sync -- now fixed). Fail loudly but don't crash the whole run.
            if (prev_newlocs.empty()) {
                std::cout << "rank: " << rank << " reverse_move: prev_newlocs empty (parallel_transfer: "
                          << parallel_transfer << ") -- ring buffer desync, skipping this reversal\n";
                return;
            }

            new_vac = prev_newlocs.at((prev_newlocs.size() - 1));
            new_vac_unmod = prev_newlocs_unmod.at((prev_newlocs_unmod.size() - 1));
            old_vac = prev_oldlocs.at((prev_oldlocs.size() - 1));
            move = prev_moves.at((prev_moves.size() - 1));
            lattice = prev_lattice.at((prev_lattice.size() - 1));
            vacs_idx  = prev_idxs.at((prev_idxs.size() - 1));

            if (lattice == 5) { std::cout << "rank: " << rank << " lattice: " << lattice << "\n"; }
            else {
                if (parallel_transfer != -1) {
                    
                    // std::cout << "rank: " << rank << " lattice: " << lattice << "\n";
                    // std::cout << "rank: " << rank << " pre add vacs_idx: " << vacs_idx << "\n";

                    i1 = old_vac[0];
                    i2 = old_vac[1];
                    i3 = old_vac[2];
                    i4 = old_vac[3];
                    
                    Matrix<int> only_vacancies = vacancies.nonzero(rank); // configuration of vacancies at current timestep
                    /*
                    std::cout << "rank: " << rank << " post rollback only_vacancies.rows(): " << only_vacancies.rows() << "  post rollback vacancies_pos.rows(): " << vacancies_pos.rows() << "\n"; 
                
                    std::cout << "rank: " << rank << " pre reverse vacancies_pos(vacs_idx,0): " << vacancies_pos(vacs_idx,0) << " vacancies_pos(vacs_idx,1): " << vacancies_pos(vacs_idx,1) 
                        << " vacancies_pos(vacs_idx,2): " << vacancies_pos(vacs_idx,2)  << " vacancies_pos(vacs_idx,3): " << vacancies_pos(vacs_idx,3) << "\n";
                    */
                    vacancies_pos.add_row(vacs_idx, rank);
                    /*
                    std::cout << "rank: " << rank << " (added) pre reverse vacancies_pos(vacs_idx,0): " << vacancies_pos(vacs_idx,0) << " vacancies_pos(vacs_idx,1): " << vacancies_pos(vacs_idx,1) 
                        << " vacancies_pos(vacs_idx,2): " << vacancies_pos(vacs_idx,2)  << " vacancies_pos(vacs_idx,3): " << vacancies_pos(vacs_idx,3) << "\n";
                        std::cout << "rank: " << rank << " (added) pre reverse vacancies_pos(vacs_idx+1,0): " << vacancies_pos(vacs_idx+1,0) << " vacancies_pos(vacs_idx+1,1): " << vacancies_pos(vacs_idx+1,1) 
                        << " vacancies_pos(vacs_idx+1,2): " << vacancies_pos(vacs_idx+1,2)  << " vacancies_pos(vacs_idx+1,3): " << vacancies_pos(vacs_idx+1,3) << "\n";
                    */
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

                }
                else {

                    i1 = old_vac[0];
                    i2 = old_vac[1];
                    i3 = old_vac[2];
                    i4 = old_vac[3];
                    
                    // moving vacancy from bc site to bc site 
                    if (lattice == 3) {
                        /*---
                        BC EDGE MOVES
                        ---*/
                                   
                        // switching occupancy for old and new site in lattice array ###
                        bc_sites((size_t)0, (size_t)new_vac[1], (size_t)new_vac[2], (size_t)new_vac[3]) = 1;
                        bc_sites((size_t)0, (size_t)old_vac[1], (size_t)old_vac[2], (size_t)old_vac[3]) = 0;

                        // switching occupancy for old and new site in vacancy and mobileion arrays ###
                        vacancies((size_t)new_vac[0], (size_t)new_vac[1], (size_t)new_vac[2], (size_t)new_vac[3]) = 0;
                        vacancies((size_t)1, (size_t)old_vac[1], (size_t)old_vac[2], (size_t)old_vac[3]) = 1;
                        
                        vacancies_pos(vacs_idx,0) = 1;
                    }
                                
                    // moving vacancy from vertex site to vertex site 
                    else if (lattice == 2) {
                        /*---
                        VERTEX EDGE MOVES
                        ---*/
                        
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

                    // std::cout << "rank: " << rank << " post reverse reverse vacancies_pos(vacs_idx,0): " << vacancies_pos(vacs_idx,0) << " vacancies_pos(vacs_idx,1): " << vacancies_pos(vacs_idx,1) 
                    //    << " vacancies_pos(vacs_idx,2): " << vacancies_pos(vacs_idx,2)  << " vacancies_pos(vacs_idx,3): " << vacancies_pos(vacs_idx,3) << "\n";

                    
                }
                
                if ( ((new_vac[1] < 2) || (new_vac[1] > (sublattice_dim[0] - 3))) || ((new_vac[2] < 2) || (new_vac[2] > (sublattice_dim[1] - 3))) ||
                ((old_vac[1] < 2) || (old_vac[1] > (sublattice_dim[0] - 3))) || ((old_vac[2] < 2) || (old_vac[2] > (sublattice_dim[1] - 3))) )
                { /*checking to see if neighbor ghost sites need to be updated */
                    // ghost_site_send's first four args double as the dispatch-direction
                    // coordinates and the "old position" fields in the transmitted buffer, so
                    // they need the raw (pre-wrap) value here -- the same reason new_loc_unmod
                    // (not new_loc) is required at the original forward-move call site.
                    ghost_site_send(new_vac_unmod[0],new_vac_unmod[1],new_vac_unmod[2],new_vac_unmod[3],old_vac,move,parallel_transfer,true);
                    ghost_site_self_reference(new_vac[0],new_vac[1],new_vac[2],new_vac[3],old_vac,parallel_transfer,"GSSR_reverse_move");
                }
            }

            prev_moves.pop_back();
            prev_newlocs.pop_back();
            prev_newlocs_unmod.pop_back();
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

            // Defensive: see the matching guard in reverse_move(). par_prev_idx/par_prev_newlocs/
            // par_prev_oldlocs/par_move_ticks are meant to stay in lockstep (one push, one pop,
            // per store_parallel_info()/reverse_move_parallel() call); a stray duplicate
            // par_prev_idx.pop_back() below used to drain it twice as fast as its siblings.
            if (par_prev_idx.empty()) {
                std::cout << "rank: " << rank << " reverse_move_parallel: par_prev_idx empty -- ring buffer desync, skipping this reversal\n";
                return;
            }

            vacs_idx = par_prev_idx.at(par_prev_idx.size()-1);

            i1 = (par_prev_newlocs.at((par_prev_newlocs.size()-1)).at(0));
            i2 = (par_prev_newlocs.at((par_prev_newlocs.size()-1)).at(1));
            i3 = (par_prev_newlocs.at((par_prev_newlocs.size()-1)).at(2));
            i4 = (par_prev_newlocs.at((par_prev_newlocs.size()-1)).at(3));
            std::vector<int> shift{0,0,0};

            std::cout << "rank: " << rank << " vacs_idx: " << vacs_idx << "\n";
            std::cout << "rank: " << rank << " vacancies_pos(vacs_idx,0): " << vacancies_pos(vacs_idx,0) << " vacancies_pos(vacs_idx,1): " << vacancies_pos(vacs_idx,1) 
                << " vacancies_pos(vacs_idx,2): " << vacancies_pos(vacs_idx,2)  << " vacancies_pos(vacs_idx,3): " << vacancies_pos(vacs_idx,3) << "\n";
            

            if (i1 == 1) {
                // moving vacancy from bc site to bc site 
                bc_sites((size_t)0, (size_t)i2, (size_t)i3, (size_t)i4) = 1;
            }
            else if (i1 == 0) {
                // moving vacancy from vertex site to vertex site 
                vertex_sites((size_t)0, (size_t)i2, (size_t)i3, (size_t)i4) = 1;
            }

            vacancies(i1,i2,i3,i4) = 0;
            
            vacancies_pos.remove_row(vacs_idx, rank);
            num_of_vacs --;
            
            /*
                for (int i=0; i<(int)par_prev_idx.size(); i++)  {
                    if ((vacs_idx <= par_prev_idx[i]) && (par_prev_idx[i] != -1)) { par_prev_idx[i] --; }
                }
                for (int i=0; i<(int)prev_idxs.size(); i++)  {
                    if ((vacs_idx <= prev_idxs[i]) && (prev_idxs[i] != -1)) { prev_idxs[i] --; }
                }
            */

            std::vector<int> old_loc(4);
            old_loc[0] = (par_prev_oldlocs.at((par_prev_oldlocs.size()-1)).at(0));
            old_loc[1] = (par_prev_oldlocs.at((par_prev_oldlocs.size()-1)).at(1));
            old_loc[2] = (par_prev_oldlocs.at((par_prev_oldlocs.size()-1)).at(2));
            old_loc[3] = (par_prev_oldlocs.at((par_prev_oldlocs.size()-1)).at(3));

            std::cout << "rank: " << rank << " paralllel new_i1: " << i1 << " new_i2: " << i2 << " new_i3: " << i3 << " new_i4: " << i4 << "\n";
            std::cout << "rank: " << rank << " parallel old_i1: " << old_loc[0] << " old_i2: " << old_loc[1] << " old_i3: " << old_loc[2] << " old_i4: " << old_loc[3] << "\n";

            par_move_ticks.pop_back();
            par_prev_newlocs.pop_back();
            par_prev_oldlocs.pop_back();
            par_prev_idx.pop_back();

            // Only i1..i4 (this rank's own position for the vacancy being deleted) is meaningful
            // here -- old_loc's coordinates belong to whichever remote rank originally owned this
            // vacancy, not to this rank, so they can't be checked against this rank's own
            // sublattice_dim boundaries.
            // (i1 is the sublattice type (0/1), not a coordinate -- an "i1 < 1" clause here was
            // testing vertex-vs-bc instead of edge proximity and has been dropped.)
            // i2/i3 are this position's x/y coordinates (i4 is z) -- the low-x edge (i2 < 1) was
            // previously missing entirely, and the high-y edge check compared i4 (z) against
            // sublattice_dim[1] (the y bound) instead of i3 (y), so a reversed parallel-transfer
            // vacancy near either of those edges silently skipped notifying its ghost neighbors.
            if ( (i2 < 1) || (i2 > (sublattice_dim[0] - 2)) || (i3 < 1) || (i3 > (sublattice_dim[1] - 2)) )
                { /*checking to see if neighbor ghost sites need to be updated */
                    // Pass {i1,i2,i3,i4} as both old and new position: this is a pure removal (the
                    // vacancy is leaving this rank entirely, back to its original remote owner,
                    // which independently restores its own ghost tracking via its own
                    // reverse_move()) -- see ghost_site_recieve's clear_only handling. Using
                    // old_loc here instead (a remote rank's local coordinates, re-sent under THIS
                    // rank's old_proc) previously caused ghost_site_recieve to decode a bogus
                    // "arrival" site using the wrong rank's chunk offset, permanently setting a
                    // ghost bit nothing would ever clear.
                    std::vector<int> self_pos = {i1,i2,i3,i4};
                    ghost_site_send(i1,i2,i3,i4,self_pos,shift,0,true);
                    ghost_site_self_reference(i1,i2,i3,i4,self_pos,true,"GSSR_reverse_move_parallel");
                }
        }

        /**
        * @brief Decrements the simulation clock based on previously stored time values.
        *
        * This function rolls back the simulation clock by subtracting the recorded
        * time intervals from the previous moves. It also clears the list of previous
        * time values after adjusting the clock.
        */
        void deincrement_time() {
            //print_1Dvector(prev_times);
            for (int i=0; i < (int)prev_times.size(); i++) {
                t = t - prev_times[i];}
            
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
        * @param new_loc_unmod The raw (pre-wrap) new location of the vacancy (default is {0,0,0,0}),
        *        needed so a later reverse_move() can re-derive which neighbor(s) the original
        *        ghost_site_send dispatched to -- see reverse_move()'s use of prev_newlocs_unmod.
        */
        void store_move_info(int idx,  int parallel_transfer, int vac_idx, int lattice = -1, std::vector<int> shift = {0,0,0}, const std::vector<int>& new_loc = {0,0,0,0}, const std::vector<int>& old_loc = {0,0,0,0}, const std::vector<int>& new_loc_unmod = {0,0,0,0}) {

            std::vector<int> old_vac(4);
            std::vector<int> new_vac(4);
            std::vector<int> new_vac_unmod(4);

            // keeping track of new location if parallel transfer
            if (parallel_transfer != -1) {

                for (int i=0; i<(int)new_vac.size(); i++)  {
                    new_vac[i] = new_loc[i];
                    old_vac[i] = old_loc[i];
                    new_vac_unmod[i] = new_loc_unmod[i];
                }
                /*
                    for (int i=0; i<(int)par_prev_idx.size(); i++)  {
                        if ((vac_idx <= par_prev_idx[i]) && (par_prev_idx[i] != -1)) { par_prev_idx[i] --; }
                    }
                    for (int i=0; i<(int)prev_idxs.size(); i++)  {
                        if ((vac_idx <= prev_idxs[i]) && (prev_idxs[i] != -1)) { prev_idxs[i] --; }
                    }
                */
            }
            // entering filler values if null move
            else if (vac_idx == -1) {
                for (int i=0; i<(int)shift.size(); i++)  {
                    shift[i] = -1;
                }
                for (int i=0; i<(int)new_vac.size(); i++)  {
                    new_vac[i] = new_loc[i];
                    old_vac[i] = old_loc[i];
                    new_vac_unmod[i] = new_loc_unmod[i];
                }
            }
            // standard move
            else {
                for (int i=0; i<(int)new_vac.size(); i++)  {
                    new_vac[i] = new_loc[i];
                    old_vac[i] = old_loc[i];
                    new_vac_unmod[i] = new_loc_unmod[i];
                }
            }

            if ((int)prev_newlocs.size() >= 2) {

                prev_moves[0] = prev_moves[1];
                prev_moves[1] = shift;

                prev_newlocs[0] = prev_newlocs[1];
                prev_newlocs[1] = new_vac;

                prev_newlocs_unmod[0] = prev_newlocs_unmod[1];
                prev_newlocs_unmod[1] = new_vac_unmod;

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
                prev_newlocs_unmod.push_back(new_vac_unmod);
                prev_oldlocs.push_back(old_vac);
                prev_idxs.push_back(vac_idx);
                prev_lattice.push_back(lattice);
            }
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
            std::vector<int> vac_new(4);
            vac_new[0] = parallel_buffer[0];
            vac_new[1] = parallel_buffer[1];
            vac_new[2] = parallel_buffer[2];
            vac_new[3] = parallel_buffer[3];

            std::vector<int> vac_old(4);
            int old_proc = (parallel_buffer[4]);
            vac_old[0] = (parallel_buffer[5]);
            vac_old[1] = (parallel_buffer[6]);
            vac_old[2] = (parallel_buffer[7]);
            vac_old[3] = (parallel_buffer[8]);
            std::vector<int> shift = {0,0,0};

            /*
            for (int i=0; i<(int)par_prev_idx.size(); i++)  {
                if ((vac_idx <= par_prev_idx[i]) && (par_prev_idx[i] != -1)) { par_prev_idx[i] ++; }
            }
            for (int i=0; i<(int)prev_idxs.size(); i++)  {
                if ((vac_idx <= prev_idxs[i]) && (prev_idxs[i] != -1)) { prev_idxs[i] ++; }
            }
            */
            

            par_prev_oldlocs.push_back(vac_old);
            prev_moves.push_back(shift);
            par_prev_newlocs.push_back(vac_new);
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

            int parallel_transfer = -1;

            // prev_newlocs/prev_oldlocs/prev_newlocs_unmod/etc. (consumed by reverse_move()) are a
            // hard-capped 2-slot ring buffer -- they only ever remember this rank's 2 most recent
            // local (type 0/1) moves, matching comm_boundary_conflict()'s "roll back two steps"
            // design intent. prev_move_type itself is pruned only by tick age, so it can still hold
            // more than 2 type-0/1 entries; calling reverse_move() for a 3rd+ one finds that history
            // already gone (reverse_move().at() on an empty vector), so we stop once both slots are
            // used. Any older type-0/1 entry is left un-reversed -- there's no data left to reverse
            // it with.
            int local_reversals_done = 0;

            for (int i=(prev_move_type.size()-1); i>=0; i--) {

                std::cout << "prev_move_type: " << prev_move_type[i] << "\n";

                if (prev_move_type[i] == 0) {
                    if (local_reversals_done >= 2) continue;
                    parallel_transfer = -1;
                    reverse_move(parallel_transfer); // reversing in-lattice move with no parallel transfer
                    local_reversals_done++;
                }
                else if (prev_move_type[i] == 1) {
                    if (local_reversals_done >= 2) continue;
                    parallel_transfer = 1; // need to determine previous proc
                    reverse_move(parallel_transfer); // reversing move with sending parallel transfer
                    local_reversals_done++;
                }
                else if (prev_move_type[i] == 2) reverse_move_parallel(); // reversing move with recieving parallel transfer
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
            
            if (vacancies((size_t)i, (size_t)j, (size_t)k, (size_t)l) == 1) {
                std::cout << "RMP_CONFLICT rank: " << rank << " move_ticks: " << move_ticks
                          << " incoming(i,j,k,l): (" << i << "," << j << "," << k << "," << l << ")"
                          << " old_proc: " << old_proc << " old(i,j,k,l): (" << i_old << "," << j_old << "," << k_old << "," << l_old << ")\n";

                Matrix<int> only_vacancies = vacancies.nonzero(rank); // configuration of vacancies at current timestep
                Matrix<int> unequal_elems_mat1 = comparison(only_vacancies, vacancies_pos);
                Matrix<int> unequal_elems_mat2 = comparison(vacancies_pos, only_vacancies);
                /*
                std::cout << "rank: " << rank << " move_ticks: " << move_ticks << " unequal_elems_mat1: \n";
                unequal_elems_mat1.print();
                std::cout << "rank: " << rank << " move_ticks: " << move_ticks << " unequal_elems_mat2: \n";
                unequal_elems_mat2.print();
                */

                if ((unequal_elems_mat1.rows() != 0) || (unequal_elems_mat2.rows() != 0)) {
                    std::cout << "rank: " << rank << " vacancies_pos: \n";
                    vacancies_pos.print();
                    std::cout << "rank: " << rank << " only_vacancies: \n";
                    only_vacancies.print();
                    exit(0);
                }

                interboundary_conflict = true;
                comm_boundary_conflict();

                return true;
            }
            else {
                std::cout << "RMP_ACCEPT rank: " << rank << " move_ticks: " << move_ticks
                          << " incoming(i,j,k,l): (" << i << "," << j << "," << k << "," << l << ")"
                          << " old_proc: " << old_proc << " old(i,j,k,l): (" << i_old << "," << j_old << "," << k_old << "," << l_old << ")\n";
                /* store previous two moves */
                size_t rows = vacancies_pos.rows();
                size_t cols = vacancies_pos.cols();
                vacancies_pos.reshape(rows+1, cols, rank);

                store_parallel_info(new_loc_buffer, move_ticks, rows);
                
                vacancies((size_t)i, (size_t)j, (size_t)k, (size_t)l) = 1;
                if (i == 0) vertex_sites((size_t)0, (size_t)j, (size_t)k, (size_t)l) = 0;
                else if (i == 1) bc_sites((size_t)0, (size_t)j, (size_t)k, (size_t)l) = 0;

                vacancies_pos(rows,0) = i;
                vacancies_pos(rows,1) = j;
                vacancies_pos(rows,2) = k;
                vacancies_pos(rows,3) = l;

                num_of_vacs ++;
                prev_move_type.push_back(2);
                prev_move_type_ticks.push_back(move_ticks);

                return false;
            }

            return interboundary_conflict;
        }

        /**
        * @brief Handles the reversal of ghost sites during a simulation.
        *
        * This function coordinates communication between MPI processes to manage
        * ghost sites and resolve conflicts in a parallel simulation environment.
        * It utilizes MPI for communication and ensures all processes reach
        * consensus on the ghost site states before proceeding.
        *
        * @param move_ticks The current simulation tick used for synchronization.
        */
        void reverse_ghost_sites(int move_ticks) {
            ghost_trace_tick = move_ticks; // diagnostic only, see GHOST_WRITE tracing
            MPI_Request request1; ///< MPI request handle for non-blocking send operations.
            MPI_Status status1;   ///< MPI status to capture information about received messages.
            std::vector<int> new_loc_buffer1(10); ///< Buffer for receiving ghost site data from tag 4.
            std::vector<int> new_loc_buffer2(10); ///< Buffer for receiving ghost site data from tag 5.
            std::vector<bool> stop_conflict_ghost_array(num_procs, 0); ///< Tracks whether each process is done with conflict resolution.
            stop_conflict_ghost_array[rank] = 1; ///< Mark the current process as done initially.
            bool stop_conflict_ghost = 0; ///< Flag indicating whether all processes have resolved conflicts.

            MPI_Barrier(MPI_COMM_WORLD); ///< Synchronize all processes before starting communication.

            // Notify all other processes of conflict resolution readiness
            for (int new_proc = 0; new_proc < num_procs; new_proc++) { 
                if (new_proc != rank) {
                    MPI_Isend(NULL, 0, MPI_CHAR, new_proc, conflict_ghost_done_tag, MPI_COMM_WORLD, &request1); 
                    MPI_Wait(&request1, MPI_STATUS_IGNORE);
                }
            }

            int loop_count = 0; ///< Counter to prevent infinite loops in case of communication issues.

            // Main loop to handle incoming messages and resolve conflicts
            while (!stop_conflict_ghost) {
                if (loop_count > 25) exit(0); ///< Safety mechanism to prevent infinite looping.

                stop_conflict_ghost = 1; ///< Assume all conflicts are resolved initially.

                status1.MPI_TAG = 0;

                MPI_Probe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status1);

                // Handle message with tag 4: Receive and process ghost site data
                if (status1.MPI_TAG == 4) {
                    MPI_Recv(new_loc_buffer1.data(), 10, MPI_INT, status1.MPI_SOURCE, 4, MPI_COMM_WORLD, &status1);
                    std::cout << "rank: " << rank << " GSR_TRACE tag: 4 from: " << status1.MPI_SOURCE
                              << " old_proc: " << new_loc_buffer1[4]
                              << " old(i,j,k,l): (" << new_loc_buffer1[5] << "," << new_loc_buffer1[6] << "," << new_loc_buffer1[7] << "," << new_loc_buffer1[8] << ")"
                              << " new(i,j,k,l): (" << new_loc_buffer1[0] << "," << new_loc_buffer1[1] << "," << new_loc_buffer1[2] << "," << new_loc_buffer1[3] << ")\n";
                    ghost_site_recieve(new_loc_buffer1, true, true);
                }

                // Handle message with tag 5: Receive and process ghost site data
                if (status1.MPI_TAG == 5) {
                    MPI_Recv(new_loc_buffer2.data(), 10, MPI_INT, status1.MPI_SOURCE, 5, MPI_COMM_WORLD, &status1);
                    std::cout << "rank: " << rank << " GSR_TRACE tag: 5 from: " << status1.MPI_SOURCE
                              << " old_proc: " << new_loc_buffer2[4]
                              << " old(i,j,k,l): (" << new_loc_buffer2[5] << "," << new_loc_buffer2[6] << "," << new_loc_buffer2[7] << "," << new_loc_buffer2[8] << ")"
                              << " new(i,j,k,l): (" << new_loc_buffer2[0] << "," << new_loc_buffer2[1] << "," << new_loc_buffer2[2] << "," << new_loc_buffer2[3] << ")\n";
                    ghost_site_recieve(new_loc_buffer2, false, true);
                }

                // Handle conflict resolution completion notification
                if (status1.MPI_TAG == conflict_ghost_done_tag) {
                    MPI_Recv(NULL, 0, MPI_CHAR, status1.MPI_SOURCE, conflict_ghost_done_tag, MPI_COMM_WORLD, &status1);
                    stop_conflict_ghost_array[status1.MPI_SOURCE] = 1;
                }

                // Check if all processes have resolved their conflicts
                for (int j = 0; j < (int)stop_conflict_ghost_array.size(); j++) {
                    stop_conflict_ghost &= stop_conflict_ghost_array[j];
                }
                loop_count++;
            }
        }

        /**
        * @brief Receive communication (ghost site update, lattice site update, conflict) using blocking communication.
        *
        * This function handles the reception of messages from other processes,
        * including updates to ghost sites, lattice sites, and conflicts.
        *
        * @param move_ticks The current iteration for moves in the simulation.`
        */
        void receive_parallel_comm_helper(int move_ticks) {
            ghost_trace_tick = move_ticks; // diagnostic only, see GHOST_WRITE tracing

            int test_flag1 = 0;
            MPI_Status status1;
            std::vector<int> new_loc_buffer1(9);
            std::vector<int> new_loc_buffer2(10);
            std::vector<int> new_loc_buffer3(10);

            std::vector<bool> stop_ghost_Array(num_procs, 0);
            std::vector<bool> stop_par_Array(num_procs, 0);
            std::vector<bool> stop_conflict_array(num_procs, 0);

            stop_ghost_Array[rank] = 1;
            stop_par_Array[rank] = 1;
            stop_conflict_array[rank] = 1;

            bool stop_ghost = 0;
            bool stop_par = 0;
            bool need_reverse = 0;
            bool stop_conflict = 0;
            bool interboundary_conflict = false;
            int curr_num_vacs = 0;
            
            while ((!stop_par) || (!stop_ghost)) {

                stop_par = 1;
                stop_ghost = 1;

                status1.MPI_TAG = 0;
                
                MPI_Probe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status1);
                //std::cout << "rank: " << rank <<  " move_ticks: " << move_ticks <<" stop par status1.MPI_TAG: " << status1.MPI_TAG << " status1.MPI_SOURCE: " << status1.MPI_SOURCE << "\n";
                

                // test for both message 1 and message 2 simulatneously in if statement
                if (status1.MPI_TAG == 1) {
                    MPI_Recv(new_loc_buffer1.data(), 9, MPI_INT, status1.MPI_SOURCE, 1, MPI_COMM_WORLD, &status1);
                    interboundary_conflict = (recieve_move_parallel(new_loc_buffer1, move_ticks) || interboundary_conflict);
                }
                
                if (status1.MPI_TAG == 2) {
                    MPI_Recv(new_loc_buffer2.data(), 10, MPI_INT, status1.MPI_SOURCE, 2, MPI_COMM_WORLD, &status1);
                    std::cout << "rank: " << rank << " GSR_TRACE tag: 2 from: " << status1.MPI_SOURCE
                              << " old_proc: " << new_loc_buffer2[4]
                              << " old(i,j,k,l): (" << new_loc_buffer2[5] << "," << new_loc_buffer2[6] << "," << new_loc_buffer2[7] << "," << new_loc_buffer2[8] << ")"
                              << " new(i,j,k,l): (" << new_loc_buffer2[0] << "," << new_loc_buffer2[1] << "," << new_loc_buffer2[2] << "," << new_loc_buffer2[3] << ")\n";
                    ghost_site_recieve(new_loc_buffer2, true, false);
                }

                if (status1.MPI_TAG == 3) {
                    MPI_Recv(new_loc_buffer3.data(), 10, MPI_INT, status1.MPI_SOURCE, 3, MPI_COMM_WORLD, &status1);
                    std::cout << "rank: " << rank << " GSR_TRACE tag: 3 from: " << status1.MPI_SOURCE
                              << " old_proc: " << new_loc_buffer3[4]
                              << " old(i,j,k,l): (" << new_loc_buffer3[5] << "," << new_loc_buffer3[6] << "," << new_loc_buffer3[7] << "," << new_loc_buffer3[8] << ")"
                              << " new(i,j,k,l): (" << new_loc_buffer3[0] << "," << new_loc_buffer3[1] << "," << new_loc_buffer3[2] << "," << new_loc_buffer3[3] << ")\n";
                    ghost_site_recieve(new_loc_buffer3, false, false);
                }
                
                if (status1.MPI_TAG == par_done_tag) {
                    MPI_Recv(NULL, 0, MPI_CHAR, status1.MPI_SOURCE, par_done_tag, MPI_COMM_WORLD, &status1);
                    stop_par_Array[status1.MPI_SOURCE] = 1;
                }
                
                for (int j = 0; j < (int)stop_par_Array.size(); j++) {
                    stop_par &= stop_par_Array[j];
                }
                
                if (status1.MPI_TAG == ghost_done_tag) {
                    MPI_Recv(NULL, 0, MPI_CHAR, status1.MPI_SOURCE, ghost_done_tag, MPI_COMM_WORLD, &status1);
                    stop_ghost_Array[status1.MPI_SOURCE] = 1;
                }
                
                for (int j = 0; j < (int)stop_ghost_Array.size(); j++) {
                    stop_ghost = (stop_ghost_Array[j] && stop_ghost);
                }
            }
            MPI_Barrier(MPI_COMM_WORLD);

            MPI_Request request1;
            for (int new_proc=0; new_proc<num_procs; new_proc++) { 
                if (new_proc != rank) {
                    //std::cout << "rank: " << rank << " move_ticks: " << move_ticks << " 1 new_proc: " << new_proc << "\n";
                    MPI_Isend( NULL, 0, MPI_CHAR, new_proc, conflict_done_tag, MPI_COMM_WORLD, &request1); 
                    MPI_Wait(&request1, MPI_STATUS_IGNORE);
                }
            }

            //std::cout << "rank: " << rank << " move_ticks: " << move_ticks << " sent conflict_done_tag \n";

            remove_old_par_moves(move_ticks);

            int loop_count = 0;
            while ((!stop_conflict)) {
                if (loop_count > 3000) { std::cout << "rank: " << rank << " exceeded receive_parallel_comm_helper loop count\n"; exit(0); }
                stop_conflict = 1;
                status1.MPI_TAG = 0;
                MPI_Probe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status1);
                //std::cout << "rank: " << rank <<  " move_ticks: " << move_ticks <<" conflict 1 status1.MPI_TAG: " << status1.MPI_TAG << " status1.MPI_SOURCE: " << status1.MPI_SOURCE << "\n";
                //std::cout << "rank: " << rank << " status1.MPI_TAG: " << status1.MPI_TAG << "\n";
                
                if (status1.MPI_TAG == 6) {
                    MPI_Recv(NULL, 0, MPI_CHAR, status1.MPI_SOURCE, 6, MPI_COMM_WORLD, &status1);
                    need_reverse = 1;
                }
                if (status1.MPI_TAG == conflict_done_tag) {
                    //std::cout << "rank: " << rank << " move_ticks: " << move_ticks << " pre stop_conflict_array: \n";
                    //print_1Dvector(stop_conflict_array);
                    MPI_Recv(NULL, 0, MPI_CHAR, status1.MPI_SOURCE, conflict_done_tag, MPI_COMM_WORLD, &status1);
                    stop_conflict_array[status1.MPI_SOURCE] = 1;

                    //std::cout << "rank: " << rank <<  " move_ticks: " << move_ticks <<" post stop_conflict_array: \n";
                    //print_1Dvector(stop_conflict_array);
                }
                for (int j = 0; j < (int)stop_conflict_array.size(); j++) {
                    stop_conflict &= stop_conflict_array[j];
                }
                loop_count ++;
            }
            //std::cout << "rank: " << rank << " move_ticks: " << move_ticks << " recv conflict_done_tag \n";
            MPI_Barrier(MPI_COMM_WORLD);

            if ((need_reverse) || (interboundary_conflict)) {
                // issue with rollback introducing new vacancies
                Matrix<int> only_vacancies = vacancies.nonzero(rank); // configuration of vacancies at current timestep
                // std::cout << "rank: " << rank << " pre rollback only_vacancies.rows(): " << only_vacancies.rows() << "  post rollback vacancies_pos.rows(): " << vacancies_pos.rows() << "\n";
                deincrement_time();
                reverse_moves_wrapper();
                only_vacancies = vacancies.nonzero(rank); // configuration of vacancies at current timestep
                //std::cout << "rank: " << rank << " post rollback only_vacancies.rows(): " << only_vacancies.rows() << "  post rollback vacancies_pos.rows(): " << vacancies_pos.rows() << "\n"; 
                
                Matrix<int> unequal_elems_mat1 = comparison(only_vacancies, vacancies_pos);
                Matrix<int> unequal_elems_mat2 = comparison(vacancies_pos, only_vacancies);

                /*
                std::cout << "rank: " << rank << " move_ticks: " << move_ticks << " unequal_elems_mat1: \n";
                unequal_elems_mat1.print();
                std::cout << "rank: " << rank << " move_ticks: " << move_ticks << " unequal_elems_mat2: \n";
                unequal_elems_mat2.print();
                */
                
                if ((unequal_elems_mat1.rows() != 0) || (unequal_elems_mat2.rows() != 0)) {
                    std::cout << "rank: " << rank << " vacancies_pos: \n";
                    vacancies_pos.print();
                    std::cout << "rank: " << rank << " only_vacancies: \n";
                    only_vacancies.print();
                    exit(0);
                }
                    
            }
            
            reverse_ghost_sites(move_ticks);
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
        // TO CONNECT TO VNC: ssh -f -N -L xxxx:localhost:yyyy cfrech@ls6.tacc.utexas.edu
        lattice_return_struct new_kmc_iterator(double time_lim, std::chrono::system_clock::time_point start, std::string folder, int iteration, int last_tick, double last_time) {
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


            std::cout << "rank: " << rank << " pre get_actions() \n";
            parallel_get_actions(); // updating list of moves in system
            std::cout << "rank: " << rank << " post get_actions() \n";
            if (rank == 2) {
                std::cout << "rank: " << rank << " POST_GET_ACTIONS_READBACK proc_neg_y_neighbors(1,49,1,124): "
                          << proc_neg_y_neighbors(1,49,1,124) << "\n";
            }
            
            t = last_time;             
            int move_ticks = last_tick;
            double old_time;
            bool restart = false;
            bool reconstruct = true;
            double total_E_allprocs = 0;
            double total_cost_allprocs = 0;

            if (move_ticks != 0) { restart = true; }

            int prev_num_vacs = 0;
            int curr_num_vacs = 0;
            bool get_rank = true;
            int write_frames = 1;

            std::ostringstream ss;

            Matrix<int> only_vacancies = vacancies.nonzero(rank);
            std::cout << "rank: " << rank << " move_ticks: " << move_ticks << " vacancies_pos.rows(): " << vacancies_pos.rows() << " vs only_vacancies.rows(): " << only_vacancies.rows() << "\n";
            
            MPI_Barrier(MPI_COMM_WORLD);
            
            double last_rate_plus1;
            double last_rate;
            double last_rate_minus1;

            while (t  < time_lim) {
                //std::cout << "rank: " << rank << " move_ticks: " << move_ticks << "\n";
                if ((rank == 2) && (move_ticks <= 25)) {
                    std::cout << "rank: " << rank << " LOOP_TOP_READBACK move_ticks: " << move_ticks
                              << " proc_neg_y_neighbors(1,49,1,124): " << proc_neg_y_neighbors(1,49,1,124) << "\n";
                }
                if ((rank == 2) && (move_ticks <= 20)) {
                    std::cout << "rank: " << rank << " LOOP_TOP_VACREADBACK move_ticks: " << move_ticks
                              << " vacancies(1,54,0,119): " << vacancies(1,54,0,119)
                              << " vacancies(0,54,1,119): " << vacancies(0,54,1,119) << "\n";
                }
                if ((rank == 0) && (move_ticks <= 20)) {
                    std::cout << "rank: " << rank << " LOOP_TOP_VACREADBACK move_ticks: " << move_ticks
                              << " proc_pos_y_neighbors(1,56,1,119): " << proc_pos_y_neighbors(1,56,1,119) << "\n";
                }
                // Investigating the 10-tick (1421-1430) mismatch on rank3's
                // proc_neg_x_neighbors(1,0,7,119) -- decodes to global (w=1,x=54,y=60,z=119),
                // which is rank2's own local (1,54,5,119). No traced SET write ever touches the
                // rank3 cell in the whole run, so this finds when the vacancy at that rank2 site
                // actually arrived, to check whether its arrival's forward dispatch should have
                // (but didn't) notify rank3.
                if ((rank == 2) && (move_ticks >= 1390) && (move_ticks <= 1435)) {
                    std::cout << "rank: " << rank << " LOOP_TOP_TARGETREADBACK move_ticks: " << move_ticks
                              << " vacancies(1,54,5,119): " << vacancies(1,54,5,119) << "\n";
                }
                // One-shot: catches the arrival even if it happened long before the 1390 window
                // above (this vacancy could have been sitting there for a while).
                {
                    static bool target_vac_seen_once = false;
                    if ((rank == 2) && (!target_vac_seen_once) && (vacancies(1,54,5,119) == 1)) {
                        target_vac_seen_once = true;
                        std::cout << "rank: " << rank << " TARGET_ARRIVAL_LATCH move_ticks: " << move_ticks
                                  << " vacancies(1,54,5,119) first seen == 1\n";
                    }
                }
                if ((rank == 3) && (move_ticks >= 1390) && (move_ticks <= 1435)) {
                    std::cout << "rank: " << rank << " LOOP_TOP_TARGETREADBACK move_ticks: " << move_ticks
                              << " proc_neg_x_neighbors(1,0,7,119): " << proc_neg_x_neighbors(1,0,7,119) << "\n";
                }
                total_cost_allprocs = sum_values_allprocs(total_cost);
                std::cout << " rank: " << rank << " total_cost: " << total_cost << "\n";
                //std::cout << " rank: " << rank << "total_cost: " << total_cost << "\n";
                total_E_allprocs = sum_values_allprocs(system_energy);
                //std::cout << "rank: " << rank << " test entering loop \n";
                
                if (rank == 0) {
                    std::cout << " rank: " << rank << " total_cost: " << total_cost << "\n";
                    std::cout << " rank: " << rank << " total_cost_allprocs: " << total_cost_allprocs << "\n";
                    std::cout << " rank: " << rank << " total_E_allprocs: " << total_E_allprocs << "\n";
                    std::cout  <<" rank: " << rank << " total E + cost: " << total_E_allprocs - total_cost_allprocs << "\n";
                }

                if (move_ticks > 10000) {exit(0);}


                if ( rate_cumsum.size() != 0) {
                    end = std::chrono::system_clock::now(); 
                    rand_idx = get_idx();
                    
                    
                    // std::cout << "rank: " << rank << " pre new_update_lattice\n";
                    new_update_lattice(rand_idx, move_ticks);
                    // std::cout << "rank: " << rank << " post new_update_lattice\n";
                    if ((rank == 2) && (move_ticks <= 25)) {
                        std::cout << "rank: " << rank << " POST_NUL_READBACK move_ticks: " << move_ticks
                                  << " proc_neg_y_neighbors(1,49,1,124): " << proc_neg_y_neighbors(1,49,1,124) << "\n";
                    }
                    if ((rank == 2) && (move_ticks <= 20)) {
                        std::cout << "rank: " << rank << " POST_NUL_VACREADBACK move_ticks: " << move_ticks
                                  << " vacancies(1,54,0,119): " << vacancies(1,54,0,119)
                                  << " vacancies(0,54,1,119): " << vacancies(0,54,1,119) << "\n";
                    }
                    if ((rank == 0) && (move_ticks <= 20)) {
                        std::cout << "rank: " << rank << " POST_NUL_VACREADBACK move_ticks: " << move_ticks
                                  << " proc_pos_y_neighbors(1,56,1,119): " << proc_pos_y_neighbors(1,56,1,119) << "\n";
                    }
                    move_counts[moves_lattice[rand_idx][0]] ++;
                    timestep = new_random_times();
                    time_count[moves_lattice[rand_idx][0]] += timestep;
                    store_time_incr(timestep);
                }
                else {
                    // std::cout << "rank: " << rank << " NO RATES\n";
                    communicate_rates(); 
                    timestep = new_random_times();
                }

                MPI_Request request1;
                for (int new_proc=0; new_proc<num_procs; new_proc++) { 
                    if (new_proc != rank) {
                        MPI_Isend( NULL, 0, MPI_CHAR, new_proc, par_done_tag, MPI_COMM_WORLD, &request1); 
                        MPI_Wait(&request1, MPI_STATUS_IGNORE);
                    }
                }
                MPI_Request request2;
                for (int new_proc=0; new_proc<num_procs; new_proc++) { 
                    if (new_proc != rank) {
                        MPI_Isend( NULL, 0, MPI_CHAR, new_proc, ghost_done_tag, MPI_COMM_WORLD, &request2); 
                        MPI_Wait(&request2, MPI_STATUS_IGNORE);
                    }
                }
                
                MPI_Barrier(MPI_COMM_WORLD);
                
                receive_parallel_comm_helper(move_ticks);
                if ((rank == 2) && (move_ticks <= 25)) {
                    std::cout << "rank: " << rank << " POST_RPCH_READBACK move_ticks: " << move_ticks
                              << " proc_neg_y_neighbors(1,49,1,124): " << proc_neg_y_neighbors(1,49,1,124) << "\n";
                }
                if ((rank == 2) && (move_ticks <= 20)) {
                    std::cout << "rank: " << rank << " POST_RPCH_VACREADBACK move_ticks: " << move_ticks
                              << " vacancies(1,54,0,119): " << vacancies(1,54,0,119)
                              << " vacancies(0,54,1,119): " << vacancies(0,54,1,119) << "\n";
                }
                if ((rank == 0) && (move_ticks <= 20)) {
                    std::cout << "rank: " << rank << " POST_RPCH_VACREADBACK move_ticks: " << move_ticks
                              << " proc_pos_y_neighbors(1,56,1,119): " << proc_pos_y_neighbors(1,56,1,119) << "\n";
                }

                MPI_Barrier(MPI_COMM_WORLD);

                if (move_ticks % write_frames == 0) {
                    
                    if (rank == 0) { 
                        std::cout << "rank: " << rank << " move_ticks: " << move_ticks << "\n"; 
                        std::cout << "rank: " << rank << " curr_num_vacs: " << curr_num_vacs << "\n";
                        std::cout << "rank: " << rank << " prev_num_vacs: " << prev_num_vacs << "\n";
                    }

                    int reconstruct_interval = std::max(1, (int)(0.1*write_frames));
                    if (move_ticks % reconstruct_interval == 0) {reconstruct = true;}
                    else {reconstruct = false;}

                    // Reassigns the outer only_vacancies (declared once before the while loop) --
                    // NOT a fresh "Matrix<int> only_vacancies = ...", which would instead declare a
                    // same-named local that shadows the outer one for this block only, leaving the
                    // outer variable frozen at its pre-loop (near-initial-configuration) value.
                    // check_ghost_sites_correctness(only_vacancies, ...) below runs after this
                    // block closes, so it was always comparing live ghost state against that frozen
                    // snapshot instead of the current tick's real vacancy configuration.
                    only_vacancies = vacancies.nonzero(rank); // configuration of vacancies at current timestep
                    std::cout << "rank: " << rank << " vac_nonzero.rows(): " << only_vacancies.rows() << " vacancies_pos.rows(): " << vacancies_pos.rows() << "\n";
                    
                    if ((only_vacancies.rows() != vacancies_pos.rows()) && (rank == 0)) { std::cout << "rank: " << rank << " ERROR: mismatch in number of vacancies in vacancy list and nonzero elements of array"; }
                    
                    
                    curr_num_vacs = sum_vacs_allprocs(only_vacancies, proc_dims[0], proc_dims[1]);

                    if (rank == 0) std::cout << "rank: " << rank << " filewrite post sum_vacs_allprocs curr_num_vacs: " << curr_num_vacs << "\n";

                    total_cost_allprocs = sum_values_allprocs(total_cost);
                    total_E_allprocs = sum_values_allprocs(system_energy);

                    if (rank == 0) {
                        std::cout << " rank: " << rank << "total_cost_allprocs: " << total_cost_allprocs << "\n";
                        std::cout << " rank: " << rank << "total_E_allprocs: " << total_E_allprocs << "\n";
                        std::cout  <<" rank: " << rank << "total E + cost: " << total_E_allprocs + total_cost_allprocs << "\n";
                    }

                    int rates_i = 0;
                    // write_output_parallel(only_vacancies, folder, proc_dims[0], proc_dims[1], iteration, rates_i, move_ticks, curr_num_vacs, reconstruct, get_rank);
                    
                }  
                
                last_rate_plus1 = rate_cumsum[(rand_idx+1)] - rate_cumsum[(rand_idx)];
                last_rate = rate_cumsum[(rand_idx)] - rate_cumsum[(rand_idx-1)];
                last_rate_minus1 = rate_cumsum[(rand_idx-1)] - rate_cumsum[(rand_idx-2)]; 


                std::cout << "rank: " << rank << " last_rate + 1: " << last_rate_plus1 << "\n";
                std::cout << "rank: " << rank << " last_rate: " << last_rate << "\n";
                std::cout << "rank: " << rank << " last_rate - 1: " << last_rate_minus1 << "\n";
                
                parallel_get_actions();
                
                fflush(stdout);
                
                std::cout << "checking: num vacs \n"; 
                if ((rank == 0) && (prev_num_vacs != curr_num_vacs)) {
                    std::cout << "Change in total number of vacancies in simulation -- curr_num_vacs: " << curr_num_vacs << " prev_num_vacs: " << prev_num_vacs << " move_ticks: " << move_ticks << "\n"; 
                    if ((move_ticks != 0) && (restart == false)) {
                        MPI_Barrier(MPI_COMM_WORLD);
                        exit(0);
                    }
                }
                
                /*
                check here for ghost sites equivalence with reconstructor
                */
                std::cout << "rank: " << rank << " move_ticks: " << move_ticks << "\n";
                check_ghost_sites_correctness(only_vacancies, proc_dims[0], proc_dims[1]);

                if (restart == true) restart = false;


                if (rank == 0 ) {
                    std::cout << "pre prev_num_vacs: " << prev_num_vacs << "\n";
                    std::cout << "pre curr_num_vacs: " << curr_num_vacs << "\n";
                }    

                elapsed_seconds = end-start;
                t += timestep;
                prev_num_vacs = curr_num_vacs;

                if (rank == 0 ) {
                    std::cout << "post prev_num_vacs: " << prev_num_vacs << "\n";
                    std::cout << "post curr_num_vacs: " << curr_num_vacs << "\n";
                }    
                move_ticks ++;
                old_time = t;

            }

            std::cout << "rank: " << rank << " elapsed_seconds: " << elapsed_seconds.count() << "\n";
            std::cout << "rank: " << rank << " exiting\n";
            std::cout << "rank: " << rank << " t: " << t << "\n";
            std::cout << "rank: " << rank << " move_ticks: " << move_ticks << "\n";

            std::vector<int> move_counts_sum = sum_vectors_allprocs(move_counts, num_procs, rank);
            std::vector<double> time_count_sum = sum_vectors_allprocs(time_count, num_procs, rank);
            

            if (rank == 0) {
                std::cout << "rank: " << rank << " move_counts_sum: \n";
                print_1Dvector(move_counts_sum);
                std::cout << "rank: " << rank << " time_count: \n";
                print_1Dvector(time_count);
            }

            t = time_lim + 1;
            fflush(stdout);
            MPI_Barrier(MPI_COMM_WORLD);
            all_times.push_back(t);

            std::cout << "rank: " << rank << " move_ticks: " << move_ticks << "\n";
            lattice_return_struct output_vals(move_counts, time_count, all_times);

            return output_vals;
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
        const void write_output_parallel(const Matrix<int>& only_vacancies, std::string folder, int xprocs, int yprocs, int i, int l, int k, int curr_num_vacs,  bool reconstruct, bool get_rank = false) {
            
            // Initialize variables
            int size = (int)( only_vacancies.rows() * only_vacancies.cols() );
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
            std::vector<int> receive_counts(num_procs);
            std::vector<int> receive_displacements(num_procs, 0);
            std::vector<int> sum_vec(num_procs);

            int num_elems = (int)(only_vacancies.rows() * 4);
            nums_and_proc[0] = num_elems;
            nums_and_proc[1] = rank;
            std::cout << "rank: " << rank << " only_vacancies.rows(): " << only_vacancies.rows() << "\n";
            MPI_Barrier(MPI_COMM_WORLD);
            int sum_of_elems = 0;

            reconstruct = false;
            
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
                for (int rec_proc = 0; rec_proc < num_procs; rec_proc++) {
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
            MPI_Gatherv(only_vacancies.data(), (only_vacancies.rows() * 4), MPI_INT, vacs_in.data(), receive_counts.data(), receive_displacements.data(), MPI_INT, 0, MPI_COMM_WORLD);

            MPI_Barrier(MPI_COMM_WORLD);

            int shift_idx = 0;
            int coords_size = (int)(vacs_in.size() / 4);
            
            if (rank == 0) {
                std::cout << "writing move_ticks: " << k << "\n";

                // Initialize output matrices for vacancies and process ranks
                FourDBoolArr vacancies_out((size_t)2, (size_t)total_dims[0], (size_t)total_dims[1], (size_t)total_dims[2]);
                FourDArr proc_rank((size_t)2, (size_t)total_dims[0], (size_t)total_dims[1], (size_t)total_dims[2]);
                vacancies_out.zero();
                proc_rank.zero();
                
                // Print the sum vector
                // std::cout << "sum_vec\n";
                // print_1Dvector(sum_vec);

                // Process data and assign to output arrays
                for (int idx = 0; idx < coords_size; idx++) {
                    while ((shift_idx < num_procs) && (idx * 4 >= sum_vec[shift_idx])) { 
                        // std::cout << "rank: " << rank << " pre shift_idx: " << shift_idx << " idx: " << idx << "\n";
                        shift_idx++; 
                        // std::cout << "rank: " << rank << " post shift_idx: " << shift_idx << " idx: " << idx << "\n\n";
                    }
                
                    x_idx = shift_idx % xprocs;
                    y_idx = floor(shift_idx / xprocs);
                    
                    x_chunk_start = (int)(total_dims[0] / xprocs * x_idx);
                    y_chunk_start = (int)(total_dims[1] / yprocs * y_idx); 

                    i1 = vacs_in[4 * idx];
                    i2 = vacs_in[4 * idx + 1] + x_chunk_start;
                    i3 = vacs_in[4 * idx + 2] + y_chunk_start;
                    i4 = vacs_in[4 * idx + 3];
                    // std::cout << "idx: " << idx << " x_idx: " << x_idx << " y_idx: " << y_idx << "\n\n";
                    if (!get_rank) {
                        vacancies_out(i1, i2, i3, i4) = 1;
                    } else {
                        proc_rank(i1, i2, i3, i4) = shift_idx + 1;
                    }
                }

                std::string output_filename = ss.str();
                ss << folder << "/vacs/vacancies_output_" << i << "_" << l << "_" << k << "_" << t << "_moves.txt";
                //std::cout << "filename: " << ss.str() << "\n\n";

                // Write the vacancies or vacancies and rank to file
                if (!get_rank) {
                    Matrix<int> all_vacancies = vacancies_out.nonzero(rank);
                    write_to_file(ss.str(), all_vacancies);
                    //std::cout << "all_vacancies.rows(): " << all_vacancies.rows() << " curr_num_vacs: " << curr_num_vacs << "\n"; 
                    
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
                    //std::cout << "vacancies_and_rank.rows(): " << vacancies_and_rank.rows() << " curr_num_vacs: " << curr_num_vacs << "\n"; 
                    
                    if (vacancies_and_rank.rows() != curr_num_vacs) {
                        std::cout << "Change in total number of vacancies in simulation -- vacancies_and_rank.rows(): " << vacancies_and_rank.rows() << " curr_num_vacs: " << curr_num_vacs << "\n"; 
                        if (k != 0) {
                            exit(0);
                        }
                    }
                }                
                if (reconstruct == true) {
                    reconstruct_ghost_sites(only_vacancies, xprocs, yprocs);
                }                

                ss.str("");
                ss.clear();               
            }

            MPI_Barrier(MPI_COMM_WORLD);  // Synchronize all processes
        }

        /* \brief Sums up the vacancies across all processes and returns the total number of vacancies.
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
        const int sum_vacs_allprocs(const Matrix<int>& only_vacancies, int xprocs, int yprocs) {
            // Initialize variables
            int size = (int)( only_vacancies.rows() * only_vacancies.cols() );  ///< Size of the vacancy matrix
            std::vector<int> output(size);
            int idx = 0;
            int size2 = 0, size3 = 0, size_in;
            int i1, i2, i3, i4, x_idx, y_idx, x_chunk_start, y_chunk_start;
            std::vector<int> vacs_in;

            std::cout << "rank: " << rank << " pre sum_vacs_allprocs only_vacancies.rows(): " << only_vacancies.rows() << "\n";

            MPI_Status status;
            MPI_Request request;
            std::ostringstream ss;

            // MPI communication buffers
            std::vector<int> nums_and_proc(2);  ///< Number of elements and process rank
            std::vector<int> num_proc_buffer(2);  ///< Buffer for receiving data from other processes
            std::vector<int> receive_counts(num_procs);  ///< Number of elements to receive from each process
            std::vector<int> receive_displacements(num_procs, 0);  ///< Displacements of received data
            std::vector<int> sum_vec(num_procs);  ///< Running sum of received data sizes

            int num_elems = (int)(only_vacancies.rows() * 4);  ///< Total number of elements in the vacancy matrix
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
                for (int rec_proc = 0; rec_proc < num_procs; rec_proc++) {
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
            MPI_Gatherv(only_vacancies.data(), (only_vacancies.rows() * 4), MPI_INT, vacs_in.data(), receive_counts.data(), receive_displacements.data(), MPI_INT, 0, MPI_COMM_WORLD);

            MPI_Barrier(MPI_COMM_WORLD);  ///< Synchronize all processes

            int shift_idx = 0;
            int coords_size = (int)(vacs_in.size() / 4);

            // Root process (rank 0) processes the gathered data
            if (rank == 0) {
                FourDBoolArr vacancies_out((size_t)2, (size_t)total_dims[0], (size_t)total_dims[1], (size_t)total_dims[2]);  ///< Output array for vacancies
                vacancies_out.zero();  ///< Initialize the output array

                //std::cout << "sum_vec\n";
                //print_1Dvector(sum_vec);

                // Process each vacancy and update the output array
                for (int idx = 0; idx < coords_size; idx++) {
                    while ((shift_idx < num_procs) && (idx * 4 >= sum_vec[shift_idx])) { 
                        //std::cout << "rank: " << rank << " sum_vacs_all_procs() pre shift_idx: " << shift_idx << " idx: " << idx << "\n";
                        shift_idx++; 
                        //std::cout << "rank: " << rank << " sum_vacs_all_procs() post shift_idx: " << shift_idx << " idx: " << idx << "\n\n";
                    }
                
                    x_idx = shift_idx % xprocs;
                    y_idx = floor(shift_idx / xprocs);
                    
                    x_chunk_start = (int)(total_dims[0] / xprocs * x_idx);
                    y_chunk_start = (int)(total_dims[1] / yprocs * y_idx); 

                    i1 = vacs_in[4 * idx];
                    i2 = vacs_in[4 * idx + 1] + x_chunk_start;
                    i3 = vacs_in[4 * idx + 2] + y_chunk_start;
                    i4 = vacs_in[4 * idx + 3];
                    
                    vacancies_out(i1, i2, i3, i4) = 1;
                }

                // Return the total number of non-zero vacancies
                Matrix<int> all_vacancies = vacancies_out.nonzero(rank);
                std::cout << "rank: " << rank << " post sum_vacs_allprocs all_vacancies.rows(): " << all_vacancies.rows() << "\n";
                return (int)all_vacancies.rows();
            }

            return 0; // In case of other ranks, no output is generated.
        }


        /**
        * @brief Reconstructs ghost sites for a lattice in a parallel simulation.
        *
        * This function identifies and updates ghost sites based on vacancies within 
        * the lattice, ensuring periodic boundary conditions and communication 
        * between neighboring processes in a distributed simulation.
        *
        * @param only_vacancies A matrix containing vacancy information. Each row specifies 
        *                       a vacancy's lattice position and coordinates (x, y, z).
        * @param xprocs         The number of processes in the x-direction.
        * @param yprocs         The number of processes in the y-direction.
        */
        const void reconstruct_ghost_sites(const Matrix<int>& only_vacancies, int xprocs, int yprocs) {
            int lattice_pos; ///< Lattice position of the vacancy.
            int w, x, y, z; ///< Coordinates of the vacancy in the lattice.
            int x_idx, y_idx, z_idx; ///< Indices for the vacancy position within process-local arrays.
            int x_unmod, y_unmod, z_unmod;

            // Compute periodic boundary indices for x and y directions.
            int xlo_edge = (((chunk_bounds[0][0] - 1) % proc_dims[0] + proc_dims[0]) % proc_dims[0]);
            int xhi_edge = (((chunk_bounds[0][1]) % proc_dims[0] + proc_dims[0]) % proc_dims[0]);
            int ylo_edge = (((chunk_bounds[1][0] - 1) % proc_dims[1] + proc_dims[1]) % proc_dims[1]);
            int yhi_edge = (((chunk_bounds[1][1]) % proc_dims[1] + proc_dims[1]) % proc_dims[1]);

            // Size vectors for neighboring processes in positive and negative directions.
            std::vector<size_t> x_dims = proc_neg_x_neighbors.size_vec;
            std::vector<size_t> y_dims = proc_neg_y_neighbors.size_vec;
            std::vector<size_t> x_dims_pos = proc_pos_x_neighbors.size_vec;
            std::vector<size_t> y_dims_pos = proc_pos_y_neighbors.size_vec;

            // Loop through all vacancies and process ghost site conditions.
            /*for (int vac_idx = 0; vac_idx < (int)only_vacancies.rows(); vac_idx++) {
                lattice_pos = only_vacancies[vac_idx][0];
                x = only_vacancies[vac_idx][1];
                y = only_vacancies[vac_idx][2];
                z = only_vacancies[vac_idx][3];

                if ((xprocs != 1) || (yprocs != 1)) {
                    // Handle ghost sites for first layer (111 and 100 moves).
                    if ((y >= (chunk_bounds[1][0] - 1)) && (y < (chunk_bounds[1][1]))) {
                        if ((x == xlo_edge) && (lattice_pos == 1)) {
                            y_idx = mod_with_bounds((y - chunk_bounds[1][0] + 1), proc_dims[1]);
                            proc_neg_x_neighbors(0, 1, y_idx, z_idx) = 1;

                            // Preserve corner periodic boundary conditions.
                            if ((y_idx == x_dims[2] - 1) && (proc_neighbors(rank, 2) == rank)) {
                                proc_neg_x_neighbors(0, 1, 0, z_idx) = 1;
                            }

                            if (y == ylo_edge) {
                                x_idx = mod_with_bounds((x - chunk_bounds[0][0] + 1), proc_dims[0]);
                                proc_neg_y_neighbors(0, 1, x_idx, z_idx) = 1;
                            }
                        }
                    }
                    // Further handling for the x and y edges...
                    if ((x >= (chunk_bounds[0][0])) && (x < (chunk_bounds[0][1] + 1))) {
                        if ((y == yhi_edge) && (lattice_pos == 1)) {
                            x_idx = mod_with_bounds((x - chunk_bounds[0][0]), proc_dims[0]);
                            proc_pos_y_neighbors(0, 1, x_idx, z_idx) = 1;

                            // Preserve corner periodic boundary conditions.
                            if ((x_idx == 0) && (proc_neighbors(rank, 0) == rank)) {
                                proc_pos_y_neighbors(0, 1, y_dims_pos[2] - 1, z_idx) = 1;
                            }
                        }
                    }

                    // Handle ghost sites for the second layer (100 moves).
                    if ((y >= (chunk_bounds[1][0] - 1)) && (y < (chunk_bounds[1][1]))) {
                        if ((x == xlo_edge) && (lattice_pos == 0)) {
                            y_idx = mod_with_bounds((y - chunk_bounds[1][0] + 1), proc_dims[1]);
                            proc_neg_x_neighbors(0, 0, y_idx, z_idx) = 1;

                            // Further handling for boundary conditions...
                        }
                    }
                }
            }*/

            for (int vac_idx = 0; vac_idx < (int)only_vacancies.rows(); vac_idx++) {
                lattice_pos = only_vacancies[vac_idx][0];
                x_unmod = only_vacancies[vac_idx][1];
                y_unmod = only_vacancies[vac_idx][2];
                z_unmod = only_vacancies[vac_idx][3];
                                
                w = lattice_pos;
                x = mod_with_bounds((x_unmod), total_dims[0]);
                y = mod_with_bounds((y_unmod), total_dims[1]);
                z = mod_with_bounds((z_unmod), total_dims[2]);
                                
                /* first layer ghost sites for (111) and (100) moves */
                if ((ylo_edge > yhi_edge)) {
                    //if (rank == 3) std::cout << "rank: " << rank << " ylo_edge > yhi_edge\n";
                    if (((y_unmod >= 0) && (y_unmod < yhi_edge)) || 
                        ((y_unmod < total_dims[1]) && (y_unmod >= ylo_edge))) {

                        //if (rank == 3) std::cout << "rank: " << rank << " in y bounds\n";

                        if (w == 1) { y_idx = mod_with_bounds((y - chunk_bounds[1][0] + 2), (total_dims[1])); }
                        else { y_idx = mod_with_bounds((y - chunk_bounds[1][0] + 2), (total_dims[1])); }
                        x_idx = x_unmod % 2;
                        //if (rank == 3) std::cout << "rank: " << rank << " w: " << w << " x_idx: " << x_idx << " y_idx: " << y_idx << " z: " << z << "\n"; 

                        if ((x == xlo_edge) || (x == (xlo_edge +1))) {   
                            //if (rank == 3) std::cout << "rank: " << rank << " case 1 neg x neigh\n";
                            proc_neg_x_neighbors(w,x_idx,y_idx,z) = 1;
                        }
                        else if ((x == (xhi_edge-1)) || (x == (xhi_edge -2))) {    
                            //if (rank == 3) std::cout << "rank: " << rank << " case 1 pos x neigh\n";                        
                            proc_pos_x_neighbors(w,x_idx,y_idx,z) = 1;
                        }
                        
                    }
                }
                else {
                    //if (rank == 3) std::cout << "rank: " << rank << " ylo_edge <= yhi_edge\n";
                    if ((y_unmod < yhi_edge)  && (y_unmod >= ylo_edge)) {
                        //if (rank == 3) std::cout << "rank: " << rank << " in y bounds\n";

                        if (w == 1) { y_idx = mod_with_bounds((y - chunk_bounds[1][0] + 2), (total_dims[1])); }
                        else { y_idx = mod_with_bounds((y - chunk_bounds[1][0] + 2), (total_dims[1])); }
                        x_idx = x_unmod % 2;
                        //if (rank == 3) std::cout << "rank: " << rank << " w: " << w << " x_idx: " << x_idx << " y_idx: " << y_idx << " z: " << z << "\n";

                        if ((x == xlo_edge) || (x == (xlo_edge +1))) {    
                            //if (rank == 3) std::cout << "rank: " << rank << " case 2 neg x neigh\n";
                            
                            proc_neg_x_neighbors(w,x_idx,y_idx,z) = 1;
                        }
                        else if ((x == (xhi_edge-1)) || (x == (xhi_edge -2))) {  
                            //if (rank == 3) std::cout << "rank: " << rank << " case 2 pos x neigh\n";
                            
                            proc_pos_x_neighbors(w,x_idx,y_idx,z) = 1;
                        }
                        
                    }
                }
                
                if ((xlo_edge > xhi_edge)) {
                    //if (rank == 3) std::cout << "rank: " << rank << " xlo_edge > xhi_edge\n";
                    if (((x_unmod >= 0) && (x_unmod < xhi_edge)) || 
                        ((x_unmod < total_dims[0]) && (x_unmod >= xlo_edge))) {

                        //if (rank == 3) std::cout << "rank: " << rank << " in x bounds\n";

                        if (w == 1) { x_idx = mod_with_bounds((x - chunk_bounds[0][0] + 2), (total_dims[0])); }
                        else { x_idx = mod_with_bounds((x - chunk_bounds[0][0] + 2), (total_dims[0])); }
                        y_idx = y_unmod % 2;
                        //if (rank == 3) std::cout << "rank: " << rank << " w: " << w << " x_idx: " << x_idx << " y_idx: " << y_idx << " z: " << z << "\n";

                        if ((y == ylo_edge) || (y == (ylo_edge + 1))) {    
                            //if (rank == 3) std::cout << "rank: " << rank << " case 1 neg y neigh\n";
                            proc_neg_y_neighbors(w,x_idx,y_idx,z) = 1;
                        }
                        else if ((y == (yhi_edge-1)) || (y == (yhi_edge -2))) {    
                            //if (rank == 3) {
                            //    std::cout << "rank: " << rank << " case 1 pos y neigh\n";
                            //   std::cout << "rank: " << rank << " [ " << w << " " << x_idx << " " << y_idx <<  " " << z << " ]\n";
                            //}
                            proc_pos_y_neighbors(w,x_idx,y_idx,z) = 1;
                        }
                        
                    }
                }
                else {
                    //if (rank == 3) std::cout << "rank: " << rank << " xlo_edge <= xhi_edge\n";
                    if ((x_unmod < xhi_edge)  && (x_unmod >= xlo_edge)) {
                        
                        //if (rank == 0) std::cout << "rank: " << rank << " in x bounds\n";

                        if (w == 1) { x_idx = mod_with_bounds((x - chunk_bounds[0][0] + 2), (total_dims[0])); }
                        else { x_idx = mod_with_bounds((x - chunk_bounds[0][0] + 2), (total_dims[0])); }
                        y_idx = y_unmod % 2;
                        //if (rank == 3) std::cout << "rank: " << rank << " w: " << w << " x_idx: " << x_idx << " y_idx: " << y_idx << " z: " << z << "\n"; 

                        if ((y == ylo_edge) || (y == (ylo_edge +1))) {   
                            //if (rank == 3) std::cout << "rank: " << rank << " case 2 neg y neigh\n";
                            proc_neg_y_neighbors(w,x_idx,y_idx,z) = 1;
                        }
                        else if ((y == (yhi_edge-1)) || (y == (yhi_edge -2))) {
                            //if (rank == 3) std::cout << "rank: " << rank << " case 2 pos y neigh\n";
                            proc_pos_y_neighbors(w,x_idx,y_idx,z) = 1;
                        }                    
                    }
                } 
            
                
            }  
        }



        /**
        * @brief Determines which ghost-tracking array (and index) a lattice position maps to,
        * purely from its coordinates relative to this process's chunk boundaries, and writes
        * `value` into it (1 to mark a position as a ghost site, 0 to clear it).
        *
        * IMPORTANT: x_unmod/y_unmod/z_unmod must be GLOBAL lattice coordinates, not local
        * per-processor ones -- this mirrors populate_lattice's own edge-detection convention
        * (chunk_bounds +/-2, wrapped against the global total_dims), which is the canonical rule
        * for what counts as being within the 2-layer ghost margin around this process's chunk.
        * Callers holding a local coordinate must add chunk_bounds[dim][0] before calling this.
        *
        * This is the single source of truth for the coordinate-to-ghost-array mapping: it does
        * not look up any rank or direction index, so it cannot suffer from the ambiguity that
        * occurs when a small process grid (e.g. 2x2 with periodic boundaries) causes the same
        * neighboring rank to occupy more than one of the 8 compass-direction slots. Used both to
        * reconstruct ground truth (check_ghost_sites_correctness) and to directly update the live
        * ghost-tracking arrays (ghost_site_self_reference), matching what populate_lattice does
        * for the initial seeding of these same arrays.
        */
        // trace_tag/trace_remote_proc: diagnostic-only source tracking for the live ghost arrays.
        // Pass a non-empty trace_tag to have every write this call makes logged as GHOST_WRITE
        // (rank, ghost_trace_tick, which array/cell, old/new value, and who wrote it). Leave
        // trace_tag empty (the default) for non-diagnostic callers -- e.g. check_ghost_sites_correctness's
        // ground-truth rebuild, which writes into separate "_new" arrays, not these live ones.
        // trace_remote_proc should be -1 for a self-originated write (this rank's own move) or the
        // sending rank for a write triggered by a received ghost update.
        void set_ghost_position(int w, int x_unmod, int y_unmod, int z_unmod, int value,
                                 FourDArr& neg_x_arr, FourDArr& pos_x_arr,
                                 FourDArr& neg_y_arr, FourDArr& pos_y_arr,
                                 const std::string& trace_tag = "", int trace_remote_proc = -1) {
            set_ghost_position_impl(w, x_unmod, y_unmod, z_unmod, value, total_dims, chunk_bounds,
                                     neg_x_arr, pos_x_arr, neg_y_arr, pos_y_arr,
                                     rank, ghost_trace_tick, trace_remote_proc, trace_tag);
        }

        /**
        * @brief Ground-truth reconstruction of the ghost-tracking arrays.
        *
        * proc_neg_x_neighbors/proc_pos_x_neighbors/proc_neg_y_neighbors/proc_pos_y_neighbors mirror
        * ADJACENT processes' vacancies, not this rank's own -- so a valid reconstruction can't rely
        * on only_vacancies (this rank's own local vacancy list) alone. This builds a dense global
        * (w,x,y,z) vacancy array by having each rank fill in its own chunk (converted to global
        * coordinates) and combining every rank's contribution via MPI_Allreduce with MPI_SUM -- safe
        * since chunks don't overlap, so at most one rank ever contributes a 1 for a given site. This
        * is a collective call that all ranks must reach together (true here, since every rank calls
        * this once per iteration at the same point). It then loops over every site in the combined
        * global array and rebuilds the ghost arrays with set_ghost_position -- which, given any
        * global coordinate, already correctly determines (using this rank's own chunk_bounds
        * internally) whether and where it belongs in this rank's own ghost arrays, regardless of
        * which rank the vacancy actually came from.
        */
        const void check_ghost_sites_correctness(const Matrix<int>& only_vacancies, int xprocs, int yprocs) {
            size_t global_size = (size_t)2 * total_dims[0] * total_dims[1] * total_dims[2];
            std::vector<int> my_global_vacancies(global_size, 0);

            for (int i = 0; i < (int)only_vacancies.rows(); i++) {
                int w = only_vacancies[i][0];
                int x = only_vacancies[i][1] + chunk_bounds[0][0];
                int y = only_vacancies[i][2] + chunk_bounds[1][0];
                int z = only_vacancies[i][3] + chunk_bounds[2][0];
                if ((w < 0) || (w >= 2) || (x < 0) || (x >= total_dims[0]) ||
                    (y < 0) || (y >= total_dims[1]) || (z < 0) || (z >= total_dims[2])) {
                    std::cout << "CGSC OUT_OF_RANGE rank: " << rank << " local: (" << only_vacancies[i][0] << "," << only_vacancies[i][1] << "," << only_vacancies[i][2] << "," << only_vacancies[i][3] << ") global: (" << w << "," << x << "," << y << "," << z << ") chunk_bounds[0][0]: " << chunk_bounds[0][0] << " chunk_bounds[1][0]: " << chunk_bounds[1][0] << "\n";
                    continue;
                }
                size_t idx = (((size_t)w * total_dims[0] + x) * total_dims[1] + y) * total_dims[2] + z;
                my_global_vacancies[idx] = 1;
            }

            std::vector<int> global_vacancies(global_size);
            MPI_Allreduce(my_global_vacancies.data(), global_vacancies.data(), (int)global_size,
                          MPI_INT, MPI_SUM, MPI_COMM_WORLD);

            // Size vectors for neighboring processes in positive and negative directions.
            std::vector<size_t> x_dims = proc_neg_x_neighbors.size_vec;
            std::vector<size_t> y_dims = proc_neg_y_neighbors.size_vec;
            std::vector<size_t> x_dims_pos = proc_pos_x_neighbors.size_vec;
            std::vector<size_t> y_dims_pos = proc_pos_y_neighbors.size_vec;

            // FourDArr's backing store is malloc'd (not zero-initialized like the old
            // FourDBoolArr), so each must be explicitly zeroed before accumulating into it below.
            FourDArr proc_neg_x_neighbors_new(x_dims[0], x_dims[1], x_dims[2], x_dims[3]);
            FourDArr proc_neg_y_neighbors_new(y_dims[0], y_dims[1], y_dims[2], y_dims[3]);
            FourDArr proc_pos_x_neighbors_new(x_dims_pos[0], x_dims_pos[1], x_dims_pos[2], x_dims_pos[3]);
            FourDArr proc_pos_y_neighbors_new(y_dims_pos[0], y_dims_pos[1], y_dims_pos[2], y_dims_pos[3]);
            proc_neg_x_neighbors_new.zero();
            proc_neg_y_neighbors_new.zero();
            proc_pos_x_neighbors_new.zero();
            proc_pos_y_neighbors_new.zero();

            for (int w = 0; w < 2; w++) {
                for (int x = 0; x < total_dims[0]; x++) {
                    for (int y = 0; y < total_dims[1]; y++) {
                        for (int z = 0; z < total_dims[2]; z++) {
                            size_t idx = (((size_t)w * total_dims[0] + x) * total_dims[1] + y) * total_dims[2] + z;
                            if ((w==0) && (x==55) && (y==44) && (z==127)) {
                                std::cout << "CGSC READ rank: " << rank << " idx: " << idx << " global_vacancies[idx]: " << global_vacancies[idx] << " my_global_vacancies[idx]: " << my_global_vacancies[idx] << "\n";
                            }
                            // Investigating a 10-tick (1421-1430) ghost mismatch on rank3's
                            // proc_neg_x_neighbors(1,0,7,119) -- decodes to global (w=1,x=54,y=60,
                            // z=119). No traced write ever touches that live cell during the whole
                            // run, yet ground truth (global_vacancies[idx] here) says 1 for those
                            // ticks. Printing my_global_vacancies[idx] on every rank tells us
                            // whether two ranks are simultaneously counting the same physical
                            // vacancy into this global index (a real double-count) rather than the
                            // live ghost array being wrong.
                            if ((w==1) && (x==54) && (y==60) && (z==119)) {
                                std::cout << "CGSC_TARGET rank: " << rank << " idx: " << idx << " global_vacancies[idx]: " << global_vacancies[idx] << " my_global_vacancies[idx]: " << my_global_vacancies[idx] << "\n";
                            }
                            if (global_vacancies[idx] != 0) {
                                set_ghost_position(w, x, y, z, 1,
                                                    proc_neg_x_neighbors_new, proc_pos_x_neighbors_new,
                                                    proc_neg_y_neighbors_new, proc_pos_y_neighbors_new);
                            }
                        }
                    }
                }
            }

            std::cout << "rank: " << rank << " checking proc_neg_x_neighbors:\n";
            proc_neg_x_neighbors.check_equal(proc_neg_x_neighbors_new, rank);
            std::cout << "rank: " << rank << " checking proc_neg_y_neighbors:\n";
            proc_neg_y_neighbors.check_equal(proc_neg_y_neighbors_new, rank);
            std::cout << "rank: " << rank << " checking proc_pos_x_neighbors:\n";
            proc_pos_x_neighbors.check_equal(proc_pos_x_neighbors_new, rank);
            std::cout << "rank: " << rank << " checking proc_pos_y_neighbors:\n";
            proc_pos_y_neighbors.check_equal(proc_pos_y_neighbors_new, rank);

        }

};

/**
 * @brief Determines which ghost-tracking array (and index) a lattice position maps to,
 * purely from its coordinates relative to a process's chunk boundaries, and writes
 * `value` into it (1 to mark a position as a ghost site, 0 to clear it).
 *
 * IMPORTANT: x_unmod/y_unmod/z_unmod must be GLOBAL lattice coordinates, not local
 * per-processor ones -- this mirrors the edge-detection convention used everywhere
 * this is called from (chunk_bounds +/-2, wrapped against the global total_dims),
 * which is the canonical rule for what counts as being within the 2-layer ghost
 * margin around a process's chunk. Callers holding a local coordinate must add
 * chunk_bounds[dim][0] before calling this.
 *
 * This is the single source of truth for the coordinate-to-ghost-array mapping: it
 * does not look up any rank or direction index, so it cannot suffer from the
 * ambiguity that occurs when a small process grid (e.g. 2x2 with periodic
 * boundaries) causes the same neighboring rank to occupy more than one of the 8
 * compass-direction slots. Shared by Lattice::set_ghost_position (used by
 * check_ghost_sites_correctness and ghost_site_self_reference/ghost_site_recieve)
 * and by populate_lattice's initial seeding of these same arrays, so the two can't
 * silently drift apart.
 */
void set_ghost_position_impl(int w, int x_unmod, int y_unmod, int z_unmod, int value,
                              const std::vector<int>& total_dims,
                              const std::vector<std::vector<int>>& chunk_bounds,
                              FourDArr& neg_x_arr, FourDArr& pos_x_arr,
                              FourDArr& neg_y_arr, FourDArr& pos_y_arr,
                              int trace_rank, int trace_tick, int trace_remote_proc,
                              const std::string& trace_tag) {
    // Diagnostic only (source-tracking for the ghost-array UNEQUAL investigation): when
    // trace_tag is non-empty, log every actual write this call makes to the live ghost arrays,
    // including the value the cell held immediately beforehand and who's writing it (this rank,
    // the tick, and -1/remote-rank for self-originated vs received-from-another-rank writes).
    // This lets us grep a specific (arr,w,x_idx,y_idx,z) cell's write history across ranks to see
    // who's contributing to/withdrawing from a given cell over time.
    bool do_trace = !trace_tag.empty();
    int x = mod_with_bounds(x_unmod, total_dims[0]);
    int y = mod_with_bounds(y_unmod, total_dims[1]);
    int z = mod_with_bounds(z_unmod, total_dims[2]);

    int xlo_edge = (((chunk_bounds[0][0]-2) % total_dims[0] + total_dims[0]) % total_dims[0]);
    int xhi_edge = (((chunk_bounds[0][1]+2) % total_dims[0] + total_dims[0]) % total_dims[0]);
    int ylo_edge = (((chunk_bounds[1][0]-2) % total_dims[1] + total_dims[1]) % total_dims[1]);
    int yhi_edge = (((chunk_bounds[1][1]+2) % total_dims[1] + total_dims[1]) % total_dims[1]);

    int x_idx, y_idx;

    // value==1 adds this call's contribution to the cell's reference count; value==0 removes it
    // (clamped at 0, so an errant extra "remove" -- e.g. the duplicate-message pattern already
    // observed in these logs -- can't push a cell negative and corrupt every later count at that
    // cell). Logs old vs. actual new count when tracing is on. See the ghost-array members'
    // comment for why this is a count rather than a single bit.
    auto apply_write = [&](const char* arr_name, FourDArr& arr, int xi, int yi) {
        int old_val = arr(w, xi, yi, z);
        int new_val = (value == 1) ? (old_val + 1) : std::max(old_val - 1, 0);
        if (do_trace) {
            std::cout << "GHOST_WRITE rank: " << trace_rank << " tick: " << trace_tick
                      << " arr: " << arr_name << " w: " << w << " x_idx: " << xi << " y_idx: " << yi
                      << " z: " << z << " old: " << old_val << " new: " << new_val
                      << " remote_proc: " << trace_remote_proc << " source: " << trace_tag << "\n";
        }
        arr(w, xi, yi, z) = new_val;
    };

    /* first layer ghost sites for (111) and (100) moves */
    if ((ylo_edge > yhi_edge)) {
        if (((y_unmod >= 0) && (y_unmod < yhi_edge)) ||
            ((y_unmod < total_dims[1]) && (y_unmod >= ylo_edge))) {
            y_idx = mod_with_bounds((y - chunk_bounds[1][0] + 2), (total_dims[1]));
            x_idx = x_unmod % 2;

            if ((x == xlo_edge) || (x == (xlo_edge +1))) {
                apply_write("neg_x", neg_x_arr, x_idx, y_idx);
            }
            else if ((x == (xhi_edge-1)) || (x == (xhi_edge -2))) {
                apply_write("pos_x", pos_x_arr, x_idx, y_idx);
            }
        }
    }
    else {
        if ((y_unmod < yhi_edge)  && (y_unmod >= ylo_edge)) {
            y_idx = mod_with_bounds((y - chunk_bounds[1][0] + 2), (total_dims[1]));
            x_idx = x_unmod % 2;

            if ((x == xlo_edge) || (x == (xlo_edge +1))) {
                apply_write("neg_x", neg_x_arr, x_idx, y_idx);
            }
            else if ((x == (xhi_edge-1)) || (x == (xhi_edge -2))) {
                apply_write("pos_x", pos_x_arr, x_idx, y_idx);
            }
        }
    }

    if ((xlo_edge > xhi_edge)) {
        if (((x_unmod >= 0) && (x_unmod < xhi_edge)) ||
            ((x_unmod < total_dims[0]) && (x_unmod >= xlo_edge))) {
            x_idx = mod_with_bounds((x - chunk_bounds[0][0] + 2), (total_dims[0]));
            y_idx = y_unmod % 2;

            if ((y == ylo_edge) || (y == (ylo_edge + 1))) {
                apply_write("neg_y", neg_y_arr, x_idx, y_idx);
            }
            else if ((y == (yhi_edge-1)) || (y == (yhi_edge -2))) {
                apply_write("pos_y", pos_y_arr, x_idx, y_idx);
            }
        }
    }
    else {
        if ((x_unmod < xhi_edge)  && (x_unmod >= xlo_edge)) {
            x_idx = mod_with_bounds((x - chunk_bounds[0][0] + 2), (total_dims[0]));
            y_idx = y_unmod % 2;

            if ((y == ylo_edge) || (y == (ylo_edge +1))) {
                apply_write("neg_y", neg_y_arr, x_idx, y_idx);
            }
            else if ((y == (yhi_edge-1)) || (y == (yhi_edge -2))) {
                apply_write("pos_y", pos_y_arr, x_idx, y_idx);
            }
        }
    }
}

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
/*Region* add_region(std::vector<std::string> info, int rank) {
    std::cout << "adding region \n";
    int id = std::stoi(tokenizer(info[0], ":")[0]); // region id number
    std::vector< std::vector<int> > params = vect_create_2D(2,3);
    std::vector<double> rates(2);
    std::vector<double> distribution(4,1);
    std::vector<double> interface(4);
    std::string reg_type = info.at(1); // region type
    std::string bias = info.at(2); //bias direction of region
    bool random = false; 
    double interface_terrace_rate = 0;

    if (info[1] == "GB") {
        // case of grain boundary region
        params[0][0] = std::stoi(tokenizer(info[3], ":")[1]);
        params[0][1] = std::stoi(tokenizer(info[4], ":")[1]);
        params[0][2] = std::stoi(tokenizer(info[5], ":")[1]);
        
        params[1][0] = std::stoi(tokenizer(info[6], ":")[1]);
        params[1][1] = std::stoi(tokenizer(info[7], ":")[1]);
        params[1][2] = std::stoi(tokenizer(info[8], ":")[1]);
    }
        
    if (info[1] == "BLOCK") {
        // case of region defined as rectangular prism (block)
        params[0][0] = std::stoi(tokenizer(info[3], ":")[1]);
        params[1][0] = std::stoi(tokenizer(info[4], ":")[1]);
        
        params[0][1] = std::stoi(tokenizer(info[5], ":")[1]);
        params[1][1] = std::stoi(tokenizer(info[6], ":")[1]);

        params[0][2] = std::stoi(tokenizer(info[7], ":")[1]);
        params[1][2] = std::stoi(tokenizer(info[8], ":")[1]);
        
        if (info.at(9) == "rate_neg") {
            rates[0] = std::stod(info.at(10));
            if (info.at(11) == "rate_pos") {
                rates[1] = std::stod(info.at(12));
            }
        }
        else if (info.at(9) == "rate_pos") {
            rates[1] = std::stod(info.at(10));
            if (info.at(11) == "rate_neg") {
                rates[0] = std::stod(info.at(12));
            }
        }
    }
    std::cout << "post block \n";

    if (info.size() > 13) { 
        if (info.at(13) == "RANDOM") {
            std::cout << "RANDOM\n";
            distribution[0] = std::stod(info.at(14));
            distribution[1] = std::stod(info.at(15));
            distribution[2] = std::stod(info.at(16));
            distribution[3] = std::stod(info.at(17));
            
            random = true;
        } 
        else if (info.at(13) == "INTERFACE") {
            std::cout << "INTERFACE\n";
            interface[0] = 1;
            interface[1] = std::stod(info.at(14));
            interface[2] = std::stod(info.at(15));
            interface[3] = std::stod(info.at(16));
            interface_terrace_rate = std::stod(info.at(17));
        }

        
    } 
    std::cout << "pre region \n";
    
    if (rank == 0) std::cout << "rank: " << rank << " id: " << id << " reg_type: "  << reg_type << " bias: "  << bias 
    << " params: [[" << params[0][0] << " " << params[0][1] << " " << params[0][2] << "] [" << params[1][0] << " " << params[1][1] << " " << params[1][2] << " ]] " 
    << " distribution: [" << distribution[0] << " " << distribution[1] << " " << distribution[2] << " " << distribution[3] << "] " 
    << " random: " << random << " rates: [ " << rates[0] << " " << rates[1] << " ] " 
    << " interface: [ " << interface[0] << " " << interface[1] << " " << interface[2] << " " << interface[3] << " ] " 
    << " interface_terrace_rate " << interface_terrace_rate << "\n";
    //Region* new_region = new Region(id, reg_type, bias, params, distribution, random, rates, interface, interface_terrace_rate);

    // generating random barriers according to bounds if RANDOM tag
    //included in region description
    if (info.size() > 9) { 
        if (info.at(9) == "RANDOM") {
            //new_region->random_blocking();
            new_region->random_barrier_assigner(rates);
        }
    }

    return new_region;
}
*/

Region* add_region_Elandscape(std::vector<std::string> info, int rank) {
    std::cout << "adding region \n";
    int id = std::stoi(tokenizer(info[0], ":")[0]); // region id number
    std::vector< std::vector<int> > params = vect_create_2D(2,3);
    double energy_below_bulk = 0;
    std::vector<double> distribution(4,1);
    std::vector<int> interface(4);
    std::string reg_type = info.at(1); // region type
    std::string bias = info.at(2); //bias direction of region
    bool random = false; 
    double interface_terrace_energy = 0;
    std::vector<double> rates(2);

    if (info[1] == "GB") {
        // case of grain boundary region
        params[0][0] = std::stoi(tokenizer(info[3], ":")[1]);
        params[0][1] = std::stoi(tokenizer(info[4], ":")[1]);
        params[0][2] = std::stoi(tokenizer(info[5], ":")[1]);
        
        params[1][0] = std::stoi(tokenizer(info[6], ":")[1]);
        params[1][1] = std::stoi(tokenizer(info[7], ":")[1]);
        params[1][2] = std::stoi(tokenizer(info[8], ":")[1]);
    }
        
    if (info[1] == "BLOCK") {
        // case of region defined as rectangular prism (block)
        params[0][0] = std::stoi(tokenizer(info[3], ":")[1]);
        params[1][0] = std::stoi(tokenizer(info[4], ":")[1]);        
        params[0][1] = std::stoi(tokenizer(info[5], ":")[1]);
        params[1][1] = std::stoi(tokenizer(info[6], ":")[1]);
        params[0][2] = std::stoi(tokenizer(info[7], ":")[1]);
        params[1][2] = std::stoi(tokenizer(info[8], ":")[1]);
        
        if (info.at(9) == "E_below_bulk") { energy_below_bulk = std::stod(info.at(10)); }
        else if (info.at(9) == "rate_neg") {
            rates[0] = std::stod(info.at(10));
            if (info.at(11) == "rate_pos") {
                rates[1] = std::stod(info.at(12));
            }
        }
        else if (info.at(9) == "rate_pos") {
            rates[1] = std::stod(info.at(10));
            if (info.at(11) == "rate_neg") {
                rates[0] = std::stod(info.at(12));
            }
        }
        else { 
            std::cout << "ERROR: issue in ordering or missing E_below_bulk\n";
            exit(0);
        }
    }
    std::cout << "post block \n";
    if (info.size() > 13) {
        if (info.at(13) == "E_below_bulk") { 
            energy_below_bulk = std::stod(info.at(14)); 
            if (info.at(15) == "INTERFACE") {
                std::cout << "INTERFACE\n";
                interface[0] = 1;
                interface[1] = std::stod(info.at(16));
                interface[2] = std::stod(info.at(17));
                interface[3] = 0;
                interface_terrace_energy = std::stod(info.at(18));
            }
        }
    }

    if (info.size() > 11) { 
        if (info.at(11) == "RANDOM") {
            std::cout << "RANDOM\n";
            distribution[0] = std::stod(info.at(12));
            distribution[1] = std::stod(info.at(13));
            distribution[2] = std::stod(info.at(14));
            distribution[3] = std::stod(info.at(15));
            
            random = true;
        } 
        else if (info.at(11) == "INTERFACE") {
            std::cout << "INTERFACE\n";
            interface[0] = 1;
            interface[1] = std::stod(info.at(12));
            interface[2] = std::stod(info.at(13));
            interface[3] = 0;
            interface_terrace_energy = std::stod(info.at(14));
        }
        else if (info.at(11) == "GB") {
            std::cout << "Gb\n";
            interface[0] = 0;
            interface[1] = std::stod(info.at(12));
            interface[2] = std::stod(info.at(13));
            interface[3] = 1;
            interface_terrace_energy = std::stod(info.at(14));
        }        
    } 

    std::cout << "pre region \n";

    Region* new_region = new Region(id, reg_type, bias, params, distribution, random, rates, energy_below_bulk, interface, interface_terrace_energy);

    // generating random barriers according to bounds if RANDOM tag
    //included in region description
    /*
    if (info.size() > 11) { 
        if (info.at(11) == "RANDOM") {
            //new_region->random_blocking();
            new_region->random_barrier_assigner(rates);
        }
    }
    */

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
std::vector<int> dims, int num_regions, std::string region_infile, int rank) {

    int read_idx = 0;
    std::vector<Region*> regions;
    std::vector<std::string> region_info;
    std::string curr_line = lines[read_idx];
    //std::cout << "curr_line: " << lines[read_idx] << "\n";
    //std::cout << "num_regions: " << num_regions << "\n";
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
        Region* region = add_region_Elandscape(region_info, rank);
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

            Region* region = add_region_Elandscape(region_info, rank);
            regions.push_back(region);
        }
    } 

    // entering region id in sites correspoding to pre-defined regions
    //std::cout << "regions.size(): " << regions.size() << "\n";
    //std::cout << "dims: [" << dims[0] << " " << dims[1] << " " << dims[2] << "]\n";
    // std::cout << "region_infile: " << region_infile << "\n";

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
 * @brief Reads miscellaneous rates from input lines and extracts relevant rate values.
 *
 * This function parses an input file's lines to extract rate values associated with specific labels. 
 * It processes numeric values and maps them to predefined rate indices.
 *
 * @param read_idx The starting index in the lines vector from where reading begins.
 * @param lines A vector of strings representing lines from an input file.
 * @return A tuple containing the updated read index and a vector of extracted rates.
 */
std::tuple< int, std::vector<double> > read_misc_rates(int read_idx, std::vector<std::string> lines) {
    
    std::cout << "read_misc_rates() \n";
    std::vector<double> rates(7);
    std::vector<std::string> rate_info;
    bool idx_found = false;
    int rate_idx = 0; 
    std::string curr_line = lines[read_idx];

    while (curr_line.find("rates end") == std::string::npos) {
        
        rate_info = tokenizer(lines[read_idx], " ");

        for (int i=0; i<(int)rate_info.size(); i++) {

            if ((is_numeric_or_scinotation(rate_info[i])) && (idx_found)) {
                rates[rate_idx] = std::stod(rate_info[i]);
                idx_found = false;
            }
            else if ((!is_numeric_or_scinotation(rate_info[i])) && (!idx_found)) {
                
                if (rate_info[i] == "diag") { rate_idx = 0; }
                else if (rate_info[i] == "lateral") { rate_idx = 1; }
                else if (rate_info[i] == "void_threshold") { rate_idx = 2; }
                else if (rate_info[i] == "void_rate") { rate_idx = 3; }
                else if (rate_info[i] == "terrace_rate_111") { rate_idx = 4; }
                else if (rate_info[i] == "terrace_rate_100") { rate_idx = 5; }
                else if (rate_info[i] == "void_gb_diss_rate") { rate_idx = 6; }
                idx_found = true;
            }
            else {
                std::cout << "ERROR: mismatch in order of rates and labels \n" << "\n";
                exit(0);
            }
        }  

        read_idx ++;
        curr_line = lines[read_idx];     
    }

    std::tuple< int, std::vector<double> > tuple_out(read_idx, rates);
    return tuple_out;
}


std::tuple< int, std::vector<double> > read_misc_rates_Elandscape(int read_idx, std::vector<std::string> lines) {
    
    //std::cout << "read_misc_rates() \n";
    std::vector<double> rates(11);
    std::vector<std::string> rate_info;
    bool idx_found = false;
    int rate_idx = 0; 
    std::string curr_line = lines[read_idx];

    while (curr_line.find("rates end") == std::string::npos) {
        
        rate_info = tokenizer(lines[read_idx], " ");

        for (int i=0; i<(int)rate_info.size(); i++) {
            //std::cout << "i: " << i << " rate_info[i]: " << rate_info[i] << "\n";
            //std::cout << "rate_idx: " << rate_idx << " \n";
            if ((is_numeric_or_scinotation(rate_info[i])) && (idx_found)) {
                rates[rate_idx] = std::stod(rate_info[i]);
                idx_found = false;
            }
            else if ((!is_numeric_or_scinotation(rate_info[i])) && (!idx_found)) {
                
                if (rate_info[i] == "diag") { rate_idx = 0; }
                else if (rate_info[i] == "lateral") { rate_idx = 1; }
                else if (rate_info[i] == "void_threshold") { rate_idx = 2; }
                else if (rate_info[i] == "void_E") { rate_idx = 3; }
                else if (rate_info[i] == "voidsurface_E_below_bulk") { rate_idx = 4; }
                else if (rate_info[i] == "terrace_E_111") { rate_idx = 5; }
                else if (rate_info[i] == "terrace_E_100") { rate_idx = 6; }
                else if (rate_info[i] == "void_gb_diss_E") { rate_idx = 7; }
                else if (rate_info[i] == "temp") { rate_idx = 8; }
                else if (rate_info[i] == "INTERFACE_E_below_bulk") { rate_idx = 9; }
                else if (rate_info[i] == "INTERFACE_barrier") { rate_idx = 10; }
                idx_found = true;
            }
            else {
                std::cout << "ERROR: mismatch in order of rates and labels \n" << "\n";
                exit(0);
            }
        }  

        read_idx ++;
        curr_line = lines[read_idx];     
    }
    //std::cout << "misc rates: \n";
    //print_1Dvector(rates);

    std::tuple< int, std::vector<double> > tuple_out(read_idx, rates);
    return tuple_out;
}


/*
 * @brief Populates a Lattice object by reading input files and initializing necessary data structures.
 *
 * This function performs the following tasks:
 * - Reads the input file to extract lattice dimensions, atomic types, and other necessary information.
 * - Parses and initializes vacancy, boundary condition sites, vertex sites, and region site FourDArr data structures.
 * - Reads and processes region-related data, initializing region objects.
 * - Reads and assigns rate catalogs for bulk and predefined regions.
 * - Initializes a Lattice object with extracted data and assigns region-specific rates.
 *
 * @param[in] infile_name Path to the input file containing lattice configuration.
 * @param[in] catalogfile_name Path to the rate catalog file for bulk and predefined regions.
 * @param[in] region_infile Path to the region file containing region-specific information.
 * @return A pointer to the populated Lattice object.
 */
Lattice* populate_lattice(std::string infile_name, std::string catalogfile_name, std::string region_infile, 
std::vector<int> total_dims, std::vector<std::vector<int>> chunk_bounds, int rank, std::vector<int> procs) {

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
        add_reg_struct regions_tuple = init_regions(slice_1Dvec_str(lines, read_idx, (int)lines.size()), dims_int, num_regions, region_infile, rank); 
        std::cout << "region!\n";
        incriment = regions_tuple.get_idx(); temp_regions = regions_tuple.get_regions(); temp_region_sites = regions_tuple.get_region_sites();
        read_idx += incriment;  
    }
    else {
        printf("ERROR: regions section mis-formatted in geometry file (check for extra newlines)");
        throw std::exception();
    }
    read_idx ++;

    // getting values of miscellaneous rates
    std::cout << "pre read_misc_rates()\n";
    //std::cout << "lines[read_idx]: " << lines[read_idx] << "\n";
    std::tuple< int, std::vector<double> > misc_rates_tuple;
    std::string rates_substring = "rates begin";
    std::vector<double> misc_rates;

    if (lines[read_idx].find(rates_substring) != std::string::npos) {
        read_idx ++;     
        misc_rates_tuple = read_misc_rates_Elandscape(read_idx, lines);
        read_idx = std::get<0>(misc_rates_tuple); misc_rates = std::get<1>(misc_rates_tuple); 
    }
    read_idx ++;     

    // reading in atoms, along with their type and coordinate //
    std::tuple<std::string, double, double, double, int> tuple_out;
    std::string lattice_pos;
    double x_raw; double y_raw; double z_raw; int x_unmod; int y_unmod; int z_unmod;
    int w; int x; int y; int z; int x_idx; int y_idx; int z_idx;

    int atomtype;
    int vacancies_count = 0;

    dims_int[0] = chunk_bounds[0][1] - chunk_bounds[0][0];
    dims_int[1] = chunk_bounds[1][1] - chunk_bounds[1][0];
    dims_int[2] = chunk_bounds[2][1] - chunk_bounds[2][0];

    FourDBoolArr temp_vacancies(2, (size_t)dims_int[0], (size_t)dims_int[1], (size_t)dims_int[2]);
    FourDBoolArr temp_vertex_sites(1, (size_t)dims_int[0], (size_t)dims_int[1], (size_t)dims_int[2]);
    FourDBoolArr temp_bc_sites(1, (size_t)dims_int[0], (size_t)dims_int[1], (size_t)dims_int[2]);

    // FourDArr's backing store is malloc'd (not zero-initialized like the old FourDBoolArr these
    // replaced), so each must be explicitly zeroed before the initial-seed scan below accumulates
    // into it -- multiple distinct vacancies near a process-grid corner can each contribute +1 to
    // the same cell, same as the live reference-counted ghost arrays these get copied into.
    FourDArr temp_proc_neg_x_neighbors(2, (size_t)(2), (size_t)(dims_int[1]+4), (size_t)(dims_int[2]));
    FourDArr temp_proc_neg_y_neighbors(2, (size_t)(dims_int[0]+4), (size_t)(2), (size_t)(dims_int[2]));
    FourDArr temp_proc_pos_x_neighbors(2, (size_t)(2), (size_t)(dims_int[1]+4), (size_t)(dims_int[2]));
    FourDArr temp_proc_pos_y_neighbors(2, (size_t)(dims_int[0]+4), (size_t)(2), (size_t)(dims_int[2]));
    temp_proc_neg_x_neighbors.zero();
    temp_proc_neg_y_neighbors.zero();
    temp_proc_pos_x_neighbors.zero();
    temp_proc_pos_y_neighbors.zero();

    std::vector<size_t> vacs_size_tuple = temp_vacancies.size_vec;

    for (size_t i=0; i<2; i++) {
        for (size_t j=0; j<(size_t)dims_int[0]; j++) {
            for (size_t k=0; k<(size_t)dims_int[1]; k++) {
                for (size_t l=0; l<(size_t)dims_int[2]; l++) {
                    if (i == 0) {
                        temp_vertex_sites(0,j,k,l) = 1;
                    }
                    else if (i == 1) {
                        temp_bc_sites(0,j,k,l) = 1;
                    }
                    temp_vacancies(i,j,k,l) = 0;
                }
            }
        }
    }

    int xhi_edge; int xlo_edge;
    int yhi_edge; int ylo_edge; 
    int xhi_edge_raw; int xlo_edge_raw;
    int yhi_edge_raw; int ylo_edge_raw; 

    xlo_edge = (((chunk_bounds[0][0]-2) % total_dims[0] + total_dims[0]) % total_dims[0]);
    xhi_edge = (((chunk_bounds[0][1]+2) % total_dims[0] + total_dims[0]) % total_dims[0]);
    ylo_edge = (((chunk_bounds[1][0]-2) % total_dims[1] + total_dims[1]) % total_dims[1]);
    yhi_edge = (((chunk_bounds[1][1]+2) % total_dims[1] + total_dims[1]) % total_dims[1]);
    
    xlo_edge_raw = (chunk_bounds[0][0]-2);
    xhi_edge_raw = (chunk_bounds[0][1]+2);
    ylo_edge_raw = (chunk_bounds[1][0]-2);
    yhi_edge_raw = (chunk_bounds[1][1]+2);

    if (rank == 1) std::cout << "rank: " << rank << " xlo_edge: " << xlo_edge << " xhi_edge: " << xhi_edge << "\n";
    if (rank == 1) std::cout << "rank: " << rank << " ylo_edge: " << ylo_edge << " yhi_edge: " << yhi_edge << "\n";
    if (rank == 1) std::cout << "rank: " << rank << " xlo_edge_raw: " << xlo_edge_raw << " xhi_edge_raw: " << xhi_edge_raw << "\n";
    if (rank == 1) std::cout << "rank: " << rank << " ylo_edge_raw: " << ylo_edge_raw << " yhi_edge_raw: " << yhi_edge_raw << "\n";

    std::vector<size_t> vertex_size_tup = temp_vertex_sites.size_vec;
    std::vector<size_t> bc_size_tup = temp_bc_sites.size_vec;
    std::vector<size_t> vacancies_size_tup = temp_vacancies.size_vec;
    
    //std::vector< std::vector<size_t> > coords; 
    //std::vector< std::vector<size_t> > out_coords;

    std::vector<size_t> x_dims = temp_proc_neg_x_neighbors.size_vec;
    std::vector<size_t> y_dims = temp_proc_neg_y_neighbors.size_vec;
    std::vector<size_t> x_dims_pos = temp_proc_pos_x_neighbors.size_vec;
    std::vector<size_t> y_dims_pos = temp_proc_pos_y_neighbors.size_vec;
    if (rank == 3) std::cout << "rank: " << rank << " x_dims: [ " << x_dims[0] << " " << x_dims[1] << " " << x_dims[2] << " " << x_dims[3] << " ]\n";
    if (rank == 3) std::cout << "rank: " << rank << " y_dims: [ " << y_dims[0] << " " << y_dims[1] << " " << y_dims[2] << " " << y_dims[3] << " ]\n";

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
    }
    
    std::cout << "rank: " << rank << " all_procs: \n";
    all_procs.print();

    for (int rank_i=0; rank_i<total_procs; rank_i++) {
        curr_x = rank_i % procs[0];
        curr_y = floor(rank_i / procs[0]);

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
                w = 0;
                x = mod_with_bounds(x_raw, dims_int[0]);
                y = mod_with_bounds(y_raw, dims_int[1]);
                z = mod_with_bounds(z_raw, dims_int[2]);
                z_idx = z;
            }
            else {
                w = 1;
                x = mod_with_bounds((x_raw - 0.5), dims_int[0]);
                y = mod_with_bounds((y_raw - 0.5), dims_int[1]);
                z = mod_with_bounds((z_raw - 0.5), dims_int[2]);
                //CHANGE WHEN TO USE + 1 IN IDX (ONLY FOR NEG NEIGHBOR ARRAYS)
                z_idx = z;
            }
            
            if (rank == 3) std::cout << "rank: " << rank << " lattice_pos: " << lattice_pos << " x_unmod: " << x_unmod << " y_unmod: " << y_unmod << " z_unmod: " << z_unmod << "\n";
            if (rank == 3) std::cout << "rank: " << rank << " lattice_pos: " << lattice_pos << " x: " << x << " y: " << y << " z: " << z << "\n";
            
            //ADD IF STATEMENT TO CHECK IF CHUNK_BOUNDS[i] == TOTAL_DIMS[i]

            // Diagnostic only: reuses the GHOST_WRITE tracing machinery (see set_ghost_position_impl)
            // to show exactly which raw input-file coordinate is responsible for each initial-seed
            // write, so a wrong seed can be traced back to the specific line/line-parsing that
            // produced it instead of inferring it from replicating this logic externally.
            set_ghost_position_impl(w, x_unmod, y_unmod, z_unmod, 1, total_dims, chunk_bounds,
                                     temp_proc_neg_x_neighbors, temp_proc_pos_x_neighbors,
                                     temp_proc_neg_y_neighbors, temp_proc_pos_y_neighbors,
                                     rank, -1, -1,
                                     "POPULATE_SEED(lattice_pos=" + lattice_pos + ",x_raw=" + std::to_string(x_raw) +
                                     ",y_raw=" + std::to_string(y_raw) + ",z_raw=" + std::to_string(z_raw) +
                                     ",x_unmod=" + std::to_string(x_unmod) + ",y_unmod=" + std::to_string(y_unmod) +
                                     ",z_unmod=" + std::to_string(z_unmod) + ")");
        }
        
        if ((x_unmod >= chunk_bounds[0][0]) && (x_unmod < chunk_bounds[0][1]) && 
        (y_unmod >= chunk_bounds[1][0]) && (y_unmod < chunk_bounds[1][1]) &&
        (z_unmod >= chunk_bounds[2][0]) && (z_unmod < chunk_bounds[2][1]) ) {
            if (lattice_pos == "v") {
                x = mod_with_bounds(x_unmod, dims_int[0]);
                y = mod_with_bounds(y_unmod, dims_int[1]);
                z = mod_with_bounds(z_unmod, dims_int[2]);

                if (atomtype == 0) {
                    temp_vertex_sites(0,(size_t)x,(size_t)y,(size_t)z) = 0;
                    temp_vacancies(0,(size_t)x,(size_t)y,(size_t)z) = 1;
                    /*if ( is_in(coords, {0,(size_t)x,(size_t)y,(size_t)z})) {
                        out_coords.push_back({0,(size_t)x,(size_t)y,(size_t)z});
                    }
                    else { coords.push_back({0,(size_t)x,(size_t)y,(size_t)z}); }
                    */
                    vacancies_count ++;
                }
                else if (atomtype == 1) {
                    temp_vertex_sites(0,x,y,z) = 1;
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
                    temp_bc_sites(0,(size_t)x,(size_t)y,(size_t)z) = 0;
                    temp_vacancies(1,(size_t)x,(size_t)y,(size_t)z) = 1;
                    /*
                    if (is_in(coords, {1,(size_t)x,(size_t)y,(size_t)z})) {
                        out_coords.push_back({1,(size_t)x,(size_t)y,(size_t)z});
                    }
                    else { coords.push_back({1,(size_t)x,(size_t)y,(size_t)z}); }
                    */
                    vacancies_count ++;
                }
                else if (atomtype == 1) {
                    temp_bc_sites(0,x,y,z) = 1;
                }
                else {
                    printf("Unrecognized atom type");
                    throw std::exception();
                }
            }
        }
    }  


    // intialzing lattice, basis vectors, vacancies, mobile ions, and fixed //
    // atoms based upon dimensions //
    
    int num_x_neigh = dims_int[0] + 1;
    int num_y_neigh = dims_int[1] + 1;

    Lattice* new_lattice = new Lattice(dims_int[0], dims_int[1], dims_int[2], vacancies_count, num_regions, nprocs, total_dims[0], total_dims[1], total_dims[2], temp_regions, rank);


    std::vector < std::vector <int> > temp_diag = {{0,0,0}, {1,0,0}, {0,1,0}, {1,1,0}, {0,0,1}, {1,0,1}, {0,1,1}, {1,1,1}};
    std::vector < std::vector <int> > temp_edge = {{0,0,1}, {-1,0,0}, {0,-1,0}, {0,1,0}, {1,0,0}, {0,0,-1}};
    
    for (int i=0; i < temp_diag.size(); i++) {
        for (int j=0; j < temp_diag[0].size(); j++) {
            new_lattice->diag_directions(i,j) =temp_diag[i][j];
        }
    } 

    for (int i=0; i < temp_edge.size(); i++) {
        for (int j=0; j < temp_edge[0].size(); j++) {
            new_lattice->edge_directions(i,j) = temp_edge[i][j];
        }
    } 

    for (size_t i=0; i<2; i++) {
        for (size_t j=0; j<(size_t)dims_int[0]; j++) {
            for (size_t k=0; k<(size_t)dims_int[1]; k++) {
                for (size_t l=0; l<(size_t)dims_int[2]; l++) {
                    if (i == 0) {
                        new_lattice->vertex_sites(0,j,k,l) = temp_vertex_sites(0,j,k,l);
                    }
                    else if (i == 1) {
                        new_lattice->bc_sites(0,j,k,l) = temp_bc_sites(0,j,k,l);
                    }

                    new_lattice->region_sites(i,j,k,l) = (*temp_region_sites)(i,j,k,l); // invalid read
                    new_lattice->vacancies(i,j,k,l) = temp_vacancies(i,j,k,l);
                }
            }
        }
    }

    size_t wlen; size_t xlen; size_t ylen; size_t zlen;
    
    x_dims = new_lattice->proc_pos_x_neighbors.size_vec; y_dims = new_lattice->proc_pos_y_neighbors.size_vec;
    
    wlen = x_dims[0]; xlen = x_dims[1]; ylen = x_dims[2]; zlen = x_dims[3];
    for (size_t w=0; w<wlen; w++) {
        for (size_t i=0; i<xlen; i++) {
            for (size_t j=0; j<ylen; j++) {
                for (size_t k=0; k<zlen; k++) {
                    new_lattice->proc_neg_x_neighbors(w,i,j,k) = temp_proc_neg_x_neighbors(w,i,j,k);
                    new_lattice->proc_pos_x_neighbors(w,i,j,k) = temp_proc_pos_x_neighbors(w,i,j,k);
                }
            }
        }
    }

    wlen = y_dims[0]; xlen = y_dims[1]; ylen = y_dims[2]; zlen = y_dims[3];
    // if (rank ==3) std::cout << "rank: " << rank << " y_dims: " << "[ " << wlen << " " << xlen << " " << ylen << " " << zlen << " ]\n";
    // Diagnostic only: brackets the temp->new_lattice copy loop to see whether temp_proc_neg_y's
    // own nonzero count (from the traced seed loop above) already includes unaccounted-for cells
    // before the copy even runs, or whether the copy loop itself is what introduces them.
    std::cout << "rank: " << rank << " PRE_COPY_CHECK temp_neg_y nonzero: " << temp_proc_neg_y_neighbors.nonzero_elems().rows()
              << " temp_pos_y nonzero: " << temp_proc_pos_y_neighbors.nonzero_elems().rows()
              << " y_dims: [ " << wlen << " " << xlen << " " << ylen << " " << zlen << " ]\n";
    for (size_t w=0; w<wlen; w++) {
        for (size_t i=0; i<xlen; i++) {
            for (size_t j=0; j<ylen; j++) {
                for (size_t k=0; k<zlen; k++) {
                    new_lattice->proc_neg_y_neighbors(w,i,j,k) = temp_proc_neg_y_neighbors(w,i,j,k);
                    new_lattice->proc_pos_y_neighbors(w,i,j,k) = temp_proc_pos_y_neighbors(w,i,j,k);
                }
            }
        }
    }
    std::cout << "rank: " << rank << " POST_COPY_CHECK new_neg_y nonzero: " << new_lattice->proc_neg_y_neighbors.nonzero_elems().rows()
              << " new_pos_y nonzero: " << new_lattice->proc_pos_y_neighbors.nonzero_elems().rows() << "\n";
    // Diagnostic only: direct read-back of the specific cell under investigation, immediately
    // after construction+seeding+copy and before any simulation ticks run, to settle whether it's
    // already wrong at this point (pointing at an aliasing/indexing bug somewhere in this function)
    // or only goes wrong later (pointing at the live update path instead).
    if (rank == 2) {
        std::cout << "rank: " << rank << " DIRECT_READBACK proc_neg_y_neighbors(1,49,1,124): "
                  << new_lattice->proc_neg_y_neighbors(1,49,1,124) << "\n";
    }

    if (rank == 0) {
        std::cout << "rank: " << rank << " proc_pos_x_neighbors\n";
        new_lattice->proc_pos_x_neighbors.print();
        std::cout << "rank: " << rank << " proc_pos_y_neighbors\n";
        new_lattice->proc_pos_y_neighbors.print();
        std::cout << "rank: " << rank << " proc_neg_x_neighbors\n";
        new_lattice->proc_neg_x_neighbors.print();
        std::cout << "rank: " << rank << " proc_neg_y_neighbors\n";
        new_lattice->proc_neg_y_neighbors.print();
    }
    MPI_Barrier(MPI_COMM_WORLD);
    if (rank == 1) {
        std::cout << "rank: " << rank << " proc_pos_x_neighbors\n";
        new_lattice->proc_pos_x_neighbors.print();
        std::cout << "rank: " << rank << " proc_pos_y_neighbors\n";
        new_lattice->proc_pos_y_neighbors.print();
        std::cout << "rank: " << rank << " proc_neg_x_neighbors\n";
        new_lattice->proc_neg_x_neighbors.print();
        std::cout << "rank: " << rank << " proc_neg_y_neighbors\n";
        new_lattice->proc_neg_y_neighbors.print();
    }
    MPI_Barrier(MPI_COMM_WORLD);
    if (rank == 2) {
        std::cout << "rank: " << rank << " proc_pos_x_neighbors\n";
        new_lattice->proc_pos_x_neighbors.print();
        std::cout << "rank: " << rank << " proc_pos_y_neighbors\n";
        new_lattice->proc_pos_y_neighbors.print();
        std::cout << "rank: " << rank << " proc_neg_x_neighbors\n";
        new_lattice->proc_neg_x_neighbors.print();
        std::cout << "rank: " << rank << " proc_neg_y_neighbors\n";
        new_lattice->proc_neg_y_neighbors.print();
    }
    MPI_Barrier(MPI_COMM_WORLD);
    if (rank == 3) {
        std::cout << "rank: " << rank << " proc_pos_x_neighbors\n";
        new_lattice->proc_pos_x_neighbors.print();
        std::cout << "rank: " << rank << " proc_pos_y_neighbors\n";
        new_lattice->proc_pos_y_neighbors.print();
        std::cout << "rank: " << rank << " proc_neg_x_neighbors\n";
        new_lattice->proc_neg_x_neighbors.print();
        std::cout << "rank: " << rank << " proc_neg_y_neighbors\n";
        new_lattice->proc_neg_y_neighbors.print();
    }
    MPI_Barrier(MPI_COMM_WORLD);
    
    

    // assigning rates to region-specific rate catalogs 
    std::cout << "rank: " << rank << " pre assign_region_rates_wrapper()\n";
    new_lattice->assign_region_rates_wrapper(temp_regions, misc_rates); // memory error here
    std::cout << "rank: " << rank << " post assign_region_rates_wrapper()\n";
    std::cout << "rank: " << rank << "test: \n";
    

    std::cout << "rank: " << rank << "pre a_types: \n";
    for (int i=0; i<(int)a_type_values.size(); i++) {new_lattice->a_types[a_type_keys[i]] = a_type_values[i];}

    std::cout << "rank: " << rank << " pre var assign: \n";
    new_lattice->chunk_bounds = chunk_bounds;
    new_lattice->proc_dims = procs;

    std::cout << "rank: " << rank << " pre get nonzero_vacs: \n";

    Matrix<int> nonzero_vacs = new_lattice->vacancies.nonzero(rank);

    std::cout << "rank: " << rank << " pre assign nonzero_vacs: \n";
    for (int i=0; i<(int)nonzero_vacs.rows(); i++) {
        for (int j=0; j<(int)nonzero_vacs.cols(); j++) {
            new_lattice->vacancies_pos(i,j) = nonzero_vacs[i][j];
        }
    }

    
    for (int i=0; i < temp_proc_neighbors.rows(); i++ ) {
        for (int j=0; j < temp_proc_neighbors.cols(); j++ ) {
            new_lattice->proc_neighbors(i,j) = temp_proc_neighbors(i,j);
        }
    }
    
    std::cout << "rank: " << rank << " new_lattice->proc_neighbors: \n";
    new_lattice->proc_neighbors.print();
    

    new_lattice->bulk_migration_111 = misc_rates[0];
    new_lattice->bulk_migration_100 = misc_rates[1];
    new_lattice->void_threshold = misc_rates[2];
    new_lattice->void_E = misc_rates[3];
    new_lattice->voidsurface_E_below_bulk = misc_rates[4];
    new_lattice->terrace_barrier_111 = misc_rates[5];
    new_lattice->terrace_barrier_100 = misc_rates[6];
    new_lattice->void_gb_diss_barrier = misc_rates[7];
    new_lattice->temperature = misc_rates[8];
    new_lattice->interface_E = misc_rates[9];
    new_lattice->interface_barrier = misc_rates[10];

    std::cout << "rank: " << rank << " pre temp catalog \n";

    delete temp_region_sites;

    std::cout << "rank: " << rank << " pre temp_vacancies  \n";


    std::cout << "rank: " << rank << " leaving populate_lattice()\n";
    MPI_Barrier(MPI_COMM_WORLD);
    
    new_lattice->rate_cumsum.resize(14*nonzero_vacs.rows());

    return new_lattice;
}
/*---------------------------------------------------------------------------*/

