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

#include "hpp_files/math_func.hpp"
#include "hpp_files/vec_func.hpp"
#include "hpp_files/str_func.hpp"
#include "hpp_files/print_func.hpp"


#define CEILING(x,y) ((x + y - 1) / y)

 /*! \brief A class for storing data (rate constants) corresponding to a subset of coordinates in a FourDArr object*/
class Region {
    public:
        int id;
        std::string type;
        std::vector<int> slopes;
        std::vector<int> shifts;
        std::vector<int> lowerbound;
        std::vector<int> upperbound;
        std::vector< std::vector<int> > params;
        std::vector< std::vector< std::vector<double> > > energies; 

        Region(int id_in, std::string reg_type, std::vector< std::vector<int> > params_in):
        id(id_in), 
        type(reg_type),
        params(params_in)

        {
            if (type == "GB") {
                slopes = params[0];
                shifts = params[1];
            }
            else if (type == "BLOCK") {
                lowerbound = params[0];
                upperbound = params[1];
            }
        }
};

class FourDBoolArr {
    size_t len1_, len2_, len3_, len4_;
    std::vector<uint8_t> data_;
    
    class BoolReference
    {
    private:
        uint8_t & value_;
        uint8_t mask_;
        
        void zero(void) noexcept { value_ &= ~(mask_); }
        
        void one(void) noexcept { value_ |= (mask_); }
        
        bool get() const noexcept { return !!(value_ & mask_); }
        
        void set(bool b) noexcept
        {
            if(b)
                one();
            else
                zero();
        }
        
    public:
        BoolReference(uint8_t & value, uint8_t nbit)
        : value_(value), mask_(uint8_t(0x1) << nbit)
        { }

        BoolReference(const BoolReference &ref): value_(ref.value_), mask_(ref.mask_) {}
        
        BoolReference & operator=(bool b) noexcept { set(b); return *this; }
        
        BoolReference & operator=(const BoolReference & br) noexcept { return *this = bool(br); }
        
        operator bool() const noexcept { return get(); }
        };
    
public:
    const std::tuple<size_t, size_t, size_t, size_t> size_tuple;

    FourDBoolArr(size_t len1, size_t len2, size_t len3, size_t len4) :
    len1_(len1), len2_(len2), len3_(len3), len4_(len4), data_(CEILING(len1 * len2 * len3 * len4, 8)), size_tuple(len1, len2, len3, len4) { }
    
    BoolReference operator() (size_t i1, size_t i2, size_t i3, size_t i4) {
        size_t flat_idx = i1 * len2_ * len3_ * len4_ + i2 * len3_ * len4_ + i3 * len4_ + i4;
        size_t coarse_idx = flat_idx / 8;
        size_t fine_idx = flat_idx % 8;
        BoolReference ref(data_[coarse_idx], fine_idx);
        return ref;
    }

    /*
    print elements of data field
    */
    void print_4Dvector() { 
        std::cout << "\n[ ";
        for (int i1=0; i1<(int)len1_; i1++) {
            std::cout << "[ ";
            for (int i2=0; i2<(int)len2_; i2++) {
                std::cout << "[ ";
                for (int i3=0; i3<(int)len3_; i3++) {
                    std::cout << "[ ";
                    for (int i4=0; i4<(int)len4_; i4++) {
                        std::cout << (*this)((size_t)i1, (size_t)i2, (size_t)i3, (size_t)i4)  << " ";
                    }
                    std::cout << "] ";
                    std::cout << "\n  ";
                }
                std::cout << "] ";
                std::cout << "\n  ";
            }
            std::cout << "] ";
            std::cout << "\n  ";
        } 
        std::cout << "] \n\n";
    }
    
    /*
    retrieve nonzero elements
    */
    std::vector< std::vector<int> > nonzero() {
        int vec_size = (len1_*len2_*len3_*len4_)/8;
        std::vector<std::vector<int>> vec_out(vec_size, std::vector<int>(4));

        int elem = 0; 
        for (int i1=0; i1<(int)len1_; i1++) {
            for (int i2=0; i2<(int)len2_; i2++) {
                for (int i3=0; i3<(int)len3_; i3++) {
                    for (int i4=0; i4<(int)len4_; i4++) {
                        if ( (*this)((size_t)i1, (size_t)i2, (size_t)i3, (size_t)i4) != false ) {
                            if (elem >= (vec_size-1)) {
                                vec_out.resize(vec_size*2, std::vector<int>(4));
                                vec_size = vec_size * 2;
                            }
                            vec_out[elem][0] = i1;
                            vec_out[elem][1] = i2;
                            vec_out[elem][2] = i3;
                            vec_out[elem][3] = i4;
                            elem ++;
                        }
                    }
                }
            }
        } 
        vec_out.resize(elem);
        
        return vec_out;
    }
};

/*! \brief A class for storing 4-D arrays of ints*/
class FourDArr {
public:
    /*! \brief Constructor
     * \param [in] len1 len2 len3 len4    Lengths of the 4 dimensions of the array
     */

    const std::tuple<size_t, size_t, size_t, size_t> size_tuple;

    FourDArr(size_t len1, size_t len2, size_t len3, size_t len4)
    : size_tuple(len1, len2, len3, len4)
    , len1_(len1)
    , len2_(len2)
    , len3_(len3)
    , len4_(len4)
    {
        data_ = (int *)malloc(sizeof(int) * len1 * len2 * len3 * len4);
    }
    
    /*! \brief Access an element of the 4-D array
     * \param [in] i1 First index
     * \param [in] i2 Second index
     * \param [in] i3 Third index
     * \param [in] i4 Fourth index
     * \returns Reference to array element
     */
    int& operator() (size_t i1, size_t i2, size_t i3, size_t i4) {
        return data_[i1 * len2_ * len3_ * len4_ + i2 * len3_ * len4_ + i3 * len4_ + i4];
    }
    
    int  operator() (size_t i1, size_t i2, size_t i3, size_t i4) const {
        return data_[i1 * len2_ * len3_ * len4_ + i2 * len3_ * len4_ + i3 * len4_ + i4];
    }

    /*! \brief Destructor*/
    ~FourDArr() {
        free(data_);
    }
    
    FourDArr(const FourDArr& m) = delete;
    
    FourDArr& operator= (const FourDArr& m) = delete;
    
    /*! \returns pointer to 0th element in the array */
    int *data() {
        return data_;
    }
    /*
    retrieve nonzero elements
    */
    std::vector< std::vector<int> > nonzero() {
        int vec_size = (len1_*len2_*len3_*len4_)/8;
        std::vector<std::vector<int>> vec_out(vec_size, std::vector<int>(4));

        int elem = 0; 
        for (int i1=0; i1<(int)len1_; i1++) {
            for (int i2=0; i2<(int)len2_; i2++) {
                for (int i3=0; i3<(int)len3_; i3++) {
                    for (int i4=0; i4<(int)len4_; i4++) {
                        if ( (*this)((size_t)i1, (size_t)i2, (size_t)i3, (size_t)i4) !=0 ) {

                            if (elem == (vec_size-2)) {
                                vec_out.resize(vec_size*2, std::vector<int>(4));
                                vec_size = vec_size * 2;
                            }
                            vec_out[elem][0] = i1;
                            vec_out[elem][1] = i2;
                            vec_out[elem][2] = i3;
                            vec_out[elem][3] = i4;
                            
                            elem ++;
                        }
                    }
                }
            }
        } 
        vec_out.resize(elem);
        
        return vec_out;
    }

    /*
    print data field of a FourDArr object
    */
    void print_4Dvector() { 
        std::cout << "\n[ ";
        for (int i1=0; i1<(int)len1_; i1++) {
            std::cout << "[ ";
            for (int i2=0; i2<(int)len2_; i2++) {
                std::cout << "[ ";
                for (int i3=0; i3<(int)len3_; i3++) {
                    std::cout << "[ ";
                    for (int i4=0; i4<(int)len4_; i4++) {
                        std::cout << (*this)((size_t)i1, (size_t)i2, (size_t)i3, (size_t)i4)  << " ";
                    }
                    std::cout << "] ";
                    std::cout << "\n  ";
                }
                std::cout << "] ";
                std::cout << "\n  ";
            }
            std::cout << "] ";
            std::cout << "\n  ";
        } 
        std::cout << "] \n\n";
    }

    /*
    retrieve elements to a 1D slice of a FourDArr object corresponding to coordinates input
    */
    std::vector<int> grab_idxs(std::vector< std::vector<int> > coords) {
        std::vector<int> output;

        for (int i=0; i<(int)coords.size(); i++) {
            std::vector<int> coord = coords[i];
            int w = coord[0];
            int x = coord[1];
            int y = coord[2];
            int z = coord[3];

            int value = (int)(*this)(w,x,y,z);
            output.push_back(value);
        }
        return output;
    }

    /*
    assign elements to a 1D slice of a FourDArr object corresponding to values input
    */
    void assign_idxs(std::vector< std::vector<int> > coords, std::vector<int> values) {
        size_t w;
        size_t x;
        size_t y;
        size_t z;
        
        

        for (int i=0; i<(int)coords.size(); i++) {
            std::vector<int> coord = coords[i];
            w = (size_t)coord[0];
            x = (size_t)coord[1];
            y = (size_t)coord[2];
            z = (size_t)coord[3];
            (*this)(w,x,y,z) = values[i];
        }
    }

private:
    size_t len1_, len2_, len3_, len4_; ///< Dimensions of the array
    int* data_; ///< The data stored in the array
};

/*! \brief A class for representing and manipulating matrices with variable dimension sizes
 * \tparam mat_type The type of elements to be stored in the matrix
 */
template <class mat_type>
class Matrix {
public:
    class RowReference {
        private:
        size_t row_idx_;
        Matrix<mat_type> &mat_;
        
        public:
        RowReference(Matrix<mat_type> &mat, size_t row) : row_idx_(row), mat_(mat) {}
        
        mat_type& operator[] (size_t idx) {
            return mat_(row_idx_, idx);
        }
    };
    
    /*! \brief Constructor for Matrix class
        * \param [in] rows     The number of rows the matrix should have initially
        * \param [in] cols     The number of columns the matrix should have initially
        */
    Matrix(size_t rows, size_t cols) : rows_(rows), cols_(cols), tot_size_(rows * cols), data_(rows * cols) {}
    
    /*! \brief Access matrix element
     * \param [in] row      Row index of element
     * \param [in] col      Column index of element
     * \return Reference to matrix element
     */
    mat_type& operator() (size_t row, size_t col) {
        return data_[cols_*row + col];
    }
    
    /*! \brief Access matrix element
     * \param [in] row      Row index of element
     * \param [in] col      Column index of element
     * \return matrix element
     */
    mat_type  operator() (size_t row, size_t col) const {
        return data_[cols_ * row + col];
    }
    
    /*! \brief Zero all matrix elements */
    void zero() {
        std::fill(data_.begin(), data_.end(), 0);
    }
    
    /*! \brief Access matrix row
     * \param [in] row      Row index
     * \return pointer to 0th element in a row of a matrix
     */
    mat_type *operator[] (size_t row) {
        return &data_[cols_ * row];
    }
    
    /*! \brief Access matrix row
     * \param [in] row      Row index
     * \return pointer to 0th element in a row of a matrix
     */
    const mat_type *operator[] (size_t row) const {
        return &data_[cols_ * row];
    }
    
    /*! \brief Increase number of columns in the matrix
     * Data are copied such that the first n[i] elements in each row remain the same before and after this operation
     * \param [in] new_col      Desired number of columns in the enlarged matrix
     * \param [in] n_keep       Number of elements to preserve in each row of the matrix (should have \p rows_ elements)
     */
    void enlarge_cols(size_t new_col, int *n_keep) {
        if (new_col > cols_) {
            size_t old_cols = cols_;
            reshape(rows_, new_col);
            
            size_t row_idx;
            for (row_idx = rows_; row_idx > 0; row_idx--) {
                auto begin = data_.begin();
                std::copy_backward(begin + (row_idx - 1) * old_cols, begin + (row_idx - 1) * old_cols + n_keep[row_idx - 1], begin + (row_idx - 1) * new_col + n_keep[row_idx - 1]);
            }
        }
    }
    
    /*! \brief Increase number of columns in the matrix
     * Data are copied such that the first n elements in each row remain the same before and after this operation
     * \param [in] new_col      Desired number of columns in the enlarged matrix
     * \param [in] n_keep       Number of elements to preserve in all rows of the matrix
     */
    void enlarge_cols(size_t new_col, int n_keep) {
        if (new_col > cols_) {
            size_t old_cols = cols_;
            reshape(rows_, new_col);
            
            size_t row_idx;
            for (row_idx = rows_; row_idx > 0; row_idx--) {
                auto begin = data_.begin();
                std::copy_backward(begin + (row_idx - 1) * old_cols, begin + (row_idx - 1) * old_cols + n_keep, begin + (row_idx - 1) * new_col + n_keep);
            }
        }
    }

    
    /*! \brief Change the dimensions without moving any of the data
     * \param [in] new_rows     Desired number of rows in the reshaped matrix
     * \param [in] new_cols     Desired number of columns in the reshaped matrix
     */
    void reshape(size_t new_rows, size_t new_cols) {
        size_t new_size = new_rows * new_cols;
        if (new_size > tot_size_) {
            tot_size_ = new_size;
            data_.resize(tot_size_);
        }
        rows_ = new_rows;
        cols_ = new_cols;
    }
    
    /*! \return Current number of rows in matrix */
    size_t rows() const {
        return rows_;
    }

    /*! \return Current number of columns in matrix*/
    size_t cols() const {
        return cols_;
    }
    
    /*! \return Pointer to the  data in the matrix*/
    mat_type *data() const {
        return (mat_type *) data_.data();
    }
    
    void copy_from(Matrix<mat_type> &mat) {
        std::copy(mat.data_.begin(), mat.data_.end(), data_.begin());
    }
    
    void print() {
        std::cout << "[ ";
        for (int m=0; m<(int)rows_; m++) {
            std::cout << "[ ";
            for (int n=0; n<(int)cols_; n++) {
                std::cout << (*this)((size_t)m, (size_t)n)  << " ";
            }
            std::cout << "] ";
            std::cout << "\n  ";
        }
        std::cout << "] \n\n";
    }
    
private:
    size_t rows_, cols_, tot_size_;
    std::vector<mat_type> data_;
};

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
        std::vector<int> lattice_dim;
        std::map<int, std::string> a_types;
        FourDArr vertex_sites;
        FourDBoolArr vacancies;
        FourDArr bc_sites;
        FourDArr region_sites;
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

        Lattice(int xdim, int ydim, int zdim, int num_vacancies, int num_regions):
            vertex_sites((size_t)1, (size_t)xdim, (size_t)ydim, (size_t)zdim),
            vacancies((size_t)2, (size_t)xdim, (size_t)ydim, (size_t)zdim),
            bc_sites((size_t)1, (size_t)xdim, (size_t)ydim, (size_t)zdim),
            region_sites((size_t)2, (size_t)xdim, (size_t)ydim, (size_t)zdim),
            moves_coords((14 * num_vacancies), 4),
            moves_shifts((14 * num_vacancies), 3),
            moves_lattice((14 * num_vacancies), 1),
            moves_vacs((14 * num_vacancies), 1),
            ratecatalog_111(2, exp_int(2,8)),
            ratecatalog_100(2, exp_int(2,14)),
            regionrates_111_L(num_regions, exp_int(2,8)),
            regionrates_111_R(num_regions, exp_int(2,8)),
            regionrates_100_L(num_regions, exp_int(2,14)),
            regionrates_100_R(num_regions, exp_int(2,14)),
            configs_111(1, exp_int(2,8)),
            configs_100(1, exp_int(2,14)),
            vacancies_pos(num_vacancies, 4),
            mt_obj((unsigned int)(std::chrono::high_resolution_clock::now().time_since_epoch().count())),
            x_rand(0, xdim),
            y_rand(0, ydim)

            {
                diag_directions = {{0,0,0}, {1,0,0}, {0,1,0}, {1,1,0}, {0,0,1}, {1,0,1}, {0,1,1}, {1,1,1}};
                edge_directions = {{0,0,1}, {-1,0,0}, {0,-1,0}, {0,1,0}, {1,0,0}, {0,0,-1}};
                lattice_dim = {xdim,ydim,zdim};
                num_of_moves = 0;
                num_of_vacs = num_vacancies;
            }

        void new_get_actions() {
            int curr_move_num = 0; // total number of moves at this current timestep  
            double rate; 
            int vacs_on_interface = 0; // vacancies at last z-index of lattice (used to calculate rate for stripping)
            int num_interface_sites = lattice_dim[0] * lattice_dim[1]; // total number of sites at last z-index 

            std::vector<int> moves(2);

            rates.resize((int)moves_coords.rows());
            int i=0; int j=0; int k=0; int l=0;

            // looping over all vacancies in system
            for (int idx=0; idx < (int)vacancies_pos.rows(); idx++) {
                // position in lattice of vacancy 
                i = vacancies_pos[idx][0];
                j = vacancies_pos[idx][1];
                k = vacancies_pos[idx][2];
                l = vacancies_pos[idx][3];

                if ((curr_move_num + (num_interface_sites - vacs_on_interface)) >= ((int)moves_shifts.rows() - 20)) {
                    // resizing data structures to accommodate all moves 

                    int newsize = 2 * moves_shifts.rows();
                    rate_cumsum.resize(newsize);
                    moves_coords.reshape(newsize, 4);
                    moves_shifts.reshape(newsize, 3);
                    moves_lattice.reshape(newsize, 1);
                    moves_vacs.reshape(newsize, 1);
                }

                // finding all moves along the {111} family of vectors
                for (int s=0; s < (int)diag_directions.size(); s++) {
                    if ((l == 0) && (i == 0) && (diag_directions[s][2] == 1)) {/* checking for leftmost non-periodic boundary along z-axis*/}
                    else if ((l == (int)(lattice_dim[2]-1)) && (i == 1) && (diag_directions[s][2] == 1)) {/* checking for rightmost non-periodic boundary along z-axis*/}
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
            }

            // finding all moves coresponding to anode stripping
            /*
            double strip_rate = 7.06;
            for (int i1=0; i1<lattice_dim[0]; i1++) {
                for (int i2=0; i2<lattice_dim[1]; i2++) {
                    if (vacancies(1, i1, i2, (lattice_dim[2]-1)) == 0) {
                        rate_cumsum[curr_move_num] = strip_rate + rate_cumsum[curr_move_num-1];
                        moves_shifts[curr_move_num][0] = 0;
                        moves_shifts[curr_move_num][1] = 0;
                        moves_shifts[curr_move_num][2] = 0;
                        moves_coords[curr_move_num][0] = 1;
                        moves_coords[curr_move_num][1] = i1;
                        moves_coords[curr_move_num][2] = i2;
                        moves_coords[curr_move_num][3] = (lattice_dim[2]-1);
                        moves_lattice[curr_move_num][0] = 4;
                        moves_vacs[curr_move_num][0] = -1; 
                        curr_move_num ++;
                    }
                }
            }
            */

            // UPDATING SIZE OF DATA STRUCTURES CONTIANING COORDINATES AND RATES OF MOVES
            num_of_moves = curr_move_num;
            rate_cumsum.resize(num_of_moves);
            moves_vacs.reshape(num_of_moves, 1);
            moves_coords.reshape(num_of_moves, 4);
            moves_shifts.reshape(num_of_moves, 3);
            moves_lattice.reshape(num_of_moves, 1);
        }
        
        // function for finding rate corresponding to NN-encoding and move type
        double new_get_rateconstants(int* coord, int* shift, int lattice) {  
            
            double rate = -1; 
            int LR_idx;  // index corresponding to direction of movement in lattice (left/right)
            int idx = 0;

            // DETERMINING DIRECTION OF MOVE //

            // moving vacancy from vertex site to bc site
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

            // return -1 if move not allowed
            if (idx == -1) {
                return -1;
            }
    
            int reg_id = region_sites(coord[0], coord[1], coord[2], coord[3]);
            
            if (reg_id != 0) {
                // checking to see if move occurs in region with specially-defined rate constants
                //std::cout << "Region: coord[0]: " << coord[0] << " coord[1]: " << coord[1] << "coord[2]: " << coord[2] << "coord[3]: " << coord[3] << "\n";

                if ((lattice == 0) || (lattice == 1)) {
                    if (LR_idx == 1) {rate = regionrates_111_L[(reg_id-1)][idx];}
                    else if (LR_idx == 0) {rate = regionrates_111_R[(reg_id-1)][idx];}
                    
                }
                else if ((lattice == 2) || (lattice == 3)) {
                    if (LR_idx == 1) {rate = regionrates_100_L[(reg_id-1)][idx];}
                    else if (LR_idx == 0) {rate = regionrates_100_R[(reg_id-1)][idx];}
                }
                else {
                    std::string str_output = "Error: invalid lattice type in search_catalog()";
                    printf("%s", str_output.c_str());
                    throw std::exception();
                }

            }
            
            else {
                // in case of no pre-defined region, use bulk rate constants
                if ((lattice == 1) || (lattice == 0)) {rate = ratecatalog_111[LR_idx][idx];}
                else if ((lattice == 2) || (lattice == 3)) {rate = ratecatalog_100[LR_idx][idx];}
            }

            return rate;
        }

        /*
        deprecated method corresponding to the case of only one rate constant
        */
        double search_catalog_oneconst(int* coord, int* shift, int lattice) {
            int LR_idx;
            int idx = 0;
            double rate = 0;

            // getting index of direction //
            if (lattice == 0) {
                if (shift[2] == 0) {LR_idx = 0;}
                else {LR_idx = 1;}
            }
                
            else if (lattice == 1) {
                if (shift[2] == 1) {LR_idx = 0;}
                else {LR_idx = 1;}
            }

            else if ((lattice == 2) || (lattice == 3)) {
                if ((shift[2] == 0) && (shift[2] == 1)) {LR_idx = 0;}
                else {LR_idx = 1;}
            }

            // getting indices of configurations in rate catalog //
            if ((lattice == 1) || (lattice == 0)) {idx = configs_111[0][0];}
            else if ((lattice == 2) || (lattice == 3)) {idx = configs_100[0][0];};

            if (idx == -1) {
                return -1;
            }

            int reg_id = region_sites(coord[0], coord[1], coord[2], coord[3]);
            if (reg_id != 0) {

                if (lattice == 0) {
                    int atom_id = bc_sites(coord[0], coord[1], coord[2], coord[3]);
                    rate = region_energies[reg_id][atom_id][0][0];
                }
                else if (lattice == 1) {
                    int atom_id = vertex_sites(coord[0], coord[1], coord[2], coord[3]);
                    rate = region_energies[reg_id][atom_id][1][0];
                }
                else if (lattice == 2) {
                    int atom_id = vertex_sites(coord[0], coord[1], coord[2], coord[3]);
                    rate = region_energies[reg_id][atom_id][2][0];
                }
                else if (lattice == 3) {
                    int atom_id = bc_sites(coord[0], coord[1], coord[2], coord[3]);
                    rate = region_energies[reg_id][atom_id][3][0];
                }
                else {
                    std::string str_output = "Error: invalid lattice type in search_catalog()";
                    printf("%s", str_output.c_str());
                    throw std::exception();
                }
            }

            else {
                if ((lattice == 1) || (lattice == 0)) {rate = ratecatalog_111[LR_idx][0];}
                else if ((lattice == 2) || (lattice == 3)) {rate = ratecatalog_100[LR_idx][0];}
            }
            return rate;
        }

        /*
        deprecated method for finding energy corresponding to NN-encoding and move type
        */
        double new_search_catalog(int encoding, int* coord, int* shift, int lattice) {
            int LR_idx; // index corresponding to direction of movement in lattice (left/right)
            int idx = 0;
            double energy = 0;

            // DETERMINING DIRECTION OF MOVE //

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
                if ((shift[2] == 0) && (shift[2] == 1)) {LR_idx = 0;}
                else {LR_idx = 1;}
            }

            // DETERMINING INDEX OF MOVE
            if ((lattice == 1) || (lattice == 0)) {idx = configs_111[0][encoding];}
            else if ((lattice == 2) || (lattice == 3)) {idx = configs_100[0][encoding];};

            // return -1 if move not allowed
            if (idx == -1) {
                return -1;
            }

            
            int reg_id = region_sites(coord[0], coord[1], coord[2], coord[3]);
            if (reg_id != 0) {
                // checking to see if move occurs in region with specially-defined rate constants

                if (lattice == 0) {
                    int atom_id = bc_sites(coord[0], coord[1], coord[2], coord[3]);
                    energy = region_energies[reg_id][atom_id][0][idx];
                }
                else if (lattice == 1) {
                    int atom_id = vertex_sites(coord[0], coord[1], coord[2], coord[3]);
                    energy = region_energies[reg_id][atom_id][1][idx];
                }
                else if (lattice == 2) {
                    int atom_id = vertex_sites(coord[0], coord[1], coord[2], coord[3]);
                    energy = region_energies[reg_id][atom_id][2][idx];
                }
                else if (lattice == 3) {
                    int atom_id = bc_sites(coord[0], coord[1], coord[2], coord[3]);
                    energy = region_energies[reg_id][atom_id][3][idx];
                }
                else {
                    std::string str_output = "Error: invalid lattice type in search_catalog()";
                    printf("%s", str_output.c_str());
                    throw std::exception();
                }
            }
            else {
                // in case of no pre-defined region, use bulk rate constants
                if ((lattice == 1) || (lattice == 0)) {energy = ratecatalog_111[LR_idx][idx];}
                else if ((lattice == 2) || (lattice == 3)) {energy = ratecatalog_100[LR_idx][idx];}
            }
            return energy;
        }

        /*
        generating encoding corresponding to configuration of nearest neighbors for a vacancy
        */
        int new_get_neighbors(int* vac, int lattice) {
            int sum = 0;
            int m = (int)a_types.size(); // number of atom types in system

            // moving vacancy from bc site to bc site 
            if (lattice == 3) {
                for (int i=0; i<(int)diag_directions.size(); i++) {
                    sum += exp_int(m,i) * vertex_sites(0, (((vac[0] - diag_directions[i][0]) % lattice_dim[0] + lattice_dim[0]) % lattice_dim[0]), (((vac[1] - diag_directions[i][1]) % lattice_dim[1] + lattice_dim[1]) % lattice_dim[1]), (((vac[2] - diag_directions[i][2]) % lattice_dim[2] + lattice_dim[2]) % lattice_dim[2]));
                }
                for (int i=0; i<(int)edge_directions.size(); i++) {
                    sum += exp_int(m,i) *  bc_sites(0, (((vac[0] + edge_directions[i][0]) % lattice_dim[0] + lattice_dim[0]) % lattice_dim[0]), (((vac[1] + edge_directions[i][1]) % lattice_dim[1] + lattice_dim[1]) % lattice_dim[1]), (((vac[2] + edge_directions[i][2]) % lattice_dim[2] + lattice_dim[2]) % lattice_dim[2]));
                }
            }

            // moving vacancy from vertex site to vertex site 
            else if (lattice == 2) {
                for (int i=0; i<(int)diag_directions.size(); i++) {
                    sum +=  exp_int(m,i) * bc_sites(0, (((vac[0] + diag_directions[i][0])  % lattice_dim[0] + lattice_dim[0]) % lattice_dim[0]), (((vac[1] + diag_directions[i][1]) % lattice_dim[1] + lattice_dim[1]) % lattice_dim[1]), (((vac[2] + diag_directions[i][2]) % lattice_dim[2] + lattice_dim[2]) % lattice_dim[2]));
                }
                for (int i=0; i<(int)edge_directions.size(); i++) {
                    sum += exp_int(m,i) *  vertex_sites(0, (((vac[0] + edge_directions[i][0]) % lattice_dim[0] + lattice_dim[0]) % lattice_dim[0]), (((vac[1] + edge_directions[i][1]) % lattice_dim[1] + lattice_dim[1]) % lattice_dim[1]), (((vac[2] + edge_directions[i][2]) % lattice_dim[2] + lattice_dim[2]) % lattice_dim[2]));
                }
            }

            // moving vacancy from bc site to vertex site
            else if (lattice == 1) {
                for (int i=0; i<(int)diag_directions.size(); i++) {
                    sum += exp_int(m,i) * vertex_sites(0, (((vac[0] - diag_directions[i][0]) % lattice_dim[0] + lattice_dim[0]) % lattice_dim[0]), (((vac[1] - diag_directions[i][1]) % lattice_dim[1] + lattice_dim[1]) % lattice_dim[1]), (((vac[2] - diag_directions[i][2]) % lattice_dim[2] + lattice_dim[2]) % lattice_dim[2]));
                }
            }

            // moving vacancy from vertex site to bc site     
            else if (lattice == 0) {
                for (int i=0; i<(int)diag_directions.size(); i++) {
                    sum += exp_int(m,i) * bc_sites(0, (((vac[0] + diag_directions[i][0]) % lattice_dim[0] + lattice_dim[0]) % lattice_dim[0]), (((vac[1] + diag_directions[i][1]) % lattice_dim[1] + lattice_dim[1]) % lattice_dim[1]), (((vac[2] + diag_directions[i][2])% lattice_dim[2] + lattice_dim[2]) % lattice_dim[2]));
                }
            }

            return sum;
        }
 
        /*
        updating the positions of atoms on lattice according to selected move
        */
        void new_update_lattice(int idx) {
            int* new_loc = moves_coords[idx]; // new location of vacancy
            std::vector<int> old_loc(3); // new location of vacancy
            int vacs_idx = *moves_vacs[idx]; // index of vacancy in master vector
            
            
            for (int i=0; i<3; i++) {
                old_loc[i] = (((new_loc[i+1] - moves_shifts[idx][i]) % lattice_dim[i]) + lattice_dim[i]) % lattice_dim[i];
            }

            // adding vacancy corresponding to stripping move 
            if (*moves_lattice[idx] == 4) {
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
                vacancies_pos.reshape(num_of_vacs, 4);
                for (int i=0; i<4; i++) { vacancies_pos[(num_of_vacs-1)][i] = new_loc[i]; }
                
            }

            // moving vacancy from bc site to bc site 
            if (*moves_lattice[idx] == 3) {
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
                
                vacancies_pos[vacs_idx][0] = 1;    
            }
                        
            // moving vacancy from vertex site to vertex site 
            else if (*moves_lattice[idx] == 2) {
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
                
                vacancies_pos[vacs_idx][0] = 0;
            }
                        
            // moving vacancy from bc site to vertex site 
            else if (*moves_lattice[idx] == 1) {
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

                vacancies_pos[vacs_idx][0] = 0; 
                    
            }
                                        
            // moving vacancy from vertex site to bc site 
            else if (*moves_lattice[idx] == 0) {
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

                vacancies_pos[vacs_idx][0] = 1;
                
            }

            // updating vector of positions of all vacacies
            vacancies_pos[vacs_idx][1] = new_loc[1];
            vacancies_pos[vacs_idx][2] = new_loc[2];
            vacancies_pos[vacs_idx][3] = new_loc[3];
            
        }

        /*
        calculating time elapsed for a move
        */
        double new_random_times(std::vector<double> cumsum) {
            // creating a random double between 0 and 1
            unsigned int random = mt_obj();
            double random_double = ((random / (1.+ UINT32_MAX)) +  (1 / (1.+ UINT32_MAX)));
            
            //calculating the time elapsed
            int last_idx = (int)cumsum.size() - 1;
            double time = ((-1/ cumsum[last_idx]) * log(random_double));

            return time;
        }

        /*
        selecting random move in vector of moves, with selection probability proportional
        to rate constant corresponding to move
        */
        int get_idx(std::vector<double> cumsum) {
            // creating a random double between 0 and 1
            unsigned int random = mt_obj();
            double random_double = ((random / (1.+ UINT32_MAX)) +  (1 / (1.+ UINT32_MAX)));

            // accessing random element into cumulative sum array, 
            // probability of access proportional to the value at that point
            int last_idx = (int)cumsum.size() - 1;
            double rand_pos = cumsum[last_idx] * random_double;
            int min_idx = searchsorted_recursive(&cumsum, rand_pos, 0, last_idx);

            return min_idx;
        }

        /*
        wrapper method containing initialization of all varables in system, timer, and 
        calls to update state of system and list of moves
        */
        std::tuple< std::vector< std::vector< std::vector<int> > >, std::vector<int>, 
        std::vector<double>, std::vector<double> > new_kmc_iterator(double time_lim, std::chrono::system_clock::time_point start, std::string folder, int iteration) {
            fprintf(stdout, "%s", "beginning kmc iterations \n\n");
            double t = 0;

            // INITIALIZING VARIABLES PRIOR TO BEGINNING FIRST KMC STEP //

            std::chrono::system_clock::time_point end; // current (real) clock time
            std::chrono::duration<double> elapsed_seconds; // elapesed (simulated) time in simulation

            std::vector<double> timesteps; // time elapsed at each step
            std::vector<int> move_counts(4, 0); // each type of move propogated by simulation 
            std::vector<double> time_count(4, 0); // time elapsed by each type of move
            int min_idx; // index of move selected
            double timestep;
            std::vector< std::vector<int> > only_vacancies; // configuration of vacancies at current timestep
            std::vector< std::vector< std::vector<int> > > all_vacancies; // vector containing trajectory of vacancies
            std::vector<double> all_times; // vector containing trajectory of time elapsed by each type of move
            new_get_actions(); // updating list of moves in system

            bool ten_k_partition = false;
            bool first_partition = false;
            bool sec_partition = false;

            only_vacancies = vacancies.nonzero();
            all_vacancies.push_back(only_vacancies); 

            int move_ticks = 0;
            double old_time;

            // output files 
            std::ofstream out_file;
            std::ostringstream ss;
            std::string output_filename;
            std::string times_filename;
            std::string count_filename;

            while (t < time_lim) {

                end = std::chrono::system_clock::now(); 
                elapsed_seconds = end-start; 
                min_idx = get_idx(rate_cumsum);
                new_update_lattice(min_idx);
                move_counts[*moves_lattice[min_idx]] ++; 
                timestep = new_random_times(rate_cumsum);
                t += timestep; 
                time_count[*moves_lattice[min_idx]] += timestep; 
                timesteps.push_back(t);
                move_ticks ++; 
                old_time = t; 
                
                // terminating simulation after real-time limit reached
                if (elapsed_seconds.count() >= 172800) {
                    only_vacancies = vacancies.nonzero();
                    all_vacancies.push_back(only_vacancies);
                    t = time_lim + 1;
                    std::cout << "t: " << t << "\n";
                    print_1Dvector(move_counts);
                    break;
                }
                
                if ((t > (int)floor(time_lim/((double)3))) && (first_partition == false)) {
                    only_vacancies = vacancies.nonzero();
                    all_vacancies.push_back(only_vacancies);
                    all_times.push_back(t);
                    first_partition = true;
                    std::cout << "t: " << t << "\n";
                    print_1Dvector(move_counts);
                }

                else if ((t > 10000) && (ten_k_partition == false)) {
                    only_vacancies = vacancies.nonzero();
                    all_vacancies.push_back(only_vacancies);
                    all_times.push_back(t);
                    ten_k_partition = true;
                    std::cout << "t: " << t << "\n";
                    print_1Dvector(move_counts);
                }
                
                else if (((t > (int)floor((double)2*time_lim/((double)3))) && (sec_partition == false)) || (t >= time_lim)) {
                    only_vacancies = vacancies.nonzero();
                    all_vacancies.push_back(only_vacancies);
                    all_times.push_back(t);
                    sec_partition = true;
                    std::cout << "t: " << t << "\n";
                    print_1Dvector(move_counts);
                }
                
                new_get_actions();

                // writing output files every 500 timesteps
                if (move_ticks % 250 == 0) {
                    only_vacancies = vacancies.nonzero();
                    ss << folder << "/vacs/vacancies_output_" << iteration << "_" << move_ticks << "_" << t << "_moves.txt";
                    output_filename = ss.str();
                    write_to_file(output_filename, only_vacancies);
                    ss.str("");
                    ss.clear();
                
                    for (int l=0; l<(int)move_counts.size(); l++) {
                        //ss << "yueqi/yueqi_LiO2_moves_stripping/counts/counts_output_" << t << "_moves.txt_" << l << ".txt";
                        ss << folder << "/counts/counts_output_" << iteration << "_" << move_ticks << "_" << t << "_moves.txt_" << l << ".txt";
                        count_filename = ss.str();
                        out_file.open(count_filename);
                        out_file << move_counts[l] << "\n";
                        ss.str("");
                        ss.clear();
                        out_file.close();

                        //ss << "yueqi/yueqi_LiO2_moves_stripping/times/times_output_" << t << "_moves.txt_" << l << ".txt";
                        ss << folder << "/times/times_output_" << iteration << "_" << move_ticks << "_" << t << "_moves.txt_" << l << ".txt";
                        times_filename = ss.str();
                        out_file.open(times_filename);
                        out_file << time_count[l] << "\n";
                        ss.str("");
                        ss.clear();
                        out_file.close();
                    }
                }
            }

            std::cout << "t: " << t << "\n";
            print_1Dvector(move_counts);

            only_vacancies = vacancies.nonzero();
            all_vacancies.push_back(only_vacancies);
            all_times.push_back(t);

            std::tuple< std::vector< std::vector< std::vector<int> > >, 
            std::vector<int>, std::vector<double>, std::vector<double> > return_tuple(all_vacancies, move_counts, time_count, all_times);
                
            return return_tuple;
        }
};

/*---------------------------------------------------------------------------*/

/*
creating entries in rate catalog corresponding to all configurations of
"len" sites with m vacancies, corresponding to all the binary strings of 
length "len" with m zeros.
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

/*
converting a inputted base-m string (type std::string) to an int
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

/*
reads in file corresponding to rate catalog and creates catalog of rates with allowed
moves, catalog of allowed configurations for moves, and rate catalog corresponding to regions
*/
std::tuple< std::vector< std::vector<int> >, 
std::vector< std::vector< std::vector<double> > >, 
std::vector< std::vector< std::vector< std::vector<double> > > > ,
int > 
updated_create_ratecatalog(std::string catalogfile, std::vector<int> atype_list) {
    
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
            configs = reorder_ints_inp(configs, unsorted_idxs);

            dft_energies[0] = (reorder_floats_inp(dft_energies[0], unsorted_idxs));
            dft_energies[1] = (reorder_floats_inp(dft_energies[1], unsorted_idxs));

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
        std::cout << "configline: " << configline << "\n";

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

    std::tuple< std::vector< std::vector<int> >, 
    std::vector< std::vector< std::vector<double> > >, 
    std::vector< std::vector< std::vector< std::vector<double> > > >,
    int > tuple_out(all_configs, all_energies, regions_catalog, region_num);

    return tuple_out;
}


/*
creating region objects corresponding to input file
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

/*
populating region_sites FourDArr with values corresponding to custom regions input file
*/
void custom_draw_regions(FourDArr* sites, std::vector<Region*> regions, int custom_reg_idx, std::vector<int> dim, std::string infile_name) {
    std::cout << "drawing regions \n";
    Region* region;

    std::fstream in_file;
    std::vector<std::string> lines;
    std::string line;
    std::string output;
    int read_idx = 0;
    std::tuple<std::string, int, int, int> tuple_out;
    std::string lattice_pos; int x; int y; int z;

    for (int i=0; i<(int)regions.size(); i++) {  
        region = regions[i];
        std::cout << "opening custom region file \n";
        in_file.open(infile_name);

        if (in_file.is_open()) {
            std::cout << "custom region file open\n";
            while ( getline (in_file,line) )
            {
                tuple_out = parse_reg_line(line);
                lattice_pos = std::get<0>(tuple_out); x = std::get<1>(tuple_out); y = std::get<2>(tuple_out); 
                z = std::get<3>(tuple_out);
                
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

/*
populating region_sites FourDArr with values corresponding to input file
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

/*
wrapper function reading input file and creating region objects / populating region_sites FourDArr
*/
std::tuple< int, std::vector<Region*>, FourDArr* > init_regions(std::vector<std::string> lines, std::vector<int> dims, int num_regions, std::string region_infile) {

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

        for (int i=(int)regions.size(); i < num_regions; i++) {
            
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

    std::tuple< int, std::vector<Region*>, FourDArr* > tuple_out(read_idx, regions, temp_region_sites);

    return tuple_out;
}

/*
wrapper function for:
- reading input files
- populating vacancy, bc_site, vertex_site, and region_site FourDArr data structures
corresponding to initial configuration of simulation
-  creating region objects 
- creating rate catalog for bulk & pre-defined regions
*/
Lattice* populate_lattice(std::string infile_name, std::string catalogfile_name, std::string region_infile, double vertex_rate, double edge_rate, std::vector<std::vector<double>> reg_rates) {
    std::fstream in_file;
    in_file.open(infile_name);
    std::vector<std::string> lines;
    std::string line;
    std::string output;
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

    std::tuple< int, std::vector<Region*>, FourDArr*> regions_tuple;

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

    if (lines[read_idx].find(substring) != std::string::npos) {
        read_idx ++;
        regions_tuple = init_regions(slice_1Dvec_str(lines, read_idx, (int)lines.size()), dims_int, num_regions, region_infile); 
        std::cout << "region!\n";
        incriment = std::get<0>(regions_tuple); temp_regions = std::get<1>(regions_tuple); temp_region_sites = std::get<2>(regions_tuple);
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
    double x;
    double y;
    double z;
    int atomtype;
    int vacancies_count = 0;

    FourDBoolArr* temp_vacancies = new FourDBoolArr(2, (size_t)dims_int[0], (size_t)dims_int[1], (size_t)dims_int[2]);
    FourDBoolArr* temp_vertex_sites = new FourDBoolArr(1, (size_t)dims_int[0], (size_t)dims_int[1], (size_t)dims_int[2]);
    FourDBoolArr* temp_bc_sites = new FourDBoolArr(1, (size_t)dims_int[0], (size_t)dims_int[1], (size_t)dims_int[2]);

    std::tuple<size_t, size_t, size_t, size_t> vacs_size_tuple = (*temp_vacancies).size_tuple;


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
    
    for (int i=read_idx; i<(int)lines.size(); i++){
        //std::cout << lines[i] << "\n";
        tuple_out = parse_line(lines[i]);
        lattice_pos = std::get<0>(tuple_out); x = std::get<1>(tuple_out); y = std::get<2>(tuple_out); 
        z = std::get<3>(tuple_out); atomtype = std::get<4>(tuple_out);
        atomtype = (int)(atomtype);

        if (lattice_pos == "v") {
            //std::cout << "lattice position v " << "\n";
            x = (int)x;
            y = (int)y;
            z = (int)z;
            if (atomtype == 0) {
                (*temp_vertex_sites)(0,x,y,z) = 0;
                (*temp_vacancies)(0,x,y,z) = 1;
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
            x = (int)(x - 0.5);
            y = (int)(y - 0.5);
            z = (int)(z - 0.5);

            if (atomtype == 0) {
                (*temp_bc_sites)(0,x,y,z) = 0;
                (*temp_vacancies)(1,x,y,z) = 1; 
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

    std::tuple< std::vector< std::vector<int> >, 
    std::vector< std::vector< std::vector<double> > >, 
    std::vector< std::vector< std::vector< std::vector<double> > > >, 
    int > 
    cat_tuple;
    if (catalogfile_name != "None") {
        cat_tuple = updated_create_ratecatalog(catalogfile_name, a_type_keys);
    }

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

    reg_energies = std::get<2>(cat_tuple);
    reg_num = std::get<3>(cat_tuple);

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
    new_configs[0] = reorder_ints_inp(new_configs[0], unsorted_idxs_8bit);
    new_energies[0] = reorder_dbl_inp_2D(new_energies[0], unsorted_idxs_8bit);

    std::vector<int> unsorted_idxs_14bit = arange(0,(int)new_configs[1].size(),1);
    auto comparator_14bit = [new_configs](int idx1, int idx2) {
                return new_configs[1][idx1] < new_configs[1][idx2];
        };
    std::sort(unsorted_idxs_14bit.begin(), unsorted_idxs_14bit.end(), comparator_14bit);
    new_configs[1] = reorder_ints_inp(new_configs[1], unsorted_idxs_14bit);
    new_energies[1] = reorder_dbl_inp_2D(new_energies[1], unsorted_idxs_14bit);

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
    std::cout << "reg_num: " << reg_num << "\n";
    Lattice* new_lattice = new Lattice(dims_int[0], dims_int[1], dims_int[2], vacancies_count, reg_num);
    //print_2Dvector(temp_region_sites->nonzero());
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

    std::cout << reg_num << "\n";

    rows = new_lattice->regionrates_111_L.rows();
    cols = new_lattice->regionrates_111_L.cols();
    for (size_t i=0; i<rows; i++) {
        std::cout << i << "\n";
        std::cout << reg_rates[i][0] << "\n";
        std::cout << reg_rates[i][1] << "\n";
        for (size_t j=0; j<cols; j++) {
            std::cout << "j:" << j << "\n";
            new_lattice->regionrates_111_L(i,j) = reg_rates[i][0];
            new_lattice->regionrates_111_R(i,j) = reg_rates[i][1];
        }
    }

    new_lattice->probs = create_vec_1D_float(vacancies_count);
    new_lattice->rates = create_vec_1D_float(vacancies_count);

    for (int i=0; i<(int)a_type_values.size(); i++) {new_lattice->a_types[a_type_keys[i]] = a_type_values[i];}

    new_lattice->region_energies = std::get<2>(cat_tuple); 

    std::vector< std::vector<int> > nonzero_vacs = new_lattice->vacancies.nonzero();

    for (int i=0; i<nonzero_vacs.size(); i++) {
        for (int j=0; j<nonzero_vacs[0].size(); j++) {
            new_lattice->vacancies_pos[i][j] = nonzero_vacs[i][j];
        }
    }

    delete temp_111_catalog;
    delete temp_100_catalog;
    delete temp_region_sites;
    delete temp_vacancies;
    delete temp_vertex_sites;
    delete temp_bc_sites;

    new_lattice->rate_cumsum.resize(14*nonzero_vacs.size());


    return new_lattice;
}

/*---------------------------------------------------------------------------*/

