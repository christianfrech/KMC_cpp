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

/*! \brief A class for storing 4-D arrays of ints*/
class FourDDoubleArr {
public:
    /*! \brief Constructor
     * \param [in] len1 len2 len3 len4    Lengths of the 4 dimensions of the array
     */

    std::vector<size_t> size_vec;
    std::vector<double> data_; ///< The data stored in the array

    FourDDoubleArr(size_t len1, size_t len2, size_t len3, size_t len4)
    : len1_(len1)
    , len2_(len2)
    , len3_(len3)
    , len4_(len4)
    , data_(len1 * len2 * len3 * len4)
    {
        //data_ = (double *)malloc(sizeof(double) * len1 * len2 * len3 * len4);
        size_vec.push_back(len1);
        size_vec.push_back(len2);
        size_vec.push_back(len3);
        size_vec.push_back(len4);
    }
    
    /*! \brief Access an element of the 4-D array
     * \param [in] i1 First index
     * \param [in] i2 Second index
     * \param [in] i3 Third index
     * \param [in] i4 Fourth index
     * \returns Reference to array element
     */
    double& operator() (size_t i1, size_t i2, size_t i3, size_t i4) {
        return data_[i1 * len2_ * len3_ * len4_ + i2 * len3_ * len4_ + i3 * len4_ + i4];
    }
    
    double  operator() (size_t i1, size_t i2, size_t i3, size_t i4) const {
        return data_[i1 * len2_ * len3_ * len4_ + i2 * len3_ * len4_ + i3 * len4_ + i4];
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

    Matrix<int> nonzero_elems() {
        int vec_size = (int)(len1_*len2_*len3_*len4_)/8;
        Matrix<double> mat_out(vec_size, 5);
        mat_out.zero();

        int elem = 0; 
        
        for (int i1=0; i1<(int)len1_; i1++) {
            for (int i2=0; i2<(int)len2_; i2++) {
                for (int i3=0; i3<(int)len3_; i3++) {
                    for (int i4=0; i4<(int)len4_; i4++) {
                        if ( (*this)((size_t)i1, (size_t)i2, (size_t)i3, (size_t)i4) != 0 ) {

                            if (elem >= (vec_size-1)) {
                                mat_out.reshape(vec_size*2, 5);
                                vec_size = vec_size * 2;
                            }
                            
                            mat_out[elem][0] = i1;
                            mat_out[elem][1] = i2;
                            mat_out[elem][2] = i3;
                            mat_out[elem][3] = i4;
                            mat_out[elem][4] = (*this)((size_t)i1, (size_t)i2, (size_t)i3, (size_t)i4);
                            elem ++;
                        }
                    }
                }
            }
        } 
        mat_out.reshape(elem, 5);
        return mat_out;
    }
    
    /*
    set all elements in data field to zero
    */
    void zero() {
        for (int i1=0; i1<len1_; i1++) {
            for (int i2=0; i2<len2_; i2++) {
                for (int i3=0; i3<len3_; i3++) {
                    for (int i4=0; i4<len4_; i4++) {
                        (*this)((size_t)i1, (size_t)i2, (size_t)i3, (size_t)i4) = 0;
                    }
                }
            }    
        }
    }

private:
    size_t len1_, len2_, len3_, len4_; ///< Dimensions of the array

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
    retrieve nonzero elements and values
    */
    std::vector< std::vector<int> > nonzero_with_values() {
        int vec_size = (len1_*len2_*len3_*len4_)/8;
        std::vector<std::vector<int>> vec_out(vec_size, std::vector<int>(5));

        int elem = 0; 
        for (int i1=0; i1<(int)len1_; i1++) {
            for (int i2=0; i2<(int)len2_; i2++) {
                for (int i3=0; i3<(int)len3_; i3++) {
                    for (int i4=0; i4<(int)len4_; i4++) {
                        if ( (*this)((size_t)i1, (size_t)i2, (size_t)i3, (size_t)i4) != 0 ) {
                            if (elem >= (vec_size-1)) {
                                vec_out.resize(vec_size*2, std::vector<int>(4));
                                vec_size = vec_size * 2;
                            }

                            vec_out[elem][0] = i1;
                            vec_out[elem][1] = i2;
                            vec_out[elem][2] = i3;
                            vec_out[elem][3] = i4;
                            vec_out[elem][4] = (*this)((size_t)i1, (size_t)i2, (size_t)i3, (size_t)i4);
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


/**
 * @class Region
 * @brief A class for storing data (rate constants) corresponding to a subset of coordinates in a FourDArr object.
 */
class Region {
public:
    int id;  ///< Unique identifier for the region.
    std::string bias;  ///< Bias direction of the region.
    std::string type;  ///< Type of region (e.g., "GB" or "BLOCK").
    std::vector<int> slopes;  ///< Slopes of the region (used for "GB" type).
    std::vector<int> shifts;  ///< Shifts of the region (used for "GB" type).
    std::vector<int> lowerbound;  ///< Lower bounds for the region (used for "BLOCK" type).
    std::vector<int> upperbound;  ///< Upper bounds for the region (used for "BLOCK" type).
    std::vector<std::vector<int>> params;  ///< Parameters defining the region.
    std::vector<std::vector<std::vector<double>>> energies;  ///< Energy values for the region.
    FourDDoubleArr region_rates_L;  ///< Left transition rates in the region.
    FourDDoubleArr region_rates_R;  ///< Right transition rates in the region.
    std::mt19937 rand_num;  ///< Random number generator for the region.
    double mean;  ///< Mean value of the energy barrier distribution.
    double stddev;  ///< Standard deviation of the energy barrier distribution.
    double barrier_min;  ///< Minimum allowed energy barrier.
    double barrier_max;  ///< Maximum allowed energy barrier.
    std::vector<int> dimensions;  ///< Dimensions of the region.
    bool random;  ///< Indicates whether rates are randomly assigned.
    std::normal_distribution<double> distribution;  ///< Normal distribution for random barrier generation.
    std::vector<double> rates;  ///< Predefined rates for the region.
    bool interface;  ///< Whether the region has an interface.
    int interface_i;  ///< Index of the interface.
    int interface_dim;  ///< Dimension of the interface.
    double interface_100_rate;  ///< Rate for the 100 interface.

    /**
     * @brief Constructs a Region object with given parameters.
     * 
     * @param id_in Unique identifier for the region.
     * @param reg_type Type of region ("GB" or "BLOCK").
     * @param bias_direc Bias direction for the region.
     * @param params_in Parameters defining the region.
     * @param distribution Distribution parameters for random barriers.
     * @param random_val Boolean indicating whether rates are randomized.
     * @param rates_in Vector of predefined rates.
     * @param interface_params Parameters for interface configuration.
     * @param interface_100_rate_in Rate for the 100 interface.
     */
    Region(int id_in, std::string reg_type, std::string bias_direc, std::vector<std::vector<int>> params_in,
           std::vector<double> distribution, bool random_val, std::vector<double> rates_in,
           std::vector<int> interface_params, double interface_100_rate_in);

    /**
     * @brief Computes the average connectivity of sites within the region.
     * 
     * @return The average connectivity of sites.
     */
    double avg_connectivity() {
        std::vector< std::vector<int> > diag_directions = {{0,0,0}, {1,0,0}, {0,1,0}, {1,1,0}, {0,0,1}, {1,0,1}, {0,1,1}, {1,1,1}};
        int degree_sum = 0;
        int site_count = 0;
        
        for (int i=0; i<(int)dimensions[0]; i++) {
            for (int j=0; j<(int)dimensions[1]; j++) {
                for (int k=0; k<(int)dimensions[2]; k++) {
                    for (int l=0; l<(int)dimensions[3]; l++) {
                        
                        if ( region_rates_L(i,j,k,l) != -1) {

                            //std::cout << "i: " << i << " j: " << j << " k: " << k << " l: " << l << "\n";
                            int connectivity = 0;
                            
                            // moving vacancy from bc site to vertex site
                            if (i == 1) {
                                for (int i=0; i<(int)diag_directions.size(); i++) {
                                    //if ((vac[3] == 0)) {/* checking for leftmost non-periodic boundary along z-axis*/}
                                    //else if ((vac[3] == (int)(dimensions[3]-1))) {/* checking for rightmost non-periodic boundary along z-axis*/}
                                    if ( region_rates_L( 0, (j + diag_directions[i][0]), (k + diag_directions[i][1]), 
                                            (l + diag_directions[i][2])) != -1) { connectivity ++; }
                                }
                            }

                            // moving vacancy from vertex site to bc site     
                            else if (i == 0) {
                                for (int i=0; i<(int)diag_directions.size(); i++) {
                                    //if ((vac[3] == 0)) {/* checking for leftmost non-periodic boundary along z-axis*/}
                                    //else if ((vac[3] == (int)(lattice_dim[2]-1))) {/* checking for rightmost non-periodic boundary along z-axis*/}
                                    if ( region_rates_L( 1, (j - diag_directions[i][0]), (k - diag_directions[i][1]), 
                                            (l - diag_directions[i][2])) != -1) { connectivity ++; }
                                }
                            }

                            degree_sum += connectivity;
                            site_count ++;
                        }
                    }
                }
            }
        }
        
        double avg_degree = (double)degree_sum / (double)site_count;
        
        return avg_degree;
    }

    /**
     * @brief Generates a random energy barrier within specified bounds.
     * 
     * @return A randomly generated barrier value.
     */
    double get_rand_barrier() {

        double barrier = distribution(rand_num);
        
        while ((barrier > barrier_max) || (barrier < barrier_min)) { barrier = distribution(rand_num); }
        
        return barrier;

    }

    /**
     * @brief Retrieves the transition rate for a given site in the region.
     * 
     * @param i First index (layer).
     * @param j Second index (x-coordinate).
     * @param k Third index (y-coordinate).
     * @param l Fourth index (z-coordinate).
     * @param LR_idx Indicates whether to retrieve from Left (1) or Right (0) rate array.
     * @return The transition rate at the specified site.
     */
    double get_rate(int i, int j, int k, int l, int LR_idx) {

        int i_new = i;
        int j_new = j - lowerbound[0];
        int k_new = k - lowerbound[1];
        int l_new = l - lowerbound[2];

        if (LR_idx) return region_rates_L(i_new, j_new, k_new, l_new);
        else return region_rates_R(i_new, j_new, k_new, l_new);

    }

    /**
     * @brief Randomly blocks sites within the region until the connectivity is within a desired range.
     */
    void random_blocking() {
        int i_new; int j_new; int k_new; int l_new;
        std::random_device dev;
        std::mt19937 rng(dev());
        std::uniform_int_distribution<std::mt19937::result_type> i_rand(0,1); 
        std::uniform_int_distribution<std::mt19937::result_type> j_rand(0,(int)(dimensions[1]));
        std::uniform_int_distribution<std::mt19937::result_type> k_rand(0,(int)(dimensions[2])); 
        std::uniform_int_distribution<std::mt19937::result_type> l_rand(0,(int)(dimensions[3])); 

        double connectivity = avg_connectivity();
        std::cout << "connectivity: " << connectivity << "\n";
        
        while ((connectivity < 5) || (7 < connectivity)) {

            i_new = i_rand(rng); 
            j_new = j_rand(rng); 
            k_new = k_rand(rng); 
            l_new = l_rand(rng); 
            //std::cout << "i_new: " << i_new << " j_new: " << j_new << " k_new: " << k_new << " l_new: " << l_new << "\n";

            region_rates_L(i_new, j_new, k_new, l_new) = -1;
            region_rates_R(i_new, j_new, k_new, l_new) = -1;

            connectivity = avg_connectivity();
        }
    }

    /**
     * @brief Generates random barriers according to a distribution, checks if they are within bounds,
     * and assigns them to sites in the region rate map.
     * 
     * @param reg_rates A vector of rates to be assigned within the region.
     */
    void random_barrier_assigner(std::vector<double> reg_rates) {

        int x_range = dimensions[1];
        int y_range = dimensions[2];
        int z_range = dimensions[3];

        std::cout << "reg_rates: \n";
        print_1Dvector(reg_rates);

        double rand_barrier;
        double rate;
        double EULER = 2.71828182845904523536;

        int direc;
        if (reg_rates[0] > reg_rates[1]) direc = 0;
        else direc = 1;

        std::cout << "x_range: " << x_range << " y_range " << y_range << " z_range " << z_range << "\n";

        // iterating over region sites
        for (int i=0; i<2; i++) {
            for (int j=0; j<x_range; j++) {
                for (int k=0; k<y_range; k++) {
                    for (int l=0; l<z_range; l++) {
                        if (region_rates_L(i,j,k,l) != -1) {
                            
                            rand_barrier = get_rand_barrier();
                            rate = 5e12 * std::pow(EULER, (-rand_barrier / (300* 8.6173e-5)));
                            
                            if (direc == 0) {

                                if (rand_barrier < 0) {
                                    region_rates_L(i,j,k,l) =  rate;
                                    region_rates_R(i,j,k,l) =  5e12;
                                }
                                else {
                                    region_rates_L(i,j,k,l) =  5e12;
                                    region_rates_R(i,j,k,l) =  rate;
                                }
                            }
                            else {

                                if (rand_barrier < 0) {
                                    region_rates_L(i,j,k,l) =  5e12;
                                    region_rates_R(i,j,k,l) =  rate;
                                }
                                else {
                                    region_rates_L(i,j,k,l) =  rate;
                                    region_rates_R(i,j,k,l) =  5e12;
                                }
                            }
                        }
                    }
                }
            }
        }
    }
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
        std::vector<std::vector<int>> onevac_vec;
        std::vector<std::vector<int>> non_onevac;
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
        int void_threshold;
        double void_barrier;
        double terrace_barrier_100;
        double terrace_barrier_111;
        double void_gb_diss_barrier;

        Lattice(int xdim, int ydim, int zdim, int num_vacancies, int num_regions, std::vector<Region*> regs_in):
            regions(regs_in),
            vertex_sites((size_t)1, (size_t)xdim, (size_t)ydim, (size_t)zdim),
            vacancies((size_t)2, (size_t)xdim, (size_t)ydim, (size_t)zdim),
            bc_sites((size_t)1, (size_t)xdim, (size_t)ydim, (size_t)zdim),
            region_sites((size_t)2, (size_t)xdim, (size_t)ydim, (size_t)zdim),
            probs(num_vacancies),
            rates(num_vacancies),
            moves_coords((14 * num_vacancies), 4),
            moves_shifts((14 * num_vacancies), 3),
            moves_lattice((14 * num_vacancies), 1),
            moves_vacs((14 * num_vacancies), 1),
            ratecatalog_111(2, exp<int>(2,8)),
            
            ratecatalog_100(2, exp<int>(2,14)),
            regionrates_111_L(num_regions, exp<int>(2,8)),
            regionrates_111_R(num_regions, exp<int>(2,8)),
            regionrates_100_L(num_regions, exp<int>(2,14)),
            regionrates_100_R(num_regions, exp<int>(2,14)),
            configs_111(1, exp<int>(2,8)),
            configs_100(1, exp<int>(2,14)),
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

        
        /**
        * @brief Wrapper function to assign rates to region-specific rate catalogs.
        * 
        * @param temp_regions A vector of pointers to Region objects representing different regions.
        * @param misc_rates A constant reference to a vector of double values representing miscellaneous rates.
        */
        void assign_region_rates_wrapper(std::vector<Region*>& temp_regions, const std::vector<double>& misc_rates) {
            for (int i=0; i<(int)temp_regions.size(); i++) {
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
            
            if ((lattice == 2) || (lattice == 3)) {
                if (lattice == 2) { direc_sign_NN = 1; }
                else if (lattice == 3) { direc_sign_NN = 1; }
                i1 = i;
                i2 = (((j + edge_directions[s][0]) % lattice_dim[0] + lattice_dim[0]) % lattice_dim[0]);
                i3 = (((k + edge_directions[s][1]) % lattice_dim[1] + lattice_dim[1]) % lattice_dim[1]);
                i4 = (((l + edge_directions[s][2]) % lattice_dim[2] + lattice_dim[2]) % lattice_dim[2]); 
            }
            else if ((lattice == 0) || (lattice == 1)) {
                if (lattice == 0) { i1 = 1; direc_sign_NN =  1; }
                else if (lattice == 1) { i1 = 0; direc_sign_NN = -1; }
                i2 = (((j + direc_sign * diag_directions[s][0]) % lattice_dim[0] + lattice_dim[0]) % lattice_dim[0]);
                i3 = (((k + direc_sign * diag_directions[s][1]) % lattice_dim[1] + lattice_dim[1]) % lattice_dim[1]);
                i4 = (((l + direc_sign * diag_directions[s][2]) % lattice_dim[2] + lattice_dim[2]) % lattice_dim[2]);
            }

            //std::cout << " i: " << i << " j: " << j << " k: " << k << " l: " << l << "\n";
            //std::cout << " i1: " << i1 << " i2: " << i2 << " i3: " << i3 << " i4: " << i4 << "\n";
            int NN_count = 0;
            for (int s2=0; s2 < (int)diag_directions.size(); s2++) {
                i1_NN = !i1;
                i2_NN = (((i2 + direc_sign_NN * diag_directions[s2][0]) % lattice_dim[0] + lattice_dim[0]) % lattice_dim[0]);
                i3_NN = (((i3 + direc_sign_NN * diag_directions[s2][1]) % lattice_dim[1] + lattice_dim[1]) % lattice_dim[1]);
                i4_NN = (((i4 + direc_sign_NN * diag_directions[s2][2]) % lattice_dim[2] + lattice_dim[2]) % lattice_dim[2]);

                if (((lattice == 0) || (lattice == 1)) && (i4 == 0) && (i1 == 0) && (diag_directions[s2][2] == 1)) {/* checking for leftmost non-periodic boundary along z-axis*/}
                else if (((lattice == 0) || (lattice == 1)) && (i4 == (int)(lattice_dim[2]-1)) && (i1 == 1) && (diag_directions[s2][2] == 1)) {/* checking for rightmost non-periodic boundary along z-axis*/}
                else if (((lattice == 2) || (lattice == 3)) && (i4 == 0) && (edge_directions[s][2] == -1)) {}  
                else if (((lattice == 2) || (lattice == 3)) && (i4 == (int)(lattice_dim[2]-1)) && (edge_directions[s][2] == 1)) {}
                else {
                    //std::cout << " i1_NN: " << i1_NN << " i2_NN: " << i2_NN << " i3_NN: " << i3_NN << " i4_NN: " << i4_NN << "\n";
                    if ((i1_NN == i) && (i2_NN == j) && (i3_NN == k) && (i4_NN == l)) {}
                    else if (vacancies(i1_NN,i2_NN,i3_NN,i4_NN)) {NN_count++;} // std::cout << "incriment \n";}
                }
            }

            return NN_count;
        }
        

        /**
        * @brief Computes all possible atomic moves in the lattice.
        *
        * This function identifies all potential atomic moves within the system,
        * calculates their corresponding rates, and stores them for further processing.
        */
        void new_get_actions() {
            int curr_move_num = 0; // total number of moves at this current timestep  
            double rate; 
            int vacs_on_interface = 0; // vacancies at last z-index of lattice (used to calculate rate for stripping)
            int num_interface_sites = lattice_dim[0] * lattice_dim[1]; // total number of sites at last z-index 

            std::vector<int> moves(2);

            rates.resize((int)moves_coords.rows());
            int i=0; int j=0; int k=0; int l=0;
            int bulk_rate_count = 0;
            int NN_vac = 0;
            int NN_newsite = 0;

            // looping over all vacancies in system
            for (int idx=0; idx < (int)vacancies_pos.rows(); idx++) {
                // position in lattice of vacancy 
                i = vacancies_pos[idx][0];
                j = vacancies_pos[idx][1];
                k = vacancies_pos[idx][2];
                l = vacancies_pos[idx][3];
                //std::cout << "idx: " << idx << " i: " << i << " j: " << j << " k: " << k << " l: " << l << "\n";

                if ((curr_move_num + (num_interface_sites - vacs_on_interface)) >= ((int)moves_shifts.rows() - 20)) {
                    // resizing data structures to accommodate all moves 

                    int newsize = 2 * moves_shifts.rows();
                    rate_cumsum.resize(newsize);
                    moves_coords.reshape(newsize, 4);
                    moves_shifts.reshape(newsize, 3);
                    moves_lattice.reshape(newsize, 1);
                    moves_vacs.reshape(newsize, 1);
                }

                NN_vac = get_NN_count(vacancies_pos[idx], i); 

                /*
                if (NN_vac < 3) {
                    onevac_vec.push_back({i,j,k,l});
                }
                else {
                    non_onevac.push_back({i,j,k,l});
                }
                */
                
            
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

                            NN_newsite = get_NNcountofNN(i, j, k, l, -1, s, 0);

                            // getting rate corresponding to move
                            //if ((NN_vac >= void_threshold) && (NN_newsite < void_threshold)) { rate = void_barrier; }  
                            if ((NN_vac >= void_threshold) && (NN_newsite >= void_threshold)) { rate = terrace_barrier_111; }
                            else { 
                                rate = new_get_rateconstants(moves_coords[curr_move_num], moves_shifts[curr_move_num], moves_lattice[curr_move_num][0], NN_vac, NN_newsite); 
                                bulk_rate_count++;
                            }

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

                            //NN_vac = get_NNcountofNN(i, j, k, l, 1, s, 1);
                            //std::cout << "diag_directions[" << s << "][i] \n";
                            //std::cout << "NN_vac: " << NN_vac << "\n";

                            
                            // getting rate corresponding to move
                            NN_newsite = get_NNcountofNN(i, j, k, l, 1, s, 1);

                            // getting rate corresponding to move
                            //if ((NN_vac >= void_threshold) && (NN_newsite < void_threshold)) { rate = void_barrier; }  
                            if ((NN_vac >= void_threshold) && (NN_newsite >= void_threshold)) { rate = terrace_barrier_111; }
                            else { 
                                rate = new_get_rateconstants(moves_coords[curr_move_num], moves_shifts[curr_move_num], moves_lattice[curr_move_num][0], NN_vac, NN_newsite); 
                                bulk_rate_count++;
                            }

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
                            NN_newsite = get_NNcountofNN(i, j, k, l, 1, s, 2);
                            //std::cout << "edge_directions[" << s << "][i] \n";
                            //std::cout << "NN_vac: " << NN_vac << "\n";

                        }
                        else if (i == 1) {
                            moves_lattice[curr_move_num][0] = 3;
                            NN_newsite = get_NNcountofNN(i, j, k, l, 1, s, 3);
                            //std::cout << "edge_directions[" << s << "][i] \n";
                            //std::cout << "NN_vac: " << NN_vac << "\n";
                        }
                        
                        // getting rate corresponding to move
                        //NN_newsite = get_NNcountofNN(i, j, k, l, -1, s, 0);

                        // getting rate corresponding to move
                        //if ((NN_vac >= void_threshold) && (NN_newsite < void_threshold)) { rate = void_barrier; }  
                        if ((NN_vac >= void_threshold) && (NN_newsite >= void_threshold)) { rate = terrace_barrier_100; }
                        else { 
                            rate = new_get_rateconstants(moves_coords[curr_move_num], moves_shifts[curr_move_num], moves_lattice[curr_move_num][0], NN_vac, NN_newsite); 
                            bulk_rate_count++;
                        }

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
            //exit(0);
        }
        
        /**
        * @brief Function for finding the rate corresponding to NN-encoding and move type.
        *
        * @param coord Pointer to an array representing the coordinates.
        * @param shift Pointer to an array representing the shift in coordinates.
        * @param lattice The lattice type.
        * @param curr_NN The current nearest neighbor count.
        * @param new_NN The new nearest neighbor count.
        * @return The rate constant for the move.
        */
        double new_get_rateconstants(int* coord, int* shift, int lattice, int curr_NN, int new_NN) {  
            
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

                if ((lattice == regions[(reg_id-1)]->interface_i) && 
                        (regions[(reg_id-1)]->interface) && (shift[regions[(reg_id-1)]->interface_dim] != 0)) {
                    rate = regions[(reg_id-1)]->interface_100_rate;
                }
                else if (regions[(reg_id-1)]->random) {
                    if ( regions[(reg_id-1)]->get_rate(coord[0], coord[1], coord[2], coord[3], LR_idx) == -1) rate = -1;
                    
                    else if ((curr_NN >= void_threshold) && (new_NN < void_threshold)) { 
                        rate = void_gb_diss_barrier; 
                    }  
                    else {
                        if (LR_idx == 1) rate = regions[(reg_id-1)]->get_rate(coord[0], coord[1], coord[2], coord[3], LR_idx);
                        if (LR_idx == 0) rate = regions[(reg_id-1)]->get_rate(coord[0], coord[1], coord[2], coord[3], LR_idx);
                    }
                }
                else {
                    if ((curr_NN >= void_threshold) && (new_NN < void_threshold)) { 
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

                if ((curr_NN >= void_threshold) && (new_NN < void_threshold)) { rate = void_barrier; }  
                else {
                    // in case of no pre-defined region, use bulk rate constants
                    if ((lattice == 1) || (lattice == 0)) {rate = ratecatalog_111[LR_idx][idx];}
                    else if ((lattice == 2) || (lattice == 3)) {rate = ratecatalog_100[0][idx];}
                }
            }

            return rate;
        }

        /**
        * @brief Method for obtaining the number of vacancies in the nearest neighbor shell of a vacancy.
        *
        * @param vac Pointer to an array representing the vacancy coordinates.
        * @param lattice The lattice type.
        * @return The number of vacancies in the nearest neighbor shell.
        */
        int get_NN_count(int* vac, int lattice) {
            int count = 0;
            
            
            // moving vacancy from bc site to vertex site
            if (lattice == 1) {
                for (int i=0; i<(int)diag_directions.size(); i++) {
                    if ((vac[3] == 0) && (diag_directions[i][2] == 1)) {/* checking for leftmost non-periodic boundary along z-axis*/}
                    else if ((vac[3] == (int)(lattice_dim[2]-1)) && (diag_directions[i][2] == 1)) {/* checking for rightmost non-periodic boundary along z-axis*/}
                    else { count += vacancies(0, (((vac[1] + diag_directions[i][0]) % lattice_dim[0] + lattice_dim[0]) % lattice_dim[0]), (((vac[2] + diag_directions[i][1]) % lattice_dim[1] + lattice_dim[1]) % lattice_dim[1]), (((vac[3] + diag_directions[i][2]) % lattice_dim[2] + lattice_dim[2]) % lattice_dim[2])); }
                    
                }
            }

            // moving vacancy from vertex site to bc site     
            else if (lattice == 0) {
                for (int i=0; i<(int)diag_directions.size(); i++) {
                    if ((vac[3] == 0) && (diag_directions[i][2] == 1)) {/* checking for leftmost non-periodic boundary along z-axis*/}
                    else if ((vac[3] == (int)(lattice_dim[2]-1)) && (diag_directions[i][2] == 1)) {/* checking for rightmost non-periodic boundary along z-axis*/}
                    else { count += vacancies(1, (((vac[1] - diag_directions[i][0]) % lattice_dim[0] + lattice_dim[0]) % lattice_dim[0]), (((vac[2] - diag_directions[i][1]) % lattice_dim[1] + lattice_dim[1]) % lattice_dim[1]), (((vac[3] - diag_directions[i][2])% lattice_dim[2] + lattice_dim[2]) % lattice_dim[2])); }
                    
                }
            }
            
            return count;
        }

        /**
        * @brief Deprecated method corresponding to the case of only one rate constant.
        *
        * @param coord Pointer to an array representing the coordinates.
        * @param shift Pointer to an array representing the shift in coordinates.
        * @param lattice The lattice type.
        * @return The rate constant for the move.
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

        /**
        * @brief Deprecated method for finding energy corresponding to NN-encoding and move type.
        *
        * @param encoding The encoding value.
        * @param coord Pointer to an array representing the coordinates.
        * @param shift Pointer to an array representing the shift in coordinates.
        * @param lattice The lattice type.
        * @return The energy corresponding to the move.
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

        /**
        * @brief Computes the encoding corresponding to the configuration of nearest neighbors for a vacancy.
        * 
        * This function calculates a unique encoding by summing the contributions of neighboring sites
        * based on their occupancy, using exponential weights. The encoding depends on whether the vacancy
        * is in a boundary condition (BC) site or a vertex site.
        * 
        * @param vac Pointer to an array of size 3 representing the coordinates of the vacancy.
        * @param lattice Integer indicating the type of lattice site:
        *               - 3: Vacancy moves from BC site to BC site
        *               - 2: Vacancy moves from vertex site to vertex site
        *               - 1: Vacancy moves from BC site to vertex site
        *               - 0: Vacancy moves from vertex site to BC site
        * 
        * @return Integer encoding of the local configuration around the vacancy.
        */
        int new_get_neighbors(int* vac, int lattice) {
            int sum = 0;
            int m = (int)a_types.size(); // number of atom types in system

            // moving vacancy from bc site to bc site 
            if (lattice == 3) {
                for (int i=0; i<(int)diag_directions.size(); i++) {
                    sum += exp<int>(m,i) * vertex_sites(0, (((vac[0] - diag_directions[i][0]) % lattice_dim[0] + lattice_dim[0]) % lattice_dim[0]), (((vac[1] - diag_directions[i][1]) % lattice_dim[1] + lattice_dim[1]) % lattice_dim[1]), (((vac[2] - diag_directions[i][2]) % lattice_dim[2] + lattice_dim[2]) % lattice_dim[2]));
                }
                for (int i=0; i<(int)edge_directions.size(); i++) {
                    sum += exp<int>(m,i) *  bc_sites(0, (((vac[0] + edge_directions[i][0]) % lattice_dim[0] + lattice_dim[0]) % lattice_dim[0]), (((vac[1] + edge_directions[i][1]) % lattice_dim[1] + lattice_dim[1]) % lattice_dim[1]), (((vac[2] + edge_directions[i][2]) % lattice_dim[2] + lattice_dim[2]) % lattice_dim[2]));
                }
            }

            // moving vacancy from vertex site to vertex site 
            else if (lattice == 2) {
                for (int i=0; i<(int)diag_directions.size(); i++) {
                    sum +=  exp<int>(m,i) * bc_sites(0, (((vac[0] + diag_directions[i][0])  % lattice_dim[0] + lattice_dim[0]) % lattice_dim[0]), (((vac[1] + diag_directions[i][1]) % lattice_dim[1] + lattice_dim[1]) % lattice_dim[1]), (((vac[2] + diag_directions[i][2]) % lattice_dim[2] + lattice_dim[2]) % lattice_dim[2]));
                }
                for (int i=0; i<(int)edge_directions.size(); i++) {
                    sum += exp<int>(m,i) *  vertex_sites(0, (((vac[0] + edge_directions[i][0]) % lattice_dim[0] + lattice_dim[0]) % lattice_dim[0]), (((vac[1] + edge_directions[i][1]) % lattice_dim[1] + lattice_dim[1]) % lattice_dim[1]), (((vac[2] + edge_directions[i][2]) % lattice_dim[2] + lattice_dim[2]) % lattice_dim[2]));
                }
            }

            // moving vacancy from bc site to vertex site
            else if (lattice == 1) {
                for (int i=0; i<(int)diag_directions.size(); i++) {
                    sum += exp<int>(m,i) * vertex_sites(0, (((vac[0] - diag_directions[i][0]) % lattice_dim[0] + lattice_dim[0]) % lattice_dim[0]), (((vac[1] - diag_directions[i][1]) % lattice_dim[1] + lattice_dim[1]) % lattice_dim[1]), (((vac[2] - diag_directions[i][2]) % lattice_dim[2] + lattice_dim[2]) % lattice_dim[2]));
                }
            }

            // moving vacancy from vertex site to bc site     
            else if (lattice == 0) {
                for (int i=0; i<(int)diag_directions.size(); i++) {
                    sum += exp<int>(m,i) * bc_sites(0, (((vac[0] + diag_directions[i][0]) % lattice_dim[0] + lattice_dim[0]) % lattice_dim[0]), (((vac[1] + diag_directions[i][1]) % lattice_dim[1] + lattice_dim[1]) % lattice_dim[1]), (((vac[2] + diag_directions[i][2])% lattice_dim[2] + lattice_dim[2]) % lattice_dim[2]));
                }
            }

            return sum;
        }
 
        /**
        * @brief Updates the positions of atoms on the lattice according to the selected move.
        * 
        * This function updates the lattice by moving vacancies according to predefined move types.
        * It modifies the lattice occupancy, updates vacancy positions, and manages boundary conditions.
        * 
        * @param idx Index of the selected move in the list of possible moves.
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

        /**
        * @brief Calculates the time elapsed for a move in the KMC simulation.
        * 
        * This function generates a random number between 0 and 1 and computes the 
        * time elapsed for a move based on the cumulative sum of move rates.
        * The time is determined using the exponential distribution formula.
        * 
        * @param cumsum A vector representing the cumulative sum of move rates.
        * @return The calculated time elapsed for the move.
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

        /**
        * @brief Selects a random move index from a cumulative sum vector, 
        *        with probability proportional to the corresponding rate constant.
        * 
        * This function generates a random number and uses it to access an index 
        * in the cumulative sum vector (`cumsum`), ensuring that selection probability 
        * is proportional to the value at that position.
        * 
        * @param cumsum A vector of cumulative sum values representing move probabilities.
        * @return The index of the selected move.
        */
        int get_idx(std::vector<double> cumsum) {
            // creating a random double between 0 and 1
            unsigned int random = mt_obj();
            double random_double = ((random / (1.+ UINT32_MAX)) +  (1 / (1.+ UINT32_MAX)));

            // accessing random element into cumulative sum array, 
            // probability of access proportional to the value at that point
            int last_idx = (int)rate_cumsum.size() - 1;
            double rand_pos = rate_cumsum[last_idx] * random_double;
            int min_idx = searchsorted_recursive(&rate_cumsum, rand_pos, 0, last_idx);

            return min_idx;
        }


        /**
        * @brief Runs a Kinetic Monte Carlo (KMC) simulation for a given time limit.
        * 
        * This function initializes the necessary system variables, updates the state of 
        * the system at each step, records vacancy configurations, and logs simulation progress.
        * It executes a loop until the simulated time exceeds `time_lim`.
        * 
        * @param time_lim The maximum simulation time.
        * @param start The start time of the simulation (real clock time).
        * @param folder The output folder for saving simulation results.
        * @param iteration The current iteration number of the simulation.
        * @param last_tick The last recorded tick of the simulation.
        * @param last_time The last recorded simulation time.
        * @return A tuple containing:
        *         - A trajectory of all vacancy configurations.
        *         - A vector with move counts of different move types.
        *         - A vector with time elapsed for different move types.
        *         - A vector of all recorded time steps.
        */
        std::tuple< std::vector< std::vector< std::vector<int> > >, std::vector<int>, 
        std::vector<double>, std::vector<double> > new_kmc_iterator(double time_lim, std::chrono::system_clock::time_point start, std::string folder, int iteration, int last_tick, double last_time) {
            fprintf(stdout, "%s", "beginning kmc iterations \n\n");
            
            std::cout << "last_time: " << last_time << "\n";
            double t = last_time;
            int old_time;

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
  

            only_vacancies = vacancies.nonzero();
            all_vacancies.push_back(only_vacancies);
            std::cout << "all_vacancies done \n"; 

            int move_ticks = last_tick;

            // output files 
            std::ofstream out_file;
            std::ostringstream ss;
            std::string output_filename;
            std::string times_filename;
            std::string count_filename;


            
            int idx = 0;

            while (t < time_lim) {
                //std::cout << "move_ticks: " << move_ticks << "\n";
                //std::cout << "rate_cumsum.size: " << rate_cumsum.size() << "\n";
                //std::cout << "rate_cumsum[-1]: " << rate_cumsum[(rate_cumsum.size()-1)] << "\n";
                //print_1Dvector(rate_cumsum);
                //std::cout << "onevac_vec: \n";
                //print_2Dvector(onevac_vec);
                /*
                if (move_ticks >= 0) {
                    std::cout << "onevac_vec: \n";
                    print_2Dvector(onevac_vec);
                    ss << folder << "/onevoid_1NN_test.txt";
                    output_filename = ss.str();
                    write_to_file(output_filename, onevac_vec);
                    ss.str("");
                    ss.clear();
                    ss << folder << "/non_onevoid.txt";
                    output_filename = ss.str();
                    write_to_file(output_filename, non_onevac);
                    ss.str("");
                    ss.clear();
                    exit(0);
                }
                */
                
                // writing output files every 100000 timesteps 
                if (move_ticks % 20000 == 0) {
                    std::cout << "move_ticks: " << move_ticks << "\n";

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

                /*
                // terminating simulation after real-time limit reached
                if (elapsed_seconds.count() >= 900) {
                    only_vacancies = vacancies.nonzero();
                    all_vacancies.push_back(only_vacancies);
                    std::cout << "t: " << t << "\n";
                    std::cout << "move_counts: " << "\n";
                    print_1Dvector(move_counts);
                    std::cout << "time_count: " << "\n";
                    print_1Dvector(time_count);
                    t = time_lim + 1;
                    break;
                    exit(0);
                }
                */
                 
                new_get_actions();
                
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

/**
 * @brief Generates all binary configurations of length `len` with `m` zeros.
 * 
 * This function creates entries in a rate catalog corresponding to all possible 
 * binary strings of length `len` containing exactly `m` zeros. The function 
 * computes all such combinations and stores their integer representations in a vector.
 * 
 * @param len The length of the binary string.
 * @param m The number of zeros in the binary string.
 * @param size The base value used for exponentiation.
 * @return A vector containing integer representations of all valid binary configurations.
 */
std::vector<int> bin_m_zeros(int len, int m, int size) {
    int num_combos = NCR(len, m);
    std::vector<int> vec(num_combos);
    int smallest_bin = 0;
    unsigned int t=0;
    unsigned int v=0;
    unsigned int w=0;

    for (int i=0; i<(len-m); i++) {
        smallest_bin += exp<int>((int)(size), i);
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
 * @brief Converts a base-m encoded string to an integer.
 * 
 * This function takes a string representation of a number in base-m 
 * and converts it to an integer. Each character in the string represents 
 * a digit in base-m, and the function calculates its decimal equivalent.
 * 
 * @param config The input string representing the number in base-m.
 * @param size The base value (m) used for conversion.
 * @return The integer equivalent of the base-m number.
 */
int base_m_to_int(std::string config, int size) {

    int result = 0;
    std::vector<char> toks = split_by_char(config);
    char tok;


    for (int i=0; i<(int)toks.size(); i++) {
        tok = toks[i];
        result += (int)((tok) - '0') * exp<int>(size, i);
    }
        
    return result;
}

/**
 * @brief Reads a rate catalog file and creates data structures for allowed moves, configurations, and region-specific rates.
 * 
 * This function parses a file containing rate catalog information, extracting configurations, 
 * associated energies, and region-based rate assignments. The output consists of:
 * - A catalog of allowed configurations.
 * - A catalog of configuration energies.
 * - A region-specific rate catalog.
 * - The total number of regions processed.
 * 
 * @param catalogfile The file path to the rate catalog.
 * @param atype_list A vector containing atom type identifiers.
 * @return A tuple containing:
 *         - A 2D vector of integers representing allowed configurations.
 *         - A 3D vector of doubles representing energies for each configuration.
 *         - A 4D vector of doubles representing region-specific energy configurations.
 *         - An integer indicating the number of regions.
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
    types = slice_1Dvec(types, 1, (int)types.size());

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
    
    lines = slice_1Dvec(lines, 1, (int)lines.size());

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
            toks = slice_1Dvec(toks, 1, (int)toks.size());
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
            toks = slice_1Dvec(toks, 1, (int)toks.size());
            k = j;

            for (int l=0; l<(int)toks.size(); l++) {
                tok = toks[l];
                dft_energies[0][k] = std::stof(tok);
                k ++;
            }
        }

        else if (configline.find("R:") != std::string::npos) {
            toks = tokenizer(configline, " ");
            toks = slice_1Dvec(toks, 1, (int)toks.size());

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

    lines = slice_1Dvec(lines, lines_read, (int)lines.size());

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

            toks = slice_1Dvec(toks, 1, (int)toks.size());

            for (j=0; j<(int)toks.size(); j++) {catalog_out[j] = std::stof(toks[j]);}

            catalog_out = slice_1Dvec(catalog_out, 1, (int)catalog_out.size());
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

/**
 * @brief Creates a Region object based on the provided input information.
 * 
 * This function parses the given input vector to extract details about a region,
 * including its type, bias direction, position parameters, rates, and optional
 * properties such as random distribution or interface properties.
 * 
 * @param info A vector of strings containing region details from the input file.
 * @return A pointer to the newly created Region object.
 */
Region* add_region(std::vector<std::string> info) {
    std::cout << "adding region \n";
    int id = std::stoi(tokenizer(info[0], ":")[0]); // region id number
    std::vector< std::vector<int> > params = vect_create_2D(2,3);
    std::vector<double> rates(2);
    std::vector<double> distribution(4,1);
    std::vector<int> interface(4);
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

    Region* new_region = new Region(id, reg_type, bias, params, distribution, random, rates, interface, interface_terrace_rate);

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

/**
 * @brief Populates the `region_sites` FourDArr with values corresponding to a custom regions input file.
 *
 * This function reads a custom region file and assigns region IDs to the appropriate locations 
 * within the `FourDArr` structure. It ensures that assigned sites do not exceed the simulation cell bounds.
 *
 * @param sites Pointer to the `FourDArr` structure where region IDs will be assigned.
 * @param regions A vector of pointers to `Region` objects representing different regions.
 * @param custom_reg_idx Index of the custom region to be processed.
 * @param dim A vector containing the dimensions [x, y, z] of the region space.
 * @param infile_name The name of the input file containing custom region definitions.
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

/**
 * @brief Populates the `region_sites` FourDArr with values corresponding to the input file.
 *
 * This function iterates over region objects and assigns region IDs to the appropriate locations 
 * within the `FourDArr` structure based on the type of region. It supports both grain boundary (GB) 
 * and block (rectangular prism) regions.
 *
 * @param sites Pointer to the `FourDArr` structure where region IDs will be assigned.
 * @param regions A vector of pointers to `Region` objects that define different regions.
 * @param dim A vector containing the dimensions [x, y, z] of the region space.
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
        std::cout << "region number: " << i << "\n";

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
            
            std::cout << "lowerbound[0]: " << region->lowerbound[0] << " lowerbound[1]: " << region->lowerbound[1] << "lowerbound[2]: " << region->lowerbound[2] << "\n";
            std::cout << "upperbound[0]: " << region->upperbound[0] << " upperbound[1]: " << region->upperbound[1] << "upperbound[2]: " << region->upperbound[2] << "\n";

            start = {0, lo[0], lo[1], lo[2]}; end = {0, hi[0], hi[1], hi[2]};
            coords = FourD_idxs(start, end);
            values = create_vec_1D((int)coords.size(),region->id);
            sites->assign_idxs(coords, values);
            start = {1, lo[0], lo[1], lo[2]}; end = {1, hi[0], hi[1], hi[2]};
            coords = FourD_idxs(start, end);
            sites->assign_idxs(coords, values);
        }
    }
    //exit(0);
}

/**
 * @brief Initializes regions by reading input file data and populating a 4D array.
 *
 * This function processes an input file to create region objects and populate a `FourDArr` structure. 
 * It first initializes a 4D array with zeros, then reads the input file to define regions, adding extra 
 * regions if necessary. The function also assigns region IDs to appropriate locations in the array.
 *
 * @param lines A vector of strings representing lines from an input file.
 * @param dims A vector containing the dimensions of the region space [x, y, z].
 * @param num_regions The expected number of regions to initialize.
 * @param region_infile Path to an optional file specifying custom region definitions.
 * @return A tuple containing the updated read index, a vector of region pointers, and a pointer to a FourDArr structure.
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
    int region_idx = 0;
    while (curr_line.find("regions end") == std::string::npos) {
        region_info = tokenizer(curr_line, " ");
        
        std::cout << "region_info\n";
        print_1Dvector(region_info);
        
        Region* region = add_region(region_info);
        regions.push_back(region);
        
        read_idx ++;
        curr_line = lines[read_idx];
        region_idx ++;
        std::cout << "End of loop\n";
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

/**
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
Lattice* populate_lattice(std::string infile_name, std::string catalogfile_name, std::string region_infile) {
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

    std::cout << "file read!\n";
    // parsing first line to grab dimensions of lattice ###
    std::string dims = lines[read_idx]; //getting dimension line
    read_idx ++;

    std::vector<std::string> dims_str = tokenizer(dims," "); 
    std::vector<int> dims_int(3);

    for (int i=0; i<(int)dims_int.size(); i++) {
        dims_int[i] = std::stoi(dims_str[i+1]); 

    }

    std::cout << "dims int\n";
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
        regions_tuple = init_regions(slice_1Dvec(lines, read_idx, (int)lines.size()), dims_int, num_regions, region_infile); 
        std::cout << "region!\n";
        incriment = std::get<0>(regions_tuple); temp_regions = std::get<1>(regions_tuple); temp_region_sites = std::get<2>(regions_tuple);
        read_idx += incriment;  
    }
    else {
        printf("ERROR: regions section mis-formatted in geometry file (check for extra newlines)");
        throw std::exception();
    }
     
    //read_idx = read_idx + incriment;
    read_idx ++;

    // getting values of miscellaneous rates
    std::cout << "pre read_misc_rates()\n";
    std::cout << "lines[read_idx]: " << lines[read_idx] << "\n";
    std::tuple< int, std::vector<double> > misc_rates_tuple;
    std::string rates_substring = "rates begin";
    std::vector<double> misc_rates;

    if (lines[read_idx].find(rates_substring) != std::string::npos) {
        read_idx ++;     
        misc_rates_tuple = read_misc_rates(read_idx, lines);
        read_idx = std::get<0>(misc_rates_tuple); misc_rates = std::get<1>(misc_rates_tuple); 
    }
    read_idx ++;     

    std::cout << "lines[read_idx]: " << lines[read_idx] << "\n"; 

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
    
    for (int i=read_idx; i<(int)lines.size(); i++) {
        //std::cout << lines[i] << "\n";
        tuple_out = parse_line(lines[i]);
        lattice_pos = std::get<0>(tuple_out); x = std::get<1>(tuple_out); y = std::get<2>(tuple_out); 
        z = std::get<3>(tuple_out); atomtype = std::get<4>(tuple_out); 
        atomtype = (int)(atomtype);

        if (lattice_pos == "v") {
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

    Matrix<double>* temp_111_catalog = new Matrix<double>((size_t)2, (size_t)exp<int>(2, 8));
    Matrix<double>* temp_100_catalog = new Matrix<double>((size_t)2, (size_t)exp<int>(2, 14));
    Matrix<int> configs_111((size_t)1, (size_t)exp<int>(2, 8));
    Matrix<int> configs_100((size_t)1, (size_t)exp<int>(2, 14)); 
    
    for (int i=0; i<9; i++) {
        temp_vec = bin_m_zeros(8, i, (int)a_type_values.size());
        diag_configs.insert(diag_configs.end(), temp_vec.begin(), temp_vec.end() );
        if (i==0) {
            temp_vec_2D = vect_create_2D_float(2, NCR(8,i), misc_rates[0]);
            diag_E = temp_vec_2D;
        }
        else {
            temp_vec_2D = vect_create_2D_float(2, NCR(8,i), misc_rates[0]);
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
            temp_vec_2D = vect_create_2D_float(2, NCR(14,i),  misc_rates[1]);
            lateral_E = temp_vec_2D;
        }
        else {
            temp_vec_2D = vect_create_2D_float(2, NCR(14,i),  misc_rates[1]);
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
    Lattice* new_lattice = new Lattice(dims_int[0], dims_int[1], dims_int[2], vacancies_count, num_regions, temp_regions);

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

    // assigning rates to region-specific rate catalogs 
    new_lattice->assign_region_rates_wrapper(temp_regions, misc_rates);

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
            new_lattice->ratecatalog_100(i,j) = misc_rates[1];
        }
    }

    for (int i=0; i<(int)a_type_values.size(); i++) {new_lattice->a_types[a_type_keys[i]] = a_type_values[i];}

    new_lattice->region_energies = std::get<2>(cat_tuple); 

    std::vector< std::vector<int> > nonzero_vacs = new_lattice->vacancies.nonzero();

    for (int i=0; i<nonzero_vacs.size(); i++) {
        for (int j=0; j<nonzero_vacs[0].size(); j++) {
            new_lattice->vacancies_pos[i][j] = nonzero_vacs[i][j];
        }
    }

    new_lattice->void_threshold = misc_rates[2];
    new_lattice->void_barrier = misc_rates[3];
    new_lattice->terrace_barrier_111 = misc_rates[4];
    new_lattice->terrace_barrier_100 = misc_rates[5];
    new_lattice->void_gb_diss_barrier = misc_rates[6];


    /*
    new_lattice->void_threshold = 3;
    new_lattice->void_barrier = 6.84e8; //3.04e6; //6.84e8; //2.88e3; // 5.06e7; //2761667.608;
    
    if (new_lattice->void_barrier == 6.84e8) {
        new_lattice->terrace_barrier_111 = 1.43e9; //1.43e9; //7.13e6; //1.67e9;
        new_lattice->terrace_barrier_100 = 1.43e9; //1.43e9; //7.13e6; //1.67e9;
    }
    else {
        new_lattice->terrace_barrier_111 = 1.67e9; //1.43e9; //7.13e6; //1.67e9;
        new_lattice->terrace_barrier_100 = 7.13e6; //1.43e9; //7.13e6; //1.67e9;
    }
    
    new_lattice->void_gb_diss_barrier = 6.84e8; //1.43e9; //7.13e6; //1.67e9;
    */

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

