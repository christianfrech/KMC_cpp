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
#include <fstream>
#include <map>
#include <cmath>
#include <math.h>
#include <chrono>
#include <execution>
#include <cstdint>
#include <filesystem>

#define CEILING(x,y) ((x + y - 1) / y)


/*! \brief A class for representing and manipulating matrices with variable dimension sizes
 * \tparam mat_type The type of elements to be stored in the matrix
 */
template <class mat_type>
class Matrix {
public:
    /*! \brief A nested class for referencing rows within the Matrix */
    class RowReference {
        private:
            /*! \brief Index of the row */
            size_t row_idx_;
            /*! \brief Reference to the parent matrix */
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

    /*! \brief Remove row from matrix
     * Data are copied such that the first n elements in each row remain the same before and after this operation
     * \param [in] new_col      Desired number of columns in the enlarged matrix
     * \param [in] n_keep       Number of elements to preserve in all rows of the matrix
     */
    void remove_row(size_t row, int rank) { 
        //std::cout << "rows_: " << rows_ << " row: " << row << "\n";

        if ((row < rows_) && (row >= 0)) {
            for (size_t row_idx = 0; row_idx < row; row_idx++) {
                for (size_t col_idx = 0; col_idx < cols_; col_idx++) {
                    (*this)(row_idx, col_idx) = (*this)(row_idx, col_idx);
                }
            }
            for (size_t row_idx = (size_t)(row+1); row_idx < rows_; row_idx++) {
                for (size_t col_idx = 0; col_idx < cols_; col_idx++) {
                    (*this)((size_t)(row_idx-1), col_idx) = (*this)(row_idx, col_idx);
                }
            }
            reshape((rows_ - 1), cols_);
        }
        else {
            std::cout << "ERROR: Attempted to remove row out of bounds for the matrix - row " << row << " on rank " << rank << "\n";
            exit(0);
        }
    }

    /*! \brief Add row to matrix
     * Data are copied such that the first n elements in each row remain the same before and after this operation
     * \param [in] new_col      Desired number of columns in the enlarged matrix
     * \param [in] n_keep       Number of elements to preserve in all rows of the matrix
     */
    void add_row(size_t row, int rank) { 
        /*std::cout << "rank: " << rank << " rows_: " << rows_ << " row: " << row << "\n";
        std::cout << "rank: " << rank << " row[0]" << (*this)(row, 0) << " row[1]: " << (*this)(row, 1) << " row[2]: " << (*this)(row, 2) << " row[3]: " << (*this)(row, 3) << "\n";
        
        std::cout << "rank: " << rank << " pre rows_: " << rows_ << " row+1: " << (row+1) << "\n";
        std::cout << "rank: " << rank << " pre row+1[0]: " << (*this)((row+1), 0) << " row+1[1]: " << (*this)((row+1), 1) << " row+1[2]: " << (*this)((row+1), 2) << " row+1[3]: " << (*this)((row+1), 3) << "\n";
        */
        reshape((rows_ + 1), cols_);

        if ((row < rows_) && (row >= 0)) {
            for (int row_idx = 0; row_idx < row; row_idx++) {
                for (int col_idx = 0; col_idx < cols_; col_idx++) {
                    (*this)(row_idx, col_idx) = (*this)(row_idx, col_idx);
                }
            }

            for (int row_idx = (int)(rows_-2); row_idx > (int)(row-1); row_idx--) {
                for (int col_idx = 0; col_idx < cols_; col_idx++) {
                    
                    /*if ((row_idx < (row + 10)) && (col_idx == 0)) {
                        std::cout << "rank: " << rank << " pre shift row_idx: " << row_idx << "\n";    
                        std::cout << "rank: " << rank << " pre shift row_idx[0]: " << (*this)(row_idx, 0) << " row_idx[1]: " << (*this)(row_idx, 1) << " row_idx[2]: " << (*this)(row_idx, 2) << " row_idx[3]: " << (*this)(row_idx, 3) << "\n";
                        std::cout << "rank: " << rank << " pre shift row_idx+1[0]: " << (*this)((row_idx+1), 0) << " row_idx+1[1]: " << (*this)((row_idx+1), 1) << " row_idx+1[2]: " << (*this)((row_idx+1), 2) << " row_idx+1[3]: " << (*this)((row_idx+1), 3) << "\n";
                    }*/
                    
                    (*this)((size_t)(row_idx+1), col_idx) = (*this)((size_t)(row_idx), col_idx); // THIS ISN'T WORKING, NOT SHIFTING

                    /*if ((row_idx < (row + 10)) && (col_idx == 0)) {
                        std::cout << "rank: " << rank << " post shift row_idx: " << row_idx << "\n";    
                        std::cout << "rank: " << rank << " post shift row_idx[0]: " << (*this)(row_idx, 0) << " row_idx[1]: " << (*this)(row_idx, 1) << " row_idx[2]: " << (*this)(row_idx, 2) << " row_idx[3]: " << (*this)(row_idx, 3) << "\n";
                        std::cout << "rank: " << rank << " post shift row_idx+1[0]: " << (*this)((row_idx+1), 0) << " row_idx+1[1]: " << (*this)((row_idx+1), 1) << " row_idx+1[2]: " << (*this)((row_idx+1), 2) << " row_idx+1[3]: " << (*this)((row_idx+1), 3) << "\n";
                    }*/   
                }
            }
        }
        else {
            std::cout << "ERROR: Attempted to add row out of bounds for the matrix - row " << row << " on rank " << rank << "\n";
        }
        /*std::cout << "rank: " << rank << " post rows_: " << rows_ << " row: " << row << "\n";
        std::cout << "rank: " << rank << " post row[0]: " << (*this)(row, 0) << " row[1]: " << (*this)(row, 1) << " row[2]: " << (*this)(row, 2) << " row[3]: " << (*this)(row, 3) << "\n";
        std::cout << "rank: " << rank << " post rows_: " << rows_ << " row+1: " << (row+1) << "\n";
        std::cout << "rank: " << rank << " post row+1[0]" << (*this)((row+1), 0) << " row+1[1]: " << (*this)((row+1), 1) << " row+1[2]: " << (*this)((row+1), 2) << " row+1[3]: " << (*this)((row+1), 3) << "\n";
        */
    }

    /*! \brief Remove row from col
     * Data are copied such that the first n elements in each row remain the same before and after this operation
     * \param [in] new_col      Desired number of columns in the enlarged matrix
     * \param [in] n_keep       Number of elements to preserve in all rows of the matrix
     */
    void remove_col(size_t col, int rank) {
        std::cout << "cols_: " << cols_ << " col: " << col << "\n";

        if ((col < cols_) && (col >= 0)) {
            for (size_t row_idx = 0; row_idx < rows_; row_idx++) {
                for (size_t col_idx = 0; col_idx < col; col_idx++) {
                    (*this)(row_idx, col_idx) = (*this)(row_idx, col_idx);
                }
            }
            for (size_t row_idx = 0; row_idx < rows_; row_idx++) {
                for (size_t col_idx = (size_t)(col+1); col_idx < cols_; col_idx++) {
                    (*this)(row_idx, (size_t)(col_idx-1)) = (*this)(row_idx, col_idx);
                }
            }
            reshape(rows_, (cols_ - 1));
        }
        else {
            std::cout << "ERROR: Attempted to remove col out of bounds for the matrix - col " << col << " on rank " << rank << "\n";
            exit(0);
        }
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
    void reshape(size_t new_rows, size_t new_cols, int rank=-1) {
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
    
    /*! \return Pointer to the data in the matrix*/
    mat_type *data() const {
        return (mat_type *) data_.data();
    }
    
    /*!
    * \brief Copies data from another matrix.
    * \param mat Matrix to copy from.
    */
    void copy_from(Matrix<mat_type> &mat) {
        std::copy(mat.data_.begin(), mat.data_.end(), data_.begin());
    }

    /*!
    * \brief Prints the matrix to standard output.
    */
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
    size_t rows_; /*!< Number of rows in the matrix */
    size_t cols_; /*!< Number of columns in the matrix */
    size_t tot_size_; /*!< Total size of the matrix */
    std::vector<mat_type> data_; /*!< Internal storage for matrix data */
};

/*! \brief A class for storing 4-D arrays of ints*/
class FourDArr {
public:
    /*! \brief Constructor
     * \param [in] len1 len2 len3 len4    Lengths of the 4 dimensions of the array
     */

    std::vector<size_t> size_vec;

    FourDArr(size_t len1, size_t len2, size_t len3, size_t len4)
    : len1_(len1)
    , len2_(len2)
    , len3_(len3)
    , len4_(len4)
    {
        data_ = (int *)malloc(sizeof(int) * len1 * len2 * len3 * len4);
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
    Matrix<int> nonzero() {
        int vec_size = (int)(len1_*len2_*len3_*len4_)/8;
        //std::cout << "nonzero \n"; 
        Matrix<int> mat_out(vec_size, 4);

        int elem = 0; 
        //std::cout << "len1_: " << len1_ << " len2_: " << len2_ << " len3_: " << len3_ << " len4_: " << len4_ << "\n"; 
        
        for (int i1=0; i1<(int)len1_; i1++) {
            for (int i2=0; i2<(int)len2_; i2++) {
                for (int i3=0; i3<(int)len3_; i3++) {
                    for (int i4=0; i4<(int)len4_; i4++) {
                        if ( (*this)((size_t)i1, (size_t)i2, (size_t)i3, (size_t)i4) != false ) {
                            
                            if (elem >= (vec_size-1)) {
                                mat_out.reshape(vec_size*2, 4);
                                vec_size = vec_size * 2;
                            }
                            
                            mat_out[elem][0] = i1;
                            mat_out[elem][1] = i2;
                            mat_out[elem][2] = i3;
                            mat_out[elem][3] = i4;

                            elem ++;
                        }
                    }
                }
            }
        } 

        mat_out.reshape(elem, 4);
        return mat_out;
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
        Matrix<int> mat_out(vec_size, 5);
        mat_out.zero();

        int elem = 0; 
        //std::cout << "nonzero_elems() \n";
        //std::cout << "len1_: " << len1_ << " len2_: " << len2_ << " len3_: " << len3_ << " len4_: " << len4_ << "\n";
        
        for (int i1=0; i1<(int)len1_; i1++) {
            for (int i2=0; i2<(int)len2_; i2++) {
                for (int i3=0; i3<(int)len3_; i3++) {
                    for (int i4=0; i4<(int)len4_; i4++) {
                        if ( (*this)((size_t)i1, (size_t)i2, (size_t)i3, (size_t)i4) != 0 ) {
                            //std::cout << "i1: " << i1 << " i2: " << i2 << " i3: " << i3 << " i4: " << i4 << " elem: " << (*this)((size_t)i1, (size_t)i2, (size_t)i3, (size_t)i4) << "\n";
                            
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
        //std::cout << "nonzero elem: " << elem << "\n";
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

    void check_elems() {

        //std::cout << "check_elems()\n";
        for (int i1=0; i1<len1_; i1++) {
            for (int i2=0; i2<len2_; i2++) {
                for (int i3=0; i3<len3_; i3++) {
                    for (int i4=0; i4<len4_; i4++) {
                        if ((*this)((size_t)i1, (size_t)i2, (size_t)i3, (size_t)i4) != 0) {
                            //td::cout << "i1: " << i1 << " i2: " << i2 << " i3: " << i3 << " i4: " << i4 << " elem: " << (*this)((size_t)i1, (size_t)i2, (size_t)i3, (size_t)i4) << "\n";
                        }
                    }
                }
            }    
        }
    }

private:
    size_t len1_, len2_, len3_, len4_; ///< Dimensions of the array
    int* data_; ///< The data stored in the array
};

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

    Matrix<double> nonzero_elems() {
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

/*!
 * \class FourDBoolArr
 * \brief Represents a 4D boolean array implemented with bit-packing for memory efficiency.
 */
class FourDBoolArr {
    size_t len1_, len2_, len3_, len4_;
    std::vector<uint8_t> data_;

    /*!
     * \class BoolReference
     * \brief A helper class to manage individual bits in the 4D array.
     */
    class BoolReference {
    private:
        uint8_t &value_;
        uint8_t mask_;

        /*! \brief Sets the referenced bit to 0. */
        void zero() noexcept { value_ &= ~mask_; }

        /*! \brief Sets the referenced bit to 1. */
        void one() noexcept { value_ |= mask_; }

        /*! \brief Retrieves the current value of the referenced bit. */
        bool get() const noexcept { return !!(value_ & mask_); }

        /*! 
         * \brief Sets the referenced bit to the specified value.
         * \param b The boolean value to set.
         */
        void set(bool b) noexcept {
            if (b)
                one();
            else
                zero();
        }

    public:
        /*!
         * \brief Constructor for BoolReference.
         * \param value The byte containing the bit.
         * \param nbit The position of the bit within the byte.
         */
        BoolReference(uint8_t &value, uint8_t nbit)
            : value_(value), mask_(uint8_t(0x1) << nbit) {}

        /*! \brief Copy constructor for BoolReference. */
        BoolReference(const BoolReference &ref) : value_(ref.value_), mask_(ref.mask_) {}

        /*! \brief Assigns a boolean value to the referenced bit. */
        BoolReference &operator=(bool b) noexcept {
            set(b);
            return *this;
        }

        /*! \brief Copy assignment operator. */
        BoolReference &operator=(const BoolReference &br) noexcept { return *this = bool(br); }

        /*! \brief Implicit conversion to boolean. */
        operator bool() const noexcept { return get(); }
    };

public:
    std::vector<size_t> size_vec;

    /*!
     * \brief Constructor for FourDBoolArr.
     * \param len1 Length of the first dimension.
     * \param len2 Length of the second dimension.
     * \param len3 Length of the third dimension.
     * \param len4 Length of the fourth dimension.
     */
    FourDBoolArr(size_t len1, size_t len2, size_t len3, size_t len4)
        : len1_(len1), len2_(len2), len3_(len3), len4_(len4),
          data_(CEILING(len1 * len2 * len3 * len4, 8)) {
        size_vec.push_back(len1);
        size_vec.push_back(len2);
        size_vec.push_back(len3);
        size_vec.push_back(len4);
    }

    /*!
     * \brief Accesses an element of the 4D array.
     * \param i1 Index along the first dimension.
     * \param i2 Index along the second dimension.
     * \param i3 Index along the third dimension.
     * \param i4 Index along the fourth dimension.
     * \return A BoolReference to the specified bit.
     */
    BoolReference operator()(size_t i1, size_t i2, size_t i3, size_t i4) {
        size_t flat_idx = i1 * len2_ * len3_ * len4_ + i2 * len3_ * len4_ + i3 * len4_ + i4;
        size_t coarse_idx = flat_idx / 8;
        size_t fine_idx = flat_idx % 8;
        return BoolReference(data_[coarse_idx], fine_idx);
    }

    /*!
     * \brief Prints the elements of the 4D array to the console.
     */
    void print() { 
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
    
    /*!
     * \brief Retrieves the indices of non-zero elements.
     * \return A Matrix<int> containing the indices of non-zero elements.
     */
    Matrix<int> nonzero(int rank) {
        
        //std::cout << "rank: "<< rank  << " nonzero enter \n";
        int vec_size = (int)(len1_*len2_*len3_*len4_)/8;
        //std::cout << "rank: "<< rank  << " vec_size: " << vec_size << " \n";
        Matrix<int> mat_out(vec_size, 4);
        //std::cout << "rank: "<< rank  << " nonzero matrix initialized \n";

        int elem = 0; 
        
        for (int i1=0; i1<(int)len1_; i1++) {
            for (int i2=0; i2<(int)len2_; i2++) {
                for (int i3=0; i3<(int)len3_; i3++) {
                    for (int i4=0; i4<(int)len4_; i4++) {
                        if ( (*this)((size_t)i1, (size_t)i2, (size_t)i3, (size_t)i4) != false ) {
                            
                            if (elem >= (vec_size-2)) {

                                //std::cout << "reshaping in loop: \n";
                                mat_out.reshape(vec_size*2, 4);
                                //std::cout << "reshaped \n";
                                vec_size = vec_size * 2;
                            }
                            
                            mat_out[elem][0] = i1;
                            mat_out[elem][1] = i2;
                            mat_out[elem][2] = i3;
                            mat_out[elem][3] = i4;
                            elem ++;
                        }
                    }
                }
            }
        } 
        
        //std::cout << "rank: "<< rank  << " reshaping: \n";
        mat_out.reshape(elem, 4);
        //std::cout << "rank: "<< rank  << " nonzero elem: " << elem << "\n";
        
        return mat_out;
    }

    /*!
     * \brief Retrieves the indices and values of non-zero elements.
     * \return A Matrix<int> containing indices and values of non-zero elements.
     */
    Matrix<int> nonzero_elems() {
        int vec_size = (int)(len1_*len2_*len3_*len4_)/8;
        Matrix<int> mat_out(vec_size, 5);

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

    /*!
     * \brief Sets all elements of the array to zero.
     */
    void zero() {
        std::fill(data_.begin(), data_.end(), 0);
    }


    void check_equal(FourDBoolArr array_to_compare, int rank) {
        std::vector<size_t> orig_dims = size_vec;
        std::vector<size_t> comparison_dims = array_to_compare.size_vec;
        assert(std::equal(orig_dims.begin(), orig_dims.begin() + orig_dims.size(), comparison_dims.begin()));

        for (size_t i=0; i<orig_dims[0]; i++) {
            for (size_t j=0; j<orig_dims[1]; j++) {
                for (size_t k=0; k<orig_dims[2]; k++) {
                    for (size_t l=0; l<orig_dims[3]; l++) {
                        //std::cout << "rank: " << rank << " i: " << i << " j: " << j << " k: " << k << " l: " << l << "\n";
                        if ((*this)(i,j,k,l) 
                            != array_to_compare(i,j,k,l)) {
                            std::cout << "rank: " << rank << " orig_array_val: " << (*this)(i,j,k,l) 
                            << " comparison_arr_val: " << array_to_compare(i,j,k,l) << "\n";
                            std::cout << "rank: " << rank << " UNEQUAL i: " << i << " j: " << j << " k: " << k << " l: " << l << "\n";
                            exit(0);
                        }
                    } 
                }  
            } 
        } 
    }

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
    double rate;  ///< Maximum allowed energy barrier.
    std::vector<int> dimensions;  ///< Dimensions of the region.
    bool random;  ///< Indicates whether rates are randomly assigned.
    std::normal_distribution<double> distribution;  ///< Normal distribution for random barrier generation.
    std::vector<double> rates;  ///< Predefined rates for the region.
    bool interface;  ///< Whether the region has an interface.
    int interface_i;  ///< Index of the interface.
    int interface_dim;  ///< Dimension of the interface.
    double interface_100_e;  ///< Rate for the 100 interface.
    bool is_gb;
    double e_below_bulk;
    double interface_100_rate;
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
    /* energetic landscape based constructor */
    Region(int id_in, std::string reg_type, std::string bias_direc, std::vector<std::vector<int>> params_in,
           std::vector<double> distribution, bool random_val, double e_below_bulk_in,
           std::vector<int> interface_params, double interface_100_e_in):
            id(id_in), 
            type(reg_type),
            bias(bias_direc),
            params(params_in),
            region_rates_L(2, (params[1][0] - params[0][0]), (params[1][1] - params[0][1]), (params[1][2] - params[0][2])),
            region_rates_R(2, (params[1][0] - params[0][0]), (params[1][1] - params[0][1]), (params[1][2] - params[0][2])),
            random(random_val),
            dimensions(3),
            e_below_bulk(e_below_bulk_in),
            interface(interface_params[0])

            {
                lowerbound = params[0];
                upperbound = params[1];
                dimensions[0] = params[1][0] - params[0][0];
                dimensions[1] = params[1][1] - params[0][1];
                dimensions[2] = params[1][2] - params[0][2];
                std::cout << "interface_params: " <<  interface_params[0] << " "<<  interface_params[1]  << " "<<  interface_params[2] << " " <<  interface_params[3] << "\n";
                if ((bool)interface_params[0]) {
                    interface_i = (int)interface_params[1];
                    interface_dim = (int)interface_params[2];
                    is_gb = (bool)interface_params[3];
                    interface_100_e = interface_100_e_in;
                    std::cout << "INTERFACE id_in: " << id_in << " is_gb: " <<  is_gb << "\n";
                }
                else if ((bool)interface_params[3]) {
                    interface_i = (int)interface_params[1];
                    interface_dim = (int)interface_params[2];
                    is_gb = (bool)interface_params[3];
                    interface_100_e = interface_100_e_in;
                    std::cout << "GB id_in: " << id_in << " is_gb: " <<  is_gb << "\n";
                }
                else {
                    interface_i = 0;
                    interface_dim = 0;
                    is_gb = 0;
                    interface_100_e = 0;
                    std::cout << "ELSE id_in: " << id_in << " is_gb: " <<  is_gb << "\n";

                }
            }

    /* rate based constructor */
    Region(int id_in, std::string reg_type, std::string bias_direc, std::vector<std::vector<int>> params_in,
           std::vector<double> distribution, bool random_val, std::vector<double> rates_in, 
           double e_below_bulk_in, std::vector<int> interface_params, double interface_100_e_in):
            id(id_in), 
            type(reg_type),
            bias(bias_direc),
            params(params_in),
            region_rates_L(2, (params[1][0] - params[0][0]), (params[1][1] - params[0][1]), (params[1][2] - params[0][2])),
            region_rates_R(2, (params[1][0] - params[0][0]), (params[1][1] - params[0][1]), (params[1][2] - params[0][2])),
            random(random_val),
            dimensions(3),
            rates(rates_in),
            e_below_bulk(e_below_bulk_in),
            interface(interface_params[0])

            {
                lowerbound = params[0];
                upperbound = params[1];
                dimensions[0] = params[1][0] - params[0][0];
                dimensions[1] = params[1][1] - params[0][1];
                dimensions[2] = params[1][2] - params[0][2];
                std::cout << "interface_params: " <<  interface_params[0] << " "<<  interface_params[1]  << " "<<  interface_params[2] << " " <<  interface_params[3] << "\n";
                if ((bool)interface_params[0]) {
                    interface_i = (int)interface_params[1];
                    interface_dim = (int)interface_params[2];
                    is_gb = (bool)interface_params[3];
                    interface_100_e = interface_100_e_in;
                    std::cout << "interface_100_e: " << interface_100_e << "\n";
                    std::cout << "INTERFACE id_in: " << id_in << " is_gb: " <<  is_gb << "\n";
                }
                else if ((bool)interface_params[3]) {
                    interface_i = (int)interface_params[1];
                    interface_dim = (int)interface_params[2];
                    is_gb = (bool)interface_params[3];
                    interface_100_e = interface_100_e_in;
                    std::cout << "interface_100_e: " << interface_100_e << "\n";
                    std::cout << "GB id_in: " << id_in << " is_gb: " <<  is_gb << "\n";
                }
                else {
                    interface_i = 0;
                    interface_dim = 0;
                    is_gb = 0;
                    interface_100_e = 0;
                    std::cout << "interface_100_e: " << interface_100_e << "\n";
                    std::cout << "ELSE id_in: " << id_in << " is_gb: " <<  is_gb << "\n";

                }
            }

   
   
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

/*!
 * \struct add_reg_struct
 * \brief Represents a structure to hold region data and related operations.
 */
struct add_reg_struct {
public:
    /*!
     * \brief Constructor for add_reg_struct.
     * \param idx Index of the structure.
     * \param regions A vector of Region pointers.
     * \param region_sites Pointer to a FourDArr object.
     */
    add_reg_struct(int idx , std::vector<Region*> regions, FourDArr* region_sites):
            idx_(idx)
            {
                regions_ = regions;
                region_sites_ = region_sites;
                std::cout << "regions.size(): " << regions.size() << "\n";
                std::cout << "regions_.size(): " << regions_.size() << "\n";
            }
            
    /*! \brief Retrieves the index field of the structure. */
    int get_idx() const { return idx_; }

    /*! \brief Retrieves the regions vector. */
    std::vector<Region *> get_regions() const { return regions_; }

    /*! \brief Retrieves the region_sites pointer. */
    FourDArr *get_region_sites() const { return region_sites_; }

private:
    int idx_;
    std::vector<Region *> regions_;
    FourDArr *region_sites_;
};
 
typedef struct add_reg_struct add_reg_struct;

/*!
 * \struct lattice_return_struct
 * \brief Represents a structure to return lattice operation results.
 */
struct lattice_return_struct {
public:
    /*!
     * \brief Constructor for lattice_return_struct.
     * \param move_counts A vector of movement counts.
     * \param time_count A vector of time counts.
     * \param all_times A vector of all times.
     */
    lattice_return_struct(std::vector<int> move_counts, std::vector<double> time_count, std::vector<double> all_times): 
        move_counts_(move_counts), time_count_(time_count), all_times_(all_times) {}

    /*! \brief Retrieves the move_counts vector. */
    std::vector<int> get_move_counts() const { return move_counts_; }

    /*! \brief Retrieves the time_count vector. */
    std::vector<double> get_time_count() const { return time_count_; }

    /*! \brief Retrieves the all_times vector. */
    std::vector<double> get_all_times() const { return all_times_; }

private:
    std::vector<int> move_counts_;
    std::vector<double> time_count_;
    std::vector<double> all_times_;
};

typedef struct lattice_return_struct lattice_return_struct;