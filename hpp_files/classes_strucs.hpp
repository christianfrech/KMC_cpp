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


/*! \brief A class for storing data (rate constants) for directional bias corresponding to a region of lattice
 *
 *  The Region class holds various parameters specific to a defined region, including the region ID, shape and bounds of region,
 *  rate constants, and direction of bias for rate constants. Depending on the type of region (e.g., "GB" or "BLOCK"),
 *  different attributes (slopes and shifts or bounds) are initialized.
 */
class Region {
    public:
        /*! \brief Unique number corresponding to the region */
        int id;
        
        /*! \brief Type of the region (BLOCK or GB), defining the parameter structure */
        std::string type;
        
        /*! \brief Slope values for the region, used if type is "GB" */
        std::vector<int> slopes;
        
        /*! \brief Shift away from origin values for the region, used if type is "GB" */
        std::vector<int> shifts;
        
        /*! \brief Lower bounds for the region, used if type is "BLOCK" */
        std::vector<int> lowerbound;
        
        /*! \brief Upper bounds for the region, used if type is "BLOCK" */
        std::vector<int> upperbound;
        
        /*! \brief General parameters for the region, format varies by region type */
        std::vector< std::vector<int> > params;
        
        /*! \brief energy barriers associated with region, structured as a multi-dimensional vector */
        std::vector< std::vector< std::vector<double> > > energies;

        /*! 
         *  \brief Constructs a Region object
         *  
         *  Initializes a Region with a given ID, type, and parameter set. Depending on the specified type,
         *  initializes slopes and shift away from origin  dfor type "GB" or lower and upper bounds for type "BLOCK".
         * 
         *  \param id_in Integer identifier for the region
         *  \param reg_type String specifying the region type (e.g., "GB" or "BLOCK")
         *  \param params_in Parameter set used to initialize region attributes
         */
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
            /*! 
             * \brief Constructs a RowReference object 
             * \param mat Reference to the matrix object
             * \param row Index of the row in the matrix
             */
            RowReference(Matrix<mat_type> &mat, size_t row) : row_idx_(row), mat_(mat) {}

            /*! 
             * \brief Access an element in the referenced row
             * \param idx Column index of the element
             * \return Reference to the element at [row, idx]
             */
            mat_type& operator[] (size_t idx) {
                return mat_(row_idx_, idx);
            }
    };
    
    /*! 
     * \brief Constructor for Matrix class
     * \param [in] rows The initial number of rows in the matrix
     * \param [in] cols The initial number of columns in the matrix
     */
    Matrix(size_t rows, size_t cols) : rows_(rows), cols_(cols), tot_size_(rows * cols), data_(rows * cols) {}
    
    /*! 
     * \brief Access matrix element (non-const version)
     * \param [in] row Row index of element
     * \param [in] col Column index of element
     * \return Reference to matrix element
     */
    mat_type& operator() (size_t row, size_t col) {
        return data_[cols_ * row + col];
    }
    
    /*! 
     * \brief Access matrix element (const version)
     * \param [in] row Row index of element
     * \param [in] col Column index of element
     * \return Value of matrix element
     */
    mat_type  operator() (size_t row, size_t col) const {
        return data_[cols_ * row + col];
    }
    
    /*! \brief Zero all matrix elements */
    void zero() {
        std::fill(data_.begin(), data_.end(), 0);
    }
    
    /*! 
     * \brief Access matrix row (non-const version)
     * \param [in] row Row index
     * \return Pointer to the 0-index element in the specified row
     */
    mat_type *operator[] (size_t row) {
        return &data_[cols_ * row];
    }
    
    /*! 
     * \brief Access matrix row (const version)
     * \param [in] row Row index
     * \return Pointer to the 0-index element in the specified row
     */
    const mat_type *operator[] (size_t row) const {
        return &data_[cols_ * row];
    }

    /*! 
     * \brief Remove a row from the matrix
     * \param [in] row Row index to be removed
     * \param [in] rank Processor rank (for error handling)
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

    /*! 
     * \brief Add an empty row into the matrix at a specified position
     * \param [in] row Position where the row should be added
     * \param [in] rank Processor rank (for error handling)
     */
    void add_row(size_t row, int rank) {
        //std::cout << "rows_: " << rows_ << " row: " << row << "\n";
        reshape((rows_ + 1), cols_);

        if ((row < rows_) && (row >= 0)) {
            for (size_t row_idx = 0; row_idx < row; row_idx++) {
                for (size_t col_idx = 0; col_idx < cols_; col_idx++) {
                    (*this)(row_idx, col_idx) = (*this)(row_idx, col_idx);
                }
            }

            //std::cout << "rank: " << rank << " vacancies_pos.print(): " << "\n";
            //(*this).print();

            for (size_t row_idx = (size_t)(rows_-2); row_idx > (size_t)(row-1); row_idx--) {
                for (size_t col_idx = 0; col_idx < cols_; col_idx++) {
                    //std::cout << "rank: " << rank << " row_idx: " << row_idx << "\n";
                    //std::cout << "rank: " << rank << " (*this).rows(): " << (*this).rows() << "\n";
                    (*this)((size_t)(row_idx+1), col_idx) = (*this)((size_t)(row_idx), col_idx);
                }
            }
        }
        else {
            std::cout << "ERROR: Attempted to add row out of bounds for the matrix - row " << row << " on rank " << rank << "\n";
            exit(0);
        }
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

/*! \brief A class for storing 4-D arrays of integers. */
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
        std::cout << "nonzero \n";
        Matrix<int> mat_out(vec_size, 4);

        int elem = 0;
        std::cout << "len1_: " << len1_ << " len2_: " << len2_ << " len3_: " << len3_ << " len4_: " << len4_ << "\n";
        
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


/*! \brief A class for storing 4-D arrays of boolean values. */
class FourDBoolArr {
    size_t len1_, len2_, len3_, len4_; ///< Dimensions of the array.
    std::vector<uint8_t> data_; ///< Vector of boolean values packed into uint8_t 
                                /// (8 boolean values per uint8_t vectoir entry).

    /*! \brief A reference wrapper for a single bit in the 4-D array.
     *
     * This class provides access to a specific bit within a byte and allows modification and retrieval
     * of that bit as a boolean value.
     */
    class BoolReference {
    private:
        uint8_t& value_; ///< Reference to the byte containing the target bit.
        uint8_t mask_; ///< Mask to isolate the target bit.

        /*! \brief Sets the target bit to zero. */
        void zero() noexcept { value_ &= ~(mask_); }

        /*! \brief Sets the target bit to one. */
        void one() noexcept { value_ |= mask_; }

        /*! \brief Gets the current value of the target bit.
         * \return Boolean value of the bit.
         */
        bool get() const noexcept { return !!(value_ & mask_); }

        /*! \brief Sets the bit to the specified boolean value.
         * \param b Boolean value to set.
         */
        void set(bool b) noexcept { b ? one() : zero(); }

    public:
        /*! \brief Constructor.
         * \param value Reference to the byte containing the bit.
         * \param nbit Position of the bit within the byte.
         */
        BoolReference(uint8_t& value, uint8_t nbit)
            : value_(value), mask_(uint8_t(0x1) << nbit) {}

        /*! \brief Copy constructor.
         * \param ref BoolReference to copy.
         */
        BoolReference(const BoolReference& ref) : value_(ref.value_), mask_(ref.mask_) {}

        /*! \brief Assignment operator for boolean values.
         * \param b Boolean value to assign.
         * \return Reference to this BoolReference.
         */
        BoolReference& operator=(bool b) noexcept { set(b); return *this; }

        /*! \brief Assignment operator for BoolReference values.
         * \param br BoolReference value to assign.
         * \return Reference to this BoolReference.
         */
        BoolReference& operator=(const BoolReference& br) noexcept { return *this = bool(br); }

        /*! \brief Conversion operator to boolean.
         * \return Boolean value of the bit.
         */
        operator bool() const noexcept { return get(); }
    };
    
public:
    std::vector<size_t> size_vec; ///< Vector storing the sizes of each dimension.

    /*! \brief Constructor.
     * \param len1 Length of the first dimension.
     * \param len2 Length of the second dimension.
     * \param len3 Length of the third dimension.
     * \param len4 Length of the fourth dimension.
     */
    FourDBoolArr(size_t len1, size_t len2, size_t len3, size_t len4)
        : len1_(len1), len2_(len2), len3_(len3), len4_(len4), data_(CEILING(len1 * len2 * len3 * len4, 8)) {
        size_vec.push_back(len1);
        size_vec.push_back(len2);
        size_vec.push_back(len3);
        size_vec.push_back(len4);
    }

    /*! \brief Access an element of the 4-D array.
     * \param i1 First index.
     * \param i2 Second index.
     * \param i3 Third index.
     * \param i4 Fourth index.
     * \return Reference to BoolReference for the element.
     */
    BoolReference operator()(size_t i1, size_t i2, size_t i3, size_t i4) {
        size_t flat_idx = i1 * len2_ * len3_ * len4_ + i2 * len3_ * len4_ + i3 * len4_ + i4;
        size_t coarse_idx = flat_idx / 8;
        size_t fine_idx = flat_idx % 8;
        BoolReference ref(data_[coarse_idx], fine_idx);
        return ref;
    }
    
    /*
    print elements of data field
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
    
    /*
    retrieve nonzero elements
    */
    Matrix<int> nonzero() {
        //std::cout << "nonzero enter \n";
        int vec_size = (int)(len1_*len2_*len3_*len4_)/8;
        //std::cout << "vec_size: " << vec_size << " \n";
        Matrix<int> mat_out(vec_size, 4);
        //std::cout << "nonzero matrix initialized \n";

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
        
        //std::cout << "reshaping: \n";
        mat_out.reshape(elem, 4);
        //std::cout << "nonzero elem: " << elem << "\n";
        return mat_out;
    }

    /*
    retrieve nonzero elements
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

    /*
    set all elements in data field to zero
    */
    void zero() {
        std::fill(data_.begin(), data_.end(), 0);
    }

};


/*! \brief Structure for managing regions and site data. */
struct add_reg_struct {
    public:
        /*! \brief Constructor
        * \param [in] len1 len2 len3 len4    Lengths of the 4 dimensions of the array
        */

        int idx_;

        add_reg_struct(int idx , std::vector<Region*> regions, FourDArr* region_sites):
            idx_(idx)
            {
                regions_ = regions;
                region_sites_ = region_sites;
                std::cout << "regions.size(): " << regions.size() << "\n";
                std::cout << "regions_.size(): " << regions_.size() << "\n";
            }

        /*! \return get idx field of struc */
        int get_idx() const {
            return idx_;
        }

        /*! \return get regions field of struc */
        std::vector<Region*> get_regions() const {
            return regions_;
        }

        /*! \return get regions_sites field of struc */
        FourDArr* get_region_sites() const {
            return region_sites_;
        }
        
    private:
        std::vector<Region*> regions_; FourDArr* region_sites_;
};

typedef struct add_reg_struct add_reg_struct;

struct lattice_return_struct {
public:
    /*! \brief Constructor
     * \param [in] move_counts Vector of integer counts of moves.
     * \param [in] time_count Vector of doubles representing time intervals.
     * \param [in] all_times Vector of doubles representing cumulative times.
     */
    lattice_return_struct(std::vector<int> move_counts, std::vector<double> time_count, std::vector<double> all_times)
        : move_counts_(move_counts), time_count_(time_count), all_times_(all_times) {}

    /*! \brief Get the counts of each type of move executed throughout simulation.
     * \return Vector of integers representing move counts by each move type.
     */
    std::vector<int> get_move_counts() const { return move_counts_; }

    /*! \brief Get the time elapsed during each type of move throughout simulation.
     * \return Vector of doubles representing times elapsed by move type.
     */
    std::vector<double> get_time_count() const { return time_count_; }

    /*! \brief Get each timestep of each type of move executed throughout simulation.
     * \return Vector of doubles representing cumulative times.
     */
    std::vector<double> get_all_times() const { return all_times_; }

private:
    std::vector<int> move_counts_; ///< Vector containing counts of moves.
    std::vector<double> time_count_; ///< Vector of time intervals.
    std::vector<double> all_times_; ///< Vector of cumulative times.
};

typedef struct lattice_return_struct lattice_return_struct;
