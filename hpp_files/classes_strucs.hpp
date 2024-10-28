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
    mat_type operator() (size_t row, size_t col) const {
        return data_[cols_ * row + col];
    }
    
    /*! \brief Set all matrix elements to zero */
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
        if ((row < rows_) && (row >= 0)) {
            for (size_t row_idx = row; row_idx < rows_ - 1; ++row_idx) {
                for (size_t col_idx = 0; col_idx < cols_; ++col_idx) {
                    (*this)(row_idx, col_idx) = (*this)(row_idx + 1, col_idx);
                }
            }
            reshape(rows_ - 1, cols_);
        } else {
            std::cout << "ERROR: Attempted to remove row out of bounds - row " << row << " on rank " << rank << "\n";
            exit(0);
        }
    }

    /*! 
     * \brief Add an empty row into the matrix at a specified position
     * \param [in] row Position where the row should be added
     * \param [in] rank Processor rank (for error handling)
     */
    void add_row(size_t row, int rank) { 
        reshape(rows_ + 1, cols_);
        if ((row < rows_) && (row >= 0)) {
            for (size_t row_idx = rows_ - 1; row_idx > row; --row_idx) {
                for (size_t col_idx = 0;

            }
        }
    }
};

/*! \brief A class for storing 4-D arrays of integers. */
class FourDArr {
public:
    /*! \brief Constructor.
     * \param [in] len1 The length of the first dimension.
     * \param [in] len2 The length of the second dimension.
     * \param [in] len3 The length of the third dimension.
     * \param [in] len4 The length of the fourth dimension.
     */
    FourDArr(size_t len1, size_t len2, size_t len3, size_t len4);

    /*! \brief Access an element of the 4-D array.
     * \param [in] i1 First index.
     * \param [in] i2 Second index.
     * \param [in] i3 Third index.
     * \param [in] i4 Fourth index.
     * \return Reference to array element.
     */
    int& operator()(size_t i1, size_t i2, size_t i3, size_t i4);

    /*! \brief Access an element of the 4-D array (const).
     * \param [in] i1 First index.
     * \param [in] i2 Second index.
     * \param [in] i3 Third index.
     * \param [in] i4 Fourth index.
     * \return Array element.
     */
    int operator()(size_t i1, size_t i2, size_t i3, size_t i4) const;

    /*! \brief Destructor to free allocated memory. */
    ~FourDArr();

    /*! \brief Returns a pointer to the 0-index element in the array.
     * \return Pointer to the first element.
     */
    int* data();

    /*! \brief Retrieve nonzero elements.
     * \return A matrix containing indices of non-zero elements.
     */
    Matrix<int> nonzero();

    /*! \brief Print FourDArr object. */
    void print_4Dvector();

    /*! \brief Retrieve elements of FourDArr the based on vector of input coordinates.
     * \param [in] coords Coordinates of elements to retrieve.
     * \return Vector of values at the specified coordinates.
     */
    std::vector<int> grab_idxs(const std::vector<std::vector<int>>& coords);

    /*! \brief Assign values to 4-D array based on input coordinates.
     * \param [in] coords Coordinates of elements to assign.
     * \param [in] values Values to assign at specified coordinates.
     */
    void assign_idxs(const std::vector<std::vector<int>>& coords, const std::vector<int>& values);

    /*! \brief Retrieve nonzero elements with their values from the 4-D array.
     * \return A matrix containing indices and values of non-zero elements.
     */
    Matrix<int> nonzero_elems();

    /*! \brief Set all elements in the data field to zero. */
    void zero();

    std::vector<size_t> size_vec; ///< Vector storing the sizes of each dimension.

private:
    size_t len1_; ///< Length of the first dimension.
    size_t len2_; ///< Length of the second dimension.
    size_t len3_; ///< Length of the third dimension.
    size_t len4_; ///< Length of the fourth dimension.
    int* data_; ///< Pointer to the data stored in the array.
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
        return BoolReference(data_[coarse_idx], fine_idx);
    }

    /*! \brief Print elements of the data field. */
    void print();

    /*! \brief Retrieve nonzero elements of the array.
     * \return A matrix of indices for nonzero elements.
     */
    Matrix<int> nonzero();

    /*! \brief Retrieve nonzero elements along with their values.
     * \return A matrix containing indices and values of nonzero elements.
     */
    Matrix<int> nonzero_elems();

    /*! \brief Set all elements in the data field to zero. */
    void zero();
};


/*! \brief Structure for managing regions and site data. */
struct add_reg_struct {
public:
    /*! \brief Constructor
     * \param [in] idx Index of the structure.
     * \param [in] regions Vector of pointers to Region objects.
     * \param [in] region_sites Pointer to a FourDArr containing region sites data.
     */
    add_reg_struct(int idx, std::vector<Region*> regions, FourDArr* region_sites)
        : idx_(idx) {
        regions_ = regions;
        region_sites_ = region_sites;
        std::cout << "regions.size(): " << regions.size() << "\n";
        std::cout << "regions_.size(): " << regions_.size() << "\n";
    }

    /*! \brief Get the idx field.
     * \return Integer value of index into vector of lines read in from input file.
     */
    int get_idx() const { return idx_; }

    /*! \brief Get the regions vector.
     * \return Vector of Region pointers.
     */
    std::vector<Region*> get_regions() const { return regions_; }

    /*! \brief Get the region_sites field.
     * \return Pointer to the FourDArr containing region site data.
     */
    FourDArr* get_region_sites() const { return region_sites_; }

private:
    int idx_; ///< Index into the vector of lines from input file to be read.
    std::vector<Region*> regions_; ///< Vector containing pointers to Region objects.
    FourDArr* region_sites_; ///< Pointer to the FourDArr with sites occupied by ids 
                             /// correpsonding to bounds of regions.
};

typedef struct add_reg_struct add_reg_struct;

/*! \brief Structure for storing lattice simulation results. */
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
