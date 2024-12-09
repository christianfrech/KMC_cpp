#include <iostream>
#include <sstream>
#include <stdio.h>
#include <string>
#include <vector>
#include <cstdlib>
#include <fstream>
#include "print_func.hpp"

/*
template <typename A, typename B>
A vect_create_2D(size_t N, size_t M) {
    A vec_out(N);

    for(int i = 0; i < (int)N; i++) { 
        vec_out[i] = B(M);
    }

    return vec_out;
}

template <typename A, typename B, typename C>
A vect_create_3D_float(int L, int N, int M) {

    std::vector< std::vector< std::vector<double> > > vec_out(L, std::vector< std::vector<double> >(N, std::vector<double>(M, 0)));
    return vec_out;
}
*/

/*! 
 * \brief Creates a 2D vector of integers with specified dimensions.
 * 
 * This function creates a 2D vector of integers with the given number of rows (N) and columns (M).
 * Each element of the vector is initialized to 0.
 * 
 * \param N The number of rows in the 2D vector.
 * \param M The number of columns in the 2D vector.
 * 
 * \return A 2D vector of integers with dimensions N x M, initialized to 0.
 */
std::vector< std::vector<int> > vect_create_2D(size_t N, size_t M) {
    std::vector< std::vector<int> > vec_out(N);

    for(int i = 0; i < (int)N; i++) { 
        vec_out[i] = std::vector<int>(M);
    }

    return vec_out;
}

/*! 
 * \brief Creates a 2D vector of doubles with specified dimensions and initial value.
 * 
 * This function creates a 2D vector of doubles with the given number of rows (N) and columns (M),
 * initializing each element with a specified value (val).
 * 
 * \param N The number of rows in the 2D vector.
 * \param M The number of columns in the 2D vector.
 * \param val The value to initialize each element of the vector.
 * 
 * \return A 2D vector of doubles with dimensions N x M, initialized to val.
 */
std::vector< std::vector<double> > vect_create_2D_float(size_t N, size_t M, double val) {
    std::vector<std::vector<double>> vec_out(N, std::vector<double>(M, val));
    return vec_out;
}

/*! 
 * \brief Creates a 1D vector of integers with a specified size and initial value.
 * 
 * This function creates a 1D vector of integers with the given size, initializing each element to the specified value.
 * 
 * \param size The size of the 1D vector.
 * \param value The value to initialize each element of the vector.
 * 
 * \return A 1D vector of integers of the specified size, initialized to value.
 */
std::vector<int> create_vec_1D(int size, int value) {
    std::vector<int> vec(size, value);
    return vec;
}

/*! 
 * \brief Creates a 3D vector of doubles with specified dimensions and initialized to zero.
 * 
 * This function creates a 3D vector of doubles with the given dimensions L x N x M, initializing each element to 0.
 * 
 * \param L The number of 2D layers in the 3D vector.
 * \param N The number of rows in each 2D layer.
 * \param M The number of columns in each 2D layer.
 * 
 * \return A 3D vector of doubles with dimensions L x N x M, initialized to 0.
 */
std::vector< std::vector< std::vector<double> > > vect_create_3D_float(int L, int N, int M) {
    std::vector< std::vector< std::vector<double> > > vec_out(L, std::vector< std::vector<double> >(N, std::vector<double>(M, 0)));
    return vec_out;
}

/*! 
 * \brief Creates a 3D vector of integers with specified dimensions and initialized to zero.
 * 
 * This function creates a 3D vector of integers with the given dimensions L x N x M, initializing each element to 0.
 * 
 * \param L The number of 2D layers in the 3D vector.
 * \param N The number of rows in each 2D layer.
 * \param M The number of columns in each 2D layer.
 * 
 * \return A 3D vector of integers with dimensions L x N x M, initialized to 0.
 */
std::vector< std::vector< std::vector<int> > > vect_create_3D_int(int L, int N, int M) {
    std::vector< std::vector< std::vector<int> > > vec_out(L, std::vector< std::vector<int> >(N, std::vector<int>(M, 0)));
    return vec_out;
}

/*! 
 * \brief Creates a 1D vector of doubles with a specified size.
 * 
 * This function creates a 1D vector of doubles with the given size, leaving each element uninitialized.
 * 
 * \param size The size of the 1D vector.
 * 
 * \return A 1D vector of doubles with the specified size.
 */
std::vector<double> create_vec_1D_float(int size) {
    std::vector<double> vec(size);
    return vec;
}

/*! 
 * \brief Transposes a 2D integer vector (matrix).
 * 
 * This function performs a matrix transposition on a 2D vector of integers. The number of rows and columns
 * in the resulting transposed matrix is switched.
 * 
 * \param src The 2D vector (matrix) to transpose.
 * 
 * \return A new 2D vector (matrix) that is the transpose of the input vector.
 */
std::vector< std::vector<int> > int_transpose(std::vector< std::vector<int> > src) {
    int N = src.size();
    int M = src[0].size();
    std::vector< std::vector<int> > dst(M, std::vector<int>(N, 0));

    for(int i=0; i<(int)src.size(); i++) {
        for(int j=0; j<(int)src[0].size(); j++) {
            dst[j][i] = src[i][j];
        }
    }
    return dst;
}

/*! 
 * \brief Finds the indices of elements in a 1D vector that match a specified value.
 * 
 * This function iterates over a 1D vector and finds all indices where the elements match the given value.
 * The function returns a vector of indices where the value is found.
 * 
 * \param vec1 The 1D vector to search.
 * \param value The value to search for in the vector.
 * 
 * \return A vector containing the indices where the value is found in the input vector.
 */
std::vector<int> vec_where_1D(std::vector<int> vec1, int value) {  
    std::vector<int> coordinates;
    int curr;
    int len1 = vec1.size();

    for (int i=0; i<len1; i++) { 
        curr = vec1[i];
        if (curr == value) {
            coordinates.push_back(i);
        } 
    }
    return coordinates;
}

/*! 
 * \brief Grabs specific rows from a 2D integer vector based on given indices.
 * 
 * This function extracts the rows from a 2D vector based on the provided list of indices. 
 * The function returns a new 2D vector containing the selected rows.
 * 
 * \param vec The 2D vector to extract rows from.
 * \param coords A vector containing the indices of the rows to select.
 * 
 * \return A 2D vector containing the selected rows.
 */
std::vector< std::vector<int> > grab_idxs_2D_int_inp(std::vector< std::vector<int> > vec, std::vector<int> coords) {
    std::vector< std::vector<int> > temp_vec;
    int len = (int)coords.size();
    int coord;

    for (int i = 0; i < len; i++) { 
        coord = coords[i];
        temp_vec.push_back(vec[coord]); 
    } 
    for (int i = 0; i < len; i++) {
        vec[i] = temp_vec[i];
    }
    vec.resize(len);
    return vec;
}

/*! 
 * \brief Grabs specific rows from a 2D integer vector based on given indices.
 * 
 * This function extracts the rows from a 2D vector based on the provided list of indices. 
 * The function returns a new 2D vector containing the selected rows.
 * 
 * \param vec The 2D vector to extract rows from.
 * \param coords A vector containing the indices of the rows to select.
 * 
 * \return A 2D vector containing the selected rows.
 */
std::vector< std::vector<int> > grab_idxs_2D_int(std::vector< std::vector<int> > vec, std::vector<int> coords) {
    std::vector< std::vector<int> > output;
    int len = coords.size();

    for (int i = 0; i < len; i++) { 
        int coord = coords[i];
        output.push_back(vec[coord]);
    }

    return output;
}

/*! 
 * \brief Grabs specific rows from a 2D vector (template function).
 * 
 * This template function extracts rows from a 2D vector (or similar container) based on given indices.
 * It can work with different container types by using template parameters.
 * 
 * \param vec The 2D vector (or similar container) to extract rows from.
 * \param coords A vector containing the indices of the rows to select.
 * 
 * \return A container of the selected rows.
 */
template <typename T, typename U>
T grab_idxs_2D(T vec, U coords) {
    T output;
    int len = coords.size();

    for (int i = 0; i < len; i++) { 
        int coord = coords[i];
        output.push_back(vec[coord]);
    }

    return output;
}

/*! 
 * \brief Slices a 1D vector within the specified range.
 * 
 * This function slices a 1D vector, returning a subvector that contains elements from index i1_start to i1_end.
 * The resulting vector contains elements within the specified range.
 * 
 * \param vec The vector to slice.
 * \param i1_start The start index of the slice (inclusive).
 * \param i1_end The end index of the slice (exclusive).
 * 
 * \return A new vector containing the sliced elements.
 */
template <typename T>
T slice_1Dvec(T vec, int i1_start, int i1_end) {
    int size = i1_end - i1_start;
    for (int i = 0; i < size; i++) { 
        vec[i] = vec[i + i1_start];
    }
    vec.resize(size);
    return vec;
}
/**
 * @brief Extracts a slice from a 1D vector of doubles.
 * 
 * This function creates a new vector containing the elements from index `i1_start` to `i1_end - 1` from the input vector `vec`.
 *
 * @param vec Input 1D vector of doubles.
 * @param i1_start The starting index of the slice (inclusive).
 * @param i1_end The ending index of the slice (exclusive).
 * @return A vector containing the sliced elements.
 */
std::vector<double> slice_1Dvec(std::vector<double> vec, int i1_start, int i1_end) {
    int size1 = i1_end - i1_start;
    std::vector<double> slice(size1);

    for (int i = i1_start; i < i1_end; i++) { 
        slice[i] = vec[i];
    }
    return slice;
}

/**
 * @brief Extracts a slice from a 1D vector of strings.
 * 
 * This function creates a new vector containing the elements from index `i1_start` to `i1_end - 1` from the input vector `vec`.
 *
 * @param vec Input 1D vector of strings.
 * @param i1_start The starting index of the slice (inclusive).
 * @param i1_end The ending index of the slice (exclusive).
 * @return A vector containing the sliced elements.
 */
std::vector<std::string> slice_1Dvec(std::vector<std::string> vec, int i1_start, int i1_end) {
    int size = i1_end - i1_start;
    for (int i = 0; i < size; i++) { 
        vec[i] = vec[i + i1_start];
    }
    vec.resize(size);
    return vec;
}

/**
 * @brief Reorders a 1D vector of integers based on given indices.
 * 
 * This function rearranges the elements of the input vector `v` according to the order specified in `idxs`.
 *
 * @param v Input vector of integers.
 * @param idxs A vector of indices specifying the new order.
 * @return A vector with elements reordered as specified by `idxs`.
 */
std::vector<int> reorder_inp(std::vector<int> v, std::vector<int> idxs) {
    std::vector<int> reordered_v((int)v.size());

    for (int i = 0; i < (int)v.size(); i++) {
        reordered_v[i] = v[idxs[i]];
    }
    for (int i = 0; i < (int)v.size(); i++) {
        v[i] = reordered_v[i];
    }

    return v;
}

/**
 * @brief Reorders a 2D vector of integers based on given indices.
 * 
 * This function rearranges the rows of the input 2D vector `v` according to the order specified in `idxs`.
 *
 * @param v Input 2D vector of integers.
 * @param idxs A vector of indices specifying the new order of rows.
 * @return A 2D vector with rows reordered as specified by `idxs`.
 */
std::vector<std::vector<int>> reorder_inp(std::vector<std::vector<int>> v, std::vector<int> idxs) {
    std::vector<std::vector<int>> reordered_v = vect_create_2D((int)v.size(), (int)v[0].size());

    for (int i = 0; i < (int)v.size(); i++) {
        for (int j = 0; j < (int)v[i].size(); j++) {
            reordered_v[i][j] = v[i][idxs[j]];
        }
        for (int j = 0; j < (int)v[i].size(); j++) {
            v[i][j] = reordered_v[i][j];
        }
    }

    return v;
}

/**
 * @brief Reorders a 2D vector of doubles based on given indices.
 * 
 * This function rearranges the rows of the input 2D vector `v` according to the order specified in `idxs`.
 *
 * @param v Input 2D vector of doubles.
 * @param idxs A vector of indices specifying the new order of rows.
 * @return A 2D vector with rows reordered as specified by `idxs`.
 */
std::vector<std::vector<double>> reorder_inp(std::vector<std::vector<double>> v, std::vector<int> idxs) {
    std::vector<std::vector<double>> reordered_v = vect_create_2D_float((int)v.size(), (int)v[0].size(), 0);

    for (int i = 0; i < (int)v.size(); i++) {
        for (int j = 0; j < (int)v[i].size(); j++) {
            reordered_v[i][j] = v[i][idxs[j]];
        }
        for (int j = 0; j < (int)v[i].size(); j++) {
            v[i][j] = reordered_v[i][j];
        }
    }

    return v;
}

/**
 * @brief Reorders a 3D vector of doubles based on given indices.
 * 
 * This function rearranges the rows of the 2D slices in the input 3D vector `v` according to the order specified in `idxs`.
 *
 * @param v Input 3D vector of doubles.
 * @param idxs A vector of indices specifying the new order of rows in each 2D slice.
 * @return A 3D vector with rows reordered as specified by `idxs`.
 */
std::vector<std::vector<std::vector<double>>> reorder_inp(std::vector<std::vector<std::vector<double>>> v, std::vector<int> idxs) {
    std::vector<std::vector<std::vector<double>>> reordered_v = vect_create_3D_float((int)v.size(), (int)v[0].size(), (int)v[0][0].size());

    for (int i = 0; i < (int)v.size(); i++) {
        for (int j = 0; j < (int)v[i].size(); j++) {
            for (int k = 0; k < (int)v[i][j].size(); k++) {
                reordered_v[i][j][k] = v[i][j][idxs[k]];
            }
            for (int k = 0; k < (int)v[i][j].size(); k++) {
                v[i][j][k] = reordered_v[i][j][k];
            }
        }
    }

    return v;
}

/**
 * @brief Reorders a 1D vector of doubles based on given indices.
 * 
 * This function rearranges the elements of the input vector `v` according to the order specified in `idxs`.
 *
 * @param v Input vector of doubles.
 * @param idxs A vector of indices specifying the new order.
 * @return A vector with elements reordered as specified by `idxs`.
 */
std::vector<double> reorder_inp(std::vector<double> v, std::vector<int> idxs) {
    std::vector<double> reordered_v((int)v.size());

    for (int i = 0; i < (int)v.size(); i++) {
        reordered_v[i] = v[idxs[i]];
    }
    for (int i = 0; i < (int)v.size(); i++) {
        v[i] = reordered_v[i];
    }

    return v;
}

/**
 * @brief Generates a 4D index grid between two vectors of start and end indices.
 * 
 * This function generates all possible index combinations between the start and end vectors for each of the 4 dimensions.
 *
 * @param start_vec A vector containing the starting indices for each dimension.
 * @param end_vec A vector containing the ending indices for each dimension.
 * @return A 2D vector of indices representing the grid of possible index combinations.
 */
std::vector<std::vector<int>> FourD_idxs(std::vector<int> start_vec, std::vector<int> end_vec) {
    std::vector<int> idx;
    std::vector<std::vector<int>> idxs;

    for (int i = start_vec[0]; i < (int)(end_vec[0] + 1); i++) {
        for (int j = start_vec[1]; j < (int)(end_vec[1]); j++) {
            for (int k = start_vec[2]; k < ((int)end_vec[2]); k++) {
                for (int l = start_vec[3]; l < (int)(end_vec[3]); l++) { 
                    idx = {i, j, k, l};
                    idxs.push_back(idx);
                }
            }
        }
    }

    return idxs;
}

/**
 * @brief Generate a range of integers from start to end with a given step.
 *
 * This function generates a vector of integers starting at `start`, ending at `end`,
 * and with a specified `step` between each value.
 *
 * @param start The starting value of the range.
 * @param end The end value (exclusive) of the range.
 * @param step The step size between values.
 * @return A vector containing the generated range of integers.
 */
std::vector<int> arange(int start, int end, int step) {
    int size = end - start;
    std::vector<int> output(size);
    int value;

    for (int i = 0; i < size; i++) {
        value = start + i * step;
        output[i] = value;
    }

    return output;
}

/**
 * @brief Recursive binary search to find the position for item insertion in a sorted array.
 *
 * This function recursively searches the correct position where an item should be inserted
 * in a sorted vector such that the vector remains sorted.
 *
 * @param a A pointer to the vector in which the search is performed.
 * @param item The item to be inserted.
 * @param low The starting index of the range to search.
 * @param high The ending index of the range to search.
 * @return The index at which the item should be inserted.
 */
int searchsorted_recursive(std::vector<double>* a, double item, int low, int high) {
    if (high <= low) {
        return (item > (*a)[low]) ? (low + 1) : low;
    }

    int mid = (low + high) / 2;

    if (item == (*a)[mid]) {
        return (mid + 1);
    }

    if (item > (*a)[mid]) {
        return searchsorted_recursive(a, item, (mid + 1), high);
    }

    return searchsorted_recursive(a, item, low, (mid - 1));
}

/**
 * @brief Find the index where the element should be inserted in a sorted vector.
 *
 * This function finds the index at which the given value can be inserted into the sorted vector
 * while maintaining the sorted order. It uses the `lower_bound` function for this.
 *
 * @tparam T The type of elements in the vector.
 * @param vec The sorted vector in which to insert the element.
 * @param value The value to insert.
 * @return The index where the value should be inserted.
 */
template <typename T>
int idx_to_insert(std::vector<T> &vec, T value) {
    auto it = std::lower_bound(vec.begin(), vec.end(), value);
    int index = std::distance(vec.begin(), it);
    return index;
}

/**
 * @brief Slice a 1D vector of strings from a start to end index.
 *
 * This function extracts a portion of a 1D vector of strings from index `i1_start` to `i1_end`
 * and returns it as a new vector.
 *
 * @param vec The original vector of strings to slice.
 * @param i1_start The starting index of the slice.
 * @param i1_end The ending index (exclusive) of the slice.
 * @return A vector containing the sliced portion of the original vector.
 */
std::vector<std::string> slice_1Dvec_str_inp(std::vector<std::string> vec, int i1_start, int i1_end) {
    int size = i1_end - i1_start;
    for (int i = 0; i < size; i++) {
        vec[i] = vec[i + i1_start];
    }
    vec.resize(size);
    return vec;
}

/**
 * @brief Slice a 1D vector of doubles from a start to end index.
 *
 * This function extracts a portion of a 1D vector of doubles from index `i1_start` to `i1_end`
 * and returns it as a new vector.
 *
 * @param vec The original vector of doubles to slice.
 * @param i1_start The starting index of the slice.
 * @param i1_end The ending index (exclusive) of the slice.
 * @return A vector containing the sliced portion of the original vector.
 */
std::vector<double> slice_1Dvec_float(std::vector<double> vec, int i1_start, int i1_end) {
    int size = i1_end - i1_start;
    for (int i = 0; i < size; i++) {
        vec[i] = vec[i + i1_start];
    }
    return vec;
}

/**
 * @brief Slice a 1D vector of strings and return it as a new vector.
 *
 * This function slices a portion of a 1D vector of strings from `i1_start` to `i1_end`
 * and returns a new vector with that portion.
 *
 * @param vec The vector to slice.
 * @param i1_start The starting index of the slice.
 * @param i1_end The ending index (exclusive) of the slice.
 * @return A vector containing the sliced portion of the original vector.
 */
std::vector<std::string> slice_1Dvec_str(std::vector<std::string> vec, int i1_start, int i1_end) {
    int size1 = i1_end - i1_start;
    std::vector<std::string> slice(size1);

    for (int i = 0; i < size1; i++) {
        slice[i] = vec[i + i1_start];
    }
    return slice;
}

/**
 * @brief Find the index of the maximum element in a vector.
 *
 * This function finds the index of the maximum value in a vector of doubles.
 *
 * @param max_rates A pointer to a vector of doubles.
 * @return The index of the maximum element.
 */
int find_max_element(std::vector<double>* max_rates) {
    int max_idx;
    double max_val = 0;

    for (int i = 0; i < (int)(*max_rates).size(); i++) {
        if ((*max_rates)[i] >= max_val) {
            max_idx = i;
            max_val = (*max_rates)[i];
        }
    }
    return max_idx;
}

/**
 * @brief Find the index of the minimum element in a vector.
 *
 * This function finds the index of the minimum value in a vector of doubles.
 *
 * @param max_rates A pointer to a vector of doubles.
 * @return The index of the minimum element.
 */
int find_min_element(std::vector<double>* max_rates) {
    int max_idx;
    double max_val = 0;

    for (int i = 0; i < (int)(*max_rates).size(); i++) {
        if ((*max_rates)[i] <= max_val) {
            max_idx = i;
            max_val = (*max_rates)[i];
        }
    }
    return max_idx;
}

/**
 * @brief Check if two arrays of integers are equal.
 *
 * This function checks if two arrays of integers of the same length are equal element-wise.
 *
 * @param elem1 The first array of integers.
 * @param elem2 The second array of integers.
 * @param len The length of the arrays.
 * @return True if the arrays are equal, false otherwise.
 */
bool equal(int* elem1, int* elem2, int len) {
    for (int i = 0; i < len; i++) {
        if (elem1[i] != elem2[i]) return false;
    }
    return true;
}

/**
 * @brief Compare two matrices and return elements from the first matrix that are not present in the second.
 *
 * This function compares two matrices element-wise and returns the rows from the first matrix
 * that are not found in the second matrix.
 *
 * @param mat1 The first matrix.
 * @param mat2 The second matrix.
 * @return A new matrix containing the rows from `mat1` not found in `mat2`.
 */
Matrix<int> comparison(Matrix<int>& mat1, Matrix<int>& mat2) {
    int vec_size = std::max((int)mat1.rows(), (int)mat2.rows());
    Matrix<int> mat_out(vec_size, 4);

    int elem = 0;
    bool in_mat = false;
    int idx = 0;

    if ((int)mat1.rows() > (int)mat2.rows()) {
        in_mat = false;
        for (int i1 = 0; i1 < (int)mat1.rows(); i1++) {
            for (int i2 = 0; i2 < (int)mat2.rows(); i2++) {
                if (equal(mat1[i1], mat2[i2], mat1.cols())) {
                    in_mat = true;
                }
            }
            if (!in_mat) {
                for (int j = 0; j < (int)mat_out.cols(); j++) {
                    mat_out[idx][j] = mat1[i1][j];
                }
                idx++;
            }
        }
    }
    else {
        for (int i1 = 0; i1 < (int)mat2.rows(); i1++) {
            in_mat = false;
            for (int i2 = 0; i2 < (int)mat1.rows(); i2++) {
                if (equal(mat2[i1], mat1[i2], mat1.cols())) {
                    in_mat = true;
                }
            }
            if (!in_mat) {
                for (int j = 0; j < (int)mat_out.cols(); j++) {
                    mat_out[idx][j] = mat2[i1][j];
                }
                idx++;
            }
        }
    }

    mat_out.reshape((idx), 4);
    return mat_out;
}

/**
 * @brief Check if a vector of vectors contains a sub-vector.
 *
 * This function checks if a vector of vectors contains a specified sub-vector.
 *
 * @param v The vector of vectors.
 * @param sub_v The sub-vector to check for.
 * @return True if the sub-vector is found in the vector, false otherwise.
 */
bool is_in(std::vector<std::vector<size_t>> v, std::vector<size_t> sub_v) {
    for (int i = 0; i < (int)v.size(); i++) {
        if (v[i] == sub_v) { return true; }
    }
    return false;
}

/**
 * @brief Check if a vector of vectors contains a sub-vector.
 *
 * This function checks if a vector of vectors of integers contains a specified sub-vector.
 *
 * @param v The vector of vectors.
 * @param sub_v The sub-vector to check for.
 * @return True if the sub-vector is found in the vector, false otherwise.
 */
bool is_in(std::vector<std::vector<int>> v, std::vector<int> sub_v) {
    for (int i = 0; i < (int)v.size(); i++) {
        if (v[i] == sub_v) { return true; }
    }
    return false;
}

/**
 * @brief Check if a vector of vectors contains a sub-vector.
 *
 * This function checks if a vector of vectors of doubles contains a specified sub-vector.
 *
 * @param v The vector of vectors.
 * @param sub_v The sub-vector to check for.
 * @return True if the sub-vector is found in the vector, false otherwise.
 */
bool is_in(std::vector<std::vector<double>> v, std::vector<double> sub_v) {
    for (int i = 0; i < (int)v.size(); i++) {
        if (v[i] == sub_v) { return true; }
    }
    return false;
}
/**
 * @brief Checks if an element is present in a vector.
 * 
 * This function iterates through the given vector and checks if the specified element
 * exists in the vector. Returns true if the element is found, false otherwise.
 * 
 * @param check_vec The vector in which the element is searched.
 * @param elem The element to search for in the vector.
 * @return True if the element is found, false otherwise.
 */
bool is_in(std::vector<int> check_vec, int elem) {
    bool elem_in;

    for (int i = 0; i < (int)check_vec.size(); i++) {
        if (check_vec[i] == elem) return true;
    }
    return false;
}

/**
 * @brief Sums vectors across all MPI processes.
 * 
 * This function performs a reduction (MPI_Reduce) to sum the elements of a vector 
 * across multiple processes in an MPI parallel environment. The resulting sum 
 * for each element is stored in the output vector.
 * 
 * @param vec The input vector whose elements are summed.
 * @param nprocs The number of processes in the MPI environment.
 * @param rank The rank of the current process in the MPI environment.
 * @return A vector of summed elements after the reduction operation.
 */
std::vector<int> sum_vectors_allprocs(const std::vector<int> vec, int nprocs, int rank) {
    std::vector<int> output(vec.size());
    int value = 0;
    int reduction_result = 0;
    print_1Dvector(vec);
    
    for (int i = 0; i < vec.size(); i++) {
        if (rank == 0) {
            std::cout << "i: " << i << "\n";
            print_1Dvector(output);
        }
        reduction_result = 0;
        value = vec[i];
        MPI_Reduce(&value, &reduction_result, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
        output[i] = reduction_result;
    }

    return output;
}

/**
 * @brief Sums vectors across all MPI processes (for double-precision elements).
 * 
 * Similar to the `sum_vectors_allprocs` function, but this one operates on a vector
 * of double-precision floating point numbers. The elements of the vector are summed
 * across all processes in the MPI environment.
 * 
 * @param vec The input vector of double values to sum.
 * @param nprocs The number of processes in the MPI environment.
 * @param rank The rank of the current process in the MPI environment.
 * @return A vector of summed elements after the reduction operation.
 */
std::vector<double> sum_vectors_allprocs(const std::vector<double> vec, int nprocs, int rank) {
    std::vector<double> output(vec.size());
    double value = 0;
    double reduction_result = 0;
    print_1Dvector(vec);

    for (int i = 0; i < vec.size(); i++) {
        if (rank == 0) {
            std::cout << "i: " << i << "\n";
            print_1Dvector(output);
        }
        reduction_result = 0;
        value = vec[i];
        MPI_Reduce(&value, &reduction_result, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        output[i] = reduction_result;
    }

    return output;
}
