#include <iostream>
#include <sstream>
#include <stdio.h>
#include <string>
#include <vector>
#include <cstdlib>
#include <fstream>
#include <cmath>
#include <math.h>

/*!
 * \brief Computes the integer power of a base raised to an exponent.
 * \param base The base value.
 * \param exponent The exponent value.
 * \return The result of base raised to the power of exponent.
 */
int exp_int(int base, int exponent) {
    int res = 0;
    if (exponent == 0) {
        return 1;
    } else {
        for (int i = 0; i < exponent; i++) {
            if (i == 0) {
                res = base;
            } else {
                res = res * base;
            }
        }
    }
    return res;
}

/*!
 * \brief Computes the sum of all elements in a vector of doubles.
 * \param vec The vector of doubles.
 * \return The sum of all elements in the vector.
 */
double sum_float(std::vector<double> vec) {
    double sum = 0;
    for (int i = 0; i < (int)vec.size(); i++) {
        sum += vec[i];
    }
    return sum;
}

/*!
 * \brief Computes the binomial coefficient (n choose r).
 * \param n The total number of items.
 * \param r The number of items to choose.
 * \return The binomial coefficient for n and r.
 */
int NCR(int n, int r) {
    if (r == 0) return 1;

    /*!
     * \note Optimized for cases where r > n / 2 using the property:
     *       N choose R = N choose (N - R).
     */
    if (r > n / 2) return NCR(n, n - r);

    long res = 1;
    for (int k = 1; k <= r; ++k) {
        res *= n - k + 1;
        res /= k;
    }
    return res;
}

/*!
 * \brief Finds the greatest common factor (GCF) of two integers.
 * \param a The first integer.
 * \param b The second integer.
 * \return The greatest common factor of a and b.
 */
int find_gcf(int a, int b) {
    std::cout << "a: " << a << " b: " << b << "\n";

    // Base cases
    if (a == 0) return b;
    if (b == 0) return a;
    if (a == b) return a;

    // Recursive cases
    if (a > b)
        return find_gcf(a - b, b);

    return find_gcf(a, b - a);
}

/*!
 * \brief Creates a partition of two integers based on their parity.
 * \param a The first integer.
 * \param b The second integer.
 * \return The GCF of the integers after partitioning.
 */
int create_partition(int a, int b) {
    if ((a % 2 == 0) && (b % 2 == 0)) {
        if (a == b) {
            return find_gcf(a, (int)(b / 2));
        } else {
            return find_gcf(a, b);
        }
    } else {
        return find_gcf(a, b);
    }
}

/*!
 * \brief Computes the modulo operation for two integers within specified bounds (2D version).
 * \param i The first integer.
 * \param j The second integer.
 * \param bound_x The upper bound for the first dimension.
 * \param bound_y The upper bound for the second dimension.
 * \return A vector containing the modulo results for i and j within bounds.
 */
std::vector<size_t> mod_with_bounds(int i, int j, int bound_x, int bound_y) {
    std::vector<size_t> idxs(2);
    idxs[0] = (size_t)((i % bound_x + bound_x) % bound_x);
    idxs[1] = (size_t)((j % bound_y + bound_y) % bound_y);
    return idxs;
}

/*!
 * \brief Computes the modulo operation for a single integer within a specified bound.
 * \param i The integer to be wrapped.
 * \param bound_x The upper bound.
 * \return The modulo result for i within the specified bound.
 */
size_t mod_with_bounds(int i, int bound_x) {
    size_t idx = (size_t)((i % bound_x + bound_x) % bound_x);
    return idx;
}
