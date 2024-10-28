#include <iostream>
#include <sstream>
#include <stdio.h>
#include <string>
#include <vector>
#include <cstdlib>
#include <fstream>
#include <cmath>
#include <math.h>

#include <iostream>
#include <sstream>
#include <stdio.h>
#include <string>
#include <vector>
#include <cstdlib>
#include <fstream>

#include "classes_strucs.hpp"

/*! \brief Calculates the power of an integer base raised to an exponent.
 * \param [in] base The base integer.
 * \param [in] exponent The exponent integer.
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

/*! \brief Sums the elements in a vector of doubles.
 * \param [in] vec The vector of double values to be summed.
 * \return The sum of all elements in the vector.
 */
double sum_float(std::vector<double> vec) {
    double sum = 0;
    for (int i = 0; i < (int)vec.size(); i++) {
        sum += vec[i];
    }
    return sum;
}

/*! \brief Calculates the binomial coefficient, "n choose r."
 * \param [in] n The total number of items.
 * \param [in] r The number of items to choose.
 * \return The binomial coefficient (n choose r).
 */
int NCR(int n, int r) {
    if (r == 0) return 1;

    /*! 
     * \note Extra computation saving for large r, using property:
     * N choose R = N choose (N - R)
     */
    if (r > n / 2) return NCR(n, n - r);

    long res = 1;

    for (int k = 1; k <= r; ++k) {
        res *= n - k + 1;
        res /= k;
    }

    return res;
}

/*! \brief Finds the greatest common factor (GCF) of two integers.
 * \param [in] a The first integer.
 * \param [in] b The second integer.
 * \return The greatest common factor of a and b.
 */
int find_gcf(int a, int b) {
    std::cout << "a: " << a << " b: " << b << "\n";

    // Everything divides 0 
    if (a == 0) 
        return b; 
    if (b == 0) 
        return a; 
  
    // base case 
    if (a == b) 
        return a; 
  
    // a is greater 
    if (a > b) 
        return find_gcf(a - b, b); 
    return find_gcf(a, b - a);
}

/*! \brief Creates a partition by finding the GCF of two integers.
 * \param [in] a The first integer.
 * \param [in] b The second integer.
 * \return The GCF of a and b, with conditional adjustments.
 */
int create_partition(int a, int b) {
    if ((a % 2 == 0) && (b % 2 == 0)) {
        if (a == b) { return find_gcf(a, (int)(b / 2)); }
        else return find_gcf(a, b);
    } else return find_gcf(a, b);
}

/*! \brief reduces two integers modulo given bounds, ensuring they fall within specified bounds.
 * \param [in] i The first integer.
 * \param [in] j The second integer.
 * \param [in] bound_x The bound for the first integer.
 * \param [in] bound_y The bound for the second integer.
 * \return A vector of size_t values containing the modulo results for i and j.
 */
std::vector<size_t> mod_with_bounds(int i, int j, int bound_x, int bound_y) {
    std::vector<size_t> idxs(2);

    idxs[0] = (size_t)((i % bound_x + bound_x) % bound_x);
    idxs[1] = (size_t)((j % bound_y + bound_y) % bound_y);

    return idxs;
}

/*! \brief Reduces an integer modulo a given bound, ensuring it falls within the bound.
 * \param [in] i The integer to modify.
 * \param [in] bound_x The bound.
 * \return A size_t value containing the modulo result of i.
 */
size_t mod_with_bounds(int i, int bound_x) {
    size_t idx = (size_t)((i % bound_x + bound_x) % bound_x);
    return idx;
}
