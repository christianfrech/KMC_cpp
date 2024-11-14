#include <iostream>
#include <sstream>
#include <stdio.h>
#include <string>
#include <vector>
#include <cstdlib>
#include <fstream>
#include <cmath>
#include <math.h>


int exp_int(int base, int exponent) {
    int res = 0;
    
    if (exponent == 0) {return 1;}
    else {
        for (int i=0; i<(exponent); i++) {
            if (i == 0) {
                res = base;
            }
            else {res = res*base;}
        }
    }
    return res;
}

double sum_float(std::vector<double> vec) {
    double sum = 0;
    for (int i=0; i<(int)vec.size(); i++) {
        sum += vec[i];
        
        }

    return sum;
}

int NCR(int n, int r) {
    if (r == 0) return 1;

    /*
     Extra computation saving for large R,
     using property:
     N choose R = N choose (N-R)
    */
    if (r > n / 2) return NCR(n, n - r);

    long res = 1;

    for (int k = 1; k <= r; ++k)
    {
        res *= n - k + 1;
        res /= k;
    }

    return res;
}

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

int create_partition(int a, int b) {
    if ((a%2 == 0) && (b%2 == 0)) {
        if (a == b) { return find_gcf(a, (int)(b/2));}
        else return find_gcf(a, b);
    }
    else return find_gcf(a, b);
}

std::vector<size_t> mod_with_bounds(int i, int j, int bound_x, int bound_y) {
    std::vector<size_t> idxs(2);

    idxs[0] = (size_t)((i % bound_x + bound_x) % bound_x);
    idxs[1] = (size_t)((j % bound_y + bound_y) % bound_y);

    return idxs;
}

size_t mod_with_bounds(int i, int bound_x) {
    size_t idx = (size_t)((i % bound_x + bound_x) % bound_x);
    return idx;
}
