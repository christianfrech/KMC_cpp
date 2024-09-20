#include <iostream>
#include <sstream>
#include <stdio.h>
#include <string>
#include <vector>
#include <cstdlib>
#include <fstream>
#include <cmath>
#include <math.h>


template <typename T>
T exp(T base, int exponent) {
    T res = 0;
    
    if (exponent == 0) {return 1;}
    else {
        for (int i=0; i<(exponent); i++) {
            if (i == 0) {res = base;}
            else {res = res*base;}
        }
    }
    return res;
}

int exp(int base, int exponent) {
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