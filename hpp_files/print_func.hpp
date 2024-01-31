#include <iostream>
#include <sstream>
#include <stdio.h>
#include <string>
#include <vector>
#include <cstdlib>
#include <fstream>



void print_3Dvector_float(std::vector< std::vector< std::vector<double> > > vec) {
    int M = (int)vec.size();;
    std::cout << "[ ";
    for (int m=0; m<M; m++) {
        std::cout << "[ ";
        for (int n=0; n<(int)vec[m].size(); n++) {
            std::cout << "[ ";
            for (int o=0; o<(int)vec[m][n].size(); o++) {
                std::cout << vec[m][n][o] << " ";
            }
            std::cout << "] ";
        }
        std::cout << "] ";
        std::cout << "\n  ";
    }
    std::cout << "] \n\n";
}

void print_3Dvector(std::vector< std::vector< std::vector<int> > > vec) {
    int M = (int)vec.size();;
    std::cout << "[ ";
    for (int m=0; m<M; m++) {
        std::cout << "[ ";
        for (int n=0; n<(int)vec[m].size(); n++) {
            std::cout << "[ ";
            for (int o=0; o<(int)vec[m][n].size(); o++) {
                std::cout << vec[m][n][o] << " ";
            }
            std::cout << "] ";
        }
        std::cout << "] ";
        std::cout << "\n  ";
    }
    std::cout << "] \n\n";
}

void print_2Dvector_float(std::vector< std::vector<double> > vec) {
    int M = (int)vec.size();;
    std::cout << "[ ";
    for (int m=0; m<M; m++) {
        std::cout << "[ ";
        for (int n=0; n<(int)vec[m].size(); n++) {
            std::cout << vec[m][n] << " ";
        }
        std::cout << "] ";
        std::cout << "\n  ";
    }
    std::cout << "] \n\n";
}

void print_2Dvector(std::vector< std::vector<int> > vec) {
    int M = (int)vec.size();;
    std::cout << "[ ";
    for (int m=0; m<M; m++) {
        std::cout << "[ ";
        for (int n=0; n<(int)vec[m].size(); n++) {
            std::cout << vec[m][n] << " ";
        }
        std::cout << "] ";
        std::cout << "\n  ";
    }
    std::cout << "] \n\n";
}

void print_1Dvector(std::vector<int> vec) {
    int N = (int)vec.size();
    std::cout << "[ ";
    for (int n=0; n<N; n++) {
        std::cout << vec[n] << " ";
    }
    std::cout << "] \n\n";
}

void print_1Darr(int* arr, int size) {
    std::cout << "[ ";
    for (int n=0; n<size; n++) {
        std::cout << arr[n] << " ";
    }
    std::cout << "] \n\n";
}

void print_1Dvector_float(std::vector<double> vec) {
    int N = (int)vec.size();
    std::cout << "[ ";
    for (int n=0; n<N; n++) {
        std::cout << vec[n] << " ";
    }
    std::cout << "] \n\n";
}

void print_1Dvector_string(std::vector<std::string> vec) {
    int N = (int)vec.size();
    std::cout << "[ ";
    for (int n=0; n<N; n++) {
        std::cout << vec[n] << " ";
    }
    std::cout << "] \n\n";
}

void print_1Dvector_size_t(std::vector<size_t> &vec) {
    int N = (int)vec.size();
    std::cout << "[ ";
    for (int n=0; n<N; n++) {
        std::cout << vec[n] << " ";
    }
    std::cout << "] \n\n";
}