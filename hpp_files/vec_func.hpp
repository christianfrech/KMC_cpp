#include <iostream>
#include <sstream>
#include <stdio.h>
#include <string>
#include <vector>
#include <cstdlib>
#include <fstream>

struct lattice_return_struct {
    public:
        /*! \brief Constructor
        */

        lattice_return_struct(std::vector< std::vector< std::vector<int> > > all_vacancies,
        std::vector<int> move_counts, 
        std::vector<double> time_count, 
        std::vector<double> all_times):
            all_vacancies_(all_vacancies),
            move_counts_(move_counts),
            time_count_(time_count),
            all_times_(all_times)
            {}

        /*! \return get x field of struc */
        std::vector< std::vector< std::vector<int> > > get_all_vacancies() const {
            return all_vacancies_;
        }

        /*! \return get y field of struc */
        std::vector<int> get_move_counts() const {
            return move_counts_;
        }

        /*! \return get z field of struc */
        std::vector<double> get_time_count() const {
            return time_count_;
        }

        /*! \return get z field of struc */
        std::vector<double> get_all_times() const {
            return all_times_;
        }
        
    private:
        std::vector< std::vector< std::vector<int> > > all_vacancies_; 
        std::vector<int> move_counts_; 
        std::vector<double> time_count_; 
        std::vector<double> all_times_;
};

typedef struct lattice_return_struct lattice_return_struct;

std::vector< std::vector<int> > vect_create_2D(size_t N, size_t M) {
    std::vector< std::vector<int> > vec_out(N);

    for(int i = 0; i < (int)N; i++) { 
        vec_out[i] = std::vector<int>(M);
    }

    return vec_out;
}

std::vector< std::vector< std::vector<double> > > vect_create_3D_float(int L, int N, int M) {

    std::vector< std::vector< std::vector<double> > > vec_out(L, std::vector< std::vector<double> >(N, std::vector<double>(M, 0)));
    return vec_out;
}

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

std::vector< std::vector<double> > vect_create_2D_float(size_t N, size_t M, double val) {
    std::vector<std::vector<double>> vec_out(N, std::vector<double>(M, val));
    return vec_out;
}

std::vector<int> vec_where_1D(std::vector<int> vec1, int value) {  
    std::vector<int> coordinates;
    int curr;
    int len1 = vec1.size();

    for (int i=0;i<len1;i++) { 
            curr = vec1[i];
            if (curr == value) {
                coordinates.push_back(i);
            } 
        }
    //std::cout << "coordinates: \n";
    //print_1Dvector(coordinates);
    return coordinates;
}

std::vector< std::vector<int> > grab_idxs_2D_int_inp(std::vector< std::vector<int> > vec, std::vector<int> coords) {        
    std::vector< std::vector<int> > temp_vec;
    int len = (int)coords.size();
    int coord;

    for (int i=0;i<len;i++) { 
        coord = coords[i];
        temp_vec.push_back(vec[coord]); 
    } 
    for (int i=0;i<len;i++) {
        vec[i] = temp_vec[i];
    }
    vec.resize(len);
    return vec;
}

std::vector< std::vector<int> > grab_idxs_2D_int(std::vector< std::vector<int> > vec, std::vector<int> coords) {        
    //std::cout << "runnin that grab_2d_int(): \n";
    std::vector< std::vector<int> > output;
    int len = coords.size();
    //std::cout << "vec: " << "\n";
    //print_2Dvector(vec);

    //std::cout << "\ncoords: " << "\n";
    //print_1Dvector(coords);

    for (int i=0;i<len;i++) { 
	    int coord = coords[i];
        //std::cout << "coord: \n";
        //std::cout << coord << "\n";
        output.push_back(vec[coord]);
        //std::cout << "vec[coord]: \n";
        //print_1Dvector(vec[coord]);

    } 

    //std::cout << "output: \n";
   //s print_2Dvector(output);
    return output;
}

std::vector<int> create_vec_1D(int size, int value) {
    std::vector<int> vec(size, value);
    return vec;
}

std::vector<double> create_vec_1D_float(int size) {
    std::vector<double> vec(size);
    return vec;
}

std::vector<double> slice_1Dvec_float(std::vector<double> vec, int i1_start, int i1_end) {
    int size1 = i1_end - i1_start;
    std::vector<double> slice(size1);

    for (int i=i1_start; i<i1_end; i++) { 
        slice[i] = vec[i];
    }
    return slice;
}

std::vector<std::string> slice_1Dvec_str(std::vector<std::string> vec, int i1_start, int i1_end) {
    int size1 = i1_end - i1_start;
    std::vector<std::string> slice(size1);


    for (int i=0; i<size1; i++) {
        slice[i] = vec[i+i1_start];
    }
    return slice;
}

std::vector<std::string> slice_1Dvec_str_inp(std::vector<std::string> vec, int i1_start, int i1_end) {
    int size = i1_end - i1_start;
    for (int i=0; i<size; i++) { 
        vec[i] = vec[i+i1_start];
    }
    vec.resize(size);
    return vec;
}

std::vector<int> reorder_ints_inp(std::vector<int> v, std::vector<int> idxs) {
    
    std::vector<int> reordered_v((int)v.size());

    for (int i=0; i<(int)v.size(); i++) {
        reordered_v[i] = v[idxs[i]];
    }
    for (int i=0; i<(int)v.size(); i++) {
        v[i] = reordered_v[i];
    }

    return v;
}

std::vector< std::vector<int> > reorder_ints_inp_2D(std::vector< std::vector<int> > v, std::vector<int> idxs) {
    
    std::vector< std::vector<int> > reordered_v = vect_create_2D((int)v.size(), (int)v[0].size());

    for (int i=0; i<(int)v.size(); i++) {
        for (int j=0; j<(int)v[i].size(); j++) {
            reordered_v[i][j] = v[i][idxs[j]];
        }
        for (int j=0; j<(int)v[i].size(); j++) {
            v[i][j] = reordered_v[i][j];
        }
    }

    return v;
}

std::vector< std::vector<double> > reorder_dbl_inp_2D(std::vector< std::vector<double> > v, std::vector<int> idxs) {
    
    std::vector< std::vector<double> > reordered_v = vect_create_2D_float((int)v.size(), (int)v[0].size(), 0);

    for (int i=0; i<(int)v.size(); i++) {
        for (int j=0; j<(int)v[i].size(); j++) {
            reordered_v[i][j] = v[i][idxs[j]];
        }
        for (int j=0; j<(int)v[i].size(); j++) {
            v[i][j] = reordered_v[i][j];
        }
    }

    return v;
}

std::vector< std::vector< std::vector<double> > > reorder_dbl_inp_3D(std::vector< std::vector< std::vector<double> > > v, std::vector<int> idxs) {

    std::vector< std::vector< std::vector<double> > > reordered_v = vect_create_3D_float((int)v.size(), (int)v[0].size(), (int)v[0][0].size());
    
    for (int i=0; i<(int)v.size(); i++) {
        for (int j=0; j<(int)v[i].size(); j++) {
            for (int k=0; k<(int)v[i][j].size(); k++) {
                reordered_v[i][j][k] = v[i][j][idxs[k]];
            }
            for (int k=0; k<(int)v[i][j].size(); k++) {
                v[i][j][k] = reordered_v[i][j][k];
            }
        }
    }

    return v;
}

std::vector<double> reorder_floats_inp(std::vector<double> v, std::vector<int> idxs) {
    std::vector<double> reordered_v((int)v.size());

    for (int i=0; i<(int)v.size(); i++) {
        reordered_v[i] = v[idxs[i]];
    }
    for (int i=0; i<(int)v.size(); i++) {
        v[i] = reordered_v[i];
    }

    return v;
}

std::vector< std::vector<int> > FourD_idxs(std::vector<int> start_vec, std::vector<int> end_vec) {
    std::vector<int> idx;
    std::vector< std::vector<int> > idxs;

    for (int i=start_vec[0]; i<(int)(end_vec[0] + 1); i++) {
        for (int j=start_vec[1]; j<(int)(end_vec[1]); j++) {
            for (int k=start_vec[2]; k<((int)end_vec[2]); k++) {
                for (int l=start_vec[3]; l<(int)(end_vec[3]); l++) { 
                    //std::cout << "[ " << i << " " << j << " " << k << " " << l << " ] \n";
                    idx = {i, j, k, l};
                    idxs.push_back(idx);
                }
            }
        }
    }

    return idxs;
}

std::vector<int> arange(int start, int end, int step) {
    int size = end-start;
    std::vector<int> output(size);
    int value;

    for (int i=0; i<size; i++) {
        value = start + i*step;
        output[i] = value;
    }

    return output;
}

// A binary search based function
// to find the position
// where item should be inserted
// in a[low..high]
int searchsorted_recursive(std::vector<double>* a, double item, int low, int high) {
    if (high <= low) { return (item > (*a)[low]) ? (low + 1) : low; }
    
    int mid = (low + high) / 2;
 
    if (item == (*a)[mid]) { return (mid + 1); }
 
    if (item > (*a)[mid]) {return searchsorted_recursive(a, item, (mid + 1), high); }
    
    return searchsorted_recursive(a, item, low, (mid - 1));
}


int find_max_element(std::vector<double> * max_rates) {
    int max_idx;
    double max_val = 0;

    for (int i=0; i<(int)(*max_rates).size(); i++) {
        if ( (*max_rates)[i] >= max_val ) {
            max_idx = i;
            max_val = (*max_rates)[i];
        }
    }
    return max_idx;
}