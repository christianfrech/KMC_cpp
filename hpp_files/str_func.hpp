#include <iostream>
#include <sstream>
#include <stdio.h>
#include <string>
#include <vector>
#include <cstdlib>
#include <fstream>

#include "classes_strucs.hpp"


struct reg_line_struct {
    public:
        /*! \brief Constructor
        */

        reg_line_struct(std::string lattice_type, int x, int y, int z):
            lattice_type_(lattice_type),
            x_(x),
            y_(y),
            z_(z)
            {}

        /*! \return get lattice_type field of struc */
        std::string get_latice_type() const {
            return lattice_type_;
        }

        /*! \return get x field of struc */
        int get_x() const {
            return x_;
        }

        /*! \return get y field of struc */
        int get_y() const {
            return y_;
        }

        /*! \return get y field of struc */
        int get_z() const {
            return z_;
        }
        
    private:
        std::string lattice_type_; int x_; int y_; int z_;
};
 
typedef struct reg_line_struct reg_line_struct;


struct line_struct {
    public:
        /*! \brief Constructor
        */

        line_struct(std::string lattice_pos, double x, double y, double z, int atype):
            lattice_pos_(lattice_pos),
            x_(x),
            y_(y),
            z_(z),
            atype_(atype)
            {}

        /*! \return get lattice_type field of struc */
        std::string get_latice_pos() const {
            return lattice_pos_;
        }

        /*! \return get x field of struc */
        double get_x() const {
            return x_;
        }

        /*! \return get y field of struc */
        double get_y() const {
            return y_;
        }

        /*! \return get z field of struc */
        double get_z() const {
            return z_;
        }

        /*! \return get z field of struc */
        int get_atype() const {
            return atype_;
        }
        
    private:
        std::string lattice_pos_; double x_; double y_; double z_; int atype_;
};
 
typedef struct line_struct line_struct;


struct ratecatalog_struct {
    public:
        /*! \brief Constructor
        */

        ratecatalog_struct(std::vector< std::vector<int> > configs, 
        std::vector< std::vector< std::vector<double> > > energies, 
        std::vector< std::vector< std::vector< std::vector<double> > > > regions_cat, 
        int region_num):
            configs_(configs),
            energies_(energies),
            regions_cat_(regions_cat),
            region_num_(region_num)
            {}

        /*! \return get x field of struc */
        std::vector< std::vector<int> > get_configs() const {
            return configs_;
        }

        /*! \return get y field of struc */
        std::vector< std::vector< std::vector<double> > > get_energies() const {
            return energies_;
        }

        /*! \return get z field of struc */
        std::vector< std::vector< std::vector< std::vector<double> > > > get_regions_cat() const {
            return regions_cat_;
        }

        /*! \return get z field of struc */
        int get_region_num() const {
            return region_num_;
        }
        
    private:
        std::vector< std::vector<int> > configs_; 
        std::vector< std::vector< std::vector<double> > > energies_; 
        std::vector< std::vector< std::vector< std::vector<double> > > > regions_cat_; 
        int region_num_;
};
 
typedef struct ratecatalog_struct ratecatalog_struct;

std::vector<std::string> tokenizer(std::string s, std::string delimiter) {
    size_t pos_start = 0, pos_end, delim_len = delimiter.length();
    std::string token;
    std::vector<std::string> res;

    while ((pos_end = s.find(delimiter, pos_start)) != std::string::npos) {
        token = s.substr (pos_start, pos_end - pos_start);
        pos_start = pos_end + delim_len;
        res.push_back (token);
    }
    if (s.substr (pos_start, pos_end - pos_start) != "") {
        res.push_back(s.substr(pos_start));
    }
    else {}
    
    return res;
}

std::vector<char> split_by_char(std::string s) {
    std::vector<char> v(s.begin(), s.end());
    return v;
}

bool str_isalpha(std::string str) {
    auto it = std::find_if(str.begin(), str.end(), [](char const &c) {
        return !std::isalpha(c);
    });
    return it == str.end();
}

line_struct parse_line(std::string line) {

    std::vector<std::string> toks = tokenizer(line, " ");
    std::string lattice_pos = (toks[0]);

    double x = std::stof(toks[1]);
    double y = std::stof(toks[2]);
    double z = std::stof(toks[3]);
    int atype = std::stoi(toks[4]);

    line_struct output_vals(lattice_pos,x,y,z,atype);
    return output_vals;
}

reg_line_struct parse_reg_line(std::string line) {

    std::vector<std::string> toks = tokenizer(line, " ");
    std::string lattice_pos = (toks[0]);
    double x = std::stof(toks[1]);
    double y = std::stof(toks[2]);
    double z = std::stof(toks[3]);

    reg_line_struct output_vals(lattice_pos,x,y,z);
    return output_vals;
}


void write_to_file(std::string filename, Matrix<int>* values) {
    std::ofstream out_file;
    out_file.open(filename);
    std::string s;

    if (out_file.is_open()) {
        for (int i=0; i<(int)values->rows(); i++) {
            for (int j=0; j<(int)values->cols(); j++) {     
                out_file << (*values)[i][j] << " ";
            }
            out_file << "\n";  
        }
    }

    out_file.close();
}


void write_to_file(std::string filename, std::vector< std::vector<int> > values) {
    std::ofstream out_file;
    out_file.open(filename);
    std::string s;

    if (out_file.is_open()) {
        for (int i=0; i<(int)values.size(); i++) {
            for (int j=0; j<(int)values[i].size(); j++) { 
                out_file << values[i][j] << " ";
            }
        out_file << "\n";
        }
    }

    out_file.close();
}
