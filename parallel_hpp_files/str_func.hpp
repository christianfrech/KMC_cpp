#include <iostream>
#include <sstream>
#include <stdio.h>
#include <string>
#include <vector>
#include <cstdlib>
#include <fstream>
#include <filesystem>

#include "classes_strucs.hpp"

/*!
 * \brief Structure for storing a regular line with lattice type and coordinates.
 * 
 * This structure holds information about a lattice line including its type and three-dimensional coordinates (x, y, z).
 */
struct reg_line_struct {
    public:
        /*! \brief Constructor
        * 
        * Initializes the `reg_line_struct` with lattice type and 3D coordinates.
        * 
        * \param lattice_type A string representing the type of the lattice.
        * \param x The x-coordinate of the lattice.
        * \param y The y-coordinate of the lattice.
        * \param z The z-coordinate of the lattice.
        */
        reg_line_struct(std::string lattice_type, int x, int y, int z):
            lattice_type_(lattice_type),
            x_(x),
            y_(y),
            z_(z)
            {}

        /*! \return Returns the lattice type */
        std::string get_latice_type() const {
            return lattice_type_;
        }

        /*! \return Returns the x-coordinate */
        int get_x() const {
            return x_;
        }

        /*! \return Returns the y-coordinate */
        int get_y() const {
            return y_;
        }

        /*! \return Returns the z-coordinate */
        int get_z() const {
            return z_;
        }
        
    private:
        std::string lattice_type_; ///< Lattice type
        int x_; ///< x-coordinate
        int y_; ///< y-coordinate
        int z_; ///< z-coordinate
};
 
typedef struct reg_line_struct reg_line_struct;


/*!
 * \brief Structure for storing a line with lattice position, 3D coordinates, and atom type.
 * 
 * This structure contains a lattice position (as a string), three-dimensional coordinates (x, y, z), and an atom type identifier.
 */
struct line_struct {
    public:
        /*! \brief Constructor
        * 
        * Initializes the `line_struct` with lattice position, coordinates, and atom type.
        * 
        * \param lattice_pos A string representing the lattice position.
        * \param x The x-coordinate of the lattice.
        * \param y The y-coordinate of the lattice.
        * \param z The z-coordinate of the lattice.
        * \param atype An integer representing the atom type.
        */
        line_struct(std::string lattice_pos, double x, double y, double z, int atype):
            lattice_pos_(lattice_pos),
            x_(x),
            y_(y),
            z_(z),
            atype_(atype)
            {}

        /*! \return Returns the lattice position */
        std::string get_latice_pos() const {
            return lattice_pos_;
        }

        /*! \return Returns the x-coordinate */
        double get_x() const {
            return x_;
        }

        /*! \return Returns the y-coordinate */
        double get_y() const {
            return y_;
        }

        /*! \return Returns the z-coordinate */
        double get_z() const {
            return z_;
        }

        /*! \return Returns the atom type */
        int get_atype() const {
            return atype_;
        }
        
    private:
        std::string lattice_pos_; ///< Lattice position
        double x_; ///< x-coordinate
        double y_; ///< y-coordinate
        double z_; ///< z-coordinate
        int atype_; ///< Atom type
};
 
typedef struct line_struct line_struct;


/*!
 * \brief Structure for storing rate catalog data.
 * 
 * This structure contains configuration data, energy values, region categories, and region number.
 * It is used to store and access rate catalog data in a structured format.
 */
struct ratecatalog_struct {
    public:
        /*! \brief Constructor
        * 
        * Initializes the `ratecatalog_struct` with configurations, energies, region categories, and region number.
        * 
        * \param configs A 2D vector containing configuration data.
        * \param energies A 3D vector containing energy values.
        * \param regions_cat A 4D vector containing region categories.
        * \param region_num An integer representing the number of regions.
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

        /*! \return Returns the configuration data */
        std::vector< std::vector<int> > get_configs() const {
            return configs_;
        }

        /*! \return Returns the energy values */
        std::vector< std::vector< std::vector<double> > > get_energies() const {
            return energies_;
        }

        /*! \return Returns the region categories */
        std::vector< std::vector< std::vector< std::vector<double> > > > get_regions_cat() const {
            return regions_cat_;
        }

        /*! \return Returns the number of regions */
        int get_region_num() const {
            return region_num_;
        }
        
    private:
        std::vector< std::vector<int> > configs_; ///< Configuration data
        std::vector< std::vector< std::vector<double> > > energies_; ///< Energy values
        std::vector< std::vector< std::vector< std::vector<double> > > > regions_cat_; ///< Region categories
        int region_num_; ///< Number of regions
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


void write_to_file(std::string filename, Matrix<int>& values) {
    std::ofstream out_file;
    out_file.open(filename);
    std::string s;

    std::cout << "writing filename: " << filename << "\n";

    if (out_file.is_open()) {
        for (int i=0; i<(int)values.rows(); i++) {
            for (int j=0; j<(int)values.cols(); j++) {     
                out_file << values[i][j] << " ";
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


std::tuple<int, double> get_last_iter_time(std::string folder) {
    std::string cwd = std::filesystem::path( "." );
    std::string path = cwd + "/" + folder + "/vacs";
    std::string dir;
    int tick = 0;
    double time = 0;
    
    for (const auto& dirEntry : std::filesystem::recursive_directory_iterator(path)) {
        //std::cout << dirEntry << std::endl; 
        
        std::string sep = "/"; 
        std::vector<std::string> tokens = tokenizer(dirEntry.path().string(), sep);
        dir = tokens[(tokens.size()-1)];

        sep = "_"; 
        tokens = tokenizer(dir, sep);

        //std::cout << "tokens[3]: " << tokens[3] << "\n";
        //std::cout << "tokens[4]: " << tokens[4] << "\n";

        if (std::stoi(tokens[4]) > tick) { tick = std::stoi(tokens[4]); }
        if (std::stod(tokens[5]) > time) { time = std::stod(tokens[5]); }
    }
    
    std::tuple<int, double> tup(tick, time);
    std::cout << "tick: " << tick << "\n";
    std::cout << "time: " << time << "\n";

    return tup;
}


bool is_numeric_or_scinotation(const std::string& s)
{
    std::string::const_iterator it = s.begin();
    while (it != s.end() && (std::isdigit(*it) || (*it == 'e') || (*it == '.') || (*it == '-'))) ++it;
    return !s.empty() && it == s.end();
}