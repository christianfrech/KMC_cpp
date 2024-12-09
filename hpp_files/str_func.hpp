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
/*! 
 * \brief Tokenizes a string based on a specified delimiter.
 * 
 * This function splits a string into substrings based on a specified delimiter.
 * It returns a vector of strings containing the substrings.
 * 
 * \param s The string to tokenize.
 * \param delimiter The delimiter used to split the string.
 * 
 * \return A vector of strings, where each element is a substring from the original string.
 */
std::vector<std::string> tokenizer(std::string s, std::string delimiter) {
    size_t pos_start = 0, pos_end, delim_len = delimiter.length();
    std::string token;
    std::vector<std::string> res;

    while ((pos_end = s.find(delimiter, pos_start)) != std::string::npos) {
        token = s.substr(pos_start, pos_end - pos_start);
        pos_start = pos_end + delim_len;
        res.push_back(token);
    }
    if (s.substr(pos_start, pos_end - pos_start) != "") {
        res.push_back(s.substr(pos_start));
    }
    
    return res;
}

/*! 
 * \brief Splits a string into a vector of characters.
 * 
 * This function takes a string and returns a vector where each element is a character from the string.
 * 
 * \param s The string to split.
 * 
 * \return A vector of characters from the string.
 */
std::vector<char> split_by_char(std::string s) {
    std::vector<char> v(s.begin(), s.end());
    return v;
}

/*! 
 * \brief Checks if a string contains only alphabetic characters.
 * 
 * This function checks if all characters in the input string are alphabetic (a-z, A-Z).
 * It returns `true` if the string contains only alphabetic characters, and `false` otherwise.
 * 
 * \param str The string to check.
 * 
 * \return `true` if the string contains only alphabetic characters, `false` otherwise.
 */
bool str_isalpha(std::string str) {
    auto it = std::find_if(str.begin(), str.end(), [](char const &c) {
        return !std::isalpha(c);
    });
    return it == str.end();
}

/*! 
 * \brief Parses a line into a `line_struct` object.
 * 
 * This function takes a line in string format, splits it into components, and constructs a `line_struct` object
 * with the parsed data. It expects the line to have at least five components: lattice position, x, y, z, and atom type.
 * 
 * \param line The string to parse.
 * 
 * \return A `line_struct` object constructed from the parsed data.
 */
line_struct parse_line(std::string line) {

    std::vector<std::string> toks = tokenizer(line, " ");
    std::string lattice_pos = (toks[0]);

    double x = std::stof(toks[1]);
    double y = std::stof(toks[2]);
    double z = std::stof(toks[3]);
    int atype = std::stoi(toks[4]);

    line_struct output_vals(lattice_pos, x, y, z, atype);
    return output_vals;
}

/*! 
 * \brief Parses a line into a `reg_line_struct` object.
 * 
 * This function takes a line in string format, splits it into components, and constructs a `reg_line_struct` object
 * with the parsed data. It expects the line to have at least four components: lattice position, x, y, and z.
 * 
 * \param line The string to parse.
 * 
 * \return A `reg_line_struct` object constructed from the parsed data.
 */
reg_line_struct parse_reg_line(std::string line) {

    std::vector<std::string> toks = tokenizer(line, " ");
    std::string lattice_pos = (toks[0]);
    double x = std::stof(toks[1]);
    double y = std::stof(toks[2]);
    double z = std::stof(toks[3]);

    reg_line_struct output_vals(lattice_pos, x, y, z);
    return output_vals;
}

/*! 
 * \brief Writes matrix data to a file.
 * 
 * This function writes the data from a `Matrix<int>` object to a file. It writes the contents of the matrix row by row.
 * Each element in a row is separated by a space, and each row is written on a new line in the file.
 * 
 * \param filename The name of the file to write to.
 * \param values The matrix of integer values to write.
 */
void write_to_file(std::string filename, Matrix<int>& values) {
    std::ofstream out_file;
    out_file.open(filename);
    std::string s;

    std::cout << "writing filename: " << filename << "\n";

    if (out_file.is_open()) {
        for (int i = 0; i < (int)values.rows(); i++) {
            for (int j = 0; j < (int)values.cols(); j++) {     
                out_file << values[i][j] << " ";
            }
            out_file << "\n";  
        }
    }

    out_file.close();
}

/*! 
 * \brief Writes a 2D vector of integers to a file.
 * 
 * This function writes the data from a 2D vector of integers to a file. Each element in a row is separated by a space,
 * and each row is written on a new line in the file.
 * 
 * \param filename The name of the file to write to.
 * \param values The 2D vector of integer values to write.
 */
void write_to_file(std::string filename, std::vector< std::vector<int> > values) {
    std::ofstream out_file;
    out_file.open(filename);
    std::string s;

    if (out_file.is_open()) {
        for (int i = 0; i < (int)values.size(); i++) {
            for (int j = 0; j < (int)values[i].size(); j++) { 
                out_file << values[i][j] << " ";
            }
            out_file << "\n";
        }
    }

    out_file.close();
}

/*! 
 * \brief Gets the last iteration's time and tick from a folder.
 * 
 * This function iterates through the files in a given folder, parses the filenames, and extracts the highest tick and 
 * corresponding time values from the filenames. The filenames are expected to contain tick and time information
 * following the format "fillertext_fillertext_trajidx_rateidx_movetick_time_moves.txt".
 * 
 * \param folder The folder to search for files.
 * 
 * \return A tuple containing the highest tick and the corresponding time.
 */
std::tuple<int, double> get_last_iter_time(std::string folder) {
    std::string cwd = std::filesystem::path(".");
    std::string path = cwd + "/" + folder + "/vacs";
    std::string dir;
    int tick = 0;
    double time = 0;
    
    for (const auto& dirEntry : std::filesystem::recursive_directory_iterator(path)) {
        std::string sep = "/"; 
        std::vector<std::string> tokens = tokenizer(dirEntry.path().string(), sep);
        dir = tokens[(tokens.size()-1)];

        sep = "_"; 
        tokens = tokenizer(dir, sep);

        if (std::stoi(tokens[4]) > tick) { tick = std::stoi(tokens[4]); }
        if (std::stod(tokens[5]) > time) { time = std::stod(tokens[5]); }
    }
    
    std::tuple<int, double> tup(tick, time);
    std::cout << "tick: " << tick << "\n";
    std::cout << "time: " << time << "\n";

    return tup;
}

