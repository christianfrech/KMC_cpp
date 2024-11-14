#include <iostream>
#include <sstream>
#include <stdio.h>
#include <string>
#include <vector>
#include <cstdlib>
#include <fstream>

#include "classes_strucs.hpp"

/*! \brief Structure for portability of information of a site read in from custom region input file . */
struct reg_line_struct {
public:
    /*! \brief Constructor
     * \param [in] lattice_type Type of lattice.
     * \param [in] x X-coordinate.
     * \param [in] y Y-coordinate.
     * \param [in] z Z-coordinate.
     */
    reg_line_struct(std::string lattice_type, int x, int y, int z)
        : lattice_type_(lattice_type), x_(x), y_(y), z_(z) {}

    /*! \brief Get the lattice type.
     * \return String representing the lattice type.
     */
    std::string get_latice_type() const { return lattice_type_; }

    /*! \brief Get the X coordinate.
     * \return Integer X coordinate.
     */
    int get_x() const { return x_; }

    /*! \brief Get the Y coordinate.
     * \return Integer Y coordinate.
     */
    int get_y() const { return y_; }

    /*! \brief Get the Z coordinate.
     * \return Integer Z coordinate.
     */
    int get_z() const { return z_; }

private:
    std::string lattice_type_; ///< Lattice type.
    int x_; ///< X coordinate.
    int y_; ///< Y coordinate.
    int z_; ///< Z coordinate.
};

typedef struct reg_line_struct reg_line_struct;

/*! \brief Structure for storing information of location of atom from geometry input file. */
struct line_struct {
public:
    /*! \brief Constructor
     * \param [in] lattice_pos Lattice position identifier.
     * \param [in] x X-coordinate.
     * \param [in] y Y-coordinate.
     * \param [in] z Z-coordinate.
     * \param [in] atype Atom type.
     */
    line_struct(std::string lattice_pos, double x, double y, double z, int atype)
        : lattice_pos_(lattice_pos), x_(x), y_(y), z_(z), atype_(atype) {}

    /*! \brief Get the lattice position identifier.
     * \return String representing the lattice position.
     */
    std::string get_latice_pos() const { return lattice_pos_; }

    /*! \brief Get the X coordinate.
     * \return Double X coordinate.
     */
    double get_x() const { return x_; }

    /*! \brief Get the Y coordinate.
     * \return Double Y coordinate.
     */
    double get_y() const { return y_; }

    /*! \brief Get the Z coordinate.
     * \return Double Z coordinate.
     */
    double get_z() const { return z_; }

    /*! \brief Get the atom type.
     * \return Integer representing the atom type.
     */
    int get_atype() const { return atype_; }

private:
    std::string lattice_pos_; ///< Lattice position.
    double x_; ///< X coordinate.
    double y_; ///< Y coordinate.
    double z_; ///< Z coordinate.
    int atype_; ///< Atom type.
};

typedef struct line_struct line_struct;

/*! \brief Structure for cataloging rate data for regions. */
struct ratecatalog_struct {
public:
    /*! \brief Constructor
     * \param [in] configs Vector of binary strings representing configurations.
     * \param [in] energies Vector of energy data.
     * \param [in] regions_cat Vector of region catalog data.
     * \param [in] region_num Number of regions.
     */
    ratecatalog_struct(std::vector<std::vector<int>> configs,
                       std::vector<std::vector<std::vector<double>>> energies,
                       std::vector<std::vector<std::vector<std::vector<double>>>> regions_cat,
                       int region_num)
        : configs_(configs), energies_(energies), regions_cat_(regions_cat), region_num_(region_num) {}

    /*! \brief Get the vector of binary strings representing configurations.
     * \return Vector of configuration data.
     */
    std::vector<std::vector<int>> get_configs() const { return configs_; }

    /*! \brief Get the vector of energies corresponding to each NN-configuration.
     * \return Vector of energy data.
     */
    std::vector<std::vector<std::vector<double>>> get_energies() const { return energies_; }

    /*! \brief Get the rate catalog corresponding to regions.
     * \return Vector of rate catalogs (tripley-nested vectors) corresponding to regions.
     */
    std::vector<std::vector<std::vector<std::vector<double>>>> get_regions_cat() const { return regions_cat_; }

    /*! \brief Get the number of regions.
     * \return Integer number of regions.
     */
    int get_region_num() const { return region_num_; }

private:
    std::vector<std::vector<int>> configs_; ///< Vector of configurations.
    std::vector<std::vector<std::vector<double>>> energies_; ///< Vector of energy data.
    std::vector<std::vector<std::vector<std::vector<double>>>> regions_cat_; ///< Region catalog data.
    int region_num_; ///< Number of regions.
};

typedef struct ratecatalog_struct ratecatalog_struct;

/*! \brief Tokenizes a string based on a given delimiter.
 * \param [in] s The string to tokenize.
 * \param [in] delimiter The delimiter string.
 * \return Vector of tokens.
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

/*! \brief Splits a string into a vector of characters.
 * \param [in] s The string to split.
 * \return Vector of characters.
 */
std::vector<char> split_by_char(std::string s) {
    std::vector<char> v(s.begin(), s.end());
    return v;
}

/*! \brief Checks if a string contains only alphabetic characters.
 * \param [in] str The string to check.
 * \return True if the string is alphabetic, false otherwise.
 */
bool str_isalpha(std::string str) {
    auto it = std::find_if(str.begin(), str.end(), [](char const &c) {
        return !std::isalpha(c);
    });
    return it == str.end();
}

/*! \brief Parses a line to create a line_struct object.
 * \param [in] line The input line to parse.
 * \return A line_struct object populated with parsed data.
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

/*! \brief Parses a line from custom region input file to create a reg_line_struct object.
 * \param [in] line The input line to parse.
 * \return A reg_line_struct object populated with parsed data.
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

/*! \brief Writes matrix data to a file.
 * \param [in] filename The name of the output file.
 * \param [in] values Matrix data to write to the file.
 */
void write_to_file(std::string filename, Matrix<int> values) {
    std::ofstream out_file;
    out_file.open(filename);

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

/*! \brief Writes vector data to a file.
 * \param [in] filename The name of the output file.
 * \param [in] values Vector of integer vectors to write to the file.
 */
void write_to_file(std::string filename, std::vector<std::vector<int>> values) {
    std::ofstream out_file;
    out_file.open(filename);

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
