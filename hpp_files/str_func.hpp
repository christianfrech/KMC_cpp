#include <iostream>
#include <sstream>
#include <stdio.h>
#include <string>
#include <vector>
#include <cstdlib>
#include <fstream>


std::vector<std::string> tokenizer(std::string s, std::string delimiter) {
    size_t pos_start = 0, pos_end, delim_len = delimiter.length();
    std::string token;
    std::vector<std::string> res;
    //std::cout << "s: " << s << "\n";

    while ((pos_end = s.find(delimiter, pos_start)) != std::string::npos) {
        token = s.substr(pos_start, pos_end - pos_start);
        
        //std::cout << "pos_start: " << pos_start << " delimiter: " << delimiter << " pos_end: " << pos_end << "\n";
        //std::cout << "token: " << token << "\n";
        pos_start = pos_end + delim_len;
        res.push_back(token);
    }
    if (s.substr (pos_start, pos_end - pos_start) != "") {
        //std::cout << "s.substr(pos_start): " << s.substr(pos_start) << "\n";
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

std::tuple<std::string, double, double, double, int, int> parse_line_cluster(std::string line) {

    std::vector<std::string> toks = tokenizer(line, " ");
    std::string lattice_pos = (toks[0]);

    double x = std::stof(toks[1]);
    double y = std::stof(toks[2]);
    double z = std::stof(toks[3]);
    int atype = std::stoi(toks[4]);
    int cluster_num = std::stoi(toks[5]);

    std::tuple<std::string, double, double, double, int, int> tuple_out(lattice_pos,x,y,z,atype, cluster_num);
    return tuple_out;
}

std::tuple<std::string, double, double, double, int> parse_line(std::string line) {

    std::vector<std::string> toks = tokenizer(line, " ");
    std::string lattice_pos = (toks[0]);

    double x = std::stof(toks[1]);
    double y = std::stof(toks[2]);
    double z = std::stof(toks[3]);
    int atype = std::stoi(toks[4]);

    std::tuple<std::string, double, double, double, int> tuple_out(lattice_pos,x,y,z,atype);
    return tuple_out;
}

std::tuple<std::string, int, int, int> parse_reg_line(std::string line) {

    std::vector<std::string> toks = tokenizer(line, " ");
    std::string lattice_pos = (toks[0]);
    double x = std::stof(toks[1]);
    double y = std::stof(toks[2]);
    double z = std::stof(toks[3]);

    std::tuple<std::string, int, int, int> tuple_out(lattice_pos,x,y,z);
    return tuple_out;
}

void write_to_file(std::string filename, std::vector< std::vector<int> > values) {
    std::ofstream out_file;
    out_file.open(filename);
    std::string s;

    if (out_file.is_open()) {
        for (int i=0; i<(int)values.size(); i++) {
            for (int j=0; j<(int)values[i].size(); j++) { 
                //std::cout << values[i][j] << " ";
                out_file << values[i][j] << " ";
            }
            //std::cout << "\n";
            out_file << " \n";
        }
    }
    out_file.close();
}


template <typename T>
void write_to_file(std::string filename, std::vector< std::vector<T> > values, std::vector<T> NN_of_vacs) {
    std::ofstream out_file;
    out_file.open(filename);
    std::string s;

    if (out_file.is_open()) {
        for (int i=0; i<(int)values.size(); i++) {
            for (int j=0; j<(int)values[i].size(); j++) { 
                out_file << values[i][j] << " ";
            }
            
            out_file << NN_of_vacs[i] << " \n";
        }
    }
    out_file.close();
}

template <typename U>
void write_to_file(std::string filename, std::vector<U> values, std::vector<U> more_values) {
    std::ofstream out_file;
    out_file.open(filename);
    std::string s;

    if (out_file.is_open()) {
        for (int i=0; i<(int)values.size(); i++) {
            out_file << values[i] << " ";            
            out_file << more_values[i] << " \n";
        }
    }
    out_file.close();
}

std::tuple<int, double> get_last_iter_time(std::string folder) {
    std::string cwd = std::filesystem::current_path();
    std::string path = cwd + "/" + folder + "/vacs";
    std::string dir;
    int tick = 0;
    double time = 0;
    
    for (const auto& dirEntry : std::filesystem::recursive_directory_iterator(path)) {
        std::cout << dirEntry << std::endl; 
        
        std::string sep = "/"; 
        std::vector<std::string> tokens = tokenizer(dirEntry.path().string(), sep);
        dir = tokens[(tokens.size()-1)];
        std::cout << "dir: " << dir << "\n";

        sep = "_"; 
        tokens = tokenizer(dir, sep);

        std::cout << "tokens[3]: " << tokens[3] << "\n";
        std::cout << "tokens[4]: " << tokens[4] << "\n";

        if (std::stoi(tokens[3]) > tick) { tick = std::stoi(tokens[3]); }
        if (std::stod(tokens[4]) > time) { time = std::stod(tokens[4]); }
        std::cout << "tick: " << tick << "\n";
        std::cout << "time: " << time << "\n";
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