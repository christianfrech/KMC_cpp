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

void write_to_file(std::string filename, std::vector< std::vector<int> > values) {
    std::ofstream out_file;
    out_file.open(filename);
    std::string s;

    //std::cout << "writing \n";

    if (out_file.is_open()) {
        //std::cout << "file open \n";
        for (int i=0; i<(int)values.size(); i++) {
            for (int j=0; j<(int)values[i].size(); j++) {         
                //s = std::to_string(values[i][j]);
                out_file << values[i][j] << " ";
            }
        out_file << "\n";
        }
    }

    out_file.close();
    //std::cout << "closing \n";
}