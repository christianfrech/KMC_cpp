//#include "kmc_onecheckactions.hpp"
#include "kmc_stripping.hpp"
#include <sstream>
#include <iostream>
#include <chrono>
#include <ctime>    


int main() {
    std::ostringstream ss;
    std::string output_filename;
    std::string alltimes_filename;
    std::string times_filename;
    std::string count_filename;
    std::vector< std::string > in_file = {"singlevoid_interface_large.txt"};
    std::vector< std::string > catalog_file = {"ratecatalog_void.txt"};
    std::string region_infile;
    std::tuple< std::vector< std::vector< std::vector<int> > >, std::vector<int>, std::vector<double>, std::vector<double> > return_tuple;
    std::vector< std::vector< std::vector<int> > > vacancies;
    std::vector<double> all_times;
    std::vector<int> counts;
    std::vector<double> times;
    Lattice* lattice;

    printf("running\n\n");
    int iterations = 1;
    double time = 5e9;
    std::vector<double> temperatures = {300,500,700,900};
    std::vector<double> vertex_rates = {7.65954e11};
    std::vector<double> edge_rates = {1.23537e6};
    std::vector< std::vector< std::vector<double> > > reg_rates = {{}}; 
    // settings of kb = 8.6173e-5, T = 300, E = 0.05, prefactor = 5e12
    // interface barrier = 8.28403e10   interface rate (uncorrected) = 7.98587e8
    std::ofstream out_file;
    std::chrono::system_clock::time_point start;
    std::vector<std::string> folders = {"singlevoid_large"};
    

    for (int h=0; h<folders.size(); h++) {
        for (int i=0; i<iterations; i++) {
            for (int k=0; k<(int)edge_rates.size(); k++) {
                std::cout << "-----------------------------------------\n";
                std::cout << "i: " << i << "\n"; 
                std::cout << "temp: " << temperatures[k] << "\n"; 
                lattice = populate_lattice(in_file[h], catalog_file[h], region_infile, vertex_rates[k], edge_rates[k], reg_rates[h]);
                std::cout << "starting run \n";
                start = std::chrono::system_clock::now();
                return_tuple = lattice->new_kmc_iterator(time, start, folders[h], i);
                vacancies = std::get<0>(return_tuple); counts = std::get<1>(return_tuple); times = std::get<2>(return_tuple); all_times = std::get<3>(return_tuple);

                print_1Dvector_float(all_times);
                
                for (int l=0; l<vacancies.size(); l++) {
                    ss << folders[h] << "/vacs/vacancies_output_" << (i) << "_" << l << "_" << k << "rate.txt";
                    output_filename = ss.str();
                    write_to_file(output_filename, vacancies[l]);
                    ss.str("");
                    ss.clear();
                }

                for (int l=0; l<all_times.size(); l++) {
                    ss << folders[h] << "/all_times/alltimes_output_" << (i) << "_" << l << "_" << k << "rate.txt";
                    std::cout << ss.str() << "\n";
                    alltimes_filename = ss.str();
                    out_file.open(alltimes_filename);
                    out_file << all_times[l] << "\n";
                    ss.str("");
                    ss.clear();
                    out_file.close();
                } 

                for (int l=0; l<(int)counts.size(); l++) {
                    ss << folders[h] << "/counts/counts_output_" << (i) << "_" << l << "_" << k << "rate.txt";
                    count_filename = ss.str();
                    out_file.open(count_filename);
                    out_file << counts[l] << "\n";
                    ss.str("");
                    ss.clear();
                    out_file.close();

                    ss << folders[h] << "/times/times_output_" << (i) << "_" << l << "_" << k << "rate.txt";
                    times_filename = ss.str();
                    out_file.open(times_filename);
                    out_file << times[l] << "\n";
                    ss.str("");
                    ss.clear();
                    out_file.close();
                }

                delete lattice;
            }
        }
    }

    return 0;
}