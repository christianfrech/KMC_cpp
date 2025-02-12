//#include "kmc_onecheckactions.hpp"
#include "kmc_stripping_with_NNofNN.hpp"
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
    std::vector< std::string > in_file = {"singlevoid_interface_gb_1k_nonrandom_fullmodel_corrected.txt"};
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
    std::vector<double> edge_rates = {1.23537e6};
    
    // gb 0.075eV = 2.74802e11 s^-1, 0.1eV = 1.04481e11 s^-1, 0.125eV = 3.97241e10 s^-1, 0.15eV = 1.51033e10 s^-1, 
    // gb 0.175eV = 5.74233e9 s^-1, 0.2eV = 2.18326e9 s^-1, 0.225eV = 8.30083e8 s^-1, 0.25eV =  3.15601e8 s^-1
    // settings of kb = 8.6173e-5, T = 300, E = 0.05, prefactor = 5e12
    // vacuum barrier (110) bulk vs surface = 0.325   vacuum rate (110) bulk vs surface = 1.734559609419e7
    // interface barrier = 0.52031605eV   interface averaged rate = 9078.442295
    // gb distribution: min 
    
    std::ofstream out_file;
    std::chrono::system_clock::time_point start;
    std::vector<std::string> folders = {"test_folder"};
    
    bool restart = false;
    int last_move_tick;
    double last_time;
    std::tuple<int, double> tup_out;


    for (int h=0; h<folders.size(); h++) {
        for (int i=0; i<iterations; i++) {
            for (int k=0; k<(int)edge_rates.size(); k++) {
                std::cout << "getting times \n";
                tup_out = get_last_iter_time(folders[h]);
                last_move_tick = std::get<0>(tup_out);
                last_time = std::get<1>(tup_out);
                std::cout << "last_time: " << last_time << "\n";
                std::cout << "last_move_tick: " << last_move_tick << "\n";
                if (last_move_tick!=0) restart = true;

                std::cout << "-----------------------------------------\n";
                std::cout << "i: " << i << "\n"; 
                
                lattice = populate_lattice(in_file[h], catalog_file[h], region_infile);
                std::cout << "starting run \n";
                start = std::chrono::system_clock::now();
                return_tuple = lattice->new_kmc_iterator(time, start, folders[h], i, last_move_tick, last_time);
                vacancies = std::get<0>(return_tuple); counts = std::get<1>(return_tuple); times = std::get<2>(return_tuple); all_times = std::get<3>(return_tuple);

                
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