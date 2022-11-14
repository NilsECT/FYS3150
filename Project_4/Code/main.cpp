#include <fstream>
#include "omp.h"
#include "Grid.hpp"
#include <armadillo>
#include <iostream>
#include <string>
#include <random>
/*
void write_to_file(int thread_seed, int L, double T, std::vector<double> thread_avg){
  std::ofstream ofile;
  st
  static int print_prec = 10;

}
*/

void var_LT(int cycles, std::vector<int> L, arma::vec T, std::string filename = "var_LT", int seed = 137) {

  // start by opening the file which will gather all the data
  std::ofstream file;
  file.open(filename + ".txt", std::ofstream::trunc);
  static int print_prec = 10;

  // header of file:
  file << "Thread seed, Temperature, Epsilon, Epsilon squared, Abs(magnetization), Magnetization squared" << std::endl;

  // start looping through all L and T in parallel
  #pragma omp parallel for
  for (int i : L) {
    for (double ii : T) {

      Grid model = Grid(i, ii); // generate model (create)
      model.fill_grid(seed);  // generate model (fill)
      model.MCMC(cycles, seed); // simulate

      // the extraction has to happen based on the fixed L and T
      #pragma omp critical
      file << std::setprecision(print_prec) << seed << ", " << std::setprecision(print_prec) << ii << ", " << std::setprecision(print_prec) << model.epsilon << ", " << std::setprecision(print_prec) << model.epsilon_squared << ", " << std::setprecision(print_prec) << model.m_abs << ", " << std::setprecision(print_prec) << model.m_squared << std::endl;
    }
  }
  file.close()
}

void var_cycles(std::vector<int> num_MC_cycles, double T, bool random_config, int L = 20, std::string filename = "var_cycles", int seed = 137) {

  // start by opening the file which will gather all the data
  std::ofstream file;
  file.open(filename + ".txt", std::ofstream::trunc);
  static int print_prec = 10;

  // header of file:
  file << "Thread seed, Temperature, Epsilon, Epsilon squared, Abs(magnetization), Magnetization squared" << std::endl;

  for (int cycles : num_MC_cycles) {

    Grid model = Grid(L, T); // generate model (create)
    model.fill_grid(seed, random_config);  // generate model (fill)
    model.MCMC(cycles, seed); // simulate

    file << std::setprecision(print_prec) << seed << ", " << std::setprecision(print_prec) << cycles << ", " << std::setprecision(print_prec) << model.epsilon << ", " << std::setprecision(print_prec) << model.m_abs << std::endl;
  }
  file.close()
}

int main(int argc, char* argv[]){

  const int L = atoi(argv[1]);
  const int num_threads = atoi(argv[2]);
  const int num_MC_cycles = atoi(argv[3]);


  const int seed = 137;
  std::mt19937 MC_seed_generator (seed);

  double T_min = 1.6;
  double T_max = 2.6;
  double dT = (T_max - T_min) / num_threads;

  std::ofstream ofile;
  static int print_prec = 10;
  //std::string filename = "averages_" + L + "_" + num_MC_cycles + "_" + num_threads + ".txt";
  ofile.open("averages_" + std::to_string(L) + "_" + std::to_string(num_MC_cycles) + "_" + std::to_string(num_threads) + ".txt", std::ofstream::trunc);

  ofile << "Thread seed, Temperature, |eps|, eps squared, |m|, m squared, C_v, Chi" << std::endl;
  #pragma omp parallel
  {

    # pragma omp for
    for (int i = 0; i < num_threads; i++){
      double T = T_min + i * dT;

      Grid G = Grid(L, T);

      auto thread_seed = MC_seed_generator();

      G.fill_grid(thread_seed);
      G.MCMC(num_MC_cycles, thread_seed);

      ofile << std::setprecision(print_prec) << thread_seed << ", "
            << std::setprecision(print_prec) << T << ", "
            << std::setprecision(print_prec) << G.epsilon << ", "
            << std::setprecision(print_prec) << G.epsilon_squared << ", "
            << std::setprecision(print_prec) << G.m_abs << ", "
            << std::setprecision(print_prec) << G.m_squared << ", "
            << std::setprecision(print_prec) << G.cv << ", "
            << std::setprecision(print_prec) << G.chi << std::endl;
    }
  }
  

  // results for P8
  // scan area of L and T
  std::vector<int> L_p8 = {40, 60, 80, 100};
  arma::vec T_p8 = arma::linspace(2.1, 2.4, 100);
  // here too we can run each of these about 10 times to get a confidence interval
  //var_LT(num_MC_cycles, L_p8, T_p8);

  arma::vec cycles_p5 = {10, 100, 1000, 10000, 100000}; // make logspace(1, 6)? and do so a few, say 10 times, to get a confidence interval?
  // all of these runs should then be in the same file
  //var_cycles(cycles_p5);

  return 0;
}
