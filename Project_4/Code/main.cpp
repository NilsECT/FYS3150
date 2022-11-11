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

  ofile << "Thread seed, Temperature, Epsilon, Epsilon squared, Abs(magnetization), Magnetization squared" << std::endl;
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
            << std::setprecision(print_prec) << G.m_squared << std::endl;

    }
  }

  return 0;
}
