#include <fstream>
#include "omp.h"
#include "Grid.hpp"
#include <armadillo>
#include <iostream>
#include <string>
#include <random>

void generate_histograms(std::vector<double> temperatures, int num_MC_cycles,
                         int seed = 137,
                         int L = 20, std::string filename = "histogram"){

  // Unordered initial state:
  // Open file to which we will write the data:
  std::ofstream file;
  file.open(filename + "_unordered.txt", std::ofstream::trunc);
  static int print_prec = 10;

  file << std::setprecision(print_prec) << "Seed, "
       << std::setprecision(print_prec) << "Temperature, "
       << std::setprecision(print_prec) << "Energy, "   // Energy per spin [1/J]
       << std::setprecision(print_prec) << "Frequency"  // Frequency, normalized
       << std::endl;

  for (double T : temperatures){

    Grid model = Grid(L, T);
    model.fill_grid(seed);
    model.MCMC(num_MC_cycles, seed);

    int bin_width = model.bins.n_elem / (model.N*2);

    for (int i = 0; i < model.bins.n_elem; i ++){

      file << std::setprecision(print_prec) << seed << ", "
           << std::setprecision(print_prec) << T << ", "
           << std::setprecision(print_prec) << (float)(-model.N + i*bin_width) / (float)(model.N) << ", "
           << std::setprecision(print_prec) << (float)(model.bins.at(i)) / (float) num_MC_cycles << std::endl;

    }
  }

  // Ordered initial state:
  // Open file to which we will write the data:
  std::ofstream ordered_file;
  ordered_file.open(filename + "_ordered.txt", std::ofstream::trunc);

  ordered_file << std::setprecision(print_prec) << "Seed, "
       << std::setprecision(print_prec) << "Temperature, "
       << std::setprecision(print_prec) << "Energy, "   // Energy per spin [1/J]
       << std::setprecision(print_prec) << "Frequency"  // Frequency, normalized
       << std::endl;

  for (double T : temperatures){

    Grid model = Grid(L, T);
    model.fill_grid(seed, false);
    model.MCMC(num_MC_cycles, seed);

    int bin_width = model.bins.n_elem / (model.N*2);

    for (int i = 0; i < model.bins.n_elem; i ++){

      ordered_file << std::setprecision(print_prec) << seed << ", "
           << std::setprecision(print_prec) << T << ", "
           << std::setprecision(print_prec) << (float)(-model.N + i*bin_width) / (float)(model.N) << ", "
           << std::setprecision(print_prec) << (float)(model.bins.at(i)) / (float) num_MC_cycles << std::endl;

    }
  }
}

void analytical_comparison(std::vector<double> T, std::string filename = "analytical", int seed = 137) {

}

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
}

void var_cycles(int num_MC_cycles, double T, bool random_config, int L = 20, std::string filename = "var_cycles", int seed = 137) {

  //
  // Something is off in this function, it is not finished yet
  //

  int sampling_distance = 100;

  // start by opening the file which will gather all the data
  std::ofstream file;
  file.open(filename + ".txt", std::ofstream::trunc);
  static int print_prec = 10;

  // header of file:
  file << "Thread seed, Number of MCMC, Epsilon, Abs(magnetization)" << std::endl;

  for (int i = 1; i < num_MC_cycles; i += sampling_distance) {

    Grid model = Grid(L, T); // generate model (create)
    model.fill_grid(seed, random_config);  // generate model (fill)
    model.MCMC(i, seed); // simulate

    std::cout << "i = " << i << "    <eps> = " << model.epsilon << std::endl;

    file << std::setprecision(print_prec) << seed << ", "
         << std::setprecision(print_prec) << i << ", "
         << std::setprecision(print_prec) << model.epsilon << ", "
         << std::setprecision(print_prec) << model.m_abs << std::endl;
  }
}

int main(int argc, char* argv[]){

  const int L = atoi(argv[1]);
  const int num_threads = atoi(argv[2]);
  const int num_MC_cycles = atoi(argv[3]);


  const int seed = 137;
  std::mt19937 MC_seed_generator (seed);

  // Problem 6:
  std::vector<double> temperatures = {1, 2.4};
  generate_histograms(temperatures, 2000000);
  return 0;


  // Code for testing whether code makes sense for different temperatures:

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
  //var_LT(num_MC_cycles, L_p8, T_p8);

  arma::vec cycles_p5 = {10, 100, 1000, 10000, 100000};
  //var_cycles(cycles_p5);

  return 0;
}
