#include <fstream>
#include "omp.h"
#include "Grid.hpp"
#include <armadillo>
#include <iostream>
#include <string>
#include <random>
#include <iomanip>      // std::setprecision

std::vector<double> analytical(int L, double T){
  double analytical_Z = 0;
  double mean_E = 0;
  double mean_E_squared = 0;

  std::vector<double> E = {-8, 0, 0, 8, 0, -8};

  for (int i = 0; i < 6; i++){
    analytical_Z += exp(-E.at(i) / T);

    mean_E += E.at(i) * exp(-E.at(i) / T);
    mean_E_squared += E.at(i) * E.at(i) * exp(-E.at(i) / T);
  }

  mean_E /= analytical_Z;
  mean_E_squared /= analytical_Z;

  std::vector<double> avgs = {mean_E/(L*L), mean_E_squared/(L*L*L*L)};
  return avgs;

}

std::vector<double> monte_carlo(int N_spinflips, int num_MC_cycles, double T, int seed, int L = 2, bool random_config = true) {

  Grid model = Grid(L, T);

  double avg_epsilon = 0.0;
  double avg_epsilon_sq = 0.0;
  double avg_m_abs = 0.0;
  double avg_m_sq = 0.0;

  for (int i = 0; i < num_MC_cycles; i ++){

    model.fill_grid(seed + i, random_config);
    model.random_walk(N_spinflips, seed + i);

    avg_epsilon = avg_epsilon + model.epsilon;
    avg_epsilon_sq = avg_epsilon_sq + model.epsilon_squared;

    avg_m_abs = avg_m_abs + model.m_abs;
    avg_m_sq = avg_m_sq + model.m_squared;

  }

  avg_epsilon /= num_MC_cycles;
  avg_epsilon_sq /= num_MC_cycles;
  avg_m_abs /= num_MC_cycles;
  avg_m_sq /= num_MC_cycles;

  std::vector<double> avgs = {avg_epsilon, avg_epsilon_sq, avg_m_abs, avg_m_sq};
  return avgs;
}

void generate_histograms(std::vector<double> temperatures, int num_MC_cycles, int seed = 137, int L = 20, std::string filename = "histogram"){

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
    model.random_walk(num_MC_cycles, seed);

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
    model.random_walk(num_MC_cycles, seed);

    int bin_width = model.bins.n_elem / (model.N*2);

    for (int i = 0; i < model.bins.n_elem; i ++){

      ordered_file << std::setprecision(print_prec) << seed << ", "
           << std::setprecision(print_prec) << T << ", "
           << std::setprecision(print_prec) << (float)(-model.N + i*bin_width) / (float)(model.N) << ", "
           << std::setprecision(print_prec) << (float)(model.bins.at(i)) / (float) num_MC_cycles << std::endl;

    }
  }
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
      model.random_walk(cycles, seed); // simulate

      // the extraction has to happen based on the fixed L and T
      #pragma omp critical
      file << std::setprecision(print_prec) << seed << ", " << std::setprecision(print_prec) << ii << ", " << std::setprecision(print_prec) << model.epsilon << ", " << std::setprecision(print_prec) << model.epsilon_squared << ", " << std::setprecision(print_prec) << model.m_abs << ", " << std::setprecision(print_prec) << model.m_squared << std::endl;
    }
  }
  file.close();
}

void analytical_comparison(std::vector<double> temperatures, int N_spinflips, int N_MC_cycles, std::string filename = "analytical_comparison", const int seed = 137){

  std::mt19937 MC_seed_generator (seed);

  std::ofstream file;
  file.open(filename + "_" + std::to_string(N_spinflips) + "_" + std::to_string(N_MC_cycles) + ".txt", std::ofstream::trunc);
  static int print_prec = 10;

  file << std::setprecision(print_prec) << "Seed, "
       << std::setprecision(print_prec) << "Temperature, "
       << std::setprecision(print_prec) << "Epsilon (abs), "
       << std::setprecision(print_prec) << "Epsilon (squared), "
       << std::setprecision(print_prec) << "M (abs), "
       << std::setprecision(print_prec) << "M (squared), "
       << std::setprecision(print_prec) << "Analytical epsilon (abs), "
       << std::setprecision(print_prec) << "Analytical epsilon (squared)"
       << std::endl;

  #pragma omp parallel
  {
    #pragma omp for
    for (double &temperature : temperatures){

      auto thread_seed = MC_seed_generator();

      std::vector<double> res = monte_carlo(N_spinflips, N_MC_cycles, temperature, thread_seed);

      std::vector<double> avgs = analytical(2, temperature);

      file << std::setprecision(print_prec) << temperature << ", "
           << std::setprecision(print_prec) << res.at(0) << ", "
           << std::setprecision(print_prec) << res.at(1) << ", "
           << std::setprecision(print_prec) << res.at(2) << ", "
           << std::setprecision(print_prec) << res.at(3) << ", "
           << std::setprecision(print_prec) << avgs.at(0) << ", "
           << std::setprecision(print_prec) << avgs.at(1) << std::endl;
      }

  }
  file.close();
}

void varying_n_mc_cycles(double temperature, arma::vec n_cycles, std::string filename = "varying_cycles", int seed = 137, bool random_config = true, int L = 20, int N_spinflips = 100){

  std::mt19937 MC_seed_generator (seed);

  std::ofstream file;
  file.open(filename + "_" + std::to_string(temperature) + ".txt", std::ofstream::trunc);
  static int print_prec = 10;

  file << std::setprecision(print_prec) << "Seed, "
       << std::setprecision(print_prec) << "Temperature, "
       << std::setprecision(print_prec) << "MC cycles, "
       << std::setprecision(print_prec) << "Epsilon (abs), "
       << std::setprecision(print_prec) << "M (abs), "
       << std::endl;

  Grid model = Grid(L, temperature);

  #pragma omp parallel
  {
    #pragma omp for
    for (double &n : n_cycles){

      auto thread_seed = MC_seed_generator();
      model.fill_grid(thread_seed, random_config);

      std::vector<double> res = monte_carlo(N_spinflips, n, temperature, thread_seed, 2, random_config);

      file << std::setprecision(print_prec) << thread_seed << ", "
           << std::setprecision(print_prec) << temperature << ", "
           << std::setprecision(print_prec) << n << ", "
           << std::setprecision(print_prec) << res.at(0) << ", "
           << std::setprecision(print_prec) << res.at(2)
           << std::endl;
      }

  }
}
/*
void phase_transitions(std::vector<int> L, arma::vec temperatures, int N_spinflips, int N_MC_cycles, int seed){

  std::mt19937 MC_seed_generator (seed);

  std::ofstream file;
  file.open(filename + "_" + ".txt", std::ofstream::trunc);
  static int print_prec = 10;

  file << std::setprecision(print_prec) << "Seed, "
       << std::setprecision(print_prec) << "Temperature, "
       << std::setprecision(print_prec) << "L, "
       << std::setprecision(print_prec) << "Epsilon (abs), "
       << std::setprecision(print_prec) << "M (abs), "
       << std::endl;


  #pragram omp parallel
  {
    #pragma omp for
    for (double &temperature : temperatures){
      std::vector<double> res = monte_carlo(N_spinflips, N_MC_cycles, temperature, thread_seed, L, true);

      Cv = (res.at(1) - res.at(0) * res.at(0))

      file << std::setprecision(print_prec) << thread_seed << ", "
           << std::setprecision(print_prec) << temperature << ", "
           << std::setprecision(print_prec) << L << ", "
           << std::setprecision(print_prec) << res.at(0) << ", "
           << std::setprecision(print_prec) << res.at(2)
           << std::endl;
    }
  }
}
*/

int main(int argc, char* argv[]){


  const int seed = 137;
  std::mt19937 MC_seed_generator (seed);

  // Problem 4: compare with analytical for L = 2
  std::vector<double> temperatures = {1.0, 2.0};

  //int N_MC_cycles = 1000000;
  //int N_spinflips = 1000;

  int N_MC_cycles = 100;
  int N_spinflips = 1000;

  analytical_comparison(temperatures, N_spinflips, N_MC_cycles);


  // Problem 5: study burn-in time

  arma::vec cycles = arma::logspace(1, 7, 7);
  varying_n_mc_cycles(1.0, cycles, "varying_cycles_ordered", MC_seed_generator(), true);
  varying_n_mc_cycles(1.0, cycles, "varying_cycles_unordered", MC_seed_generator(), false);

  varying_n_mc_cycles(2.4, cycles, "varying_cycles_ordered", MC_seed_generator(), true);
  varying_n_mc_cycles(2.4, cycles, "varying_cycles_unordered", MC_seed_generator(), false);



  /*

  // Problem 6:
  std::vector<double> temperatures = {1, 2.4};
  generate_histograms(temperatures, 1000000);
  return 0;


  */


  // results for P8
  // scan area of L and T
  std::vector<int> L_p8 = {40, 60, 80, 100};
  arma::vec T_p8 = arma::linspace(2.1, 2.4, 10);
  //var_LT(num_MC_cycles, L_p8, T_p8);

  //phase_transitions(L_p8, T_p8, N_spinflips, N_MC_cycles);

  return 0;
}
