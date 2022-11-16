#include <fstream>
#include "omp.h"
#include "Grid.hpp"
#include <armadillo>
#include <iostream>
#include <string>
#include <random>
#include <iomanip>      // std::setprecision
#include "Ising.hpp"

int main(){

  const int seed = 137;
  std::mt19937 MC_seed_generator (seed);

  // Problem 4: compare with analytical for L = 2
  //std::vector<double> temperatures = {1.0, 2.0};

  int N_MC_cycles = 1000;
  int N_spinflips = 1000;

  // model instance:
  Ising model = Ising();

  // std::cout << "Starting analytical comparison" << std::endl;
  // model.analytical_comparison(temperatures, N_spinflips, N_MC_cycles);


  // // Problem 5: study burn-in time
  arma::vec cycles = arma::logspace(1, 7, 50);
  arma::vec temperatures_p5 = arma::vec {1., 1.5, 2., 2.4, 3.};
  model.varying_n_mc_cycles(temperatures_p5, cycles);

  std::vector<double> temperatures_2 = {1., 1.5, 2., 2.4, 3.};
  model.epsilon_dist(temperatures_2, 20, 100000, 10);

  // // Problem 8:
  std::vector<int> lattice_sizes_p8 = {40, 60, 80, 100};
  arma::vec temperatures_p8 = arma::linspace(2.1, 2.4, 4);  // if we get results we need to upgrade this to more values
  model.phase_transitions(lattice_sizes_p8, temperatures_p8, N_spinflips, N_MC_cycles, seed);

  return 0;
}
