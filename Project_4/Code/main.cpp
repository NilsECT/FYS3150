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

  // Problem 4: compare with analytical for L = 2
  //std::vector<double> temperatures = {1.0, 2.0};

  // int N_MC_cycles = 500000;

  // model instance:
  Ising model = Ising();

  // Working code for comparing with analytical solution (currently missing magnetization, will fix soon):
  // std::vector<double> temps = {1.0, 2.4};
  // model.analytical_comparison(temps, 10000000, "analytical_comparison", seed);

  // // sigurd:
  double T_max = 2.5;
  double T_min = 2.1;

  std::vector<int> lattices = std::vector<int> {20, 40, 60, 80, 100};

  // // Problem 8:
  arma::vec temperatures_p8 = arma::linspace(T_min, T_max, 100);  // if we get results we need to upgrade this to more values
  model.phase_transitions(lattices, temperatures_p8, 10000, 3, seed, 5000);

  // trying to find burn in time based on cycles
  arma::vec temperatures = arma::vec {1., 2.4, 3.};

  // Sigurd spesial
  model.varying_n_mc_cycles(temperatures, 100000, 5, lattices, 0, "varying_cycles_no_burn_20");
  model.varying_n_mc_cycles(temperatures, 10000, 5, lattices, 5000, "varying_cycles_yes_burn_20");


  // looking for the distribution of the energies
  model.epsilon_dist(temperatures, lattices, 5000, 5);

  return 0;
}
