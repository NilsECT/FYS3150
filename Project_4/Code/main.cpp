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

  // model instance:
  Ising model = Ising();

  // Compare with analytical solution for L = 2
  std::vector<double> temperatures = {1.0};
  // Different number of MC cycles (where the values are
  // chosen to obtain an acceptable computation time and convergence):
  arma::vec cycles = arma::logspace(4, 7.8, 38);
  model.analytical_comparison(temperatures, cycles, "analytical_comparison", seed);

  // // sigurd:
  double T_max = 2.6;
  double T_min = 2.1;

  std::vector<int> lattices = std::vector<int> {40};

  // // Problem 8:
  arma::vec temperatures_p8 = arma::linspace(T_min, T_max, 50);  // if we get results we need to upgrade this to more values
  model.phase_transitions(lattices, temperatures_p8, 15000, 3, seed, 5000);

  // trying to find burn in time based on cycles
  arma::vec temperatures = arma::vec {1., 1.5, 2., 2.4};

  // Sigurd spesial
  model.varying_n_mc_cycles(temperatures, 10000, 25, lattices, 0, "varying_cycles_no_burn_lattices_true", 250);
  model.varying_n_mc_cycles(temperatures, 15000, 25, lattices, 5000, "varying_cycles_yes_burn_lattices_true", 250);


  // looking for the distribution of the energies
  // model.epsilon_dist(temperatures, lattices, 5000, 5, 0, "energy_dist_lattices");

  return 0;
}
