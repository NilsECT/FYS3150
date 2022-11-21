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

  // trying to find burn in time based on cycles
  // arma::vec temperatures = arma::vec {1., 1.5, 2., 2.4};

  // Sigurd spesial
  // std::vector<int> lattices = std::vector<int> {20, 40, 60, 100};
  // model.varying_n_mc_cycles(temperatures, 10000, 25, 20, 0, "varying_cycles_no_burn");
  // model.varying_n_mc_cycles(temperatures, 5000, 25, 20, 1000, "varying_cycles_yes_burn");


  // // looking for the distribution of the energies
  // model.epsilon_dist(temperatures, lattices, 2000, 25);

  // Everyone else DO NOT RUN IF YOUR PC IS WEAK
  // model.varying_n_mc_cycles(temperatures, 10000, 25, 20, 0, "varying_cycles_no_burn");
  // model.varying_n_mc_cycles(temperatures, 5000, 25, 20, 1000, "varying_cycles_yes_burn");
  // looking for the distribution of the energies
  // model.epsilon_dist(temperatures, 20, 2000, 25);


  // // sigurd:
  double T_max = 2.35;
  double T_min = 2.2;


  // // Nils:
  // double T_max = 2.4;
  // double T_min = 3.;


  // // Brage:
  // double T_max = 1.;
  // double T_min = 2.1;


  // // Problem 8:
  // std::vector<int> lattice_sizes_p8 = {40, 60, 80, 100};
  std::vector<int> lattice_sizes_p8 = {40};
  arma::vec temperatures_p8 = arma::linspace(T_min, T_max, 100);  // if we get results we need to upgrade this to more values
  model.phase_transitions(lattice_sizes_p8, temperatures_p8, 2000, 5, seed, 1000);

  return 0;
}
