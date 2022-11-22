#include <fstream>
#include "omp.h"
#include "Grid.hpp"
#include <armadillo>
#include <iostream>
#include <string>
#include <random>
#include <iomanip>      
#include "Ising.hpp"

int main(){

  const int seed = 137; // Our favorite seed(based on the fine structure constant) 
  // model instance:
  Ising model = Ising();

  // Compare with analytical solution for L = 2
  std::vector<double> temperature_a = {1.0};

  // Different number of MC cycles (where the values are chosen to obtain an acceptable computation time and convergence):
  arma::vec cycles = arma::logspace(4, 7.8, 38);
  model.analytical_comparison(temperature_a, cycles, "analytical_comparison", seed);

  // Our lattice size: 
  std::vector<int> lattice= {20};
  arma::vec temperatures = arma::vec {1., 1.5, 2., 2.4};
  // looking for the distribution of the energies
  model.epsilon_dist(temperatures, lattice, 5000, 5, 0, "energy_dist_lattices");

  // trying to find burn in time based on cycles
  model.varying_n_mc_cycles(temperatures, 10000, 25, lattice, 0, "varying_cycles_no_burn_lattices_true", 250);
  model.varying_n_mc_cycles(temperatures, 15000, 25, lattice, 5000, "varying_cycles_yes_burn_lattices_true", 250);

  // more lattice sizes:
  std::vector<int> lattices = std::vector<int> {20, 40, 60, 80, 100};
  std::string filename = "phase_transition_varL_";

  // First scan:
  double T_max = 2.6;
  double T_min = 2.1;
  std::vector<int> steps = std::vector<int> {23, 24, 25};
  for (int step : steps){
    arma::vec temperatures_p8 = arma::linspace(T_min, T_max, step);  
    model.phase_transitions(lattices, temperatures_p8, 15000, 3, seed, 5000,filename+"step_"+std::to_string(step));
  }

  // Detailed scan:
  T_max = 2.35;
  T_min = 2.25;
  int detailed_step = 100;
  arma::vec temperatures_p8 = arma::linspace(T_min, T_max, detailed_step);  
  model.phase_transitions(lattices, temperatures_p8, 15000, 3, seed, 5000,filename+"step"+std::to_string(detailed_step));

  return 0;
}
