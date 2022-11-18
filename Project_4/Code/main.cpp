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

  int N_MC_cycles = 500000;

  // model instance:
  Ising model = Ising();

  // Working code for comparing with analytical solution (currently missing magnetization, will fix soon):
  std::vector<double> temps = {1.0, 2.4};
  model.analytical_comparison(temps, 10000000, "analytical_comparison", seed);

  // --------
  
  // Code that currently does not work:

  // trying to find burn in time based on walks
  arma::vec temperatures_burn = arma::vec {1., 1.5, 2., 2.4, 3.};
  // std::vector<int> n_walks = {15000};
  // int n_walks = 12000;
  // std::vector<int> lattices = {20, 40};
  // model.varying_n_walk(temperatures_burn, n_walks, lattices, 50);

  // std::cout << "Starting analytical comparison" << std::endl;
  // model.analytical_comparison(temperatures, N_spinflips, N_MC_cycles);

  // SETT N=L^2 !!!!!!!!!!!!!!!!!!!!!!!

  // // Problem 5: study burn-in time
  arma::vec cycles = arma::logspace(1, 6, 50);
  arma::vec temperatures_p5 = arma::vec {1., 1.5, 2., 2.4, 3.};
  model.varying_n_mc_cycles(temperatures_p5, cycles);

  std::vector<double> temperatures_2 = {1., 1.5, 2., 2.4, 3.};
  model.epsilon_dist(temperatures_2, 20, 1000000, 10);

  // // sigurd:
  // int T_max = 2.4;
  // int T_min = 2.3;


  // // Nils:
  int T_max = 2.2;
  int T_min = 2.1;


  // // Brage:
  // int T_max = 2.28;
  // int T_min = 2.22;


  // // // Problem 8:
  std::vector<int> lattice_sizes_p8 = {40, 60, 80, 100};
  arma::vec temperatures_p8 = arma::linspace(T_min, T_max, 5);  // if we get results we need to upgrade this to more values
  model.phase_transitions(lattice_sizes_p8, temperatures_p8, N_spinflips, N_MC_cycles, seed);

  return 0;
}
