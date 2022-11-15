#include "Ising.hpp"
#include "omp.h"
#include "Grid.hpp"
#include <armadillo>
#include <random>
#include <string>
#include <iomanip>
// #include <fstream>

/**
 * Analytical solution.
 *
 * @param L Lattice size
 * @param T Temperature.
 *
 * return avgs Vector containing the four different averages.
 */
std::vector<double> Ising::analytical(int L, double T){
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

/**
 * Perform one monte carlo cycle by creating an instance of Grid.
 *
 * @param N_spinflips
 * @param num_MC_cycles
 * @param avg_eps, avg_eps_sq, avg_m_abs, avg_m_sq References to averages that we add to
 * @param T
 * @param seed Seed for?????????????????????
 * @param L
 * @param random_config Whether initial configuration is randomized
 */
void Ising::monte_carlo(int N_spinflips, int num_MC_cycles, double &avg_eps, double &avg_eps_sq, double &avg_m_abs, double &avg_m_sq, double T, int seed, int L, bool random_config) {

  Grid model = Grid(L, T);

  for (int i = 0; i < num_MC_cycles; i ++){

    model.fill_grid(seed + i, random_config);
    model.random_walk(N_spinflips, seed + i);

    avg_eps = avg_eps + model.epsilon;
    avg_eps_sq = avg_eps_sq + model.epsilon_squared;

    avg_m_abs = avg_m_abs + model.m_abs;
    avg_m_sq = avg_m_sq + model.m_squared;

  }

  avg_eps /= num_MC_cycles;
  avg_eps_sq /= num_MC_cycles;
  avg_m_abs /= num_MC_cycles;
  avg_m_sq /= num_MC_cycles;

}

/**
 * Compare the algorithm with analytical solution at given temperature
 * for lattice size L = 2.
 *
 * @param temperatures
 * @param N_spinflips
 * @param N_MC_cycles
 * @param filename
 * @param seed
 */
void Ising::analytical_comparison(std::vector<double> temperatures, int N_spinflips, int N_MC_cycles, std::string filename, const int seed){

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

  #pragma omp parallel for
    for (double &temperature : temperatures){

        double avg_eps = 0.0;
        double avg_eps_sq = 0.0;
        double avg_m_abs = 0.0;
        double avg_m_sq = 0.0;

        auto thread_seed = MC_seed_generator();

        // Compute averages from model:
        monte_carlo(N_spinflips, N_MC_cycles, avg_eps, avg_eps_sq, avg_m_abs, avg_m_sq, temperature, thread_seed, 2, true);

        // Analytical values:
        std::vector<double> avgs = analytical(2, temperature);

        file << std::setprecision(print_prec) << temperature << ", "
            << std::setprecision(print_prec) << avg_eps << ", "
            << std::setprecision(print_prec) << avg_eps_sq << ", "
            << std::setprecision(print_prec) << avg_m_abs << ", "
            << std::setprecision(print_prec) << avg_m_sq << ", "
            << std::setprecision(print_prec) << avgs.at(0) << ", "
            << std::setprecision(print_prec) << avgs.at(1) << std::endl;
    }
  file.close();
}

/**
 * Vary the number of MC cycles for grid with given size and temperature.
 *
 * @param temperature
 * @param n_cycles Armadillo-vector containing the different MC-cycle-numbers.
 * @param filename
 * @param seed For generating thread seeds.
 * @param random_config Whether the initial config is random or not.
 * @param L Lattice size.
 * @param N_spinflips Number of spin flips to perform in each MC cycle.
 */
void Ising::varying_n_mc_cycles(double temperature, arma::vec n_cycles, std::string filename, int seed, bool random_config, int L, int N_spinflips){

  std::mt19937 MC_seed_generator (seed);

  // Open file to which we write the data:
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

      // Initialize averages:
      double avg_eps = 0.0;
      double avg_eps_sq = 0.0;
      double avg_m_abs = 0.0;
      double avg_m_sq = 0.0;

      auto thread_seed = MC_seed_generator();
      model.fill_grid(thread_seed, random_config);  // Fill grid

      monte_carlo(N_spinflips, n, avg_eps, avg_eps_sq, avg_m_abs, avg_m_sq, temperature, thread_seed, 20, random_config);

      file << std::setprecision(print_prec) << thread_seed << ", "
           << std::setprecision(print_prec) << temperature << ", "
           << std::setprecision(print_prec) << n << ", "
           << std::setprecision(print_prec) << avg_eps << ", "
           << std::setprecision(print_prec) << avg_m_abs
           << std::endl;
      }
  }
  file.close();
}

/**
 * For a given lattice size L, iterate through temperatures and write the
 * average values to a file.
 *
 * @param L
 * @param temperatures
 * @param N_spinflips
 * @param N_MC_cycles
 * @param seed
 * @param filename
 */
void Ising::phase_transitions(int L, arma::vec &temperatures, int N_spinflips, int N_MC_cycles, int seed, std::string filename){

  std::mt19937 MC_seed_generator (seed);

  std::ofstream file;
  file.open(filename + "_" + std::to_string(L) + ".txt", std::ofstream::trunc);
  static int print_prec = 10;

  file << std::setprecision(print_prec) << "Seed, "
       << std::setprecision(print_prec) << "Temperature, "
       << std::setprecision(print_prec) << "L, "
       << std::setprecision(print_prec) << "Epsilon (abs), "
       << std::setprecision(print_prec) << "Epsilon (squared), "
       << std::setprecision(print_prec) << "M (abs), "
       << std::setprecision(print_prec) << "M (squared)"
       << std::endl;

  #pragma omp parallel
  {
    auto thread_seed = MC_seed_generator();

    #pragma omp for
    for (double &temperature : temperatures){

       double avg_eps = 0.0;
       double avg_eps_sq = 0.0;
       double avg_m_abs = 0.0;
       double avg_m_sq = 0.0;

       monte_carlo(N_spinflips, N_MC_cycles, avg_eps, avg_eps_sq, avg_m_abs, avg_m_sq, temperature, thread_seed, L, true);

       file << std::setprecision(print_prec) << thread_seed << ", "
            << std::setprecision(print_prec) << temperature << ", "
            << std::setprecision(print_prec) << L << ", "
            << std::setprecision(print_prec) << avg_eps << ", "
            << std::setprecision(print_prec) << avg_eps_sq << ", "
            << std::setprecision(print_prec) << avg_m_abs << ", "
            << std::setprecision(print_prec) << avg_m_sq
            << std::endl;
      }
  }
  file.close();


}

/*

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
*/