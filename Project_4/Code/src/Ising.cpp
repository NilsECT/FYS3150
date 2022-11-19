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
 * @param burn Burn in time in number of cycles
 * @param random_config Whether initial configuration is randomized
 */
void Ising::monte_carlo(int num_MC_cycles, double &avg_eps, double &avg_eps_sq, double &avg_m_abs, double &avg_m_sq,double T, int seed, int L, int burn, bool random_config) {

  Grid model = Grid(L, T);
  model.fill_grid(seed, random_config);
  int flips = model.N;

  // Burn in
  for (int i = 0; i < burn; i ++) {
    int new_seed = seed - i;

    model.random_walk(new_seed);
  }

  for (int i = 0; i < num_MC_cycles; i ++){

    int new_seed = seed + i;

    model.random_walk(new_seed);

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
void Ising::analytical_comparison(std::vector<double> temperatures, int N_MC_cycles, std::string filename, const int seed){

  std::mt19937 MC_seed_generator (seed);

  std::ofstream file;
  file.open(filename + "_" + std::to_string(N_MC_cycles) + ".txt", std::ofstream::trunc);
  static int print_prec = 10;

  file<< std::setprecision(print_prec) << "Temperature [J/kb], "
       << std::setprecision(print_prec) << "Energy [J], "
       << std::setprecision(print_prec) << "Energy squared [J^2], "
       << std::setprecision(print_prec) << "Magnetisation, "
       << std::setprecision(print_prec) << "Magnetisation squared, "
       << std::setprecision(print_prec) << "Analytical energy [J], "
       << std::setprecision(print_prec) << "Analytical energy squared [J^2]"
       << std::endl;

  #pragma omp parallel for
    for (double &temperature : temperatures){

        double avg_eps = 0.0;
        double avg_eps_sq = 0.0;
        double avg_m_abs = 0.0;
        double avg_m_sq = 0.0;

        auto thread_seed = MC_seed_generator();

        // Compute averages from model:
        monte_carlo(N_MC_cycles, avg_eps, avg_eps_sq, avg_m_abs, avg_m_sq, temperature, thread_seed, 2);

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
void Ising::phase_transitions(int lattice, arma::vec &temperatures, int N_MC_cycles, int seed, int burn, std::string filename){

  std::mt19937 MC_seed_generator (seed);

  std::ofstream file;
  file.open(filename + "_" + std::to_string(lattice) + ".txt", std::ofstream::trunc);
  static int print_prec = 10;

  file << std::setprecision(print_prec) << "Burn in [cycles], "
       << std::setprecision(print_prec) << "Temperature [J/kb], "
       << std::setprecision(print_prec) << "Lattice, "
       << std::setprecision(print_prec) << "Energy [J], "
       << std::setprecision(print_prec) << "Energy squared [J^2], "
       << std::setprecision(print_prec) << "Magnetisation, "
       << std::setprecision(print_prec) << "Magnetisation squared, "
       << std::setprecision(print_prec) << "Heat capacity, "
       << std::setprecision(print_prec) << "Susceptibility"
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

      monte_carlo(N_MC_cycles, avg_eps, avg_eps_sq, avg_m_abs, avg_m_sq, temperature, thread_seed, lattice, burn);

      double chi = (lattice*lattice) * (avg_m_sq - avg_m_abs*avg_m_abs) / temperature;
      double cv = (lattice*lattice) * (avg_eps_sq - avg_eps*avg_eps) / (temperature*temperature);

      file << std::setprecision(print_prec) << burn << ", "
          << std::setprecision(print_prec) << temperature << ", "
          << std::setprecision(print_prec) << lattice << ", "
          << std::setprecision(print_prec) << avg_eps << ", "
          << std::setprecision(print_prec) << avg_eps_sq << ", "
          << std::setprecision(print_prec) << avg_m_abs << ", "
          << std::setprecision(print_prec) << avg_m_sq << ", "

            << std::setprecision(print_prec) << cv << ", "
            << std::setprecision(print_prec) << chi
          << std::endl;
      }
  }
  file.close();
}

void Ising::phase_transitions(std::vector<int> lattice, arma::vec temperatures,int N_MC_cycles, int n_samples, int seed, int burn, std::string filename){

  std::mt19937 MC_seed_generator (seed);

  std::ofstream file;
  file.open(filename + ".txt", std::ofstream::trunc);
  static int print_prec = 10;

  file << std::setprecision(print_prec) << "Burn in [cycles], "
       << std::setprecision(print_prec) << "Temperature [J/kb], "
       << std::setprecision(print_prec) << "Lattice, "
       << std::setprecision(print_prec) << "Energy [J], "
       << std::setprecision(print_prec) << "Energy squared [J^2], "
       << std::setprecision(print_prec) << "Magnetisation, "
       << std::setprecision(print_prec) << "Magnetisation squared, "
       << std::setprecision(print_prec) << "Sample, "
       << std::setprecision(print_prec) << "Heat capacity, "
       << std::setprecision(print_prec) << "Susceptibility"
       << std::endl;

  #pragma omp parallel for
  for (int i = 0; i < n_samples; i++) {
    for (int &L : lattice) {
      auto thread_seed = MC_seed_generator();

      #pragma omp parallel for
      for (double &temperature : temperatures){

        double avg_eps = 0.0;
        double avg_eps_sq = 0.0;
        double avg_m_abs = 0.0;
        double avg_m_sq = 0.0;

        monte_carlo(N_MC_cycles, avg_eps, avg_eps_sq, avg_m_abs, avg_m_sq, temperature, thread_seed, L, burn);

        double chi = (L*L) * (avg_m_sq - avg_m_abs*avg_m_abs) / temperature;
        double cv = (L*L) * (avg_eps_sq - avg_eps*avg_eps) / (temperature*temperature);

        #pragma omp critical
        file << std::setprecision(print_prec) << burn << ", "
              << std::setprecision(print_prec) << temperature << ", "
              << std::setprecision(print_prec) << L << ", "
              << std::setprecision(print_prec) << avg_eps << ", "
              << std::setprecision(print_prec) << avg_eps_sq << ", "
              << std::setprecision(print_prec) << avg_m_abs << ", "
              << std::setprecision(print_prec) << avg_m_sq << ", "
              << std::setprecision(print_prec) << i << ", "
              << std::setprecision(print_prec) << cv << ", "
              << std::setprecision(print_prec) << chi
              << std::endl;
      }
    }
  }
  file.close();
}

/**
 * @brief Samples all energies given by individual cycles and prints them to a file
 * 
 * @param temperature Armadillo vector containing the temperatures you want ot sample over
 * @param L int of lattice size for the model
 * @param N_cycles 
 * @param n_samples 
 * @param burn
 * @param filename 
 * @param seed 
 */
void Ising::epsilon_dist(arma::vec temperature, int L, int N_cycles, int n_samples, int burn, std::string filename, int seed) {
  std::mt19937 MC_seed_generator (seed);

  // Open file to which we write the data:
  std::ofstream file;
  file.open(filename + ".txt", std::ofstream::trunc);
  static int print_prec = 10;

  file << std::setprecision(print_prec) << "Temperature [J/kb], "
      << std::setprecision(print_prec) << "Lattice, "
      << std::setprecision(print_prec) << "Energy [J], "
      << std::setprecision(print_prec) << "Magnetisation, "
      << std::setprecision(print_prec) << "Sample, "
      << std::setprecision(print_prec) << "Cycle"
      << std::endl;

  int l_sq = L*L;

  #pragma omp prallel for
  for (double &T : temperature){
    for (int j = 0; j < n_samples; j++) {
      auto thread_seed = MC_seed_generator();

      Grid model = Grid(L, T);
      model.fill_grid(thread_seed);

      for (int i=0; i < burn; i++) {
        model.random_walk(thread_seed);
      }

      for (int i = 0; i < N_cycles; i++) {
          int seed = thread_seed + i;

          model.random_walk(seed);

          #pragma omp critical
          file << std::setprecision(print_prec) << T << ", "
              << std::setprecision(print_prec) << L << ", "
              << std::setprecision(print_prec) << model.epsilon << ", "
              << std::setprecision(print_prec) << model.m_abs << ", "
              << std::setprecision(print_prec) << j << ", "
              << std::setprecision(print_prec) << i
              << std::endl;
      }
    }
  }
  file.close();
}

/**
 * @brief Samples all energies given by individual cycles and prints them to a file
 * 
 * @param temperature Armadillo vector containing the temperatures you want to sample over
 * @param lattice standard vector<int> of all the lattice sizes you want to sample over
 * @param N_cycles 
 * @param n_samples 
 * @param burn 
 * @param filename 
 * @param seed 
 */
void Ising::epsilon_dist(arma::vec temperature, std::vector<int> lattice, int N_cycles, int n_samples, int burn, std::string filename, int seed) {
  std::mt19937 MC_seed_generator (seed);

  // Open file to which we write the data:
  std::ofstream file;
  file.open(filename + ".txt", std::ofstream::trunc);
  static int print_prec = 10;

  file << std::setprecision(print_prec) << "Temperature [J/kb], "
      << std::setprecision(print_prec) << "Lattice, "
      << std::setprecision(print_prec) << "Energy [J], "
      << std::setprecision(print_prec) << "Magnetisation, "
      << std::setprecision(print_prec) << "Sample, "
      << std::setprecision(print_prec) << "Cycle"
      << std::endl;

  #pragma omp parallel for
  for (int &L : lattice) {
    int l_sq = L*L;

    #pragma omp prallel for
    for (double &T : temperature){
      for (int j = 0; j < n_samples; j++) {
        auto thread_seed = MC_seed_generator();

        Grid model = Grid(L, T);
        model.fill_grid(thread_seed);

        for (int i=0; i < burn; i++) {
          model.random_walk(thread_seed);
        }

        for (int i = 0; i < N_cycles; i++) {
            int seed = thread_seed + i;

            model.random_walk(seed);

            #pragma omp critical
            file << std::setprecision(print_prec) << T << ", "
                << std::setprecision(print_prec) << L << ", "
                << std::setprecision(print_prec) << model.epsilon << ", "
                << std::setprecision(print_prec) << model.m_abs << ", "
                << std::setprecision(print_prec) << j << ", "
                << std::setprecision(print_prec) << i
                << std::endl;
        }
      }
    }
  }
  file.close();
}

/**
 * Vary the number of MC cycles for grid with given size and temperature.
 *
 * @param temperature Armadillo-vector containing different temperatures
 * @param n_cycles Total number of cycles which to check, has to be multiple of 25
 * @param n_samples Number of samples for each combination of temperature and cycles
 * @param L Lattice size.
 * @param filename
 * @param seed For generating thread seeds.
 */
void Ising::varying_n_mc_cycles(arma::vec temperature, int n_cycles, int n_samples, int lattice, std::string filename, int seed){

  std::mt19937 MC_seed_generator (seed);

  // Open file to which we write the data:
  std::ofstream file;
  file.open(filename + ".txt", std::ofstream::trunc);
  static int print_prec = 10;

  file << std::setprecision(print_prec) << "Order, "
       << std::setprecision(print_prec) << "Temperature [J/kb], "
       << std::setprecision(print_prec) << "Cycle, "
       << std::setprecision(print_prec) << "Energy [J], "
       << std::setprecision(print_prec) << "Magnetisation, "
       << std::setprecision(print_prec) << "Sample"
       << std::endl;

  std::string order = "Ordered";
  std::string unorder = "Unordered";

  #pragma omp parallel for
  for (int i = 0; i < n_samples; i++) {

    #pragma omp parallel for
    for (double &T : temperature) {
      for(int n = 25; n <= n_cycles; n+=25) {

        // Initialize averages:
        double avg_eps = 0.0;
        double avg_eps_sq = 0.0;
        double avg_m_abs = 0.0;
        double avg_m_sq = 0.0;

        auto thread_seed = MC_seed_generator();
        // model.fill_grid(thread_seed, random_config);  // Fill grid

        monte_carlo(n, avg_eps, avg_eps_sq, avg_m_abs, avg_m_sq, T, thread_seed, lattice, 0, false);

        #pragma omp critical
        file << std::setprecision(print_prec) << order << ", "
              << std::setprecision(print_prec) << T << ", "
              << std::setprecision(print_prec) << n << ", "
              << std::setprecision(print_prec) << avg_eps << ", "
              << std::setprecision(print_prec) << avg_m_abs << ", "
              << std::setprecision(print_prec) << i
              << std::endl;
      }
    }

    #pragma omp parallel for
    for (double &T : temperature) {
      for(int n = 25; n <= n_cycles; n+=25) {
        // Initialize averages:
        double avg_eps = 0.0;
        double avg_eps_sq = 0.0;
        double avg_m_abs = 0.0;
        double avg_m_sq = 0.0;

        auto thread_seed = MC_seed_generator();
        // model.fill_grid(thread_seed, random_config);  // Fill grid

        monte_carlo(n, avg_eps, avg_eps_sq, avg_m_abs, avg_m_sq, T, thread_seed, lattice, 0, true);

        #pragma omp critical
        file << std::setprecision(print_prec) << unorder << ", "
              << std::setprecision(print_prec) << T << ", "
              << std::setprecision(print_prec) << n << ", "
              << std::setprecision(print_prec) << avg_eps << ", "
              << std::setprecision(print_prec) << avg_m_abs << ", "
              << std::setprecision(print_prec) << i
              << std::endl;
      }
    }
  }


  file.close();
}

/**
 * Vary the number of MC cycles for grid with given size and temperature.
 *
 * @param temperature Armadillo-vector containing different temperatures
 * @param n_cycles Total number of cycles which to check, has to be multiple of 25
 * @param n_samples Number of samples for each combination of temperature and cycles
 * @param L Lattice sizes for wihch to loop over in a standard int vector
 * @param filename
 * @param seed For generating thread seeds.
 */
void Ising::varying_n_mc_cycles(arma::vec temperature, int n_cycles, int n_samples, std::vector<int> lattice, std::string filename, int seed){

  std::mt19937 MC_seed_generator (seed);

  // Open file to which we write the data:
  std::ofstream file;
  file.open(filename + ".txt", std::ofstream::trunc);
  static int print_prec = 10;

  file << std::setprecision(print_prec) << "Order, "
       << std::setprecision(print_prec) << "Temperature [J/kb], "
       << std::setprecision(print_prec) << "Cycle, "
       << std::setprecision(print_prec) << "Energy [J], "
       << std::setprecision(print_prec) << "Magnetisation, "
       << std::setprecision(print_prec) << "Lattice, "
       << std::setprecision(print_prec) << "Sample"
       << std::endl;

  std::string order = "Ordered";
  std::string unorder = "Unordered";

  #pragma omp parallel for
  for (int i = 0; i < n_samples; i++) {
    for (int &L : lattice) {

      #pragma omp parallel for
      for (double &T : temperature) {
        for(int n = 25; n <= n_cycles; n+=25) {

          // Initialize averages:
          double avg_eps = 0.0;
          double avg_eps_sq = 0.0;
          double avg_m_abs = 0.0;
          double avg_m_sq = 0.0;

          auto thread_seed = MC_seed_generator();
          // model.fill_grid(thread_seed, random_config);  // Fill grid

          monte_carlo(n, avg_eps, avg_eps_sq, avg_m_abs, avg_m_sq, T, thread_seed, L, 0, false);

          #pragma omp critical
          file << std::setprecision(print_prec) << order << ", "
                << std::setprecision(print_prec) << T << ", "
                << std::setprecision(print_prec) << n << ", "
                << std::setprecision(print_prec) << avg_eps << ", "
                << std::setprecision(print_prec) << avg_m_abs << ", "
                << std::setprecision(print_prec) << L << ", "
                << std::setprecision(print_prec) << i
                << std::endl;
        }
      }

      #pragma omp parallel for
      for (double &T : temperature) {
        for(int n = 25; n <= n_cycles; n+=25) {
          // Initialize averages:
          double avg_eps = 0.0;
          double avg_eps_sq = 0.0;
          double avg_m_abs = 0.0;
          double avg_m_sq = 0.0;

          auto thread_seed = MC_seed_generator();
          // model.fill_grid(thread_seed, random_config);  // Fill grid

          monte_carlo(n, avg_eps, avg_eps_sq, avg_m_abs, avg_m_sq, T, thread_seed, L, 0, true);

          #pragma omp critical
          file << std::setprecision(print_prec) << unorder << ", "
                << std::setprecision(print_prec) << T << ", "
                << std::setprecision(print_prec) << n << ", "
                << std::setprecision(print_prec) << avg_eps << ", "
                << std::setprecision(print_prec) << avg_m_abs << ", "
                << std::setprecision(print_prec) << L << ", "
                << std::setprecision(print_prec) << i
                << std::endl;
        }
      }
    }
  }
  file.close();
}


// THESE TWO LAST ONES STILL NEED CHANGING

// void Ising::varying_n_walk(arma::vec temperature, std::vector<int> n_walks, std::vector<int> lattice, int num_samples, std::string filename, int seed) {

//   // Open file to which we write the data:
//   std::ofstream file;
//   file.open(filename + ".txt", std::ofstream::trunc);
//   static int print_prec = 10;

//   file << std::setprecision(print_prec) << "Sample, "
//        << std::setprecision(print_prec) << "Total walks, "
//        << std::setprecision(print_prec) << "Walk, "
//        << std::setprecision(print_prec) << "Temperature [J/kb], "
//        << std::setprecision(print_prec) << "Lattice, "
//        << std::setprecision(print_prec) << "Energy [J]"
//        << std::endl;

//   #pragma omp parallel for
//   for (int &L : lattice) {
//     for (double &T : temperature) {

//       #pragma omp prallel for
//       for (int &i : n_walks) {
//         for (int ii = 0; ii < num_samples; ii++) {

//           int seed_other = seed +100 +i +ii;
//           Grid model = Grid(L, T);
//           model.fill_grid(seed_other);

//           for (int iii = 0; iii < i; iii++) {
//             int rand_seed = seed + iii + ii + i;
//             model.one_walk(rand_seed);

//             #pragma omp critical
//             file << std::setprecision(print_prec) << ii << ", "
//                   << std::setprecision(print_prec) << i << ", "
//                   << std::setprecision(print_prec) << iii << ", "
//                   << std::setprecision(print_prec) << T << ", "
//                   << std::setprecision(print_prec) << L << ", "
//                   << std::setprecision(print_prec) << model.epsilon
//                   << std::endl;
//           }
//         }
//       }
//     }
//   }
//   file.close();
// }

// void Ising::varying_n_walk(arma::vec temperature, int n_walks, std::vector<int> lattice, int num_samples, std::string filename, int seed) {

//   // Open file to which we write the data:
//   std::ofstream file;
//   file.open(filename + ".txt", std::ofstream::trunc);
//   static int print_prec = 10;

//   file << std::setprecision(print_prec) << "Sample, "
//        << std::setprecision(print_prec) << "Total walks, "
//        << std::setprecision(print_prec) << "Walk, "
//        << std::setprecision(print_prec) << "Temperature [J/kb], "
//        << std::setprecision(print_prec) << "Lattice, "
//        << std::setprecision(print_prec) << "Energy [J]"
//        << std::endl;

//   #pragma omp parallel for
//   for (int &L : lattice) {
//     for (double &T : temperature) {

//       #pragma omp prallel for
//       for (int ii = 0; ii < num_samples; ii++) {

//         int seed_other = seed +100 +ii;
//         Grid model = Grid(L, T);
//         model.fill_grid(seed_other);

//         for (int iii = 0; iii < n_walks; iii++) {
//           int rand_seed = seed + iii + ii;
//           model.one_walk(rand_seed);

//           #pragma omp critical
//           file << std::setprecision(print_prec) << ii << ", "
//                 << std::setprecision(print_prec) << n_walks << ", "
//                 << std::setprecision(print_prec) << iii << ", "
//                 << std::setprecision(print_prec) << T << ", "
//                 << std::setprecision(print_prec) << L << ", "
//                 << std::setprecision(print_prec) << model.epsilon
//                 << std::endl;
//         }
//       }
//     }
//   }
//   file.close();
// }

// /**
//  * Vary the number of MC cycles for grid with given size and temperature.
//  *
//  * @param temperature
//  * @param n_cycles Armadillo-vector containing the different MC-cycle-numbers.
//  * @param filename
//  * @param seed For generating thread seeds.
//  * @param random_config Whether the initial config is random or not.
//  * @param L Lattice size.
//  * @param N_spinflips Number of spin flips to perform in each MC cycle.
//  */
// void Ising::varying_n_mc_cycles(double temperature, arma::vec n_cycles, std::string filename, int seed, bool random_config, int L, int n_samples){

//   std::mt19937 MC_seed_generator (seed);

//   // Open file to which we write the data:
//   std::ofstream file;
//   file.open(filename + "_" + std::to_string(temperature) + ".txt", std::ofstream::trunc);
//   static int print_prec = 10;

//   file << std::setprecision(print_prec) << "Seed, "
//        << std::setprecision(print_prec) << "Temperature [J/kb], "
//        << std::setprecision(print_prec) << "MC cycles, "
//        << std::setprecision(print_prec) << "Energy [J], "
//        << std::setprecision(print_prec) << "Magnetisation, "
//        << std::setprecision(print_prec) << "Sample"
//        << std::endl;


//   #pragma omp parallel for
//   for (double &n : n_cycles){
//     auto thread_seed = MC_seed_generator();

//     #pragma omp parallel for
//     for (int i = 0; i < n_samples; i++) {

//       // Initialize averages:
//       double avg_eps = 0.0;
//       double avg_eps_sq = 0.0;
//       double avg_m_abs = 0.0;
//       double avg_m_sq = 0.0;

//       int sample_seed = thread_seed + n_samples + i;

//       monte_carlo(n, avg_eps, avg_eps_sq, avg_m_abs, avg_m_sq, temperature, sample_seed, 20, random_config);

//       #pragma omp critical
//       file << std::setprecision(print_prec) << thread_seed << ", "
//           << std::setprecision(print_prec) << temperature << ", "
//           << std::setprecision(print_prec) << n << ", "
//           << std::setprecision(print_prec) << avg_eps << ", "
//           << std::setprecision(print_prec) << avg_m_abs << ", "
//           << std::setprecision(print_prec) << i
//           << std::endl;
//     }
//   }
//   file.close();
// }