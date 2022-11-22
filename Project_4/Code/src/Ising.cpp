#include "Ising.hpp"
#include "omp.h"
#include "Grid.hpp"
#include <armadillo>
#include <random>
#include <string>
#include <iomanip>
// #include <fstream>

/**
 * Analytical solution for 2D Ising model of dimension 2x2.
 *
 * @param L Lattice size
 * @param T Temperature.
 *
 * return avgs Vector containing mean epsilon, mean absolute m, heat capacity, and chi.
 */
std::vector<double> Ising::analytical(int L, double T){
  double Z;
  double mean_E;
  double mean_E_squared;
  double chi;
  double Cv;
  double m_abs;
  double mean_m_squared;
  double N = L*L;

  Z = 2 * (6 + exp(8/T) + exp(-8/T));  // Partition function

  // Mean energy per spin and magnetisation per spin:
  mean_E = (4/Z) * (exp(-8/T) - exp(8/T));
  m_abs = (4 + 2*exp(8/T)) / Z;

  // Mean energy squared and magnetisation squared:
  mean_E_squared = N * N * (8/Z) * (exp(8/T) + exp(-8/T));
  mean_m_squared = (2/Z) * (1 + exp(8/T));

  // Heat capacity and susceptibility per spin:
  Cv = (mean_E_squared - mean_E*mean_E*N*N) / (N * T * T);
  chi = (16 / (Z*Z*T)) * (3 + 3*exp(8/T) + exp(-8/T));

  std::vector<double> avgs = {mean_E, m_abs, Cv, chi};
  return avgs;
}

/**
 * Perform one monte carlo cycle on a given grid, and add to average values.
 *
 * @param num_MC_cycles Number of Monte Carlo cycles.
 * @param avg_eps, avg_eps_sq, avg_m_abs, avg_m_sq References to averages that we add to.
 * @param model An instance of Grid.
 * @param seed
 * @param burn Burn in time in number of Monte Carlo cycles.
 */
void Ising::monte_carlo(int num_MC_cycles, double &avg_eps, double &avg_eps_sq, double &avg_m_abs, double &avg_m_sq, Grid &model, int seed,int burn) {
  int flips = model.N;

  this->burn_in(model, burn, seed-1);

  for (int i = 0; i < (num_MC_cycles); i ++){

    int new_seed = seed + i;

    // Perform one Markov Chain walk:
    model.random_walk(new_seed);

    // Add current values to averages:
    avg_eps = avg_eps + model.epsilon;
    avg_eps_sq = avg_eps_sq + model.epsilon_squared;

    avg_m_abs = avg_m_abs + model.m_abs;
    avg_m_sq = avg_m_sq + model.m_squared;
  }

  avg_eps /= (num_MC_cycles);
  avg_m_abs /= (num_MC_cycles);

  avg_eps_sq /= (num_MC_cycles);
  avg_m_sq /= (num_MC_cycles);

}

/**
 * Compare the algorithm with analytical solution at given temperature
 * for lattice size L = 2.
 *
 * @param temperatures
 * @param cycles
 * @param filename
 * @param seed
 */
void Ising::analytical_comparison(std::vector<double> temperatures, arma::vec cycles, std::string filename, const int seed){

  std::mt19937 MC_seed_generator (seed);
  int grid_L = 2;
  int grid_N = grid_L * grid_L;

  std::ofstream file;
  file.open(filename + ".txt", std::ofstream::trunc);
  static int print_prec = 10;

  file<< std::setprecision(print_prec) << "Temperature [J/kb], "
       << std::setprecision(print_prec) << "MC cycles, "
       << std::setprecision(print_prec) << "Energy per spin [J], "
       << std::setprecision(print_prec) << "Absolute magnetisation per spin, "
       << std::setprecision(print_prec) << "Heat capacity, "
       << std::setprecision(print_prec) << "Susceptibility, "
       << std::setprecision(print_prec) << "Analytical energy per spin [J], "
       << std::setprecision(print_prec) << "Analytical absolute magnetisation per spin, "
       << std::setprecision(print_prec) << "Analytical heat capacity, "
       << std::setprecision(print_prec) << "Analytical susceptibility"
       << std::endl;

    for (double &temperature : temperatures){

        #pragma omp parallel for
        for (double &MC_cycles : cycles){

          double avg_eps = 0.0;
          double avg_eps_sq = 0.0;
          double avg_m_abs = 0.0;
          double avg_m_sq = 0.0;

          const int threadID = omp_get_thread_num();
          auto thread_seed = seed + threadID;

          // Initialize grid
          Grid model = Grid(grid_L, temperature);
          model.fill_grid(thread_seed);

          // Compute averages from model:
          monte_carlo(MC_cycles, avg_eps, avg_eps_sq, avg_m_abs, avg_m_sq, model, thread_seed);

          // Heat capacity and susceptibility:
          double chi = (avg_m_sq - (avg_m_abs*avg_m_abs)) / (temperature*grid_N);
          double cv = (avg_eps_sq - (avg_eps*avg_eps)) / (temperature*temperature*grid_N);

          // Divide by number of spins:
          avg_eps /= grid_N;
          avg_m_abs /= (grid_N);

          // Analytical values:
          std::vector<double> avgs = analytical(grid_L, temperature);

          // Write to file:
          #pragma omp critical
          file << std::setprecision(print_prec) << temperature << ", "
               << std::setprecision(print_prec) << MC_cycles << ", "
               << std::setprecision(print_prec) << avg_eps << ", "
               << std::setprecision(print_prec) << avg_m_abs << ", "
               << std::setprecision(print_prec) << cv << ", "
               << std::setprecision(print_prec) << chi << ", "
               << std::setprecision(print_prec) << avgs.at(0) << ", "
               << std::setprecision(print_prec) << avgs.at(1) << ", "
               << std::setprecision(print_prec) << avgs.at(2) << ", "
               << std::setprecision(print_prec) << avgs.at(3) << std::endl;
      }
    }

  file.close();
}

void Ising::phase_transitions(std::vector<int> lattice, arma::vec temperatures,int N_MC_cycles, int n_samples, int seed, int burn, std::string filename){

  std::ofstream file;
  file.open(filename + "fuck.txt", std::ofstream::trunc);
  static int print_prec = 10;

  file << std::setprecision(print_prec) << "Burn in [cycles], "
       << std::setprecision(print_prec) << "Temperature [J/kb], "
       << std::setprecision(print_prec) << "Lattice, "
       << std::setprecision(print_prec) << "Average energy [J], "
       << std::setprecision(print_prec) << "Energy squared [J^2], "
       << std::setprecision(print_prec) << "Average magnetisation, "
       << std::setprecision(print_prec) << "Magnetisation squared, "
       << std::setprecision(print_prec) << "Sample, "
       << std::setprecision(print_prec) << "Heat capacity, "
       << std::setprecision(print_prec) << "Susceptibility"
       << std::endl;


  //omp_set_nested(1);
  // #pragma omp parallel for
  for (int i = 0; i < n_samples; i++) {
    int sample_seed = seed + i;

    // #pragma omp parallel for
    for (int &L : lattice) {

      #pragma omp parallel for
      for (double &temperature : temperatures){

        const int threadID = omp_get_thread_num();

        std::mt19937 MC_seed_generator;
        auto thread_seed = sample_seed + threadID;
        MC_seed_generator.seed(thread_seed);

        int new_seed = thread_seed + temperature + L + i;

        double avg_eps = 0.0;
        double avg_eps_sq = 0.0;
        double avg_m_abs = 0.0;
        double avg_m_sq = 0.0;

        // make grid
        Grid model = Grid(L, temperature);
        model.fill_grid(thread_seed);

        monte_carlo(N_MC_cycles, avg_eps, avg_eps_sq, avg_m_abs, avg_m_sq, model, thread_seed, burn);

        int temp_N = model.N;

        double chi = (avg_m_sq - (avg_m_abs*avg_m_abs)) / (temperature*temp_N);
        double cv = (avg_eps_sq - (avg_eps*avg_eps)) / (temperature*temperature*temp_N);

        avg_eps /= temp_N;
        avg_m_abs /= temp_N;

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
 * @param temperature Armadillo vector containing the temperatures you want to sample over
 * @param lattice standard vector<int> of all the lattice sizes you want to sample over
 * @param N_cycles
 * @param n_samples
 * @param burn
 * @param filename
 * @param seed
 */
void Ising::epsilon_dist(arma::vec temperature, std::vector<int> lattice, int N_cycles, int n_samples, int burn, std::string filename, int seed) {
  // std::mt19937 MC_seed_generator (seed);

  // Open file to which we write the data:
  std::ofstream file;
  file.open(filename + ".txt", std::ofstream::trunc);
  static int print_prec = 10;

  file << std::setprecision(print_prec) << "Temperature [J/kb], "
      << std::setprecision(print_prec) << "Lattice, "
      << std::setprecision(print_prec) << "Average energy [J], "
      << std::setprecision(print_prec) << "Average magnetisation, "
      << std::setprecision(print_prec) << "Sample, "
      << std::setprecision(print_prec) << "Cycle"
      << std::endl;

  for (int j = 0; j < n_samples; j++) {

    int sample_seed = seed + j;

    #pragma omp parallel for
    for (int &L : lattice) {
      for (double &T : temperature){

        int l_sq = L*L;

        const int threadID = omp_get_thread_num();

        std::mt19937 MC_seed_generator;
        auto thread_seed = sample_seed + threadID;
        MC_seed_generator.seed(thread_seed);

        Grid model = Grid(L, T);
        model.fill_grid(thread_seed);

        for (int i=0; i < burn; i++) {
          model.random_walk(thread_seed);
        }

        int temp_N = model.N;

        for (int i = 0; i < N_cycles; i++) {
            int seed = thread_seed + i;

            model.random_walk(seed);

            #pragma omp critical
            file << std::setprecision(print_prec) << T << ", "
                << std::setprecision(print_prec) << L << ", "
                << std::setprecision(print_prec) << model.epsilon / temp_N << ", "
                << std::setprecision(print_prec) << model.m_abs / temp_N << ", "
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
 * @param L Lattice sizes for wihch to loop over in a standard int vector
 * @param filename
 * @param seed For generating thread seeds.
 */
void Ising::varying_n_mc_cycles(arma::vec temperature, int n_cycles, int n_samples, std::vector<int> lattice, int burn, std::string filename, int step, int seed){

  // std::mt19937 MC_seed_generator (seed);

  // Open file to which we write the data:
  std::ofstream file;
  file.open(filename + ".txt", std::ofstream::trunc);
  static int print_prec = 10;

  file << std::setprecision(print_prec) << "Order, "
       << std::setprecision(print_prec) << "Temperature [J/kb], "
       << std::setprecision(print_prec) << "Cycle, "
       << std::setprecision(print_prec) << "Average energy [J], "
       << std::setprecision(print_prec) << "Average magnetisation, "
       << std::setprecision(print_prec) << "Lattice, "
       << std::setprecision(print_prec) << "Sample"
       << std::endl;

  for (int i = 0; i < n_samples; i++) {

    int sample_seed = seed + i;
    std::vector<bool> container = std::vector<bool>{true, false};

    #pragma omp parallel for
    for (double &T : temperature) {
      for (int &L : lattice) {

        std::vector<bool> container = std::vector<bool>{true, false};

        const int threadID = omp_get_thread_num();

        std::mt19937 MC_seed_generator;
        auto thread_seed = sample_seed + threadID;
        MC_seed_generator.seed(thread_seed);

        for (bool start : container) {

          Grid model = Grid(L, T);
          model.fill_grid(thread_seed, start);
          this->burn_in(model, burn, thread_seed);

          // Initialize averages:
          double avg_eps = 0.0;
          double avg_eps_sq = 0.0;
          double avg_m_abs = 0.0;
          double avg_m_sq = 0.0;

          std::string order = start ? ("Unordered") : ("Ordered");  // labels are wrong some weird reason we can't find

          int temp_N = model.N;

          for (int n=step; n < n_cycles; n+=step) {
            monte_carlo(n, avg_eps, avg_eps_sq, avg_m_abs, avg_m_sq, model, thread_seed);

            #pragma omp critical
            file << std::setprecision(print_prec) << order << ", "
                    << std::setprecision(print_prec) << T << ", "
                    << std::setprecision(print_prec) << n << ", "
                    << std::setprecision(print_prec) << avg_eps / temp_N << ", "
                    << std::setprecision(print_prec) << avg_m_abs / temp_N << ", "
                    << std::setprecision(print_prec) << L << ", "
                    << std::setprecision(print_prec) << i
                    << std::endl;
          }
        }
      }
    }
  }
  file.close();
}

void Ising::burn_in(Grid &model, int burn, int seed) {
  // Burn in
  for (int i = 0; i < burn; i ++) {
    int new_seed = seed - i;
    model.random_walk(new_seed);
  }
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
