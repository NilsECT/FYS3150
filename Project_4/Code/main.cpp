#include <fstream>
#include "omp.h"
#include "Grid.hpp"
#include <armadillo>
#include <iostream>
#include <string>
#include <random>
/*
void write_to_file(int thread_seed, int L, double T, std::vector<double> thread_avg){
  std::ofstream ofile;
  st
  static int print_prec = 10;

}
*/

std::vector<double> MCMC(Grid G, int M, int thread_seed){

  std::mt19937 generator (thread_seed);
  std::uniform_int_distribution<int> dis(0, G.L-1);
  std::uniform_real_distribution<double> spinflip(0, 1);

  double r;

  double epsilon = 0.0;
  double epsilon_squared = 0.0;
  double m_squared = 0.0;
  double m_abs = 0.0;
  double current_magnetization;

  //std::vector<double> dEs = {8, 4, 0};

  for (int i = 0; i < M; i++){
    int s_x = dis(generator);
    int s_y = dis(generator);

    // Compute change in energy:
    double dE = G.grid((s_x+1) % G.L, s_y) + G.grid((G.L + s_x - 1) % G.L, s_y)
        + G.grid(s_x, (s_y+1) % G.L) + G.grid(s_x, (G.L + s_y-1) % G.L);
    dE = 2 * G.grid(s_x, s_y) * dE;

    // Ratio between probability of new state and old state:
    dE = exp(-dE / G.T);

    // Add current energy and magnetization values to average sum:
    epsilon += G.E;
    epsilon_squared = epsilon_squared + G.E * G.E;

    current_magnetization = G.M;
    m_abs = m_abs + std::sqrt(current_magnetization * current_magnetization);
    m_squared = m_squared + current_magnetization * current_magnetization;

    r = spinflip(generator);

    if (r <= dE) {
      G.flip_spin(s_x, s_y);
    }

  }

  // Divide averages by number of Monte Carlo cycles (and number of spins):
  epsilon = epsilon / (G.L*G.L*M);
  epsilon_squared = epsilon_squared / (M * std::pow(G.L, 4));
  m_abs = m_abs / (M * G.L * G.L);
  m_squared = m_squared / (M * std::pow(G.L, 4));

  std::cout << "THREAD SEED: " << thread_seed
            << ", AVG MAGNETIZATION: " << m_abs
            << ", EPSILON_MEAN: " << epsilon << ", TEMPERATURE: " << G.T << std::endl;

  std::vector<double> averages = {epsilon, epsilon_squared, m_abs, m_squared};
  return averages;
}


int main(int argc, char* argv[]){

  const int L = atoi(argv[1]);
  const int num_threads = atoi(argv[2]);
  const int num_MC_cycles = atoi(argv[3]);

  const int seed = 137;
  std::mt19937 MC_seed_generator (seed);

  double T_min = 1.6;
  double T_max = 2.6;
  double dT = (T_max - T_min) / num_threads;

  std::ofstream ofile;
  static int print_prec = 10;
  //std::string filename = "averages_" + L + "_" + num_MC_cycles + "_" + num_threads + ".txt";
  ofile.open("averages_" + std::to_string(L) + "_" + std::to_string(num_MC_cycles) + "_" + std::to_string(num_threads) + ".txt", std::ofstream::trunc);

  ofile << "Thread seed, Temperature, Epsilon, Epsilon squared, Abs(magnetization), Magnetization squared" << std::endl;
  #pragma omp parallel
  {

    # pragma omp for
    for (int i = 0; i < num_threads; i++){
      double T = T_min + i * dT;

      Grid G = Grid(L, T);

      auto thread_seed = MC_seed_generator();

      G.fill_grid(thread_seed);

      std::vector<double> thread_avgs = MCMC(G, num_MC_cycles, thread_seed);

      /*std::cout << "THREAD SEED: " << thread_seed
          << ", EPSILON_MEAN: " << epsilon_mean << ", TEMPERATURE: " << T << std::endl;
      */

      //std::cout << my_thread << ", " << i << std::endl;
      /*ofile << std::setprecision(print_prec) << scientific << i << ", "
            << std::setprecision(print_prec) << scientific << epsilon_mean << std::endl;
      */
      ofile << std::setprecision(print_prec) << thread_seed << ", "
            << std::setprecision(print_prec) << T << ", "
            << std::setprecision(print_prec) << thread_avgs.at(0) << ", "
            << std::setprecision(print_prec) << thread_avgs.at(1) << ", "
            << std::setprecision(print_prec) << thread_avgs.at(2) << ", "
            << std::setprecision(print_prec) << thread_avgs.at(3) << std::endl;
    }
  }

  return 0;
}
