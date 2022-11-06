#include "omp.h"  // OpenMP header
#include <armadillo>
#include "Grid.hpp"
#include <iostream>
#include <math.h>
#include <random>

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
    double dE = G.grid((s_x+1) % G.L, s_y) + G.grid((G.L + s_x - 1) % G.L, s_y) + G.grid(s_x, (s_y+1) % G.L) + G.grid(s_x, (G.L + s_y-1) % G.L);
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

  //std::mt19937 generator(137);

  const int L = atoi(argv[1]);
  const int num_threads = atoi(argv[2]);
  const int num_MC_cycles = atoi(argv[3]);

  //double J = 1.0;//std::pow(10, -23);
  //double T = 1.0;

  const int seed = 137;
  std::mt19937 MC_seed_generator (seed);

  double T_min = 1.6;
  double T_max = 2.6;
  double dT = (T_max - T_min) / num_threads;

  #pragma omp parallel
  {

    # pragma omp for
    for (int i = 0; i < num_threads; i++){
      double T = T_min + i * dT;

      Grid G = Grid(L, T);

      auto thread_seed = MC_seed_generator();

      G.fill_grid(thread_seed);

      std::vector<double> thread_averages = MCMC(G, num_MC_cycles, thread_seed);
      /*std::cout << "THREAD SEED: " << thread_seed
          << ", EPSILON_MEAN: " << epsilon_mean << ", TEMPERATURE: " << T << std::endl;
      */

      //std::cout << my_thread << ", " << i << std::endl;
      /*ofile << setprecision(print_precision) << scientific << i << ", "
            << setprecision(print_precision) << scientific << epsilon_mean << std::endl;
      */
    }
  }


  return 0;
}
