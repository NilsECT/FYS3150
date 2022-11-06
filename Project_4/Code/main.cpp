#include "omp.h"  // OpenMP header
#include <armadillo>
#include "Grid.hpp"
#include <iostream>
#include <math.h>
#include <random>

double MCMC(Grid G, int M, int thread_seed){
  //int seed = 137;
  std::mt19937 generator (thread_seed);
  std::uniform_int_distribution<int> dis(0, G.L-1);
  std::uniform_real_distribution<double> spinflip(0, 1);
  double r;
  double epsilon = 0.0;
  double magnetization = 0.0;
  double current_magnetization;

  //std::vector<double> dEs = {8, 4, 0};

  for (int i = 0; i < M; i++){
    int s_x = dis(generator);
    int s_y = dis(generator);

    double dE = G.grid((s_x+1) % G.L, s_y) + G.grid((G.L + s_x - 1) % G.L, s_y) + G.grid(s_x, (s_y+1) % G.L) + G.grid(s_x, (G.L + s_y-1) % G.L);
    dE = (2) * G.grid(s_x, s_y) * dE;// * G.grid(s_x, s_y); //dE = dE * (-2);

    //std::cout << dE <<  << std::endl;

    dE = exp(-dE / G.T);// * G.J * G.beta);

    epsilon += G.E;

    //std::cout << G.get_E() << "    " << G.get_M() << std::endl;

    current_magnetization = G.M;
    current_magnetization = std::sqrt(current_magnetization * current_magnetization);
    //std::cout << "CURRENT MAGNETIZATION: " << current_magnetization/(G.L * G.L) << std::endl;
    magnetization += current_magnetization;


    r = spinflip(generator);

    if (r <= dE) {
      G.flip_spin(s_x, s_y);
    }

  }

  epsilon = epsilon / (G.L*G.L);
  epsilon = epsilon / M;

  magnetization = magnetization / (M * G.L * G.L);

  std::cout << "THREAD SEED: " << thread_seed
            << ", AVG MAGNETIZATION: " << magnetization
            << ", EPSILON_MEAN: " << epsilon << ", TEMPERATURE: " << G.T << std::endl;


  return epsilon;
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
    //const int my_thread = omp_get_thread_num();

    # pragma omp for
    for (int i = 0; i < num_threads; i++){
      double T = T_min + i * dT;

      Grid G = Grid(L, T);

      auto thread_seed = MC_seed_generator();

      G.fill_grid(thread_seed);

      double epsilon_mean = MCMC(G, num_MC_cycles, thread_seed);
      /*std::cout << "THREAD SEED: " << thread_seed
          << ", EPSILON_MEAN: " << epsilon_mean << ", TEMPERATURE: " << T << std::endl;
      */

      //std::cout << my_thread << ", " << i << std::endl;
      /*ofile << setprecision(print_precision) << scientific << i << ", "
            << setprecision(print_precision) << scientific << epsilon_mean << std::endl;
      */
    }
  }

  /////////////////////


  /*arma::mat grid2;
  grid2.eye(L, L);
  grid2 = grid2*2 - 1;

  double E = 0.0;

  for (int i = 0; i < L; i ++){
    for (int j = 0; j < L; j ++){

      E = E + grid2(i, j) * (grid2((i+1) %  L, j) + grid2(i, (j+1) %  L));

    }
  }

  std::cout << -E/(L*L) << "    " << -E <<std::endl;
  */

  //std::cout << "analytical epsilon_mean: " << G.epsilon_mean() << std::endl;

  return 0;
}
