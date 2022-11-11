#include "Grid.hpp"
#include <armadillo>

Grid::Grid(int L, double T, double J){
  this->L = L;    // Number of spins in each dimension.
  this->T = T;    // Units: [1/k_B]
  this->J = J;    // Units: [J]
  this->N = L*L;  // Total number of spins.
  this->beta = 1/T; // Thermodynamic beta.

  arma::mat grid = arma::mat(L, L, arma::fill::zeros);
  this->grid = grid;

  // Store energy and magnetization in initial state:
  this->E = this->get_E();
  this->M = this->get_M();
}

/**
 * Flip one of the spins in the grid.
 * Returns nothing, but makes changes to grid, E, and M.
 *
 * @param x First index of spin that will flip.
 * @param y Second index of spin that will flip.
 */
void Grid::flip_spin(int x, int y){
  double new_spin = (-1) * this->grid(x, y);
  this->grid(x, y) = new_spin;

  // Update energy and magnetization:
  this->M = this->M + 2*new_spin;
  this->E = this->get_E();
}

/**
 * Calculate the energy of the current state of the grid.
 *
 * @return E Energy (in units of J).
 */
double Grid::get_E(){
  double E = 0.0;

  for (int i = 0; i < this->L; i ++){
    for (int j = 0; j < this->L; j ++){
      E = E + -this->grid(i, j) * (this->grid((i+1) % this-> L, j) + this->grid(i, (j+1) % this-> L));
    }
  }
  return E;
}

/**
 * Flip one of the spins in the grid.   <----?
 *
 * @return M Magnetization of the system.
 */
double Grid::get_M(){
  double M = 0.0;

  for (int i = 0; i < this->L; i++){
    for (int j = 0; j < this->L; j++){
      M = M + this->grid(i, j);
    }
  }
  return M;
}
/**
 * Fills the 2-dimensional grid with random spins using a
 * Mersenne Twister-generator and a uniform distribution function.
 *
 * @param seed Seed for the psuedo-random number generator.
 * @param random_config Whether the configuration is randomized (with the seed).
 */
void Grid::fill_grid(int seed, bool random_config){
  std::mt19937 generator (130);  // <------- seed istedet for 130?
  std::uniform_int_distribution<int> dis(0, 1);

  int num;

  for (int i = 0; i < L; i ++){
    for (int j = 0; j < L; j ++){
      num = dis(generator);
      num = num*2 - 1;
      this->grid(i, j) = num;
    }
  }
}

double Grid::Z(){
  double Z = 0.0;

  return Z;
}

void Grid::MCMC(int num_MC_cycles, int thread_seed){
  int L = this->L;

  std::mt19937 generator (thread_seed);
  std::uniform_int_distribution<int> dis(0, L-1);
  std::uniform_real_distribution<double> spinflip(0, 1);
  double r;

  double epsilon = 0.0;
  double epsilon_squared = 0.0;
  double m_squared = 0.0;
  double m_abs = 0.0;
  double current_magnetization;

  //std::vector<double> dEs = {8, 4, 0};

  for (int i = 0; i < num_MC_cycles; i++){
    int s_x = dis(generator);
    int s_y = dis(generator);


    // Compute change in energy:
    double dE = this->grid((s_x+1) % L, s_y)
        + this->grid((L + s_x - 1) % L, s_y)
        + this->grid(s_x, (s_y+1) % L)
        + this->grid(s_x, (L + s_y-1) % L);
    dE = 2 * this->grid(s_x, s_y) * dE;

    //std::cout << dE << "    " << exp(-dE / G.T) << std::endl;

    // Ratio between probability of new state and old state:
    dE = exp(-dE / this->T);

    // Add current energy and magnetization values to average sum:
    epsilon += this->E;
    epsilon_squared = epsilon_squared + this->E * this->E;

    current_magnetization = this->M;
    m_abs = m_abs + std::sqrt(current_magnetization * current_magnetization);
    m_squared = m_squared + current_magnetization * current_magnetization;

    r = spinflip(generator);

    if (r <= dE) {
      this->flip_spin(s_x, s_y);
    }

  }

  // Divide averages by number of Monte Carlo cycles (and number of spins):
  epsilon = epsilon / (L * L * num_MC_cycles);
  epsilon_squared = epsilon_squared / (num_MC_cycles * std::pow(L, 4));
  m_abs = m_abs / (num_MC_cycles * L * L);
  m_squared = m_squared / (num_MC_cycles * std::pow(L, 4));

  this->epsilon = epsilon;
  this->epsilon_squared = epsilon_squared;
  this->m_abs = m_abs;
  this->m_squared = m_squared;

  std::cout << "THREAD SEED: " << thread_seed
            << ", AVG MAGNETIZATION: " << m_abs
            << ", EPSILON_MEAN: " << epsilon << ", TEMPERATURE: " << this->T << std::endl;

  //std::vector<double> averages = {epsilon, epsilon_squared, m_abs, m_squared};
  //return averages;
}
