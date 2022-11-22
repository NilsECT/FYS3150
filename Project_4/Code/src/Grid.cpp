#include "Grid.hpp"
#include <armadillo>
#include <string>
#include <assert.h>

Grid::Grid(int L, double T, double J){
  this->L = L;    // Number of spins in each dimension.
  this->T = T;    // Units: [1/k_B]
  this->J = J;    // Units: [J]
  this->N = L*L;  // Total number of spins.
  this->beta = 1/T; // Thermodynamic beta.

  arma::mat grid = arma::mat(L, L, arma::fill::zeros);
  this->grid = grid;

  // Ratio between probability of new state and old state, for the
  // different obtainable values:
  this->dEs = {exp(-1 * this->J * (-8) / this->T), exp(-1 * this->J * (-4) / this->T), 1, exp(-1 * this->J * (4) / this->T), exp(-1 * this->J * (8) / this->T)};

  // Store energy and magnetization in initial state:
  this->E = this->get_E();
  this->M = this->get_M();

  // Initialize averages:
  this->epsilon = 0.0;
  this->epsilon_squared = 0.0;
  this-> m_squared = 0.0;
  this->m_abs = 0.0;
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

  for (int j = 0; j < this->L; j ++){
    arma::vec temp_col = this->grid.col(j);
    arma::vec temp_col_jj = this->grid.col(((j+1) == L) ? 0 : (j + 1)); // optimise prime

    for (int i = 0; i < this->L; i ++){   // Optimize prime
      E = E + -temp_col(i)* (temp_col(((i+1) == L) ? 0 : (i + 1)) + temp_col_jj(i));
    }
  }
  return this->J * E;
}

/**
 * Calculate the magnetization of the current state.
 *
 * @return M Magnetization of the system.
 */
double Grid::get_M(){
  double M = 0.0;

  for (int j = 0; j < this->L; j++){
    arma::vec temp_col = this->grid.col(j); // optimise prime

    for (int i = 0; i < this->L; i++){ // Optimize prime
      M = M + temp_col(i);
    }
  }
  return M;
}

/**
 * Fills the 2-dimensional grid with random spins using a
 * Mersenne Twister-generator and a uniform distribution function.
 * Calculate and initialize energy and magnetization values.
 *
 * @param seed Seed for the psuedo-random number generator.
 * @param random_config Whether the configuration is randomized (with the seed).
 */
void Grid::fill_grid(int seed, bool random_config){
  std::mt19937 generator (seed);
  std::uniform_int_distribution<int> dis(0, 1);

  if (random_config){
    int num;

    for (int j = 0; j < L; j ++){
      for (int i = 0; i < L; i ++){   // Optimize prime
        num = dis(generator);
        num = num*2 - 1;
        this->grid(i, j) = num; // I think we can optimise this, but unsure, we'd have to deal with references or pointers
      }
    }
  }
  else {
    this->grid.ones(this->L, this->L);
  }

  // Initialize averages:
  this->epsilon = 0.0;
  this->epsilon_squared = 0.0;
  this->m_squared = 0.0;
  this->m_abs = 0.0;

  this->E = this->get_E();
  this->M = this->get_M();

}

/**
 * Calculate heat capacity with current energy per spin values.
 *
 * @param none
 */
void Grid::compute_cv(){
  //assert this->epsilon != NULL;

  this->cv = (this->epsilon_squared - this->epsilon*this->epsilon);
  this->cv = this->N * this->cv / ( this->T * this->T);

}

/**
 * Calculate magnetic susceptibility with current
 * magnetization per spin values.
 *
 * @param none
 */
void Grid::compute_chi(){
  this->chi = (this->m_squared - this->m_abs*this->m_abs);
  this->chi = this->N * this->chi / (this->T);
}

void Grid::random_walk(int seed){
  int L = this->L;

  // Mersenne twister algorithm for generating random numbers:
  std::mt19937 generator (seed);
  // Distribution for choosing random grid indexes:
  std::uniform_int_distribution<int> dis(0, L-1);
  // Distribution for r:
  std::uniform_real_distribution<double> spinflip(0, 1);
  double r; // Random number to which we will compare the probability rati

  double epsilon = 0.0;
  double epsilon_squared = 0.0;
  double m_squared = 0.0;
  double m_abs = 0.0;

  double current_magnetization;
  int index;

  for (int i = 0; i < L*L; i++){
    int s_x = dis(generator); // Index of "chosen" spin in x-direction
    int s_y = dis(generator); // Index in y-direction

    // ? : is a conditional checks a condition () ? and takes a value True : False
    // that is faster than modulus

    // Compute change in energy [in units of 1/J]:
    double dE = this->grid(((s_x+1) == L) ? 0 : (s_x + 1), s_y)
              + this->grid(((s_x - 1) < 0) ? (L - 1) : (s_x - 1), s_y)
              + this->grid(s_x, ((s_y+1) == L) ? 0 : (s_y + 1))
              + this->grid(s_x, ((s_y - 1) < 0) ? (L - 1) : (s_y - 1));
    dE = 2 * this->grid(s_x, s_y) * dE;

    // Index corresponding to change in energy (to avoid multiple exp()-calls):
    index = dE / 4 + 2;
    dE = this->dEs.at(index);

    // Generate random float between 0 and 1:
    r = spinflip(generator);

    if (r <= dE) {
      this->flip_spin(s_x, s_y); // Flip the chosen spin
    }

  }

  this->epsilon = this->get_E();
  this->epsilon_squared = this->get_E() * this->get_E();
  this->m_abs = std::sqrt(this->M * this->M);
  this->m_squared = this->M * this->M;
}

void Grid::one_walk(int thread_seed){
  int L = this->L;

  std::mt19937 generator (thread_seed);
  // Distribution for choosing random grid indexes:
  std::uniform_int_distribution<int> dis(0, L-1);
  // Distribution for r:
  std::uniform_real_distribution<double> spinflip(0, 1);
  double r; // Random number to which we will compare the probability rati

  double current_magnetization;
  int index;

  int s_x = dis(generator); // Index of "chosen" spin in x-direction
  int s_y = dis(generator); // Index in y-direction

  // Compute change in energy [in units of 1/J]:
    double dE = this->grid(((s_x+1) == L) ? 0 : (s_x + 1), s_y)
              + this->grid(((s_x - 1) < 0) ? (L - 1) : (s_x - 1), s_y)
              + this->grid(s_x, ((s_y+1) == L) ? 0 : (s_y + 1))
              + this->grid(s_x, ((s_y - 1) < 0) ? (L - 1) : (s_y - 1));
    dE = 2 * this->grid(s_x, s_y) * dE;

  // Index corresponding to change in energy (to avoid multiple exp()-calls):
  index = dE / 4 + 2;
  dE = this->dEs.at(index);

  // Generate random float between 0 and 1:
  r = spinflip(generator);

  if (r <= dE) {
    this->flip_spin(s_x, s_y); // Flip the chosen spin
  }

  // what we want
  epsilon = this->E/this->N;
}
