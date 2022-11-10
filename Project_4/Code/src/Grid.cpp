#include "Grid.hpp"
#include <armadillo>

Grid::Grid(int L, double T, double J){
  this->L = L;    // Number of spins in each dimension.
  this->T = T;    // Units: [k_B]                          <--------- J/k_B?
  this->J = J;
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
 */
void Grid::fill_grid(int seed){
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
