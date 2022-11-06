#include "Grid.hpp"
#include <armadillo>

Grid::Grid(int L, double T, double J){
  this->L = L;
  this->T = T;
  this->J = J;
  this->N = L*L;
  this->beta = 1/T;//*1.3806 * std::pow(10, -23));

  arma::mat grid = arma::mat(L, L, arma::fill::zeros);
  this->grid = grid;

  this->E = this->get_E();
  this->M = this->get_M();
}

void Grid::flip_spin(int x, int y){
  double new_spin = (-1) * this->grid(x, y);
  this->grid(x, y) = new_spin;

  this->M = this->M + 2*new_spin;
  this->E = this->get_E();
}

double Grid::get_E(){
  double E = 0.0;
  //arma::mat *spins = &this->grid;

  for (int i = 0; i < this->L; i ++){
    for (int j = 0; j < this->L; j ++){

      E = E + -this->grid(i, j) * (this->grid((i+1) % this-> L, j) + this->grid(i, (j+1) % this-> L));

    }
  }
  //E = this-> J * E;
  return E;
}


double Grid::get_M(){
  double M = 0.0;

  for (int i = 0; i < this->L; i++){
    for (int j = 0; j < this->L; j++){
      M = M + this->grid(i, j);
    }
  }
  return M;
}

void Grid::fill_grid(int seed){
  std::mt19937 generator (130);
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

double Grid::epsilon_mean(){
  // Supposed to calculate analytical epsilon mean, but is defined wrong
  return 0.0;

  /*double eps = 0.0;

  for (int i = -this->N; i < this->N; i ++){
    eps = eps + i * exp(- this->beta * i);
  }

  eps = eps / (this->N * this->Z());
  return eps;*/
}

double Grid::Z(){
  double Z = 0.0;

  for (int i = -this->N; i < this->N; i ++){
    Z = Z + exp(- this->beta * i);
  }

  return Z;
}
