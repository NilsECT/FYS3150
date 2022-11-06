#ifndef __Grid__
#define __Grid__
#include <armadillo>

class Grid{

public:
  int L;
  double T;
  int N;
  double J;
  double beta;

  double E;
  double M;

  Grid(int L, double T, double J = 1.0);
  arma::mat grid;

  double get_E();
  double get_M();
  double Z();
  
  void fill_grid(int seed = 137);
  void flip_spin(int x, int y);
};


#endif
