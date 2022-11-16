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

  double epsilon;
  double epsilon_squared;
  double m_abs;
  double m_squared;

  double cv;
  double chi;

  double E;
  double M;

  arma::vec dEs;

  Grid(int L, double T, double J = 1.0);
  arma::mat grid;

  double get_E();
  double get_M();

  void fill_grid(int seed = 137, bool random_config = true);
  void flip_spin(int x, int y);

  void random_walk(int num_MC_cycles, int thread_seed);
  void compute_cv();
  void compute_chi();

  void one_walk(int thread_seed = 137);
};


#endif
