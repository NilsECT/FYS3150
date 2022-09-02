#include <iostream>
#include <armadillo>

// messing around with Armadillo

int main(int argc, char* argv[])
{

  double N = 200; // number of points, which means N-1 steps to take
  int xmin = 0;
  int xmax = 1;
  double stepsize = (xmax - xmin)/(N-1); // "length" of each step

  arma::vec x = arma::vec(N);
  arma::vec v = arma::vec(N);

  for (int i = 0; i < N; i++) {
    x(i) = xmin + i*stepsize;
    v(i) = 1 - (1 - exp(-10.)) * x(i) - exp(-10.*x(i));
  }


  if (argc > 1) // if input from command line is provided, print vectors
  {
    std::cout << "x:" << std::endl << x;

    std::cout << "v:" << std::endl << v;

    std::cout << "Step size:" << stepsize << ", number of points N = " << N << std::endl;
  }

  return 0;
}
