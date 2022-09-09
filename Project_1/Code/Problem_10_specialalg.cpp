#include <armadillo>
#include <iostream>
#include <time.h>

double f(double x){
  return 100*exp(-10*x);
}

double u(double x){
  return 1 - (1-exp(-10))*x - exp(-10*x);
}

int main(int argc, char* argv[]){


  /*int n = atoi(argv[1]);  // Number of steps in g, v*/
  int N = atoi(argv[1]);  // Number of times to "clock" the algorithm
  arma::vec n_vec = arma::logspace(1, 6, 6); // Vector with n= 10, 100, ..., 10^6

  int n = 10;  // Initial number of steps

  double width = 14;
  double prec = 2;

  std::cout << "#" << std::setw(width) << "n" << std::setw(width) << "dt_mean" << std::setw(width) << "dt_stddev" << std::endl;

  for (int n : n_vec) {

    double dt_mean;
    double dt_stddev;

    std::cout << std::setw(width) << std::setprecision(prec) << std::scientific << (double) n;

    int m = n + 2;  // Length of x

    double x_min = 0.;
    double x_max = 1.;

    double h = x_max/(m-1); // Step length

    // Remember that x and v* have length m,
    // v, g have length n
    arma::vec x = arma::linspace(x_min, x_max, m);
    arma::vec v = arma::vec(n);
    arma::vec g = arma::vec(n);

    arma::vec b_tilde = arma::vec(n);
    arma::vec g_tilde = arma::vec(n);


    // Set initial value of b_tilde equal to b:
    b_tilde(0) = 2.;

    /////////
    //    std::vector<double> times(N);
    arma::vec dt_vec = arma::vec(N);

    for (int j = 0; j < N; j ++){

      // Start measuring time
      clock_t t1 = clock();

      // Forward substitution:
      for (int i = 1; i < n; i ++) {
        // RHS of matrix equation:
        g(i) = f(x(i+1))*h*h;

        g_tilde(i) = g(i) + g_tilde(i-1)/b_tilde(i-1);
        b_tilde(i) = 2 - 1/b_tilde(i-1);      // b_tilde_i
      }

      // Set last element of v:
      v(n-1) = g_tilde(n-1) / b_tilde(n-1);

      // Back substitution:
      for (int i = n-2; i > -1; --i){
        v(i) = (g_tilde(i) + v(i+1)) / b_tilde(i);
      }

      // Stop measuring time:
      clock_t t2 = clock();

      // Add to list of times:
      double dt = ((double) (t2 - t1)) / CLOCKS_PER_SEC;

      dt_vec(j) = dt;
    }

    dt_mean = arma::mean(dt_vec); // Mean
    dt_stddev = arma::stddev(dt_vec); // Standard deviation
    // Calculate the elapsed time.

    std::cout << std::setw(width) << std::setprecision(prec) << std::scientific << dt_mean;
    std::cout << std::setw(width) << std::setprecision(prec) << std::scientific << dt_stddev << std::endl;

  }

  return 0;

}
