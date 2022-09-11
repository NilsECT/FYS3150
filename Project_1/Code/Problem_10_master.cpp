#include <time.h>
#include <armadillo>
#include <iostream>
#include <fstream>
#include <iomanip>

double f(double x){
  // Calculates the source term in our setup of the Poisson equation
    return 100*exp(-10*x);
}

double u(double x){
  // Calculates exact solution of our setup of the Poisson equation
  return 1 - (1-exp(-10))*x - exp(-10*x);
}

int main (){
    double N = 1; // Max value of x
    arma::vec n_vec = arma::logspace(1, 6, 6); // Vector with n= 10, 100, ..., 10^6
    int n_times = 100;  // No. of times to clock the time

    double width = 16;
    double prec = 2;

    std::cout << "#" << std::setw(width) << "n" << std::setw(width) << "dt_mean" << std::setw(width) << "dt_stddev";
    std::cout << std::setw(width) << "dt_mean_spec" << std::setw(width) << "dt_stddev_spec" << std::endl;

    for (int n : n_vec) {

        // Variables for the general algorithm
        double dt_mean;
        double dt_stddev;

        //Variables for the special algorithm
        double dt_mean_special;
        double dt_stddev_special;

        std::cout << std::setw(width) << std::setprecision(prec) << std::scientific << (double) n;
        double h = N/n; // Step length
        arma::vec x = arma::linspace(0,N,n+1);

        // Initialize arrays
        arma::vec g = arma::vec(n+1);   // RHS of matrix equation
        arma::vec gt = arma::vec(n+1);  // g-tilde
        arma::vec v = arma::vec(n+1);   // discretized solution

        v(0) = 0;
        v(n) = 0;

        g(0) = 0;
        g(n) = 0;

        // Boundary terms of RHS equation:
        g(1) = f(x(1))*h*h + v(0);
        g(n-1) = f(x(n-1))*h*h + v(n);

        // Fill RHS of equation:
        for (int i = 2; i < n-1; i++){
            g(i) = f(x(i))*h*h;
        }

        // Initialize arrays:
        arma::vec a = arma::vec(n).fill(-1.);
        arma::vec b = arma::vec(n).fill(2.); // b
        arma::vec bt = arma::vec(n); // b-tilde
        arma::vec c = arma::vec(n).fill(-1.);

        // For storing the measured times:
        arma::vec dt_vec = arma::vec(n_times);    // General algorithm
        arma::vec dt_vec_special = arma::vec(n_times);  // Specialized algorithm

        // Measures the time for the general algorithm
        for (int j = 0; j < n_times; j++) {

            // Start measuring time
            clock_t t1 = clock();

            // Forward substitution:
            gt(1) = g(1);
            bt(1) = b(1);
            for (int i = 2; i < n; i++ ){
                bt(i) = b(i) - (a(i)/bt(i-1))*c(i-1);
                gt(i) = g(i) - (a(i)/bt(i-1))*gt(i-1);
            }

            v(n-1) = gt(n-1)/bt(n-1);

            // Backward substitution
            for (int i = n-2; i > 0; i--){
                v(i) = (gt(i) - c(i)*v(i+1))/bt(i);
            }

            // Stop measuring time
            clock_t t2 = clock();

            double dt = ((double) (t2 - t1)) / CLOCKS_PER_SEC;

            dt_vec(j) = dt;

        }

        // Initialize arrays:
        arma::vec g_tilde = arma::vec(n+1); // g-tilde
        arma::vec v_ = arma::vec(n+1);      // v, for special algorithm
        arma::vec b_tilde = arma::vec(n+1); // b-tilde

        // Boundary terms:
        v_(0) = 0;
        v_(n) = 0;

        // Now for the special algorithm:
        for (int j = 0; j < n_times; j ++){

            // Start measuring time
            clock_t t1 = clock();

            g_tilde(1) = g(1);
            b_tilde(1) = 2;

            // Forward substitution:
            for (int i = 2; i < n; i ++) {
                g_tilde(i) = g(i) + g_tilde(i-1)/b_tilde(i-1);
                b_tilde(i) = 2 - 1/b_tilde(i-1);      // b_tilde_i
            }

            // Set last element of v:
            v_(n-1) = g_tilde(n-1) / b_tilde(n-1);

            // Backward substitution:
            for (int i = n-2; i > 0; i--) {
                v_(i) = (g_tilde(i) + v_(i+1)) / b_tilde(i);
            }

            // Stop measuring time:
            clock_t t2 = clock();

            // Elapsed time:
            double dt = ((double) (t2 - t1)) / CLOCKS_PER_SEC;

            dt_vec_special(j) = dt;
        }

        dt_mean = arma::mean(dt_vec); // Mean
        dt_stddev = arma::stddev(dt_vec); // Standard deviation

        // Mean and standard deviation, but for special algorithm:
        dt_mean_special = arma::mean(dt_vec_special); // Mean
        dt_stddev_special = arma::stddev(dt_vec_special); // Standard deviation

        // Calculate the elapsed time and print to terminal:
        std::cout << std::setw(width) << std::setprecision(prec) << std::scientific << dt_mean;
        std::cout << std::setw(width) << std::setprecision(prec) << std::scientific << dt_stddev;
        std::cout << std::setw(width) << std::setprecision(prec) << std::scientific << dt_mean_special;
        std::cout << std::setw(width) << std::setprecision(prec) << std::scientific << dt_stddev_special << std::endl;

    }

    return 0;
}
