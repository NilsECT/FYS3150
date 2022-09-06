#include <armadillo>
#include <iostream>

double f(double x){
  return 100*exp(-10*x);
}

double u(double x){
  return 1 - (1-exp(-10))*x - exp(-10*x);
}

int main(){
  double x_min = 0.;
  double x_max = 1.;

  int n = 1000;   // Length of g
  int m = n + 2;  // Length of x

  double h = x_max/(m-1); // Step length

  arma::vec x = arma::linspace(x_min, x_max, m);
  arma::vec v = arma::vec(n);
  arma::vec g = arma::vec(n);

  arma::vec b_tilde = arma::vec(n);
  arma::vec g_tilde = arma::vec(n);

  // Set initial value of b_tilde equal to b:
  b_tilde(0) = 2.;

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

  // Initialize an armadillo matrix in which we put the values:
  arma::mat A = arma::mat(m, 3);

  for (int i = 0; i < n; i++){
    A(i, 0) = x(i+1);     // Insert x
    A(i, 1) = u(x(i+1));  // Insert u(x)
    A(i, 2) = v(i);
  }

  // Print the matrix
  // (which can be written to a text file for
  // example using "./main > problem9.txt")
  std::cout << A << std::endl;
  
  return 0;

}
