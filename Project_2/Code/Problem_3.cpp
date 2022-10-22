#include <armadillo>
#include <iostream>
#include "Jacobi.hpp"

using namespace std;


int main() {
    arma::mat b = arma::mat(2,2);
    Jacobi jacobi = Jacobi(b);
    jacobi.test_find_k_l();
}