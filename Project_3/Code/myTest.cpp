
#include <armadillo>
#include <iostream>

int main() {
    int N = 1000000;
    arma::vec u = arma::vec(N).randn();

    std::cout << arma::any(u > 5) << std::endl;
}