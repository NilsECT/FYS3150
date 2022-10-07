
#include <armadillo>
#include <iostream>
#include <string>

int main() {
    int N = 1000000;
    arma::vec u = arma::vec(N).randn();

    std::cout << arma::any(u > 5.6) << std::endl;

    std::string s = std::to_string(15);

    std::cout << "String to num " << s << " yey!" << std::endl;
}