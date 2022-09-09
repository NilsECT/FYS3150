
#include <iostream>
#include <armadillo>
#include <sstream>
#include <fstream>
#include <string>


int main() {

    std::string u_filename = "Problem_2.txt";
    std::string filename10 = "out10.txt";
    std::string filename100 = "out100.txt";
    std::string filename1000 = "out1000.txt";

    arma::mat u_mat, mat10, mat100, mat1000;
    u_mat.load(u_filename, arma::raw_ascii);
    mat10.load(filename10, arma::raw_ascii);
    mat100.load(filename100, arma::raw_ascii);
    mat1000.load(filename1000, arma::raw_ascii);

    arma::vec Delta10 = mat10.col(1);
    arma::vec Delta100 = mat100.col(1);
    arma::vec Delta1000 = mat1000.col(1);
    std::cout << Delta10 << std::endl;


    return 0;
}