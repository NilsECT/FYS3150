#include <time.h>
#include <armadillo>
#include <iostream>
#include <fstream>


int main (int argc, char* argv[]){

    arma::vec vecc = arma::logspace(1, 6, 6);

    std::cout << vecc << std::endl;

    return 0;
}