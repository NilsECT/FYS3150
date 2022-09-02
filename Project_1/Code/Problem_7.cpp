#include <time.h>
#include <armadillo>
#include <iostream>
#include <fstream>

int main (){
    int n = 10;
    arma::vec v = arma::vec(n); 
    arma::vec g = arma::vec(n);   
    arma::vec a = arma::vec(n).fill(-1.); 
    arma::vec b = arma::vec(n).fill(2.); 
    arma::vec c = arma::vec(n).fill(-1.); 


    // Start measuring time
    clock_t t1 = clock();
    
    for     (int i = 1; i<n-1; i++){

    }




    // Stop measuring time
    clock_t t2 = clock();

    // Calculate the elapsed time.
    double duration_seconds = ((double) (t2 - t1)) / CLOCKS_PER_SEC;

    // ...
}