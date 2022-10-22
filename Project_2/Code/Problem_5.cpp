#include <armadillo>
#include <iostream>
#include "Jacobi.hpp"

using namespace std;

arma::mat find_A(int N) {
    // Creates the NxN tridiagonal matrix A, with subdiagonal and superdiagonal
    // -1/h^2 and diagonal 2/h^2
    arma::mat A = arma::mat(N,N).fill(0);
    double h = 1/(N+1.);
    double k = 1/(h*h);

    for(int i = 0; i<N; i++){
        A(i,i) = 2*k;
    }

    for(int i = 0; i<N-1; i++){
        A(i,i+1) = -1*k;
        A(i+1,i) = -1*k;
    }

    return A;
}

int main() {
    // Defining a matrix for every N form 2 to reps and adding number of transformation 
    // per N and outing it to a txt file.
    int reps = 100;
    arma::vec num = arma::vec(reps);
    arma::vec i_values = arma::linspace(2, reps+1, reps);

    for (int i : i_values){
        int N = i;
        arma::mat A = find_A(N);  
        Jacobi jacobi = Jacobi(A);
        jacobi.solve();
        int num_trans = jacobi.get_sim_trans();
        num(i-2) = num_trans;
    }

    arma::mat output = arma::mat(reps, 2);
    output.col(0) = i_values;
    output.col(1) = num;

    ofstream MyFile("Problem_5.txt"); // "(o)fstream" writes to file
    MyFile << output << endl; //What we want in the file
    MyFile.close(); //Wlose the file
}