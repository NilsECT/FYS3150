#include <armadillo>
#include <iostream>
#include "Jacobi.hpp"

using namespace std;

int main() {
    int N = 4;
    arma::mat A = arma::mat(N, N).randn();  
    /*int N = 4;
    arma::mat A = arma::mat(N, N).eye();
    A(3, 0) = 0.5;
    A(0, 3) = A(3, 0);
    A(1, 2) = -0.7;
    A(2, 1) = A(1, 2);*/
    // Symmetrize the matrix by reflecting the upper triangle to lower triangle
    A = arma::symmatu(A);  
    cout << A << endl;
    Jacobi jacobi = Jacobi(A);
    // cout << jacobi.get_k() << endl;
    // cout << jacobi.get_l() << endl;
    // cout << A(jacobi.get_k(), jacobi.get_l()) << endl;
    jacobi.solve();
    cout << "A:" << jacobi.get_A() << endl;
    cout << "S:" << jacobi.get_eigvec() << endl;
    
    arma::vec eigval;
    arma::mat eigvec;
    arma::eig_sym(eigval, eigvec, A);
    cout << eigval << endl;
    cout << eigvec << endl;
}