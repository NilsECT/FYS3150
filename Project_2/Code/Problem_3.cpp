#include <armadillo>
#include <iostream>
#include "Jacobi.hpp"

using namespace std;

int main() {

    // Generating the test matrix:
    int N = 4;
    arma::mat A = arma::mat(N, N).eye();
    A(3, 0) = 0.5;
    A(0, 3) = A(3, 0);
    A(1, 2) = -0.7;
    A(2, 1) = A(1, 2);
    Jacobi jacobi = Jacobi(A);
   
    jacobi.solve();
    arma::mat Aj = jacobi.get_A();
    arma::mat Sj = jacobi.get_eigvec();

    arma::vec eigval_sort = arma::diagvec(Aj);
    arma::mat eigvec_sort = Sj;
    
    arma::vec eigval;
    arma::mat eigvec;
    arma::eig_sym(eigval, eigvec, A);
    

    double tol = 1e-7; //Defining a tollerance for where we are happy with the matrix.

    // True if the numerical and analytical eigenvalues are consisten with each
    // other within the tolerance tol
    bool consistent = compare(v_anal(N), eigvec, tol);

    // Eig_sym retunes the eigen vectors in an order so we sort the eigen vectors in order form smalest to largest eigenvalue.
    // sorting:
    for (int i=0; i<N-1; i++){
        double min = eigval_sort(i);
        int min_ind = i;
        for( int j = i; j<N; j++){
            if (min > eigval_sort(j)){

                min = eigval_sort(j);
                min_ind = j;
            }
        }

        eigval_sort(min_ind) = eigval_sort(i);
        eigval_sort(i) = min;
        eigvec_sort.swap_cols(i, min_ind);
    }

    double tol = 1e-7;

    // True if the numerical and analytical eigenvalues are consisten with each
    // other within the tolerance tol
    bool consistent = compare(eigvec_sort, eigvec, tol);


    // Printing a statment telling us if the matrix is correct to some tollerance:
    if (consistent) {
        std::cout << "Numerical and analytical matrices are consistent!" << std::endl;
    } else {
        std::cout << "Numerical and analytical matrices are NOT consistent!" << std::endl;
    }
}