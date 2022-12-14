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

bool compare(arma::mat A, arma::mat B, double tol) {
    // Compares two square matrices A and B to check if they are the same,
    // within a tolerance tol. Returns true if all elements in A are equal to
    // the corresponding elements in B. This method will be used for comparing
    // the analytical and the numerical eigenvectors.

    arma::SizeMat mat_size = arma::size(A);
    int N = mat_size(0);

    for (int i=0; i < N; i++) {
        for (int j=0; j < N; j++) {
            bool diff_minus = std::abs(A(i,j)-B(i,j)) > tol;
            
            // We add this part so that the method still returns true if two
            // eigenvectors are equal, but with opposite signs.
            bool diff_plus = std::abs(A(i,j)+B(i,j)) > tol;

            if (diff_minus && diff_plus) {
                return false;
            }
        }
    }
    return true;
}

int main() {
    int N = 6;
    arma::mat A = find_A(N);  
    // Symmetrize the matrix by reflecting the upper triangle to lower triangle
    A = arma::symmatu(A);  
    //cout << A << endl;
    Jacobi jacobi = Jacobi(A);
    // cout << jacobi.get_k() << endl;
    // cout << jacobi.get_l() << endl;
    // cout << A(jacobi.get_k(), jacobi.get_l()) << endl;
    jacobi.solve();

    arma::mat Aj = jacobi.get_A();
    arma::mat Sj = jacobi.get_eigvec();

    //cout << "A:" << Aj << endl;
    //cout << "S:" << Sj << endl;
    
    arma::vec eigval;
    arma::mat eigvec;
    arma::eig_sym(eigval, eigvec, A);
    //cout << eigval << endl;
    //cout << eigvec << endl;


    // Setting new vectors for sorting and plotting 
    arma::vec eigval_sort = arma::diagvec(Aj);
    arma::mat eigvec_sort = Sj;

    // sorting
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

    if (consistent) {
        std::cout << "Numerical and analytical matrices are consistent!" << std::endl;
    } else {
        std::cout << "Numerical and analytical matrices are NOT consistent!" << std::endl;
    }
}