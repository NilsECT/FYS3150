#include <armadillo>
#include <iostream>
#include <cmath>

#define pi 3.14159265

double lambda_anal(int i, int N){
    // Calculates analytical eigenvalues using the formula in the problem text
    double eig = 2. + 2. * cos((2*i*pi)/(N+1));
    return eig;
}

arma::mat v_anal(int N) {
    // Calculates analytical eigenvectors using the formula in the problem text
    arma::mat v_anal = arma::mat(N, N);

    for (int i = 1; i <= N; i++) {
        arma::vec vi = arma::vec(N);
        
        for (int j = 1; j <= N; j++) {
            vi(j-1) = sin(j*i*pi/(N+1));
        }

        v_anal.col(i-1) = arma::normalise(vi);
    }

    return v_anal;
}

arma::mat find_A(int N, double h) {
    // Creates the NxN tridiagonal matrix A, with subdiagonal and superdiagonal
    // -1/h^2 and diagonal 2/h^2
    arma::mat A = arma::mat(N,N).fill(0);
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


int main(){
    int N = 6;
    double n_steps = N+1;
    double h = 1/n_steps;
    arma::mat A = find_A(N, h);
    
    arma::vec eigval; 
    arma::mat eigvec;
    arma::eig_sym(eigval, eigvec, A, "dc");

    double tol = 1e-7;

    // True if the numerical and analytical eigenvalues are consisten with each
    // other within the tolerance tol
    bool consistent = compare(v_anal(N), eigvec, tol);

    if (consistent) {
        std::cout << "Numerical and analytical matrices are consistent!" << std::endl;
    } else {
        std::cout << "Numerical and analytical matrices are NOT consistent!" << std::endl;
    }
}