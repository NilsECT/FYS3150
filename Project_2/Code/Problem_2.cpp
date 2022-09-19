#include <armadillo>
#include <iostream>
#include <cmath>

#define pi 3.14159265

double lambda_anal(int i, int N){
    double eig = 2. + 2. * cos((2*i*pi)/(N+1));
    return eig;
}

arma::mat v_anal(int N) {
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

arma::mat find_A(int N) {
    arma::mat A = arma::mat(N,N).fill(0);

    for(int i = 0; i<N; i++){
        A(i,i) = 2;
    }

    for(int i = 0; i<N-1; i++){
        A(i,i+1) = -1;
        A(i+1,i) = -1;
    }

    return A;
}

bool compare(arma::mat A, arma::mat B, double tol) {
    arma::mat diff_mat = A - B;
    double diff = arma::mean(arma::mean(diff_mat));
    bool nice_diff = diff < tol;
    if (nice_diff) {
        return true;
    } else {
        return false;
    }
}


int main(){
    int N = 6;
    arma::mat A = find_A(N);

    std::cout << A << std::endl;
   
    
    arma::vec eigval; 
    arma::mat eigvec;
    arma::eig_sym(eigval, eigvec, A, "dc");

    std::cout << eigval << std::endl;
    std::cout << eigvec << std::endl;

    std::cout << v_anal(N) << std::endl;

    std::cout << "Numerical and analytical matrices are consistent: " << compare(v_anal(N), eigvec, 1e-5) << std::endl;
}