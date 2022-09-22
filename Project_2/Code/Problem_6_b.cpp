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
    int N = 99;
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
    
    arma::mat output = arma::mat(N,6);
    output.col(0) = eigvec_sort.col(0);
    output.col(1) = eigvec_sort.col(1);
    output.col(2) = eigvec_sort.col(2);
    
    output.col(3) = eigvec.col(0);
    output.col(4) = eigvec.col(1);
    cout << "HOla" << endl;
    output.col(5) = eigvec.col(2);
    cout << "Hola" << endl;

    
    ofstream MyFile("Problem_6_b.txt"); // "(o)fstream" writes to file
    MyFile << output << endl; //What we want in the file
    MyFile.close(); //Wlose the file
}