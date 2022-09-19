#include <armadillo>
#include <iostream>


//#define PI 3.14159265

//double lambda(int i, int N){
  //  double eig = 2.+2.*cmath::cos((2*i*PI)/(N+1));
    //return eig;
//}

int main(){
    int N = 6;
    arma::mat A = arma::mat(N,N).fill(0);
    for(int i = 0; i<N; i++){
        A(i,i) = 2;
    }
    for(int i = 0; i<N-1; i++){
        A(i,i+1) = -1;
        A(i+1,i) = -1;
    }
    std::cout << A << std::endl;
   
    
    arma::vec eigval; 
    arma::mat eigvec;
    arma::eig_sym(eigval, eigvec, A, "dc");

    //std::cout << eigval << std::endl;
    //std::cout << eigvec << std::endl;

}