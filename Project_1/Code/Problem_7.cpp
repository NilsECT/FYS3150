#include <time.h>
#include <armadillo>
#include <iostream>
#include <fstream>
using namespace std;

double u(double x){ //Function for u-values
    //Since we only have 4-desimals in the values we use floats without loosing precision
    return 1-(1-exp(-10))*x - exp(-10*x);
}

double f(double x){
    return 100*exp(-10*x);
}

int main (){
    int N = 1;
    int n = 10;
    double h = 1./n;
    arma::vec x = arma::linspace(0,N,n);
    arma::vec g = arma::vec(n);   
    for     (int i=0; i<n; i++){
        g(i) = f(x(i))*h*h;
    }

    arma::vec v = arma::vec(n); 
    arma::vec a = arma::vec(n).fill(-1.); 
    arma::vec b = arma::vec(n).fill(2.); 
    arma::vec c = arma::vec(n).fill(-1.); 


    // Start measuring time
    clock_t t1 = clock();
    
    for     (int i = 1; i<n; i++){
        b(i) = b(i) - (a(i)/b(i-1))*c(i);
        g(i) = g(i) - (a(i)/b(i-1))*g(i);

    }
    v(n-1) = g(n-1)/b(n-1);
    for     (int i=n-2; i>=0; i--){
        v(i) = (g(i) - c(i)*v(i+1))/b(i);
    }
    std::cout << g << std::endl;
    std::cout << v << std::endl;


    // Stop measuring time
    clock_t t2 = clock();

    // Calculate the elapsed time.
    double duration_seconds = ((double) (t2 - t1)) / CLOCKS_PER_SEC;
    std::cout << duration_seconds << std::endl;

    arma::mat A = arma::mat(n, 3); //Defines matrix for easyer writing to file
    for  (int i = 0; i < n; i++){ //for-loop of every element n 
        A(i,0) = x(i); //Defines x as 1. column in the matrix
        A(i,1) = u(x(i)); //Defines u valus as 2. column in the matrix
        A(i,2) = v(i);
    }

    ofstream MyFile("Problem_7.txt"); // "(o)fstream" writes to file
    MyFile << A; //What we want in the file
    MyFile.close(); //close the file
    return 0;
}