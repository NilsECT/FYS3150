#include <armadillo>
#include <iostream>
#include <fstream>
using namespace std;

float u(float x){ //Function for u-values
    //Since we only have 4-desimals in the values we use floats without loosing precision
    return 1-(1-exp(-10))*x - exp(-10*x);
}
int main(){
    //This coes creates to array. One linspace form 0 to N, and one as a function of the first.
    //It then writes to both arrays to a matrix. Last it writes the matrix to a txt-file.
    int N = 1; //lengt of x-array
    int n = 10000; //number og slices
    arma::vec x = arma::linspace(0,N,n); //Defines x as linspace form 0 to N wit n slices
    arma::mat A = arma::mat(n, 2); //Defines matrix for easyer writing to file
    for  (int i = 0; i < n; i++){ //for-loop of every element n 
        A(i,1) = u(x(i)); //Defines u valus as 2. column in the matrix
        A(i,0) = x(i); //Defines x as 1. column in the matrix

    }

    ofstream MyFile("Problem_2.txt"); // "(o)fstream" writes to file
    MyFile << A; //What we want in the file
    MyFile.close(); //Wlose the file
    
    return 0;
}