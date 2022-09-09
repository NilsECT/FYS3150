#include <time.h>
#include <armadillo>
#include <iostream>
#include <fstream>


double f(double x){
    return 100*exp(-10*x);
}

int main (){
    double N = 1;


    for (int I = 1; I <= 7; I++){

        int n = (int) std::pow(10, I);
        double h = N/n;
        arma::vec x = arma::linspace(0,N,n+1);
        arma::vec g = arma::vec(n+1); 
        arma::vec gt = arma::vec(n+1);  // g-tilde
        arma::vec v = arma::vec(n+1); 
        arma::vec vt = arma::vec(n+1); // v-tilde
        
        v(0) = 0;
        v(n) = 0;

        g(1) = f(x(1))*h*h + v(0);
        g(n-1) = f(x(n-1))*h*h + v(n);
        for     (int i = 2; i<n-1; i++){
            g(i) = f(x(i))*h*h;
        }

        arma::vec a = arma::vec(n).fill(-1.); 
        arma::vec b = arma::vec(n).fill(2.); //b
        arma::vec bt = arma::vec(n); //b-tilde
        arma::vec c = arma::vec(n).fill(-1.); 


        // Start measuring time
        clock_t t1 = clock();

        gt(1) = g(1);
        bt(1) = b(1);
        for     (int i = 2; i<n; i++){
            bt(i) = b(i) - (a(i)/bt(i-1))*c(i-1);
            gt(i) = g(i) - (a(i)/bt(i-1))*gt(i-1);
        }
        
        //std::cout << "Hola" << std::endl;
        
        v(n-1) = gt(n-1)/bt(n-1);
        
        for     (int i=n-2; i>0; i--){
            v(i) = (gt(i) - c(i)*v(i+1))/bt(i);
        }
        //std::cout << v << std::endl;


        // Stop measuring time
        clock_t t2 = clock();

        // Calculate the elapsed time.
        double duration_seconds = ((double) (t2 - t1)) / CLOCKS_PER_SEC;
        //std::cout << duration_seconds << std::endl;

        arma::mat A = arma::mat(n+1, 2); //Defines matrix for easyer writing to file
        for  (int i = 0; i < n+1; i++){ //for-loop of every element n 
            A(i,0) = x(i); //Defines x as 1. column in the matrix
            //A(i,1) = u(x(i)); //Defines u valus as 2. column in the matrix
            A(i,1) = v(i);
        }

        // need to change I into string
        std::string name = "out_e"+std::to_string(I)+".txt";

        std::ofstream MyFile(name); // "(o)fstream" writes to file
        MyFile << A; //What we want in the file
        MyFile.close(); //Wlose the file

    }
    
    return 0;
}