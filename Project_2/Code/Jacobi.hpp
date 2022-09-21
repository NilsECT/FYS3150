#ifndef __Jacobi_hpp__
#define __Jacobi_hpp__
#include <armadillo>

class Jacobi
{
    private:

        int N;
        double tau;
        double tan;
        double cos;
        double sin;
        int sim_trans = 0;

        void test_find_k_l();
        bool is_valid_double(double x);

    protected:

        int k;
        int l;
        double max_val;
        arma::mat A;
        arma::mat S;

        void compute_tau();
        void compute_tan();
        void compute_cos();
        void compute_sin();
        void compute_trig();

        void update_S();
        void update_A();

        void find_k_l();
        //void loop();



    public:

        Jacobi(arma::mat &matrix);

        int solve(double tol=1e-10);

        arma::mat get_A();
        arma::mat get_eigvec();
        arma::vec get_eigval();

        void set_A(arma::mat matrix);

        int trans_count();

        int get_N();
        int get_k();
        int get_l();
        double get_max_val();
};

/**
 * Can be tested through the following:
 * 
 * compile and run:
 * 
 * g++ test.cpp Jacobi.cpp -I ./ -o test.exe -larmadillo -llapack -lblas
 * 
 * ./test.exe
 * 
#include <armadillo>
#include <iostream>
#include "Jacobi.hpp"

using namespace std;

int main() {
    int N = 4;
    arma::mat A = arma::mat(N, N).eye();
    A(3, 0) = 0.5;
    A(0, 3) = A(3, 0);
    A(1, 2) = -0.7;
    A(2, 1) = A(1, 2);
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
 * 
 */


#endif
