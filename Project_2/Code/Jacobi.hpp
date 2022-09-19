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
        void update_A(arma::mat &Sm);

        void find_k_l();
        //void loop();
    
    public:

        Jacobi(arma::mat &matrix);

        void solve(double tol);

        arma::mat get_A();
        arma::mat get_S();
        arma::vec get_eigvec();
        arma::vec get_egival();
        arma::mat get_eig();

        arma::mat set_A();
};




#endif