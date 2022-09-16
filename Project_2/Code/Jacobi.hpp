#ifndef __Jacobi_hpp__
#define __Jacobi_hpp__
#include <armadillo>

class Jacobi
{
    protected:

        int k;
        int l;
        double max_val;
        arma::mat A;
        arma::mat S;

        void compute_tan();
        void compute_cos();
        void compute_sin();
        //void compute_trig();

        void update_S();
        void update_A();
        //void update();

        void find_k_l();
        void set_max();
        //void loop();
    
    public:

        Jacobi(arma::mat A);

        void solve();

        arma::mat get_A();
        arma::mat get_S();
        arma::vec get_eigvec();
        arma::vec get_egival();
        arma::mat get_eig();
};




#endif