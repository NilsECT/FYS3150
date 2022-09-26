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

        int sim_trans = 0;


    public:
        void test_find_k_l();

        Jacobi(arma::mat &matrix);

        void solve(double tol=1e-10);

        arma::mat get_A();
        arma::mat get_eigvec();
        arma::vec get_eigval();

        void set_A(arma::mat matrix);

        int trans_count();

        int get_N();
        int get_k();
        int get_l();
        double get_max_val();
        int get_sim_trans();
};

#endif
