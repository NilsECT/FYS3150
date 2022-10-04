#ifndef __Particle__
#define __Particle__
#include <armadillo>
 
class Particle{
    private:
    double q;
    double m;
    arma::vec r;
    arma::vec v;
    arma::vec force;
    
    double ke;

    public:
    Particle(const double q, const double m, arma::vec r, arma::vec v);
    arma::vec get_r();
    arma::vec get_v();
    void set_v(arma::vec v_new);
    void set_r(arma::vec r_new);
    double get_q();
    double get_m();
    arma::vec find_coulomb_force(std::vector<Particle> particles);
    arma::vec find_E_field(double V_0, double d);
    arma::vec find_B_field(double B_0);
    void set_force(arma::vec F);
    void print();
}; 

#endif