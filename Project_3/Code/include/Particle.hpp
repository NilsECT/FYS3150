#ifndef __Particle__
#define __Particle__
#include <armadillo>
 
class Particle{

    public:

    double q;
    double m;
    arma::vec force;
    double ke;

    arma::vec r;
    arma::vec v;

    Particle(const double q, const double m, arma::vec r, arma::vec v);

    arma::vec find_coulomb_force(std::vector<Particle> particles);
    arma::vec find_E_field(double V_0, double d);
    arma::vec find_B_field(double B_0);
    arma::vec find_Lorentz_force(double V_0, double B_0, double d);

    void print();
}; 

#endif