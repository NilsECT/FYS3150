#ifndef __Penningtrap__
#define __Penningtrap__
#include <armadillo>
#include <Particle.hpp>

class Penningtrap{
    private:

    double B_0;
    double V_0;
    double d;
    std::vector<Particle> particles;
    
    public:
    Penningtrap(double B_0, double V_0, double d, std::vector<Particle> particles);

    //arma::vec find_E_field();
    //arma::vec find_B_field();
    void find_coulomb_force();

    //arma::mat forward_euler(arma::mat r, arma::mat v);

    //arma::mat rk4(arma::mat r,arma::mat v,double m, double q);


};

#endif