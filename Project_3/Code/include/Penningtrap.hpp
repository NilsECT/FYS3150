#ifndef __Penningtrap__
#define __Penningtrap__
#include <armadillo>
#include "Particle.hpp"
 
class Penningtrap{
    private:

    double B_0;
    double V_0;
    double d;
    
    public:
    std::vector<Particle> particles;
    Penningtrap(double B_0, double V_0, double d, std::vector<Particle> particles);

    void find_force();

    //arma::mat forward_euler(arma::mat r, arma::mat v);

    //arma::mat rk4(arma::mat r,arma::mat v,double m, double q);
};

#endif