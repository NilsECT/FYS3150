#include <armadillo>
#include <iostream>
//#include "Penningtrap.hpp"
#include "Particle.hpp"


int main(){
    double q = 1;
    double m = 1;
    arma::vec r = arma::vec("0.0,1.0,5.0");
    arma::vec v = arma::vec("1.0,2.0,3.0");
    Particle particle_1 = Particle(q, m, r, v);


    std::cout << "r " << particle_1.get_r() << std::endl;
    std::cout << "v " << particle_1.get_v() << std::endl;
    std::cout << "q " << particle_1.get_q() << std::endl;
    std::cout << "m " << particle_1.get_m() << std::endl;
    r = arma::vec("2.0,3.0,10.0");
    particle_1.set_r(r);
    v = arma::vec("20.0,7.0,8.0");
    particle_1.set_v(v);
    std::cout << "r " << particle_1.get_r() << std::endl;
    std::cout << "v " << particle_1.get_v() << std::endl;

    std::cout << "c " << particle_1.get_coulomb_force() << std::endl;
}