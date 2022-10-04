#include <armadillo>
#include <iostream>
#include "Penningtrap.hpp"
#include "Particle.hpp"


int main(){
    double q = 1;
    double m = 1;
    arma::vec r_1 = arma::vec("0.0,1.0,5.0");
    arma::vec v_1 = arma::vec("1.0,2.0,3.0");
    Particle particle_1 = Particle(q, m, r_1, v_1);
    arma::vec r_2 = arma::vec("0.0,4.0,3.0");
    arma::vec v_2 = arma::vec("20.0,7.0,8.0");
    Particle particle_2 = Particle(q, m, r_2, v_2);
    std::vector<Particle> particles;
    particles.push_back(particle_1);
    particles.push_back(particle_2);

    double B_0 = 1;
    double V_0 = 1;
    double d = 1;
    Penningtrap penningtrap = Penningtrap(B_0, V_0, d, particles);

    penningtrap.find_force(true);

    for (Particle p : penningtrap.particles) {
        p.print();
    }
}