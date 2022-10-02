#include "Particle.hpp"
#include <armadillo>
#include <iostream>
#include <iomanip>
#include <unistd.h>
#include <assert.h>

Particle::Particle(double q, double m, arma::vec r, arma::vec v){
    this->q = q;
    this->m = m;
    this->r = r;
    this->v = v;
    this->ke = 1;
    this->coulomb_force = arma::vec(3).fill(0); // Zero until we find it using find_coulomb_force()
}

arma::vec Particle::get_r(){
    return this->r;
}

arma::vec Particle::get_v(){
    return this->v;
}

void Particle::set_r(arma::vec r_new){
    this->r = r_new;
}

void Particle::set_v(arma::vec v_new){
    this->v = v_new;
}

double Particle::get_q(){
    return this->q;
}

double Particle::get_m(){
    return this->m;
}

void Particle::find_coulomb_force(std::vector<Particle> particles) {
    arma::vec E_ke = arma::vec(3).fill(0.); // Total electric field over ke.

    for (Particle p : particles){
        arma::vec r_diff = this->r -p.r;
        double r_norm = arma::norm(r_diff);
        double r_3 = std::pow(r_norm,3);
        E_ke = E_ke + r_diff/r_3;
    }

    this->coulomb_force = this->q*this->ke*E_ke;
}

arma::vec Particle::get_coulomb_force(){
    return this->coulomb_force;
}