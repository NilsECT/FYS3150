#include "Particle.hpp"
#include "Penningtrap.hpp"
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
    this->force = arma::vec(3).fill(0); // Zero until we find it using find_coulomb_force()
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

arma::vec Particle::find_coulomb_force(std::vector<Particle> particles) {
    arma::vec E_ke = arma::vec(3).fill(0.); // Total electric field over ke.

    for (Particle p : particles){
        arma::vec r_diff = this->r -p.r;
        double r_norm = arma::norm(r_diff);
        double tol = 1e-6;
        if (r_norm<tol){
            continue;
        }
        double r_3 = std::pow(r_norm,3);
        E_ke = E_ke + r_diff/r_3;
    }
    arma::vec C= this->q*this->ke*E_ke;
    return C;
}

arma::vec Particle::find_E_field(double V_0, double d){
    arma::vec E = arma::vec(3).fill(0.);
    E(0) = -this->r(0)*V_0*d;
    E(1) = -this->r(1)*V_0*d;
    E(2) = 2*this->r(2)*V_0*d;
    return E;
}

arma::vec Particle::find_B_field(double B_0){
    arma::vec B = arma::vec(3).fill(0.);
    B(2) = B_0;
    return B;
}

void Particle::set_force(arma::vec F){
    //std::cout << "FFF " << F << std::endl;
    this->force = F;
    //std::cout << "this " << this->force << std::endl;
}

void Particle::print(){
    //std::cout << "q = " << this->q << std::endl;
    //std::cout << "m = " << this->m << std::endl;
    //std::cout << "ke = " << this->ke << std::endl;
    std::cout << "r = " << this->r << std::endl;
    //std::cout << "v = " << this->v << std::endl;
    //std::cout << "force = " << this->force << std::endl;
}