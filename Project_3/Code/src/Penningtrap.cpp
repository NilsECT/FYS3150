#include "Penningtrap.hpp"
#include "Particle.hpp"
#include <armadillo>
#include <iostream>
#include <iomanip>
#include <unistd.h>
#include <assert.h>


 
/**
 * @brief What it does.
 *
 * @param What it needs.
 *
 */
Penningtrap::Penningtrap(double B_0, double V_0, double d, std::vector<Particle> particles){
    this->B_0 = B_0;
    this->V_0 = V_0;
    this->d = d;
    this->particles = particles;
}


/**
 * @brief 
 *
 * @param 
 *
 */
void Penningtrap::find_force(bool particle_interactions) {
    for (Particle &particle : this->particles) {
        arma::vec C = arma::vec(3).fill(0);

        if (particle_interactions) {
            C = particle.find_coulomb_force(this->particles);
        }

        arma::vec E = particle.find_E_field(this->V_0, this->d);
        arma::vec B = particle.find_B_field(this->B_0);

        arma::vec F = C + E + B;
        std::cout << C << E << B << F << std::endl;
        particle.set_force(F);
    }
}

int Penningtrap::num_particles_inside() {
    int counter = 0;
    for (Particle p : this->particles) {
        arma::vec r = p.get_r();
        double r_norm = arma::norm(r);

        if (r_norm < this->d) {
            counter++;
        }
    }
    return counter;
}

/**
 * @brief 
 *
 * @param 
 *
 */
/* 
arma::vec Penningtrap::find_E_field(){
    double V_d = this->V_d;
    this->E = arma::vec(3).fill(0.);
    this->E(0) = -this->r(this->i,0)*V_d;
    this->E(1) = -this->r(this->i,1)*V_d;
    this->E(2) = 2*this->r(this->i,2)*V_d;
}
*/
/**
 * @brief 
 *
 * @param 
 *
 */
 /*
void Penningtrap::B_field(){
    double B_0 = this->B_0;
    this->B = arma::vec(3).fill(0.);
    this->B(2) = B_0;
}
*/
/**
 * @brief 
 *
 * @param 
 *
 */
 /*
void Penningtrap::force(){
    arma::mat r = this->r
    arma::vec this->F = arma::vec(3).fill(0);
    int k = this->k
    this->F += this->E_field(r,k);
    this->F += this->B_field();
    this->F += this->columb_force(r,k);
}
*/
/**
 * @brief 
 *
 * @param 
 *
 */
 /*
arma::vec Penningtrap::forward_euler(arma::mat r,arma::mat v){
    int N = this->N;
    double h = 1./N;
    for (int k = 0; k<N; k++){
        arma::vec force = this->force(r,k)
         v(k) = v(k) + force*h;
         r(k) = r(k) + v(k)*h;

    }
    arma::mat retur = arma::mat(2,1);
    retur(0) = r;
    retur(1) = v;
    return retur; 
}
*/
/**
 * @brief 
 *
 * @param 
 *
 */
//arma::mat Penningtrap::rk4(){}