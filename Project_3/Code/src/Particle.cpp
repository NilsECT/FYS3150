#include "Particle.hpp"
#include <armadillo>
#include <iostream>
#include <iomanip>
#include <unistd.h>
#include <assert.h>

/**
 * @brief Construct a new Particle object
 * 
 * @param q Charge of the particle.
 * @param m Mass of the particle.
 * @param r Position of the particle.
 * @param v Velocity of the particle.
 */
Particle::Particle(double q, double m, arma::vec r, arma::vec v){
    this->q = q;
    this->m = m;
    this->r = r;
    this->v = v;
    this->ke = 1.38935333 * std::pow(10, 5); //1;
    this->force = arma::vec(3).fill(0); // Zero until we find it using the find_force() method of the penningtrap this particle is in
}


/**
 * @brief Finds the coulomb force on this particle from other particles.
 * 
 * @param particles The particles we want to calculate the coulomb force from.
 * 
 * @return The sum of the coulomb force from each particle in particles as a 3D-vector.
 */
arma::vec Particle::find_coulomb_force(std::vector<Particle> particles) {
    arma::vec E_ke = arma::vec(3).fill(0.); // Total electric field over ke.

    for (Particle p : particles){
        arma::vec r_diff = this->r - p.r;
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

/**
 * @brief Finds the electric field on the particle from the electrodes of the
 * penningtrap the particle is in.
 * 
 * @param V_0 Electric potential of the electrodes of the penningtrap this
 * particle is in.
 * @param d Characteristic dimension of the penningtrap this particle is in.
 * 
 * @return The electric field on the particle as a 3D-vector.
 */
arma::vec Particle::find_E_field(double V_0, double d){
    arma::vec E = arma::vec(3).fill(0.);
    E(0) = this->r(0)*V_0/(d*d);
    E(1) = this->r(1)*V_0/(d*d);
    E(2) = -2*this->r(2)*V_0/(d*d);
    return E;
}

/**
 * @brief Finds the magnetic field on the particle from the penningtrap it is in.
 * 
 * @param B_0 The magnetic field strength of the penningtrap.
 * 
 * @return The magnetic field on the particle as a 3D-vector.
 */
arma::vec Particle::find_B_field(double B_0){
    arma::vec B = arma::vec(3).fill(0.);
    B(2) = B_0;
    return B;
}

/**
 * @brief Finds the total Lorentz force on the particle from the electric and
 * magnetic field.
 * 
 * @param E Electric field as a 3D-vector.
 * @param B Magnetic field as a 3D-vector.
 * 
 * @return The total Lorentz force from both fields as a 3D-vector.
 */
arma::vec Particle::find_Lorentz_force(arma::vec E, arma::vec B) {

    arma::vec FE = this->q * E; // Force from E-field
    arma::vec FB = this->q * arma::cross(this->v, B); // Force from B-field

    arma::vec lorentz_force = FE + FB; // Sum of forces from E-field and B-field

    return lorentz_force;
}

/**
 * @brief Prints information about the particle to the terminal, mainly used for
 * debugging and testing.
 * 
 */
void Particle::print(){
    std::cout << "q = " << this->q << std::endl;
    std::cout << "m = " << this->m << std::endl;
    std::cout << "ke = " << this->ke << std::endl;
    std::cout << "r = " << this->r << std::endl;
    std::cout << "v = " << this->v << std::endl;
    std::cout << "force = " << this->force << std::endl;
}
