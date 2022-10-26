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
    this->r_origin = r;
    this->v_origin = v;
    this->ke = 1.38935333 * std::pow(10, 5); //1;
    this->force = arma::vec(3).fill(0); // Zero until we find it using the find_force()
    // method of the penningtrap this particle is in

    // Variables for Runge Kutta 4-method:
    this->runge_kutta_k = arma::zeros(3, 2); // Contains weighted sum of k-values
    this->v_temp = arma::zeros(3);  // Stores particle's velocity before evolving
    this->r_temp = arma::zeros(3);  // Stores particle's position
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
    E(0) = this->r(0);
    E(1) = this->r(1);
    E(2) = -2*this->r(2);
    return E*V_0/(d*d);
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

/**
 * @brief checks if the particle is labeled as outside
 *
 * @return true
 * @return false
 */
bool Particle::check_outside() {
    return outside;
}

/**
 * @brief Sets the particle's outside variable to true
 *
 */
void Particle::is_outside() {
    this->outside = true;
}
/**
 * @brief Resets the particle position and velocity to the initial values.
 * Also defines the particle to be inside the Penning trap to be true.
 * Empties all the temp vectors for velocity and position 
 * 
 */
void Particle::reset() {
    this->r = this->r_origin;
    this->v = this->v_origin;
    this->outside = false;
    this->force = arma::vec(3).fill(0); // Zero until we find it using the find_force()
    // method of the penningtrap this particle is in

    // Variables for Runge Kutta 4-method:
    this->runge_kutta_k = arma::zeros(3, 2); // Contains weighted sum of k-values
    this->v_temp = arma::zeros(3);  // Stores particle's velocity before evolving
    this->r_temp = arma::zeros(3);  // Stores particle's position
}
