#include "Penningtrap.hpp"
#include "Particle.hpp"
#include <armadillo>
#include <iostream>
#include <iomanip>
#include <unistd.h>
#include <assert.h>
#include <string>


 
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
void write_to_file() {
    std::vector<std::string> name =  {"x", "y", "z", "vx", "vy", "vz"};
    for (int i = 0; i<3;i++){
        std::ofstream name[i]_outfile;
        std::ofstream name[i+3]_outfile;
        name[i]_outfile.open(name[i]".txt", std::ios_base::app); // append instead of overwrite
        name[i+3]_outfile.open(name[i+3]".txt", std::ios_base::app); // append instead of overwrite
        for (Particle p : this->particles) {  // particle number
            name[i]_outfile << p.r[i] << "   "; 
            name[i+3]_outfile << p.v[i] << "   "; 

        }
        name[i]_outfile << "\n";
        name[i+3]_outfile << "\n";
    }
}


