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
Penningtrap::Penningtrap(double B_0, double V_0, double d){
    this->B_0 = B_0;
    this->V_0 = V_0;
    this->d = d;
    std::vector<Particle> particles;
    this->particles = particles;
    this->num_particles_inside = 0; // Number of particles in the trap
}


/**
 * @brief 
 *
 * @param particle_interactions 
 *
 */
void Penningtrap::find_force(bool has_coloumb_force, bool has_E_field, bool has_B_field) {
    for (Particle &particle : this->particles) {

        // Initialize coulomb force, E-field and B-field to zero-vectors
        arma::vec C = arma::vec(3).fill(0);
        arma::vec E = arma::vec(3).fill(0);
        arma::vec B = arma::vec(3).fill(0);

        if (has_coloumb_force) {

            C = particle.find_coulomb_force(this->particles);
        }

        if (has_E_field) {
            // Calculate E-field if it is turned on
            E = particle.find_E_field(this->V_0, this->d);
        }

        if (has_B_field) {
            // Calculate B-field if it is turned on
            B = particle.find_B_field(this->B_0);
        }
        
        arma::vec L = particle.find_Lorentz_force(E, B);

        arma::vec F = C + L;
        std::cout << C << L << F << std::endl;
        particle.force = F;
    }
}

void Penningtrap::add_particle(Particle particle) {
    this->particles.push_back(particle);
    this->num_particles_inside++;
}

void Penningtrap::generate_particles(int N, double q, double m, int seed) {
    this->particles.clear();
    arma::vec r;
    arma::vec v;

    arma::arma_rng::set_seed(seed);

    for (int i = 0; i < N; i++) {
        r = arma::vec(3).randn() * 0.1 * this->d;
        v = arma::vec(3).randn() * 0.1 * this->d;
        Particle particle = Particle(q, m, r, v);
        this->particles.push_back(particle);
    }
    this->num_particles_inside += N;
}

/**
 * @brief 
 *
 * @param 
 *
*/
/*
void write_to_file() {
    std::vector<std::string> name =  {"x", "y", "z", "vx", "vy", "vz"};
    for (int i = 0; i<3;i++){
        std::ofstream name[i]_outfile;
        std::ofstream name[i+3]_outfile;
        name[i]_outfile.open(name[i] + ".txt", std::ios_base::app); // append instead of overwrite
        name[i+3]_outfile.open(name[i+3] + ".txt", std::ios_base::app); // append instead of overwrite
        for (Particle p : this->particles) {  // particle number
            name[i]_outfile << p.r[i] << "   "; 
            name[i+3]_outfile << p.v[i] << "   "; 

        }
        name[i]_outfile << "\n";
        name[i+3]_outfile << "\n";
    }
}
*/