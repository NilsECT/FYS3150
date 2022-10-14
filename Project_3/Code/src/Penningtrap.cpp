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
 * @param update_V_0
 *
 */
void Penningtrap::update_V_0(double f, double w, double t){
    this->V_0 = this->V_0*(1. + f * std::cos(w*t));
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
        v = arma::vec(3).randn() ;//* 0.01 * this->d;//arma::vec(3).randn() * 0.1 * this->d;
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

void Penningtrap::write_to_file(std::string h, std::string inter){
    std::string N = std::to_string(this->num_particles_inside);

    std::vector<std::string> names =  {"x", "y", "z", "vx", "vy", "vz"};

    for (int i = 0; i < 3; i++) {

        std::ofstream r_outfile;
        std::ofstream v_outfile;

        r_outfile.open(N+"_"+inter+"_"+h+"_"+names[i]+".txt", std::ios_base::app); // append instead of overwrite
        v_outfile.open(N+"_"+inter+"_"+h+"_"+names[i+3]+".txt", std::ios_base::app); // append instead of overwrite
        for (Particle p : this->particles) {  // particle number
            r_outfile << p.r(i) << "   ";
            v_outfile << p.v(i) << "   ";
        }
        r_outfile << "\n";
        v_outfile << "\n";
    }


}

/**
 * @brief Evolves the position and velocity of each particle one step in time using Forward Euler.
 *
 * @param Length of time step dt, booleans indicating which forces to consider.
 *
*/

void Penningtrap::evolve_forwardeuler(double dt, bool has_coulomb_force, bool has_E_field, bool has_B_field) {

  for (Particle& p : this->particles) {  // Iterate through particles

      this->find_force(has_coulomb_force, has_E_field, has_B_field);

      p.r = p.r + p.v * dt;
      p.v = p.v + dt * p.force / p.m;
  }

}


/**
 * @brief Evolves the position and velocity of each particle one step in time using Runge-Kutta 4.
 *
 * @param Length of time step dt, booleans indicating which forces to consider.
 *
*/

void Penningtrap::evolve_RK4(double dt, bool has_coulomb_force, bool has_E_field, bool has_B_field) {

  for (Particle& p : this->particles) {  // Iterate through particles

      // Store current position and velocity:
      arma::vec r_temp = p.r;
      arma::vec v_temp = p.v;

      this->find_force(has_coulomb_force, has_E_field, has_B_field);

      arma::vec k1_r = dt * p.v;
      arma::vec k1_v = dt * p.force/p.m;

      p.r = p.r + k1_r/2;
      p.v = p.v + k1_v/2;

      arma::vec k2_r = dt * p.v;
      arma::vec k2_v = dt * p.force/p.m;

      p.r = r_temp + k2_r/2;
      p.v = v_temp + k2_v/2;

      arma::vec k3_r = dt * p.v;
      arma::vec k3_v = dt * p.force/p.m;

      p.r = r_temp + k3_r;
      p.v = v_temp + k3_v;


      arma::vec k4_r = dt * p.v;
      arma::vec k4_v = dt * p.force/p.m;

      // Update position and velocity:
      p.r = r_temp + (1/6.0) * (k1_r + 2.0*k2_r + 2.0*k3_r + k4_r);
      p.v = v_temp + (1/6.0) * (k1_v + 2.0*k2_v + 2.0*k3_v + k4_v);

  }

}
