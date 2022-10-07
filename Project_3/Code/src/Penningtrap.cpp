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
}


/**
 * @brief
 *
 * @param particle_interactions
 *
 */
void Penningtrap::find_force(bool particle_interactions) {
    for (Particle &particle : this->particles) {
        arma::vec C = arma::vec(3).fill(0);

        if (particle_interactions) {
            C = particle.find_coulomb_force(this->particles);
        }

        arma::vec L = particle.find_Lorentz_force(this->V_0, this->B_0, this->d);

        arma::vec F = C + L;
        std::cout << C << L << F << std::endl;
        particle.force = F;
    }
}

int Penningtrap::num_particles_inside() {
    int counter = 0;
    for (Particle particle : this->particles) {
        arma::vec r = particle.r;
        double r_norm = arma::norm(r);

        if (r_norm < this->d) {
            counter++;
        }
    }
    return counter;
}

void Penningtrap::add_particle(Particle particle) {
    this->particles.push_back(particle);
}

void Penningtrap::generate_particles(int N, double q, double m) {
    this->particles.clear();
    arma::vec r;
    arma::vec v;

    arma_rng::set_seed();

    for (int i = 0; i < N; i++) {
        r = arma::vec(3).randn() * 0.1 * this->d;
        v = arma::vec(3).randn() * 0.1 * this->d;
        Particle particle = Particle(q, m, r, v);
        this->particles.push_back(particle);
    }
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

/**
 * @brief Evolves the position and velocity of each particle one step in time using Forward Euler.
 *
 * @param Length of time step dt, booleans indicating which forces to consider.
 *
*/

void evolve_forwardeuler(double dt, bool particle_interaction) {

  for (Particle p : this->particles) {  // Iterate through particles
      p.r = p.r + p.v * dt;
      p.v = p.v + dt * p.find_force(particle_interaction) / p.m;
  }

}


/**
 * @brief Evolves the position and velocity of each particle one step in time using Runge-Kutta 4.
 *
 * @param Length of time step dt, booleans indicating which forces to consider.
 *
*/

void evolve_forwardeuler(double dt, bool particle_interaction) {

  for (Particle p : this->particles) {  // Iterate through particles

      // Store current position and velocity:
      r_temp = p.r;
      v_temp = p.v;

      arma::vec k1_r = dt * p.v;
      arma::vec k1_v = dt * p.find_force()/p.m;

      p.r = p.r + k1_r/2;
      p.v = p.v + k1_v/2;

      arma::vec k2_r = dt*p.v;
      arma::vec k2_v = dt*p.find_force()/p.m;

      p.r = r_temp + k2_r/2;
      p.v = v_temp + k2_v/2;

      arma::vec k3_r = dt*p.v;
      arma::vec k3_v = dt*p.find_force()/p.m;

      p.r = r_temp + k3_r;
      p.v = v_temp + k3_v;

      arma::vec k4_r = dt*p.v;
      arma::vec k4_v = dt*p.find_force()/p.m;

      // Update position and velocity:
      p.r = r_temp + (1/6.0) * (k1_r + 2.0*k2_r + 2.0*k3_r + k4_r);
      p.v = v_temp + (1/6.0) * (k1_v + 2.0*k2_v + 2.0*k3_v + k4_v);

  }

}
