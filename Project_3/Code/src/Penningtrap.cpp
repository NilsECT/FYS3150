#include "Penningtrap.hpp"
#include "Particle.hpp"
#include <armadillo>
#include <iostream>
#include <iomanip>
#include <unistd.h>
#include <assert.h>
#include <string>



/**
 * @brief Constructor for creating an object of Penningtrap. Assigns B_0, V_0
 * and d to the member variables. Also creates an std::vector with all the
 * particles of the trap, and initializes a counter for the number of particles
 * in the trap to 0.
 *
 * @param B_0 Magnetic field strength of the trap
 * @param V_0 Electric potential of the electrodes in the trap
 * @param d Characteristic dimension of the penningtrap
 *
 */
Penningtrap::Penningtrap(double B_0, double V_0, double d){
  this->B_0 = B_0;
  this->V_0 = V_0;
  this->d = d;
  std::vector<Particle> particles;
  this->particles = particles;
  this->num_particles_inside = 0; // Number of particles in the trap
  this->num_particles_out = 0;
}


/**
 * @brief Goes through each particle object and computes the forces working on
 * them from coulomb interactions, the electric field and the magnetic field. If either
 * has_coulomb_force, has_E_field or has_B_field is set to false, the
 * method will not calculate the contribution from that force.
 *
 * @param has_coulomb_force Ignores the coulomb force between particles if set to false.
 * @param has_E_field Ignores the force contribution from the electric field if
 * set to false.
 * @param has_B_field Ignores the force contribution from the magnetic field if
 * set to false.
 *
 */
void Penningtrap::find_force(bool has_coulomb_force, bool has_E_field, bool has_B_field, bool func_V, double f, double w, double ti) {

    int i = 0;

    // Contains Coulomb force for each particle:
    arma::mat C = arma::zeros(3, this->num_particles_inside);

    for (Particle &particle : this->particles) {

        // skips finding the force if the particle is outside
        if (particle.outside) {
            continue;
        }

        // Initialize E-field and B-field to zero-vectors
        arma::vec E = arma::vec(3).fill(0.);
        arma::vec B = arma::vec(3).fill(0.);

        if (has_coulomb_force) {

            int j = i;

            for (Particle &particle_i : this->particles){

              // excludes the particles which are outside
              if (particle_i.outside){
                continue;
              }

              arma::vec r_diff = particle.r - particle_i.r_coulomb;
              double r_norm = std::sqrt(r_diff(0)*r_diff(0) + r_diff(1)*r_diff(1) + r_diff(2)*r_diff(2));//arma::norm(r_diff);
              double tol = 1e-3;
              if (r_norm<tol){
                continue;
              }

              double r_3 = r_norm*r_norm*r_norm;

              C.col(i) = C.col(i) + r_diff/r_3;
              C.col(j) = -C.col(j) + r_diff/r_3;

              j++;
            }
        }

        if (has_E_field) {

            if (particle.outside) {
                continue;
            }
            else if (func_V) {
                E = particle.find_E_field(this->V_0 * (1 + std::cos(w*ti)), this->d);
            }
            else {
                // Calculate E-field if has_E_field is set to true
                E = particle.find_E_field(this->V_0, this->d);
            }
        }

        if (has_B_field) {
            if (!particle.outside) {
                // Calculate B-field if has_B_field is set to true
                B = particle.find_B_field(this->B_0);
            }
            else {
                continue;
            }

        }

        arma::vec L = particle.find_Lorentz_force(E, B); // Calculate the Lorentz force from the electric and magnetic field.

        arma::vec F = C.col(i) + L; // The total force is the sum of the coulomb force and the Lorentz force
        particle.force = F; // Updates the particles force to the one we calculated here

        i = i + 1;
    }


}

/**
 * @brief Adds a single particle to the Penningtraps std::vector of particles,
 * and increases the particle-counter by one.
 *
 * @param particle The particle to be added to the trap.
 */
void Penningtrap::add_particle(Particle &particle) {
    this->particles.push_back(particle); // Add particle to particles
    this->num_particles_inside++; // Increase the number of particles in the trap
}

/**
 * @brief Empties the std::vector containing particles
 *
 */
void Penningtrap::clear_particles() {
    this->particles.clear();
    this->num_particles_inside = 0;
}

/**
 * @brief Generates N particles randomly distributed in the trap, with randomly
 * distributed velocities (both normal distributions). If there are any other
 * particles in the trap before this method is called, these will be overwritten
 * by the new particles. This method also sets the particle counter to N.
 *
 * @param N Number of particles to be generated.
 * @param q Charge of the particles.
 * @param m Mass of the particles.
 * @param seed Seed for the random number generator.
 */
void Penningtrap::generate_particles(int N, double q, double m, int seed) {
    this->particles.clear();
    arma::vec r;
    arma::vec v;

    arma::arma_rng::set_seed(seed);

    for (int i = 0; i < N; i++) {
        r = arma::vec(3).randn() * 0.1 * this->d;
        v = arma::vec(3).randn() * 0.1 * this->d; // must be times 0.1 and d according to the problem text.
        Particle particle = Particle(q, m, r, v);
        this->particles.push_back(particle);
    }
    this->num_particles_inside = N; // Updates the particle counter
}

/**
 * @brief Writes position and velocity data to .txt-files.
 *
 * @param dt_str time_step as a string (used for naming the outfile)
 * @param has_coulomb_str "y" if the simulation had used coulomb forces, "n" if
 * not. Used for naming the outfile.
 * @param count_particles If true, the method also creates a .txt-file which
 * keeps track of the number of particles inside the penningtrap.
 *
*/
void Penningtrap::write_to_file(std::string evolve, std::string dt_str, std::string has_coulomb_str){
    std::string N = std::to_string(this->num_particles_inside);

    std::vector<std::string> names =  {"x", "y", "z", "vx", "vy", "vz"};

    /*
    if (count_particles) {
        std::ofstream count_outfile;
        count_outfile.open("Particle_Counter_" + evolve + N + "_" + has_coulomb_str + "_" + dt_str, std::ios_base::app);
        count_outfile << this->num_particles_inside << std::endl;
    }
    */

    for (int i = 0; i < 3; i++) {

        std::ofstream r_outfile;
        std::ofstream v_outfile;

        r_outfile.open(evolve+N+"_"+has_coulomb_str+"_"+dt_str+"_"+names[i]+".txt", std::ios_base::app); // append instead of overwrite
        v_outfile.open(evolve+N+"_"+has_coulomb_str+"_"+dt_str+"_"+names[i+3]+".txt", std::ios_base::app); // append instead of overwrite
        for (Particle p : this->particles) {  // particle number
            r_outfile << p.r(i) << "   ";
            v_outfile << p.v(i) << "   ";
        }
        r_outfile << "\n";
        v_outfile << "\n";

        r_outfile.close();
        v_outfile.close();
    }

}

/**
 * @brief Evolves the position and velocity of each particle one step in time using Forward Euler.
 *
 * @param dt Time step of the simulation.
 * @param has_coulomb_force True if the simulation should include coulomb
 * forces between the particles, false if these forces should be ignored.
 * @param has_E_field True if the simulation should include forces from the
 * electric field, false if these forces should be ignored.
 * @param has_B_field True if the simulation should include forces from the
 * magnetic field, false if these forces should be ignored.
 *
*/
void Penningtrap::evolve_forwardeuler(double dt, bool has_coulomb_force, bool has_E_field, bool has_B_field, bool func_V, double f, double w, int i) {

  // freezing time explicitly for the calculation of the coulomb force so it's updated homogenously for all particles
  for (Particle& p : this->particles){
    p.r_coulomb = p.r;
  }

  for (Particle& p : this->particles) {  // Iterate through particles

    if (std::sqrt(p.r(0)*p.r(0) + p.r(1)*p.r(1) + p.r(2)*p.r(2)) > this->d && !p.check_outside()) {

      p.is_outside();

      std::cout << "A particle left the trap!" << std::endl;

      this->num_particles_out++;
    }

    // skips the particle if it's outside
    if (p.check_outside()) {
      continue;
    }

    this->find_force(has_coulomb_force, has_E_field, has_B_field, func_V, f, w, i*dt);

    // for updating the coulomb force on all particles before they move.
    // basically freezing time and calculating the coulombforce from all particles in their original position.

    p.v = p.v + dt * p.force / p.m;
    p.r = p.r + p.v * dt;
  }

}


/**
 * @brief Evolves the position and velocity of each particle one step in time using Runge-Kutta 4.
 *
 * @param dt Time step of the simulation.
 * @param has_coulomb_force True if the simulation should include coulomb
 * forces between the particles, false if these forces should be ignored.
 * @param has_E_field True if the simulation should include forces from the
 * electric field, false if these forces should be ignored.
 * @param has_B_field True if the simulation should include forces from the
 * magnetic field, false if these forces should be ignored.
 *
*/
void Penningtrap::evolve_RK4(double dt, bool has_coulomb_force, bool has_E_field, bool has_B_field, bool func_V, double f, double w, int i) {

  // freezing time explicitly for the calculation of the coulomb force so it's updated homogenously for all particles
  for (Particle& p : this->particles){
    p.r_coulomb = p.r;
  }

  for (Particle& p : this->particles) {  // Iterate through particles
    // Empty current particle's array of weighted sum of k-values:

    // skips the particle if it's outside
    if (p.check_outside()) {
      continue;
    }

    p.runge_kutta_k = arma::zeros(3, 2);

    // Store particle's current position and velocity:
    p.r_temp = p.r;
    p.v_temp = p.v;

    this->find_force(has_coulomb_force, has_E_field, has_B_field, func_V, f, w, i*dt);

    arma::vec k1_r = dt * p.v;
    arma::vec k1_v = dt * p.force/p.m;

    p.r = p.r + k1_r/2;
    p.v = p.v + k1_v/2;

    // Add current k-value to the weighted sum of ks:
    p.runge_kutta_k.col(0) = p.runge_kutta_k.col(0) + (1/6.0) * k1_r;
    p.runge_kutta_k.col(1) = p.runge_kutta_k.col(1) + (1/6.0) * k1_v;
  }

  // freezing time explicitly for the calculation of the coulomb force so it's updated homogenously for all particles
  for (Particle& p : this->particles){
    p.r_coulomb = p.r;
  }

  for (Particle& p : this->particles) {  // Iterate through particles

    // skips the particle if it's outside
    if (p.check_outside()) {
      continue;
    }
    this->find_force(has_coulomb_force, has_E_field, has_B_field, func_V, f, w, i*dt);

    arma::vec k2_r = dt * p.v;
    arma::vec k2_v = dt * p.force/p.m;

    p.r = p.r_temp + k2_r/2;
    p.v = p.v_temp + k2_v/2;

    p.runge_kutta_k.col(0) = p.runge_kutta_k.col(0) + (1/6.0) * 2 * k2_r;
    p.runge_kutta_k.col(1) = p.runge_kutta_k.col(1) + (1/6.0) * 2 * k2_v;
  }


  // freezing time explicitly for the calculation of the coulomb force so it's updated homogenously for all particles
  for (Particle& p : this->particles){
    p.r_coulomb = p.r;
  }

  for (Particle& p : this->particles) {  // Iterate through particles

    // skips the particle if it's outside
    if (p.check_outside()) {
      continue;
    }

    this->find_force(has_coulomb_force, has_E_field, has_B_field, func_V, f, w, i*dt);

    arma::vec k3_r = dt * p.v;
    arma::vec k3_v = dt * p.force/p.m;

    p.r = p.r_temp + k3_r;
    p.v = p.v_temp + k3_v;

    p.runge_kutta_k.col(0) = p.runge_kutta_k.col(0) + (1/6.0) * 2 * k3_r;
    p.runge_kutta_k.col(1) = p.runge_kutta_k.col(1) + (1/6.0) * 2 * k3_v;
  }

  // freezing time explicitly for the calculation of the coulomb force so it's updated homogenously for all particles
  for (Particle& p : this->particles){
    p.r_coulomb = p.r;
  }

  for (Particle& p : this->particles) {  // Iterate through particles
    if (std::sqrt(p.r(0)*p.r(0) + p.r(1)*p.r(1) + p.r(2)*p.r(2)) > this->d && !p.check_outside()) {//if (arma::norm(p.r) > this->d && !p.check_outside()) {

      p.is_outside();

      std::cout << "A particle left the trap!" << std::endl;

      this->num_particles_out++;
    }

    // skips the particle if it's outside
    if (p.check_outside()) {
      continue;
    }
    this->find_force(has_coulomb_force, has_E_field, has_B_field, func_V, f, w, i*dt);

    arma::vec k4_r = dt * p.v;
    arma::vec k4_v = dt * p.force/p.m;

    p.runge_kutta_k.col(0) = p.runge_kutta_k.col(0) + (1/6.0) * k4_r;
    p.runge_kutta_k.col(1) = p.runge_kutta_k.col(1) + (1/6.0) * k4_v;

    // Update position and velocity:
    p.r = p.r_temp + p.runge_kutta_k.col(0);
    p.v = p.v_temp + p.runge_kutta_k.col(1);
  }

}

int Penningtrap::particles_inside() {
    return this->num_particles_inside;
}

void Penningtrap::simulate(bool has_coulomb_force,int N, double dt, std::string evolve, bool func_V, double f, double w) {
  this->num_particles_out = 0;

  std::string dt_str = std::to_string(dt);
  std::string has_col;
  if (has_coulomb_force) {
      has_col = "y";
  }
  else {
      has_col = "n";
  }

  // float time = dt*N;

  // std::cout << "total time: " << time << " microseconds" << std::endl;

  if (func_V) {
    for (int i = 0; i < N; i++) {
      if (evolve=="RK4"){
        this->evolve_RK4(dt, has_coulomb_force, true, true, func_V, f, w, i);
      }
      else{
        this->evolve_forwardeuler(dt, has_coulomb_force, true, true, func_V, f, w, i);
      }
    }
  }

  else {
    for (int i = 0; i < N; i++) {
      this->write_to_file(evolve,dt_str, has_col);
      if (evolve=="RK4"){
        this->evolve_RK4(dt, has_coulomb_force, true, true);

      // std::cout << i << std::endl;
      }
      else{
        this->evolve_forwardeuler(dt, has_coulomb_force, true, true);
      }
    }
    this->write_to_file(evolve, dt_str, has_col);
  }
}