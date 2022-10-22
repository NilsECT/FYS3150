#include "Particle.hpp"
#include "Penningtrap.hpp"

#include <armadillo>
#include <iostream>

int main(){
  // Defining core values used for simulation:
  double q = 1.0;
  double V_0 = 2.41 * std::pow(10, 6);
  double B_0 = 9.65*10;
  double m = 40.078;
  double d = 500;

  int n_part = 25;

  int seed = 5764;

  int N = 5000;
  double dt = 0.02;

  bool has_coulomb_force = false;

  Penningtrap trap_100 = Penningtrap(B_0, V_0, d);
  trap_100.generate_particles(n_part, q, m, seed);

  std::string dt_str = std::to_string(dt);

  double dw = 0.02;   // [MHz]
  double w_min = 0.2; // [MHz]
  double w_max = 2.5; // [MHz]
  int M = w_max/dw;   // number of values for angular frequency w (ie. omega)

  arma::vec f = {0.1, 0.4, 0.7};

  for (double &amp : f) {
    std::cout << std::endl;
    std::cout << "Simulate for " << N*dt << " microseconds, with amplitude f = " << amp << std::endl;
    std::cout << "         ";

    for (double w = w_min; w <= w_max; w += dw) {
      std::cout << w << "    ";
      trap_100.reset_particles();
      trap_100.simulate(has_coulomb_force, N, dt, "RK4", true, amp, w);
      trap_100.write_to_file_perturbation(amp, w, has_coulomb_force, n_part);
      std::cout << "Number of particles that left the trap: " << trap_100.num_particles_out << std::endl;
    }
  }

  return 0;
}
