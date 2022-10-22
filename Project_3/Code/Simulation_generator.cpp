#include "Particle.hpp"
#include "Penningtrap.hpp"

#include <armadillo>
#include <iostream>
// #include <string>

int main(){
  // Defining core values used for simulation:

  double q = 1.0;
  double V_0 = 2.41 * std::pow(10, 6);
  double B_0 = 9.65*10;
  double m = 40.078;
  double d = 500;

  int seed = 137;

  // Single particle
  // With RK4:
  double dt;
  Penningtrap trap = Penningtrap(B_0, V_0, d);
  arma::vec r_1 = arma::vec{20.,0.,20.}; // micro m
  arma::vec v_1 = arma::vec{0.,25.,0.}; // micro m per micro sec
  arma::vec r_2 = arma::vec{25.,25.,0.}; // micro m
  arma::vec v_2 = arma::vec{0.,4.,5.};  // micro m per micro sec
  bool has_coulomb_force = true;
  bool no_coulomb_force = false;
  Particle part_1 = Particle(q, m, r_1, v_1);
  Particle part_2 = Particle(q, m, r_2, v_2);

  std::vector<std::string> methods = {"RK4", "Forward_Euler"};
  std::vector<int> n_particles = {1, 2};
  //4 Different time steps:
  std::vector<double> n = {4000.,8000., 10000., 16000.,32000.}; // micro-sec
  for (std::string evolve : methods) {

    std::cout << "At " << evolve << std::endl;

    trap.add_particle(part_1);

    for (int i = 0; i < 5; i++) {
      std::cout << "At simulation with 1 particle, using N = " << n[i] << std::endl;
      dt = 50./n[i];
      if (i == 2) {
        dt = 0.01;
      }
      trap.simulate(no_coulomb_force,int(n[i]),dt,evolve);
      trap.reset_particles();
    }

    dt = 0.01;

    std::cout << "At two particles without Coulomb" << std::endl;
    trap.add_particle(part_2);
    trap.simulate(no_coulomb_force,int(n[2]),dt,evolve);
    trap.reset_particles();
    std::cout << "At two particles with Coulomb" << std::endl;
    trap.simulate(has_coulomb_force, int(n[2]), dt, evolve);
    trap.reset_particles();

    trap.clear_particles();
  }

  std::cout << "Calculating the analytical solutions." << std::endl;
  trap.clear_particles();
  trap.add_particle(part_1);
  for (int i = 0; i < 5; i++) {
    dt = 50./n[i];
    trap.reset_particles();
    trap.analytical(dt, int(n[i]));
  }

  return 0;
}
