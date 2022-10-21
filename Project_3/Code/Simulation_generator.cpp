#include "Particle.hpp"
#include "Penningtrap.hpp"

#include <armadillo>
#include <iostream>
// #include <string>

void simulate(Penningtrap trap, bool has_coulomb_force,int N, double dt, std::string evolve, bool func_V=false, double f=0, double w=0) {
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

    std::vector<std::string> names =  {"x", "y", "z", "vx", "vy", "vz"};

    if (func_V) {
      for (int i = 0; i < N-1; i++) {
        if (evolve=="RK4"){
          trap.evolve_RK4(dt, has_coulomb_force, true, true, func_V, f, w, i);
        }
        else{
          trap.evolve_forwardeuler(dt, has_coulomb_force, true, true, func_V, f, w, i);
        }
      }
    }

    else {
      for (int i = 0; i < N-1; i++) {
        trap.write_to_file(evolve,dt_str, has_col);
        if (evolve=="RK4"){
          trap.evolve_RK4(dt, has_coulomb_force, true, true);

        // std::cout << i << std::endl;
        }
        else{
          trap.evolve_forwardeuler(dt, has_coulomb_force, true, true);
        }
      }
      trap.write_to_file(evolve, dt_str, has_col);
    }
}

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
  double dt = 0.01;
  Penningtrap trap = Penningtrap(B_0, V_0, d);
  arma::vec r_1 = arma::vec{20.,0.,20.}; // micro m
  arma::vec v_1 = arma::vec{0.,25.,0.};
  arma::vec r_2 = arma::vec{25.,25.,0.}; // micro m
  arma::vec v_2 = arma::vec{0.,4.,5.};
  bool has_coulomb_force = true;
  bool no_coulomb_force = false;
  Particle part_1 = Particle(q, m, r_1, v_1);
  Particle part_2 = Particle(q, m, r_2, v_2);

  std::vector<std::string> methods = {"RK4", "Forward_Euler"};
  std::vector<int> n_particles = {1, 2};
  //4 Different time steps:
  std::vector<double> n = {4000.,8000., 10000., 16000.,32000.}; // micro-sec
  for (std::string evolve : methods) {

    trap.add_particle(part_1);

    for (int i = 0; i < 5; i++) {
      dt = 50./n[i];
      if (i == 2) {
        dt = 0.01;
      }
      simulate(trap,no_coulomb_force,int(n[i]),dt,evolve);
      part_1.reset(r_1, v_1);
    }

    dt = 0.01;
    
    trap.add_particle(part_2);
    simulate(trap,no_coulomb_force,int(n[2]),dt,evolve);
    part_1.reset(r_1, v_1);
    part_2.reset(r_2, v_2);
    simulate(trap, has_coulomb_force, int(n[2]), dt, evolve);
    
    part_1.reset(r_1, v_1);
    part_2.reset(r_2, v_2);
    trap.clear_particles();
  }
  return 0;
}
