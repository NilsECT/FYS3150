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

    float time = dt*N;

    // std::cout << "total time: " << time << " microseconds" << std::endl;

    std::vector<std::string> names =  {"x", "y", "z", "vx", "vy", "vz"};

    if (func_V && evolve=="RK4"){
      for (int i = 0; i < N; i++) {
        trap.write_to_file(evolve + "f", dt_str, has_col);
        trap.evolve_RK4(dt, has_coulomb_force, true, true, func_V, f, w, i);
      }
    }
    else if (func_V && evolve!="RK4"){
      for (int i = 0; i < N; i++) {
        trap.write_to_file(evolve + "f", dt_str, has_col);
        trap.evolve_forwardeuler(dt, has_coulomb_force, true, true, func_V, f, w, i);
      }
    }
    else if (!func_V && evolve=="RK4"){
      for (int i = 0; i < N; i++) {
        trap.write_to_file(evolve,dt_str, has_col);
        trap.evolve_RK4(dt, has_coulomb_force, true, true);
      }
    }
    else {
      for (int i = 0; i < N; i++) {
        trap.write_to_file(evolve,dt_str, has_col);
        trap.evolve_forwardeuler(dt, has_coulomb_force, true, true);
      }
    }

    trap.write_to_file(evolve,dt_str, has_col);

}

void simulate_perturbation(Penningtrap trap, bool has_coulomb_force, int N, double dt, arma::vec f) {
  std::string dt_str = std::to_string(dt);

  double dw = 0.02;   // [MHz]
  double w_min = 0.2; // [MHz]
  double w_max = 2.5; // [MHz]

  int M = w_max/dw;   // number of values for angular frequency w (ie. omega)


  for (double &amp : f) {
    std::cout << std::endl;
    std::cout << "Simulate for " << N*dt << " microseconds, with amplitude f = " << amp << std::endl;
    std::cout << "         ";

    for (double w = w_min; w <= w_max; w += dw) {
      std::cout << w;
      simulate(trap, has_coulomb_force, N, dt, "RK4", true, amp, w);
      trap.write_to_file_perturbation(amp, w, has_coulomb_force, 100);
    }
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
  int N = 1000;
  double dt = 0.01;

  bool has_coulomb_force = false;

  Penningtrap trap_100 = Penningtrap(B_0, V_0, d);
  trap_100.generate_particles(1, q, m, seed);

  arma::vec f = {0.1, 0.4, 0.7};
  simulate_perturbation(trap_100, has_coulomb_force, N, dt, f);

  return 0;
}

int old_main(){
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
  bool no_coulomb_force = true;
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
      trap.simulate(no_coulomb_force,int(n[i]),dt,evolve);;
    }

    dt = 0.01;

    std::cout << "At two particles without Coulomb" << std::endl;
    trap.add_particle(part_2);
    trap.simulate(no_coulomb_force,int(n[2]),dt,evolve);
    std::cout << "At two particles with Coulomb" << std::endl;
    trap.simulate(has_coulomb_force, int(n[2]), dt, evolve);

    trap.clear_particles();
  }
  return 0;
}
