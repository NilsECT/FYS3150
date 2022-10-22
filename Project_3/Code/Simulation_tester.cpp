#include "Particle.hpp"
#include "Penningtrap.hpp"

#include <armadillo>
#include <iostream>
#include <string>

void simulate(Penningtrap trap, bool has_coulomb_force, int N, double dt, std::string evolve, bool func_V=false, double f=0, double w=0) {
    std::string dt_str = std::to_string(dt);
    std::string has_col;
    if (has_coulomb_force) {
        has_col = "y";
    }
    else {
        has_col = "n";
    }

    std::string evolve = "RK4";

    float time = dt*N;

    std::cout << "total time: " << time << " microseconds" << std::endl;

    std::vector<std::string> names =  {"x", "y", "z", "vx", "vy", "vz"};

    for (int i = 0; i < N; i++) {
      trap.write_to_file(evolve, dt_str, has_col);
      trap.update_V(f, w, dt*i); // Define updating parameters of the strength of the E-field.
      evolve_method(dt, has_coulomb_force, true, true, func_V, f, w, i);
      std::cout << i << std::endl;
    }
    trap.write_to_file(evolve, dt_str, has_col);
}

int main(){
    // Defining core values used for simulation:

    double q = 1.0;
    double V_0 = 2.41 * std::pow(10, 6);
    double B_0 = 9.65*10;
    double m = 40.078;
    double d = 500;

    int seed = 137;

    /*
    Penningtrap trap = Penningtrap(B_0, V_0, d);


    std::string num_part = "2";

    //trap.generate_particles(std::stod(num_part), q, m, seed);
    arma::vec r = arma::vec{20,0,20}; // micro m
    arma::vec v = arma::vec{0,25,0};
    Particle particle_1 = Particle(q,m,r,v);
    trap.add_particle(particle_1);
    r = arma::vec{25,25,0}; // micro m
    v = arma::vec{0,40,5};
    Particle particle_2 = Particle(q,m,r,v);
    trap.add_particle(particle_2);
    */
    bool check = q/m > 4.0*V_0/(B_0*B_0*d*d);

    if (!check) {
      std::cout << "Warning: q/m is not greater than that other thing" << std::endl;
    }


    int N = 100;
    double dt = 0.01;
    bool has_coulomb_force = true;
    Penningtrap trap_100 = Penningtrap(B_0, V_0, d);
    trap_100.generate_particles(1,q,m,seed);
    //has_coulomb_force = true;
    //simulate(trap_100,has_coulomb_force,N,dt);

    // Time-varying electric field, with no Coulomb interaction:


    arma::vec f = {0.1, 0.4, 0.7};

    has_coulomb_force = false;
    N = 10000;

    double dw = 0.02; // unit: MHz
    double w_min = 0.2;
    double w_max = 2.5;

    simulate(trap_100, false, N, dt, "RK4", f.at(2), w_max);
    std::cout << trap_100.num_particles_inside << std::endl;

    /*

    for (int i = 0; i < w_max/dw; i ++) {
      simulate(trap_100, false, N, dt, f.at(2), dw*i);
      //std::cout << dw*i << std::endl;
      std::cout << trap_100.num_particles_inside << std::endl;
    }*/


    return 0;
}
