#include "Particle.hpp"
#include "Penningtrap.hpp"

#include <armadillo>
#include <iostream>
#include <string>

void simulate(Penningtrap trap, bool has_coulomb_force,int N, double dt) {
    std::string dt_str = std::to_string(dt);
    std::string has_col;
    if (has_coulomb_force) {
        has_col = "y";
    }
    else {
        has_col = "n";
    }

    std::string evolve = "RK";

    float time = dt*N;

    std::cout << "total time: " << time << " microseconds" << std::endl;

    std::vector<std::string> names =  {"x", "y", "z", "vx", "vy", "vz"};

    for (int i = 0; i < N; i++) {
      trap.write_to_file(evolve, dt_str, has_col);
      //trap.update_V_0(f,w,t); // Define updating parameters of the strength of the E-field.
      trap.evolve_RK4(dt, has_coulomb_force, true, true);
      //std::cout << i << std::endl;
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

    bool check = q/m > 4.0*V_0/(B_0*B_0*d*d);

    if (!check) {
      std::cout << "Warning: q/m is not greater than that other thing" << std::endl;
    }


    int N = 1000;
    double dt = 0.01;
    bool has_coulomb_force = true;
    /*
    simulate(trap,has_coulomb_force,N,dt);

    has_coulomb_force = false;

    simulate(trap,has_coulomb_force,N,dt);
    */
    Penningtrap trap_100 = Penningtrap(B_0, V_0, d);
    trap_100.generate_particles(100,q,m,seed);
    //has_coulomb_force = false;
    //simulate(trap_100,has_coulomb_force,N,dt);
    has_coulomb_force = true;
    simulate(trap_100,has_coulomb_force,N,dt);

    return 0;
}
