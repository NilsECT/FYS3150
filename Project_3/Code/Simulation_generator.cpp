#include "Particle.hpp"
#include "Penningtrap.hpp"

#include <armadillo>
#include <iostream>
// #include <string>

void simulate(Penningtrap trap, bool has_coulomb_force,int N, double dt, std::string evolve, bool update_V, bool count_particles = false) {
    std::string dt_str = std::to_string(dt);
    std::string has_col;
    if (has_coulomb_force) {
        has_col = "y";
    }
    else {
        has_col = "n";
    }
    
    float time = dt*N;

    std::cout << "total time: " << time << " microseconds" << std::endl;

    std::vector<std::string> names =  {"x", "y", "z", "vx", "vy", "vz"};

    for (int i = 0; i < N; i++) {
      trap.write_to_file(evolve,dt_str, has_col);
      /*if (update_V){
        trap.update_V_0(f,w,t); // Define updating parameters of the strength of the E-field.
      }*/
      if (evolve=="RK4"){
        trap.evolve_RK4(dt, has_coulomb_force, true, true);

      std::cout << i << std::endl;
      }
      else{
        trap.evolve_forwardeuler(dt, has_coulomb_force, true, true);
      }
    }
    trap.write_to_file(evolve, dt_str, has_col);
}

int main(){
    // Defining core values used for simulation:

    double q = 1.0;
    double V_0 = 9.65 * std::pow(10, 8);
    double B_0 = 9.65*10;
    double m = 40.078;
    double d = 1 * std::pow(10, 4);

    int seed = 137;
    // Single particle
    // With RK4:
    int N = 10000;
    double dt = 0.01;
    Penningtrap trap_1 = Penningtrap(B_0, V_0, d);
    arma::vec r = arma::vec{20,0,20}; // micro m
    arma::vec v = arma::vec{0,25,0};
    Particle particle_1 = Particle(q,m,r,v);
    trap_1.add_particle(particle_1);
    bool has_coulomb_force = false;
    bool update_V = false;
    std::string evolve = "RK4";
    simulate(trap_1,has_coulomb_force,N,dt,evolve,update_V);

    // With Forward Euler:
    Penningtrap trap_2 = Penningtrap(B_0, V_0, d);
    r = arma::vec{20,0,20}; // micro m
    v = arma::vec{0,25,0};
    Particle particle_2 = Particle(q,m,r,v);
    trap_2.add_particle(particle_2);
    has_coulomb_force = false;
    evolve = "Forward_Euler";
    simulate(trap_2,has_coulomb_force,N,dt,evolve,update_V);

    // Two particles with interaction:
    // With RK4:
    Penningtrap trap_3 = Penningtrap(B_0, V_0, d);
    r = arma::vec{20,0,20}; // micro m
    v = arma::vec{0,25,0};
    Particle particle_3 = Particle(q,m,r,v);
    trap_3.add_particle(particle_3);
    r = arma::vec{25,25,0}; // micro m
    v = arma::vec{0,40,5};
    Particle particle_4 = Particle(q,m,r,v);
    trap_3.add_particle(particle_4);
    has_coulomb_force = true;
    evolve = "RK4";
    simulate(trap_3,has_coulomb_force,N,dt,evolve,update_V);

    // With Forward Euler:
    Penningtrap trap_4 = Penningtrap(B_0, V_0, d);
    r = arma::vec{20,0,20}; // micro m
    v = arma::vec{0,25,0};
    Particle particle_5 = Particle(q,m,r,v);
    trap_4.add_particle(particle_5);
    r = arma::vec{25,25,0}; // micro m
    v = arma::vec{0,40,5};
    Particle particle_6 = Particle(q,m,r,v);
    trap_4.add_particle(particle_6);
    has_coulomb_force = true;
    evolve = "Forward_Euler";
    simulate(trap_4,has_coulomb_force,N,dt,evolve,update_V);

    // Two particles without interaction:
    // With RK4:
    Penningtrap trap_5 = Penningtrap(B_0, V_0, d);
    r = arma::vec{20,0,20}; // micro m
    v = arma::vec{0,25,0};
    Particle particle_7 = Particle(q,m,r,v);
    trap_5.add_particle(particle_7);
    r = arma::vec{25,25,0}; // micro m
    v = arma::vec{0,40,5};
    Particle particle_8 = Particle(q,m,r,v);
    trap_5.add_particle(particle_8);
    has_coulomb_force = false;
    evolve = "RK4";
    simulate(trap_5,has_coulomb_force,N,dt,evolve,update_V);

    // With Forward Euler:
    Penningtrap trap_6 = Penningtrap(B_0, V_0, d);
    r = arma::vec{20,0,20}; // micro m
    v = arma::vec{0,25,0};
    Particle particle_9 = Particle(q,m,r,v);
    trap_6.add_particle(particle_9);
    r = arma::vec{25,25,0}; // micro m
    v = arma::vec{0,40,5};
    Particle particle_10 = Particle(q,m,r,v);
    trap_6.add_particle(particle_10);
    has_coulomb_force = false;
    evolve = "Forward_Euler";
    simulate(trap_6,has_coulomb_force,N,dt,evolve,update_V);


    //4 Different time steps:
    std::vector<int> n = {4000,8000,16000,32000};
    for (int k : n){
      int N = k;
      dt = 50./k;
      trap_1 = Penningtrap(B_0, V_0, d);
      r = arma::vec{20,0,20}; // micro m
      v = arma::vec{0,25,0};
      particle_1 = Particle(q,m,r,v);
      trap_1.add_particle(particle_1);
      bool has_coulomb_force = false;
      bool update_V = false;
      evolve = "RK";
      simulate(trap_1,has_coulomb_force,k,dt,evolve,update_V);
    }
    for (int k : n){
      int N = k;
      dt = 50./k;
      trap_1 = Penningtrap(B_0, V_0, d);
      r = arma::vec{20,0,20}; // micro m
      v = arma::vec{0,25,0};
      particle_1 = Particle(q,m,r,v);
      trap_1.add_particle(particle_1);
      bool has_coulomb_force = false;
      bool update_V = false;
      evolve = "FE";
      simulate(trap_1,has_coulomb_force,k,dt,evolve,update_V);
    }


    /*
    // 100 particles with interaction:
    Penningtrap trap_100 = Penningtrap(B_0, V_0, d);
    trap_100.generate_particles(100,q,m,seed);
    update_V = true;
    has_coulomb_force = false;
    evolve = ?;
    simulate(trap_100,has_coulomb_force,N,dt,evolve,update_V);
    */
    return 0;
}
