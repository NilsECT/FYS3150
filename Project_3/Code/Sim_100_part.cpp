#include "Particle.hpp"
#include "Penningtrap.hpp"

#include <armadillo>
#include <iostream>
#include <chrono>
//#include <string>

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

    if (func_V) {
      for (int i = 0; i < N; i++) {

        // std::cout << "Evolving RK4. Now at timestep: " << i << std::endl;

        if (evolve=="RK4"){
          trap.evolve_RK4(dt, has_coulomb_force, true, true, func_V, f, w, i);
        }
        else{
          trap.evolve_forwardeuler(dt, has_coulomb_force, true, true, func_V, f, w, i);
        }
      }
    }

    else {
      std::cout << "Not updating V" << std::endl;
    }
}

int main(){
    // Defining core values used for simulation:

    double q = 1.0;
    double V_0 = 2.41 * std::pow(10, 6);
    double B_0 = 9.65*10;
    double m = 40.078;
    double d = 500 * std::pow(10, -6);

    int seed = 137;
    int T = 500;
    double dt = 0.01;
    int N = T/dt;

    // 100 particles with interaction:
    arma::vec omega_v = arma::linspace(0.2, 2.5, 15);
    std::vector<double> f_list = {0.1, 0.4, 0.7};
    bool var_V = true;

    std::cout << "At trap with 100 particles." << std::endl;

    Penningtrap trap_100 = Penningtrap(B_0, V_0, d);
    for (double f : f_list) {

        std::string f_name = std::to_string(f);
        std::ofstream counter_file;

        std::cout << "At f = " << f_name << std::endl;

        counter_file.open("counter_" + f_name + ".txt", std::ios_base::app);

      for (double omega : omega_v) {

        std::cout << "At omega = " << omega << std::endl;

        trap_100.generate_particles(100,q,m,seed);
        bool has_coulomb_force = false;
        std::string evolve = "RK4";
        simulate(trap_100,has_coulomb_force,N,dt,evolve, var_V, f, omega);

        // write number of contained particles to a file
        // omega  conatined_particles_after 500 micro sec
        counter_file << omega << "   " << trap_100.particles_inside();
        counter_file << "\n";
      }
      counter_file.close();
    }

    return 0;
}
