#include "Particle.hpp"
#include "Penningtrap.hpp"

#include <armadillo>
#include <iostream>
//#include <string>

int main(){
    // Defining core values used for simulation:
    double q = 20;
    double m = 1;
    double B_0 = 2;
    double V_0 = 1;
    double d = 1;

    q = 1.0;
    V_0 = 9.65 * std::pow(10, 8);
    B_0 = 9.65*10;
    m = 20.0;
    d = 1 * std::pow(10, 4);

    int seed = 137;

    Penningtrap trap = Penningtrap(B_0, V_0, d);

    std::string num_part = "1";

    trap.generate_particles(std::stod(num_part), q, m, seed);

    bool check = q/m > 4.0*V_0/(B_0*B_0*d*d);

    if (!check) {
      std::cout << "Warning: q/m is not greater than that other thing" << std::endl;
    }

    int N = 10000;
    std::string inter = "y";
    std::string dt = "0.01";

    float time = std::stod(dt)*N;

    std::cout << "total time: " << time << " microseconds" << std::endl;

    std::vector<std::string> names =  {"x", "y", "z", "vx", "vy", "vz"};

    for (int i = 0; i < 6; i++) {
      std::ofstream ofs;
      ofs.open(num_part + "_y_" + dt + "_"+ names[i] + ".txt", std::ofstream::out | std::ofstream::trunc);
      ofs.close();
    }


    for (int i = 0; i < N; i++) {
      trap.write_to_file(dt, inter);
      trap.evolve_RK4(std::stod(dt), false, true, true);

    }

    return 0;
}
