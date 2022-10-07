#include <armadillo>
#include <iostream>
#include "Penningtrap.hpp"
#include "Particle.hpp"


int main(){
    // Defining core values used for simulation:
    double q = 1;
    double m = 1;
    double B_0 = 1;
    double V_0 = 1;
    double d = 1;
    
    // Initializing A Penningtrap:
    Penningtrap penningtrap = Penningtrap(B_0, V_0,d);



    // Creating random particles in the Penningtrap:
    int N = 1 // Numbers of particles:
    penningtrap.generate_particles(N,q,m);


    // Creating a single particle:
    /*
    arma::vec r; // Position of particle
    arma::vec v; // Velocity of particle
    Particle particle = Particle(q,m,r,v); // Creating a particle object
    penningtrap.add_particle(particle); // Giving Penningtrap the particle
    */


    // Defining simulation variables:
    double t = 1; // Total runtime of simulation in seconds
    double n_step = 10; // Number of time steps
    double dt = t/n_step // time step size


    // Simulation:
    for (int i = 0; i < t*n_step; i++){
        // Choose Forward Euler or Runge-Kutta-4 as simulation method:

        // Forward Euler:
        /*
        penningtrap.Forward_Euler(dt);
        */ 

        // Runge-Kutta
        /*
        penningtrap.Runge_kutta_4(dt);
        */



        // Extracting Output to file: 
        penningtrap.write_to_file();
    }
    return 0;
}

