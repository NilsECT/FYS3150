#include <iomanip>
#include "Experiment.hpp"

int main() {
    //---------------------------
    // Defining core values:
    //---------------------------
    double h = 0.005; // Step size in the matrix where x,y \in [0,1]. Number of grid points is then 1/h^2.
    double dt = 2.5e-5; // Time step in each time evolution.
    double T = 0.008; // Total time for simulation. 
    double xc = 0.25; // x-coordinate for the center of where the wave starts.
    double yc = 0.5; // y-coordinate for the center of where the wave starts.
    double px = 200; // Momentum of the wave in x-direction.
    double py = 0; // Momentum of the wave in y-direction.
    double width_x = 0.05; // Deviation in x-direction of the gauss-shaped wave.
    double width_y = 0.05; // Deviation in y-direction of the gauss-shaped wave.
    int n_slit = 0; // Number of slit. 

    //---------------------------
    // Initializing experiment:
    //---------------------------
    Experiment exp = Experiment(h, dt, T, xc, yc, px, py, width_x, width_y, n_slit);

    //---------------------------
    // Running experiment:
    //---------------------------
    exp.run();

    //---------------------------
    // Saving data as armadillo cube:
    //---------------------------
    exp.print("test");

    return 0;
}