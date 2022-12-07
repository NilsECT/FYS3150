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
    double potential = 0; // Potential of the wall
    // int n_slit = 0; // Number of slit. 

    //---------------------------
    // Initializing experiment. nr.1 (No potential):
    //---------------------------
    Experiment exp_1 = Experiment(h, dt, T, xc, yc, px, py, width_x, width_y, potential);

    //---------------------------
    // Running experiment:
    //---------------------------
    exp_1.run();

    //---------------------------
    // Saving data as armadillo cube:
    //---------------------------
    exp_1.print("Experiment_1");


    //---------------------------
    // Initializing experiment nr. 2 (Two slits):
    //---------------------------
    potential = 1e10;
    width_y = 0.1; // Deviation in y-direction of the gauss-shaped wave.
    Experiment exp_2 = Experiment(h, dt, T, xc, yc, px, py, width_x, width_y, potential);

    //---------------------------
    // Running experiment:
    //---------------------------
    exp_2.run();

    //---------------------------
    // Saving data as armadillo cube:
    //---------------------------
    exp_2.print("Experiment_2");

    //---------------------------
    // Saving potential as armadillo matrix:
    //---------------------------
    exp_2.print_potential("Experiment_2");

    //---------------------------
    // Initializing experiment nr. 3 (Two slits):
    //---------------------------
    width_y = 0.2; // Deviation in y-direction of the gauss-shaped wave.
    T = 0.002; // Total time for simulation.
    Experiment exp_3 = Experiment(h, dt, T, xc, yc, px, py, width_x, width_y, potential);

    //---------------------------
    // Running experiment:
    //---------------------------
    exp_3.run();

    //---------------------------
    // Saving data as armadillo cube:
    //---------------------------
    exp_3.print("Experiment_3");
    exp_3.save_u("Experiment_3");

    //---------------------------
    // Initializing experiment nr. 4 (single-slit):
    //--------------------------- 
    int n_slit = 1;
    Experiment exp_4 = Experiment(h, dt, T, xc, yc, px, py, width_x, width_y, potential, n_slit);
    //---------------------------
    // Running experiment:
    //---------------------------
    exp_4.run();

    //---------------------------
    // Saving data as armadillo cube:
    //---------------------------
    exp_4.print("Experiment_4");

    //---------------------------
    // Initializing experiment nr. 5 (triple-slit):
    //--------------------------- 
    int n_slit = 3;
    Experiment exp_4 = Experiment(h, dt, T, xc, yc, px, py, width_x, width_y, potential, n_slit);
    //---------------------------
    // Running experiment:
    //---------------------------
    exp_4.run();

    //---------------------------
    // Saving data as armadillo cube:
    //---------------------------
    exp_4.print("Experiment_5");
    


    


    return 0;
}