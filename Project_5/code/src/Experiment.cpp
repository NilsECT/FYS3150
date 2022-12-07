#include "Experiment.hpp"
#include <armadillo>
#include <iostream>
#include "Matrix.hpp"

/**
 * @brief Construct a new Experiment:: Experiment object Sets up the values for the experiments. Create the object, run the experiment then print to file.
 * 
 * @param h Stepsize in the spatial dimensions (x and y).
 * @param dt Stepsize in the time dimension.
 * @param T Total time over which to run the simulation.
 * @param xc The center of the wave packet in the x direction.
 * @param yc The center of the wave packet in the y direction.
 * @param px The momentum of the wave packet in the x direction.
 * @param py The momentum of the wave packet in the y direction.
 * @param widthx The standard deviation squared (sigma^2) of the wave packet in the x direction.
 * @param widthy The standard deviation squared (sigma^2) of the wave packet in the y direction.
 * @param potential Numerical value of the desired potential.
 * @param n_slit Number of slits.
 * @param thickx Thickness of the walls in the x direction.
 * @param centerx The center of the wall placement in the x direction.
 * @param slit_sep The separation distance between the slits (the small walls in between).
 * @param aperture The opening of the slits.
 */
Experiment::Experiment(double h, double dt, double T, double xc, double yc, double px, double py, double widthx, double widthy, double potential, int n_slit, double thickx, double centerx, double slit_sep, double aperture) {
    // calculate M
    M = std::round(1./h);
    len = (M-2);

    // Generate potential
    V = this->potential(potential, M, n_slit, thickx, centerx, slit_sep, aperture);

    // Initialise wave
    this->u = wave_init(xc, yc, widthx, widthy, px, py, M);

    // checking start probability
    // arma::mat temp = this->probability(u);
    // std::cout << arma::accu(temp) << std::endl;

    // Generate A and B
    Matrix matrix = Matrix();
    AB = matrix.create_AB(V, h, dt, M);

    // generate storage
    n_timesteps = std::round(T/dt);
    storage = arma::cube(len, len, n_timesteps);
    Re = arma::cube(len, len, n_timesteps);
    Im = arma::cube(len, len, n_timesteps);
    // store remaining parameters
    this->T = T;
    this->dt = dt;
}

/**
 * @brief Runs the experiment with the values and places the values into a cube. This cube can be printed once the simulation is complete.
 * 
 */
void Experiment::run() {
    Matrix matrix = Matrix();
    
    // store initial state
    storage.slice(0) = this->probability(this->u);
    
    // simulate
    for (int i=1; i < n_timesteps ; i++) {
        // calc b
        // b = matrix.mult_Bu(this->u, AB.at(1));
        b = AB.at(1) * this->u;

        // find next timestep
        this->u = arma::spsolve(AB.at(0), b);
        storage.slice(i) = this->probability(this->u);
        print_u(i);
    }
}

/**
 * @brief Prints the cube containing the simulation to a file.
 * 
 * @param filename 
 */
void Experiment::print(std::string filename) {
    // storage.save(filename + ".csv", arma::file_type::arma_ascii);
    storage.save(filename + ".bin");
}

/**
 * @brief Prints the matrix containing the potential to a file.
 * 
 * @param filename 
 */
void Experiment::print_potential(std::string filename) {
    V.save(filename + ".bin");
}

/**
 * @brief Prints the matrix containing the potential to a file.
 * 
 * @param filename 
 */
void Experiment::save_u(std::string filename) {
    Re.save(filename + "_Re.bin");
    Im.save(filename + "_Im.bin");
}

/**
 * @brief Takes the computed simulation and produces a matrix containing the probabilities of where to find the particle at each grid coordinates.
 * 
 * @param u Vector containing the next timestep data.
 * @param len Number of rows/columns in the wave packet matrix (square).
 * @return arma::mat MAtrix containing the probabilities.
 */
arma::mat Experiment::probability(arma::cx_dvec &u) {
    //
    // Why do we send in u, when we just use this->u?
    //
    arma::mat prob = arma::mat(len, len);
    for (int j=0; j < len; j++) {
        for (int i=0; i < len; i++) {
            
            prob(i, j) = std::real(std::conj(this->u(Matrix::pair_to_single(i, j, len))) * this->u(Matrix::pair_to_single(i, j, len)));
        }
    }

    return prob;
}

/**
 * @brief Takes inn the current time step and stores the real an imaginary values of i in a cube.
 * 
 * @param u Vector containing the next timestep data.
 * @param len Number of rows/columns in the wave packet matrix (square).
 */
void Experiment::print_u(int t) {
    Matrix matrix_2 = Matrix();
    arma::cx_mat temp = arma::cx_mat(len,len);
    int temp_len = len*len;
    for (int i = 0; i < temp_len; i++) {
        std::tuple<int,int> a = matrix_2.single_to_pair(i,len);
        temp(std::get<0>(a),std::get<1>(a)) = this->u(i);
        }
    Re.slice(t) = arma::real(temp);
    Im.slice(t) = arma::imag(temp);

}

/**
 * @brief Initialiszes the wave packet and normalises it, returns the pointer to the initialised vector. All values are normalised to be in [0, 1]
 * 
 * @param centerx The center of the wave packet in the x direction.
 * @param centery The center of the wave packet in the y direction.
 * @param widthx The standard deviation squared (sigma^2) of the wave packet in the x direction.
 * @param widthy The standard deviation squared (sigma^2) of the wave packet in the y direction.
 * @param px The momentum of the wave packet in the x direction.
 * @param py The momentum of the wave packet in the y direction.
 * @param M The size of the full grid, including the boundaries.
 * @return arma::cx_dvec*  Pointer to the initialised wave packet.
 */
arma::cx_dvec Experiment::wave_init(double centerx, double centery, double widthx, double widthy, double px, double py, int M) {
    int len = M-2;
    arma::cx_dvec u = arma::cx_dvec(len*len);
    double h = 1./len;
    std::complex<double> im(0., 1.);
    //std::complex<double> norm = 0.;
    
    // i is y and j is x
    for (int i=0; i<(len); i++) {
        for (int j=0; j<(len); j++) {

            double temp_x = ((h*(j+1)) - centerx);
            double temp_y = ((h*(i+1)) - centery);
            
            std::complex<double> arg = std::exp(-(temp_x*temp_x/(2*widthx*widthx)) - (temp_y*temp_y/(2*widthy*widthy)) + im*(px*temp_x) + im*(py*temp_y));
            
            u(Matrix::pair_to_single(i, j, len)) = arg;
        }
    }
    
    u = u / std::sqrt(arma::accu(arma::conj(u) % u));

    return u;
}

/**
 * @brief Creates the potential with a chosen number of slits.
 * 
 * @param potential The numerical value of the potential.
 * @param M The size of the full grid, including the boundaries.
 * @param n_slit Number of slits.
 * @param thickx Thickness of the walls in the x direction.
 * @param centerx The center of the wall placement in the x direction.
 * @param slit_sep The separation distance between the slits (the small walls in between).
 * @param aperture The opening of the slits.
 * @return arma::mat Matrix containing the potential.
 */
arma::mat Experiment::potential(double potential, int M, int n_slit, double thickx, double centerx, double slit_sep, double aperture) {
    int len = M-2;
    double h = 1./len;
    arma::mat V = arma::mat(len, len).zeros();

    int j_start = std::round((centerx - (thickx/2))/h);
    int j_end = std::round((centerx + (thickx/2))/h) - 1;
    bool center_open = n_slit%2 != 0;   // pair number of slits means the center is open

    // attempt at making a dynamic slit builder
    int counter = n_slit;
    int i_center = std::round(0.5/h) -1;
    int i_loc_up = i_center;
    int i_loc_down = i_center;
    int len_aperture = std::round(aperture/h);
    int len_slit_sep = std::round(slit_sep/h);

    if (center_open) {
        // do nothing with the center
        int i_loc_up = i_center - std::round(aperture/(2*h));
        int i_loc_down = i_center + std::round(aperture/(2*h));
        i_loc_up--;
        i_loc_down++;
        counter--;  // placed the opening

        // if we have more than one slit then we enter this while loop and place the remaining slits
        while (counter > 0) {
            
            // first we place the wall
            for (int ii=0; ii<len_slit_sep; ii++) {
                for (int j=j_start; j<j_end; j++) {
                    // upper wall
                    V(i_loc_up, j) = potential;

                    // lower wall
                    V(i_loc_down, j) = potential;
                    
                }
                i_loc_up--;
                i_loc_down++;
            }

            // then we move over the opening
            i_loc_up -= len_aperture;
            i_loc_down += len_aperture;

            // we have now made two new slits
            counter -= 2;
        }

        // place the remaining walls
        for (int i=i_loc_up; i>=0; i--) {
            for (int j=j_start; j<j_end; j++) {
                // upper wall
                V(i, j) = potential;
            }
        }
        for (int i=i_loc_down; i<len; i++) {
            for (int j=j_start; j<j_end; j++) {
                // lower wall
                V(i, j) = potential;
            }
        }

        return V;
    }
    else if (n_slit == 0) {
        // all wall
        for (int j=j_start; j<j_end; j++) {
            V.col(j).fill(potential);
        }

        return V;
    }
    else {
        // we must place a wall on the center
        int temp = std::round(len_slit_sep/2.);

        for (int ii=0; ii<temp; ii++) {
            for (int j=j_start; j<j_end; j++) {
                // upper wall
                V(i_loc_up, j) = potential;

                // lower wall
                V(i_loc_down, j) = potential;
            }
            i_loc_up--;
            i_loc_down++;
        }

        // we create the two slits
        i_loc_up -= len_aperture;
        i_loc_down += len_aperture;
        counter -= 2;

        // then we make the remaining slits if there are some
        // this is identical to the while loop above
        while (counter > 0) {
            
            // first we place the wall
            for (int ii=0; ii<len_slit_sep; ii++) {
                for (int j=j_start; j<j_end; j++) {
                    // upper wall
                    V(i_loc_up, j) = potential;

                    // lower wall
                    V(i_loc_down, j) = potential;
                    
                }
                i_loc_up--;
                i_loc_down++;
            }

            // then we move over the opening
            i_loc_up -= len_aperture;
            i_loc_down += len_aperture;

            // we have now made two new slits
            counter -= 2;
        }

        // place the remaining walls
        // place the remaining walls
        for (int i=i_loc_up; i>=0; i--) {
            for (int j=j_start; j<j_end; j++) {
                // upper wall
                V(i, j) = potential;
            }
        }
        for (int i=i_loc_down; i<len; i++) {
            for (int j=j_start; j<j_end; j++) {
                // lower wall
                V(i, j) = potential;
            }
        }

        return V;
    }
}