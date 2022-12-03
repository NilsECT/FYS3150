#include "Experiment.hpp"
#include <armadillo>
#include <iostream>
#include "Matrix.hpp"

Experiment::Experiment(double h, double dt, double T, double xc, double yc, double px, double py, double widthx, double widthy, double potential, int n_slit, double thickx, double centerx, double slit_sep, double aperture) {
    // calculate M
    M = std::round(1./h);
    len = (M-2);

    // Generate potential
    V = this->potential(potential, M, n_slit, thickx, centerx, slit_sep, aperture);

    // Initialise wave and give it to a pointer
    *u_now = this->wave_init(xc, yc, widthx, widthy, px, py, M);

    // Generate A and B
    Matrix matrix = Matrix();
    AB = matrix.create_AB(V, h, dt, M);

    // generate storage
    n_timesteps = std::round(T/dt);
    storage = arma::cx_dcube(len, len, n_timesteps);

    // store remaining parameters
    this->T = T;
    this->dt = dt;
}

void Experiment::run() {
    Matrix matrix = Matrix();

    // store initial state
    storage.slice(0) = matrix.reshape(*u_now, M);

    // simulate
    for (int i=1; i < n_timesteps ; i++) {

        // calc b
        *b = matrix.mult_Bu(*u_now, *AB.at(1));

        // cross fingers for pointer galore
        matrix.gauss_seidel(AB.at(0), b, u_now, u_next, 1e-5);
        storage.slice(i) = matrix.reshape(*u_now, M);
    }
}

void Experiment::print(std::string filename) {
    storage.save(filename + ".bin");
}


arma::cx_dvec Experiment::wave_init(double centerx, double centery, double widthx, double widthy, double px, double py, int M) {
    int len = M-2;
    arma::cx_dvec u =  arma::cx_dvec(len*len);
    double h = 1./len;
    std::complex<double> im(0., 1.);
    std::complex<double> norm = 0.;
    
    // i is y and j is x

    for (int i=0; i<(len); i++) {
        for (int j=0; j<(len); j++) {
            double temp_x = ((h*(j+1)) - centerx);
            double temp_y = ((h*(i+1)) - centery);

            std::complex<double> arg = std::exp(-(temp_x*temp_x/(2*widthx*widthx)) - (temp_y*temp_y/(2*widthy*widthy)) + im*(px*temp_x) + im*(py*temp_y));

            norm += std::conj(arg) * arg;

            u.at(Matrix::pair_to_single(i, j, len)) = arg;
        }
    }

    u /= norm;

    return u;
}

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
                    V.at(i_loc_up, j) = potential;

                    // lower wall
                    V.at(i_loc_down, j) = potential;
                    
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
                V.at(i, j) = potential;
            }
        }
        for (int i=i_loc_down; i<len; i++) {
            for (int j=j_start; j<j_end; j++) {
                // lower wall
                V.at(i, j) = potential;
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
                V.at(i_loc_up, j) = potential;

                // lower wall
                V.at(i_loc_down, j) = potential;
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
                    V.at(i_loc_up, j) = potential;

                    // lower wall
                    V.at(i_loc_down, j) = potential;
                    
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
                V.at(i, j) = potential;
            }
        }
        for (int i=i_loc_down; i<len; i++) {
            for (int j=j_start; j<j_end; j++) {
                // lower wall
                V.at(i, j) = potential;
            }
        }

        return V;
    }
}