#include <armadillo>
#include <iostream>

arma::mat potential(double potential, int M, int n_slit=2, double thickx=0.02, double centerx=0.5, double slit_sep=0.05, double aperture=0.05) {
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