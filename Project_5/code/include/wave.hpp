#include <armadillo>
#include <iostream>
#include "matrix.hpp"

// this will be a class!

arma::cx_dvec wave_init(double centerx, double centery, double widthx, double widthy, double px, double py, int M) {
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

            std::complex<double> arg = std::exp(-(temp_x*temp_x/(2*widthx*widthx)) - (temp_y*temp_y/(2*widthy*widthy))
                                                + im*(px*temp_x) + im*(py*temp_y));

            norm += std::conj(arg) * arg;

            u.at(pair_to_single(i, j, len)) = arg;
        }
    }

    u /= norm;

    return u;
}