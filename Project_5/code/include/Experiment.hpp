#ifndef __Experiment__
#define __Experiment__
#include <armadillo>
#include <string>

class Experiment{

private:
    arma::cx_dvec* u_now;
    arma::cx_dvec* u_next;
    arma::mat V;
    std::vector<arma::sp_cx_dmat> AB;
    arma::cx_dvec b;
    arma::cx_dcube storage;

    double T;
    double dt;
    int n_timesteps;
    int len;
    int M;

protected:

    arma::mat potential(double potential, int M, int n_slit, double thickx, double centerx, double slit_sep, double aperture);

    arma::cx_dvec* wave_init(double centerx, double centery, double widthx, double widthy, double px, double py, int M);

public:
    Experiment(double h, double dt, double T, double xc, double yc, double px, double py, double widthx, double widthy, double potential, int n_slit=2, double thickx=0.02, double centerx=0.5, double slit_sep=0.05, double aperture=0.05);
    void run();
    void print(std::string filename);
};

#endif