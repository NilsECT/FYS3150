#ifndef __Ising__
#define __Ising__
#include <armadillo>
#include <string>
#include "Grid.hpp"

class Ising{

    public:
    Ising() = default;

    void monte_carlo(int num_MC_cycles, double &avg_eps, double &avg_eps_sq, double &avg_m_abs, double &avg_m_sq, Grid &model, int seed, int burn=0);

    std::vector<double> analytical(int L, double T);

    void analytical_comparison(std::vector<double> temperatures, arma::vec cycles, std::string filename = "analytical_comparison", const int seed = 137);

    void phase_transitions(std::vector<int> lattice, arma::vec temperatures, int N_MC_cycles, int n_samples, int seed, int burn=0, std::string filename = "phase_transition_varL");

    void epsilon_dist(arma::vec temperature, std::vector<int> lattice, int N_cycles, int n_samples, int burn=0, std::string filename="epsilon_distribution", int seed = 137);

    void varying_n_mc_cycles(arma::vec temperature, int n_cycles, int n_samples, std::vector<int> lattice, int burn=0, std::string filename = "varying_cycles", int step=250, int seed = 137);

    void burn_in(Grid &model, int burn, int seed=137);
};

#endif
