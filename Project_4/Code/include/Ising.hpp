#ifndef __Ising__
#define __Ising__
#include <armadillo>
#include <string>

class Ising{

    public:
    Ising() = default;

    void monte_carlo(int N_spinflips, int num_MC_cycles, double &avg_eps, double &avg_eps_sq, double &avg_m_abs, double &avg_m_sq, double T, int seed, int L, bool random_config = true);

    std::vector<double> analytical(int L, double T);

    void analytical_comparison(std::vector<double> temperatures, int N_spinflips, int N_MC_cycles, std::string filename = "analytical_comparison", const int seed = 137);

    void varying_n_mc_cycles(double temperature, arma::vec n_cycles, std::string filename = "varying_cycles", int seed = 137, bool random_config = true, int L = 20, int N_spinflips = 100);

    void phase_transitions(int L, arma::vec &temperatures, int N_spinflips, int N_MC_cycles, int seed, std::string filename = "phase_transition");
};

#endif