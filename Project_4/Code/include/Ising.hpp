#ifndef __Ising__
#define __Ising__
#include <armadillo>
#include <string>

class Ising{

    public:
    Ising() = default;

    void monte_carlo(int num_MC_cycles, double &avg_eps, double &avg_eps_sq, double &avg_m_abs, double &avg_m_sq, double T, int seed, int L, int burn=0, bool random_config=true);

    std::vector<double> analytical(int L, double T);

    void analytical_comparison(std::vector<double> temperatures, int N_MC_cycles, std::string filename = "analytical_comparison", const int seed = 137);

    void varying_n_mc_cycles(double temperature, arma::vec n_cycles, std::string filename = "varying_cycles", int seed = 137, bool random_config = true, int L = 20, int N_spinflips = 100);

    void phase_transitions(int lattice, arma::vec &temperatures, int N_MC_cycles, int seed, int burn=0, std::string filename = "phase_transition");

    void phase_transitions(std::vector<int> lattice, arma::vec temperatures, int N_MC_cycles, int n_samples, int seed, int burn=0, std::string filename = "phase_transition_varL");

    void epsilon_dist(arma::vec temperature, int L, int N_cycles, int n_samples, int burn=0, std::string filename="epsilon_distribution", int seed = 137);

    void epsilon_dist(arma::vec temperature, std::vector<int> lattice, int N_cycles, int n_samples, int burn=0, std::string filename="epsilon_distribution", int seed = 137);

    void varying_n_mc_cycles(arma::vec temperature, int n_cycles, int n_samples, int lattice = 20, std::string filename = "varying_cycles", int seed = 137);

    void varying_n_mc_cycles(arma::vec temperature, int n_cycles, int n_samples, std::vector<int> lattice, std::string filename = "varying_cycles", int seed = 137);

    void varying_n_walk(arma::vec temperature, std::vector<int> n_walks, std::vector<int> lattice, int num_samples = 100, std::string filename = "varying_walk", int seed = 137);

    void varying_n_walk(arma::vec temperature, int n_walks, std::vector<int> lattice, int num_samples = 100, std::string filename = "varying_walk", int seed = 137);
};

#endif