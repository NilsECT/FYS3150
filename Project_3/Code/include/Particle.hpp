#ifndef __Particle__
#define __Particle__
#include <armadillo>


class Particle{
    private:
    double q;
    double m;
    arma::vec r;
    arma::vec v;
    arma::vec coulomb_force;
    double ke;

    protected:

    public:
    Particle(const double q, const double m, arma::vec r, arma::vec v);
    arma::vec get_r();
    arma::vec get_v();
    void set_v(arma::vec v_new);
    void set_r(arma::vec r_new);
    double get_q();
    double get_m();
    void find_coulomb_force(std::vector<Particle> particles);
    arma::vec get_coulomb_force();
}; 

#endif