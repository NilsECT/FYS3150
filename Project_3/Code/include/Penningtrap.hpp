#ifndef __Penningtrap__
#define __Penningtrap__
#include <armadillo>
#include "Particle.hpp"
 
class Penningtrap{
    public:

    double B_0;
    double V_0;
    double d;
    std::vector<Particle> particles;
    int num_particles_inside;
    
    Penningtrap(double B_0, double V_0, double d);
    void find_force(bool has_coloumb_force, bool has_E_field, bool has_B_field);
    void add_particle(Particle particle);
    void generate_particles(int N, double q, double m, int seed k);

    //arma::mat forward_euler(arma::mat r, arma::mat v);

    //arma::mat rk4(arma::mat r,arma::mat v,double m, double q);
    
    void write_to_file();
};

#endif