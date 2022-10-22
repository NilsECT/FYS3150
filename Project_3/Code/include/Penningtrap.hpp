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
    int num_particles_out;

    Penningtrap(double B_0, double V_0, double d);
    void find_force(bool has_coloumb_force, bool has_E_field, bool has_B_field, bool func_V=false, double f=0, double w=0, double ti=0);
    void add_particle(Particle &particle);
    void clear_particles();
    void generate_particles(int N, double q, double m, int seed);

    void evolve_forwardeuler(double dt, bool particle_interaction, bool has_E_field, bool has_B_field, bool func_V=false, double f=0, double w=0, int i=0);
    void evolve_RK4(double dt, bool particle_interaction, bool has_E_field, bool has_B_field, bool func_V=false, double f=0, double w=0, int i=0);
    
    void write_to_file(std::string evolve,std::string h, std::string inter);
    void write_to_file_perturbation(double f, double w, bool has_coulomb_force, int N_particles=100);

    void simulate(bool has_coulomb_force,int N, double dt, std::string evolve, bool func_V=false, double f=0, double w=0);
    int particles_inside();
};

#endif
