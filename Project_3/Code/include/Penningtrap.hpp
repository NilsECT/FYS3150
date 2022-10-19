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
    void generate_particles(int N, double q, double m, int seed);

    void evolve_forwardeuler(double dt, bool particle_interaction, bool has_E_field, bool has_B_field);
    void evolve_RK4(double dt, bool particle_interaction, bool has_E_field, bool has_B_field);
    void update_V_0(double f, double w, double t);
    void write_to_file(std::string evolve,std::string h, std::string inter, bool count_particles);
};

#endif
