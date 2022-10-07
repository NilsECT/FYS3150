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

    Penningtrap(double B_0, double V_0, double d);
    
    void find_force(bool particle_interactions);

    int num_particles_inside();
    void add_particle(Particle particle);
    void generate_particles(int N, double q, double m);

    void evolve_forwardeuler(double dt, bool particle_interaction);
    void evolve_rk4(double dt, bool particle_interaction);

    void write_to_file();
};

#endif
