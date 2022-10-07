#include "Penningtrap.hpp"
#include "Particle.hpp"
#include <armadillo>
#include <iostream>
#include <iomanip>
#include <unistd.h>
#include <assert.h>
#include <string>


 
/**
 * @brief What it does.
 *
 * @param What it needs.
 *
 */
Penningtrap::Penningtrap(double B_0, double V_0, double d){
    this->B_0 = B_0;
    this->V_0 = V_0;
    this->d = d;
    std::vector<Particle> particles;
    this->particles = particles;
}


/**
 * @brief 
 *
 * @param particle_interactions 
 *
 */
void Penningtrap::find_force(bool particle_interactions) {
    for (Particle &particle : this->particles) {
        arma::vec C = arma::vec(3).fill(0);

        if (particle_interactions) {
            C = particle.find_coulomb_force(this->particles);
        }

        arma::vec L = particle.find_Lorentz_force(this->V_0, this->B_0, this->d);

        arma::vec F = C + L;
        std::cout << C << L << F << std::endl;
        particle.force = F;
    }
}

int Penningtrap::num_particles_inside() {
    int counter = 0;
    for (Particle particle : this->particles) {
        arma::vec r = particle.r;
        double r_norm = arma::norm(r);

        if (r_norm < this->d) {
            counter++;
        }
    }
    return counter;
}

void Penningtrap::add_particle(Particle particle) {
    this->particles.push_back(particle);
}

void Penningtrap::generate_particles(int N, double q, double m) {
    this->particles.clear();
    arma::vec r;
    arma::vec v;

    arma_rng::set_seed();

    for (int i = 0; i < N; i++) {
        r = arma::vec(3).randn() * 0.1 * this->d;
        v = arma::vec(3).randn() * 0.1 * this->d;
        Particle particle = Particle(q, m, r, v);
        this->particles.push_back(particle);
    }
}

/**
 * @brief 
 *
 * @param 
 *
*/

void write_to_file(std::string h, std::string inter){

        
    std::string N = std::to_string(this->num_particles_inside);
        
    std::ofstream r_outfile;
    std::ofstream v_outfile;
    std::vector<std::string> names =  {"x", "y", "z", "vx", "vy", "vz"};
    
    for (int i = 0: i<3;i++) {
        r_outfile.open(N+"_"+inter+"_"+h+"_"+name[i]+".txt", std::ios_base::app); // append instead of overwrite
        v_outfile.open(N+"_"+inter+"_"+h+"_"+name[i+3]+".txt", std::ios_base::app); // append instead of overwrite
        for (Particle p : this->particles) {  // particle number
            r_outfile << p.r(i) << "   "; 
            v_outfile << p.v(i) << "   "; 
        }
        r_outfile << "\n";
        v_outfile << "\n";
    }
        

}


