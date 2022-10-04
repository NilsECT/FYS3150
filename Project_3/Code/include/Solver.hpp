#ifndef __Solver__
#define __Solver__
#include <armadillo>
#include "Particle.hpp"
#include "Penningtrap.hpp"

class Solver {
    private:
    Penningtrap penningtrap;

    public:
    Solver(Penningtrap penningtrap);
    void forward_euler();
    void runge_kutta();
};

#endif