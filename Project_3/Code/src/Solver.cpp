
#include "Solver.hpp"
#include "Penningtrap.hpp"
#include "Particle.hpp"
#include <armadillo>
#include <iostream>
#include <iomanip>
#include <unistd.h>
#include <assert.h>

Solver::Solver(Penningtrap penningtrap) {
    this->penningtrap = penningtrap;
}

void Solver::forward_euler() {
    std::cout << "Do forward euler!" << std::endl;
}

void Solver::runge_kutta() {
    std::cout << "Do Runge-Kutta!" << std::endl;
}