#include <iomanip>
#include "Experiment.hpp"

int main() {
    Experiment exp = Experiment(0.005, 2.5e-5, 0.008, 0.25, 0.5, 200, 0, 0.05, 0.05, 0);
    exp.run();
    exp.print("test");

    return 0;
}