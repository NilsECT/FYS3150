#include <iostream>
#include <string>
#include <fstream>
#include <iomanip>

int main(){
    // setting a filename
    std::string filename = "numbers.txt";

    std::ofstream ofile;
    ofile.open(filename);

    // setting parameters for computation
    double x_min = 0.0;
    double x_max = 1.0;
    int n_steps = 100;  // 99 points
    double h = (x_max - x_min) / n_steps;

    // parameters for the output format
    int width = 12;
    int prec = 4;

    double x = x_min;
    double y = x*x;

    // generating dataset
    for (int i = 0; i <= n_steps; i++){
        // write a line with the current x and y values
        // nicely formatted to file
        ofile << std::setw(width) << std::setprecision(prec) << std::scientific << x << std::setw(width) << std::setprecision(prec) << std::scientific << y << std::endl;
        
        x += h;
        y = x*x;
    }

    ofile.close();

    return 0;
}