#include <armadillo>
#include <iostream>
#include <iostream>

int main() {

    std::complex<double> i(0., 1.);
    double T = 50;
    double n_steps = 10000;
    double timestep = T/n_steps;

    arma::vec t = arma::linspace(0, T, n_steps);
    arma::vec x = arma::vec(n_steps);
    arma::vec y = arma::vec(n_steps);
    arma::vec z = arma::vec(n_steps);
    arma::vec r = arma::vec{20.,0.,20.}; // micro m
    arma::vec v = arma::vec{0.,25.,0.};
    x[0] = r[0];
    y[0] = 0.;
    z[0] = r[2];

    double q = 1.0;
    double V_0 = 9.65 * std::pow(10, 8);
    double B_0 = 9.65*10;
    double m = 40.078;
    double d = 1 * std::pow(10, 4);


    double omega0 = q * B_0 / m;
    double omega_z2 = 2 * q * V_0 / (m * d*d);

    double omegap = (omega0 + std::sqrt(omega0*omega0 - 2 * omega_z2) / 2);
    double omegam = (omega0 - std::sqrt(omega0*omega0 - 2 * omega_z2) / 2);

    double Ap = (v[1] + omegam * x[0]) / (omegam - omegap);
    double Am = - (v[1] + omegap * x[0]) / (omegam - omegap);

    std::complex<double> f;

    std::ofstream r_outfile;
    r_outfile.open("analytical.txt", std::ios_base::app); // append instead of overwrite
    // legg til 0
    r_outfile << x[0] << "   " << y[0] << "   " << z[0];
    r_outfile << "\n";

    for (int j = 1; j < n_steps; j++) {
        f = Ap * std::exp(- i * (omegap * (j * timestep))) + Am * std::exp(- i * (omegam * (j * timestep)));

        r_outfile << std::real(f) << "   " << std::imag(f) << "   " << z[0] * std::cos(std::sqrt(omega_z2) * j * timestep);

        r_outfile << "\n";
    }



    return 0;
}