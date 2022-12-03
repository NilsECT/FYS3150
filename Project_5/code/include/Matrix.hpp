#ifndef __Matrix__
#define __Matrix__
#include <armadillo>
#include <iostream>

class Matrix {

protected:

    arma::cx_dmat create_tri(arma::cx_dvec &a, const double r, const int len, const int i);

    arma::cx_dmat create_rdiag(const double r, const int len);

public:

    Matrix() = default;

    static int pair_to_single(const int i, const int j, const int len=0);

    std::vector<arma::cx_dmat*> create_AB(arma::mat &V, const double h, const double dt, const int M);

    std::tuple<int, int> single_to_pair(const int k, const int len);

    arma::cx_dmat* create_mat(arma::cx_dvec &a, const double r, const int len);

    arma::cx_dvec mult_Bu(arma::cx_dvec &u, arma::cx_dmat &B);

    int gauss_seidel(arma::cx_dmat* mat, arma::cx_dvec* b, arma::cx_dvec* u_old, arma::cx_dvec* u_new, double criteria);

    arma::cx_dmat reshape(arma::cx_dvec &u, int M);
};

#endif