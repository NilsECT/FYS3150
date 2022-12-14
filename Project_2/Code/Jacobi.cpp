#include "Jacobi.hpp"
#include <armadillo>
#include <iostream>
#include <iomanip>
#include <unistd.h>
#include <assert.h>

/**
 * @brief Construct a new Jacobi:: Jacobi object. Finds k and l immediately.
 * 
 * @param matrix A which you want to diagonalize.
 *
 */
Jacobi::Jacobi(arma::mat &matrix) {
    // runs test function
    this->test_find_k_l();

    this->N = matrix.n_cols;
    this->A = matrix;
    this->find_k_l();

    // initializes S as the I-matrix
    this->S = arma::mat(this->N, this->N).fill(0.0);
    this->S.diag() += 1.;
}

/**
 * @brief Diagonalizes the matrix A stored in the object using Jacobi's Rotation Method. Stops when the highest off-diagonal value is <= tol.
 *
 * @param double tol Tolerance for the highest off-diagonal value of A once diagonalized.
 *
 * @return Number of transformations used to solve eq. Stores the updated diagonalised matrix A and the matrix containing the eigenvectors in the matrix S. They can be extracted through the get_A() and get_S() methods.
 */
void Jacobi::solve(double tol) {
    while (this->max_val > tol) {
        this->compute_trig();
        this->update_S();
        this->find_k_l();
        sim_trans++;
    }
}

/**
 * @brief Given a matrix to diagonalize, find the position (k, l) of the highest off-diagonal value in the upper triagonal part of a symmetric matrix.
 *
 * @param Nothing Uses the internally stored variables of the object
 *
 * @return Updates the k and l values.
 *
 */
void Jacobi::find_k_l() {
    this->k = 0;
    this->l = 1;
    this->max_val = std::abs(this->A(this->k, this->l));

    for(int i = 0; i < this->N; i++){

        for(int j = i+1; j < this->N; j++){

            if (std::abs(this->A(i,j)) > this->max_val){

                this->max_val = std::abs(this->A(i,j));
                this->k = i;
                this->l = j;
            }
        }
    }
}

/**
 * @brief Sets the A matrix to be an arbitrary symmetric matric given in problem set 2 and checks whether find_k_l() finds correct k, l
 *
 * @param Nothing
 *
 * @return Nothing. Causes assertion error and terminates program if find_k_l() fails
 *
 */

void Jacobi::test_find_k_l() {
  arma::mat test_mat = arma::mat(4, 4, arma::fill::eye);
  test_mat(1, 2) = -0.7;
  test_mat(2, 1) = -0.7;
  test_mat(0, 3) = 0.5;
  test_mat(3, 0) = 0.5;

  this->set_A(test_mat);

  assert(this->k == 1 && this->l == 2);
}

/**
 * @brief Sets a temporary rotation matrix from the cosine and sine and rotates A locally. Then updates the S matrix that will eventually contain the eigenvectors.
 *
 * @param Nothing Uses the matrices stored within the object. k, l, cos and sin have to be up to date.
 *
 * @return Updates A and S locally.
 *
 */
void Jacobi::update_S() {

    // update A
    this->update_A();

    // calculates the total S
    for (int i = 0; i < this->N; i++) {
        double temp_ik = S(i, k);

        this->S(i, this->k) = temp_ik * this->cos - this->S(i, this->l) * this->sin;
        this->S(i, this->l) = this->S(i, this->l) * this->cos + temp_ik * this->sin;
    }

}

/**
 * @brief Rotates the matrix A locally given the rotation matrix Sm.
 * 
 * @return Updates A locally.
 */
void Jacobi::update_A() {

    // updates the places where we want to diagonalise and the diagonal entries
    double temp_kk = this->A(this->k, this->k);
    double temp_ll = this->A(this->l, this->l);

    this->A(this->k, this->k) = temp_kk * this->cos * this->cos - 2. * this->A(this->k, this->l) * this->cos * this->sin + temp_ll * this->sin * this->sin;
    this->A(this->l, this->l) = temp_ll * this->cos * this->cos + 2. * this->A(this->k, this->l) * this->cos * this->sin + temp_kk * this->sin * this->sin;

    this->A(this->k, this->l) = 0.;
    this->A(this->l, this->k) = 0.;

    // updates the rest of the matrix
    for (int i = 0; i < this->N; i++) {

        // we don't want to update kk, ll, kl, lk.
        if ((i == this->k) || (i == this->l)) {
            continue;
        }

        double temp_ik = this->A(i, this->k);
        double temp_il = this->A(i, this->l);


        this->A(i, this->k) = temp_ik * this->cos - temp_il * this->sin;

        this->A(this->k, i) = this->A(i, this->k);

        this->A(i, this->l) = temp_il * this->cos + temp_ik * this->sin;

        this->A(this->l, i) = this->A(i, this->l);
    }
}

/**
 * @brief Computes the cos and sin values to set in the Jacobi rotation matrix and stores them within the object.
 *
 * @param Nothing Uses the private variables within the object. Must have found k and l!
 *
 * @return Updates the values for tau, tan, cos and sin within the object.
 *
 */
void Jacobi::compute_trig() {
    this->compute_tau();
    this->compute_tan();
    this->compute_cos();
    this->compute_sin();
    assert(is_valid_double(this->tau));
    assert(is_valid_double(this->tan));
    assert(is_valid_double(this->cos));
    assert(is_valid_double(this->sin));
}

/**
 * @brief Having found k and l for the highest value in the off diagonal entries of a symmetric matrix A. Computes the tau needed for the jacobi rotation matrix.
 *
 * @param Nothing Uses the private variables of the object.
 *
 * @return updates the private variable tau.
 */
void Jacobi::compute_tau() {
    this->tau = ( this->A(this->l, this->l)- this->A(this->k, this->k) ) / ( 2. * this->A(this->k,this->l) );
}

/**
 * @brief Having computed tau for the jacobi rotation matrix (compute_tau()), calculates the tan value need to find cosine and sine.
 *
 * @param Nothing Takes the private tau variable within the object.
 *
 * @return Updates the private tan variable.
 *
 */
void Jacobi::compute_tan() {
    if (this->tau > 0.){
        this->tan = 1. / ( this->tau + std::sqrt(1. + this->tau * this->tau) );
    }
    else {
        this->tan = -1. / ( -this->tau + std::sqrt(1. + this->tau * this->tau) );
    }
}

/**
 * @brief With the found tan value (compute_tau() then compute_tan()) computes the cosine value for the jacobi rotation matrix.
 *
 * @param Nothing Takes the private tan variable within the object.
 *
 * @return Updates the private cos variable.
 *
 */
void Jacobi::compute_cos() {
    this->cos = 1. / std::sqrt(1. + this->tan * this->tan);
}

/**
 * @brief With the found cos value (compute_tau() then compute_tan() then compute_cos()) computes the sine value for the jacobi rotation matrix.
 *
 * @param Nothing Takes the private cos and tan variables within the object.
 *
 * @return Updates the private sin variable.
 *
 */
void Jacobi::compute_sin() {
    this->sin = this->cos * this->tan;
}

// get and set are reserved for the bottom!

/**
 * @brief Functionality to extract the matrix A.
 *
 * @return arma::mat A
 */
arma::mat Jacobi::get_A() {
    return this->A;
}

/**
 * @brief Functionality to extract the eigenvectors. Every column is an eigenvector
 *
 * @return arma::mat eigenvectors (S in the object)
 */
arma::mat Jacobi::get_eigvec() {
    return this->S;
}

/**
 * @brief Functionality to extract the eigenvalues of A once solved. It's the diagonal entries of A after having used the Jacobi Rotation Method.
 *
 * @return arma::vec of eigenvalues.
 */
arma::vec Jacobi::get_eigval() {
    return this->A.diag();
}

/**
 * @brief Replaces the matrix A contained within the object with the one provided as a parameter. Usefull if you don't want to create a new object for every matrix you want to diagonalized.
 *
 * @param arma::mat matrix to diagonalize.
 *
 */
void Jacobi::set_A(arma::mat matrix) {
    this->N = matrix.n_cols;
    this->A = matrix;
    this->find_k_l();

    // initializes S as the I-matrix
    this->S = arma::mat(this->N, this->N).fill(0.0);
    this->S.diag() += 1.;

    sim_trans = 0;
}

int Jacobi::get_N() {
    return this->N;
}

int Jacobi::get_k() {
    return this->k;
}

int Jacobi::get_l() {
    return this->l;
}

double Jacobi::get_max_val() {
    return this->max_val;
}

/**
 * @brief Gives the number of similarity transformations done when diagonalizing A.
 *
 * @return int
 */
int Jacobi::trans_count() {
    return sim_trans;
}

bool Jacobi::is_valid_double(double x) {
    return x*0.0==0.0;
}

int Jacobi::get_sim_trans(){
    return this->sim_trans;
}