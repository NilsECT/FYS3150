#include "Jacobi.hpp"
#include <armadillo>
#include <iostream>

/**
 * @brief Construct a new Jacobi:: Jacobi object.
 * 
 * @param matrix A which you want to diagonalize.
 * 
 */
Jacobi::Jacobi(arma::mat matrix) {
    this->N = matrix.n_cols;
    this->A = matrix;
    this->max_val = matrix(1, 2);

    // initializes S as the I-matrix
    this->S.eye(this->N, this->N);

    // starts by finding the first k_l, haven't written solve() yet.
    this->find_k_l();
}

/**
 * @brief Diagonalizes the matrix A stored in the object using Jacobi's Rotation Method. Stops when the highest off-diagonal value is <= tol.
 * 
 * @param tol Tolerance for the highest off-diagonal value of A once diagonalized.
 * 
 * @return Nothing. Stores the updated diagonalised matrix A and the matrix containing the eigenvectors in the matrix S. They can be extracted through the get_A() and get_S() methods.
 */
void Jacobi::solve(double tol) {
    while (this->max_val > tol) {
        this->compute_trig();
        this->update_S();
        this->find_k_l();
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

    for(int i = 0; i < this->N; i++){

        for(int j = i+1; j < this->N; j++){

            if (std::abs(A(i,j)) > this->max_val){

                this->max_val = A(i,j);
                this->k = i;
                this->l = j;
            }
        }
    }
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

    // to not call these a million times
    int temp_k = this->k;
    int temp_l = this->l;
    double cosine = this->cos;
    double sine = this->sin;

    // creates a temporary rotation matrix Sm (called temp)
    arma::mat temp;
    temp.eye(this->N, this->N);
    temp(temp_k, temp_k) = cosine;
    temp(temp_k, temp_l) = sine;
    temp(temp_l, temp_k) = - sine;
    temp(temp_l, temp_l) = cosine;

    // update A
    this->update_A(temp);

    // calculates the total S
    for (int i = 0; i < this->N; i++) {
        double temp_ik = S(i, k);

        this->S(i, temp_k) = temp_ik * cosine - this->S(i, temp_l) * sine;
        this->S(i, temp_l) = this->S(i, temp_l) * cosine + temp_ik * sine;
    }

}

/**
 * @brief Rotates the matrix A locally given the rotation matrix Sm.
 * 
 * @param Sm The rotation matrix found from the position of the highest off-diagonal value in A.
 * 
 * @return Updates A locally.
 */
void Jacobi::update_A(arma::mat &Sm) {
    int temp_k = this->k;
    int temp_l = this->l;

    // updates the places where we want to diagonalise and the diagonal entries
    double temp_kk = this->A(temp_k, temp_k);
    double temp_ll = this->A(temp_l, temp_l);

    double cosine = this->cos;
    double sine = this->sin;

    this->A(k, k) = temp_kk * cosine * cosine - 2 * this->A(k, l) * cosine * sine + temp_ll * sine * sine;
    this->A(l, l) = temp_ll * cosine * cosine + 2 * this->A(k, l) * cosine * sine + temp_kk * sine * sine;

    this->A(k, l) = 0;
    this->A(l, k) = 0;

    // updates the rest of the matrix
    for (int i = 0; i < this->N; i++) {

        // we don't want to update kk, ll, kl, lk.
        if ((i == temp_k) || (i == temp_l)) {
            continue;
        }

        double temp_ik = this->A(i, temp_k);
        double temp_il = this->A(i, temp_l);

        this->A(i, temp_k) = temp_ik * cosine - temp_il * sine;
        this->A(temp_k, i) = this->A(i, temp_k);

        this->A(i, temp_l) = temp_il * cosine + temp_ik * sine;
        this->A(temp_l, i) = this->A(i, temp_l);
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
}

/**
 * @brief Having found k and l for the highest value in the off diagonal entries of a symmetric matrix A. Computes the tau needed for the jacobi rotation matrix.
 * 
 * @param Nothing Uses the private variables of the object.
 * 
 * @return updates the private variable tau.
 */
void Jacobi::compute_tau() {
    this->tau = ( this->A(this->l, this->l)- this->A(this->k, this->k) ) / ( 2 * this->A(this->k,this->l) );
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
    if (this->tau > 0){
        this->tan = 1 / ( this->tau + std::sqrt(1 + this->tau * this->tau) );
    }
    else {
        this->tan = -1 / ( -this->tau + std::sqrt(1 + this->tau * this->tau) );
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
    this->cos = 1 / std::sqrt(1 + this->tan * this->tan);
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
    this->max_val = matrix(1, 2);

    // initializes S as the I-matrix
    this->S.eye(this->N, this->N);

    // starts by finding the first k_l, haven't written solve() yet.
    this->find_k_l();
}