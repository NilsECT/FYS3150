#include "Jacobi.hpp"
#include <armadillo>
#include <iostream>

/**
 * @brief Construct a new Jacobi:: Jacobi object.
 * 
 * @param matrix A which you want to diagonalize.
 * 
 */
Jacobi::Jacobi(arma::mat &matrix) {
    this->N = matrix.n_cols;
    this->A = matrix;
    this->max_val = matrix(1, 2);

    this->find_k_l();
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
 * @param Uses the private variables of the object.
 * 
 * @return updates the private variable tau.
 */
void Jacobi::compute_tau() {
    this->tau = ( this->A(this->l, this->l)- this->A(this->k, this->k) ) / ( 2 * this->A(this->k,this->l) );
}

/**
 * @brief Having computed tau for the jacobi rotation matrix (compute_tau()), calculates the tan value need to find cosine and sine.
 * 
 * @param Takes the private tau variable within the object.
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
 * @param Takes the private tan variable within the object.
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
 * @param Takes the private cos and tan variables within the object.
 * 
 * @return Updates the private sin variable.
 * 
 */
void Jacobi::compute_sin() {
    this->sin = this->cos * this->tan;
}