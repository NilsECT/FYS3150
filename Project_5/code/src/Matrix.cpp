#include "Matrix.hpp"
#include <armadillo>
#include <iostream>

/**
 * @brief Converts matrix coordinates to respective vector coordinates.
 * 
 * @param i Row coordinate.
 * @param j Column coordianate.
 * @param len Number of rows/columns in the matrix (square).
 * @return int Vector coordinate.
 */
int Matrix::pair_to_single(const int i, const int j, const int len) {
    // len = M - 2
    return (j) * len + (i);
}

/**
 * @brief Converts vector coordinates to respective matrix coordinates.
 * 
 * @param k vector coordinate.
 * @param len Number of rows/columns in the matrix (square).
 * @return std::tuple<int, int> Tuple containing the (row, column) coordinates.
 */
std::tuple<int, int> Matrix::single_to_pair(const int k, const int len) {
    int i = k % len;
    int j = k / len;
    
    return std::tuple<int, int>{i, j};
}

/**
 * @brief Creates a tridaiagonal submatrix to implement in the mother matrix.
 * 
 * @param a Vector containing the values to be set in the main diagonal.
 * @param r Value to be set in the lower and upper diagonals.
 * @param len Number of rows/columns in the matrix (square).
 * @param i Current column for which we are creating the subatrix.
 * @return arma::sp_cx_dmat Sparse submatrix to be used.
 */
arma::sp_cx_dmat Matrix::create_tri(arma::cx_dvec &a, const double r, const int len, const int i) {
    arma::sp_cx_dmat temp = arma::sp_cx_dmat(len, len);

    for (int ii = 0; ii<len; ii++) {
        temp(ii, ii) = a(pair_to_single(ii, i, len));
        temp.diag(1).fill(r);
        temp.diag(-1).fill(r);
    }

    return temp;
}

/**
 * @brief Creates a diagonal matrix where the diagonal contains the value r.
 * 
 * @param r Value to be set in the diagonals.
 * @param len Number of rows/columns in the matrix (square).
 * @return arma::sp_cx_dmat Sparse submatrix to be used.
 */
arma::sp_cx_dmat Matrix::create_rdiag(const double r, const int len) {
    arma::sp_cx_dmat temp = arma::sp_cx_dmat(len, len);

    temp.diag().fill(r);

    return temp;
}

/**
 * @brief Creates and returns a sparse matrix as intended for the Crank-Nicolson approach.
 * 
 * @param a Vector containing the values to be inserted in the main diagonal.
 * @param r Value to be set in the lower and upper diagonals.
 * @param len Number of rows/columns in the wave packet matrix (square).
 * @return arma::sp_cx_dmat 
 */
arma::sp_cx_dmat Matrix::create_mat(arma::cx_dvec &a, const double r, const int len) {
    int lenlen = len*len;
    arma::sp_cx_dmat A = arma::sp_cx_dmat(lenlen, lenlen); // L+U+D, will be named mother matrix

    // filling and placing main diagonal matrices in mother
    for (int i = 0; i < len; i++) {
        int start = len*i;
        int end = len*i + (len-1);
        A.submat(start, start, end, end) = create_tri(a, r, len, i);
    }

    // filling and placing the off-diagonal matrices in mother
    for (int i = 1; i<(len); i++) {
        int start = len*i;
        int end = len*i + (len-1);
        A.submat(start-len, start, end-len, end) = create_rdiag(r, len);
        A.submat(start, start-len, end, end-len) = create_rdiag(r, len);
    }

    return A;
}

/**
 * @brief Creates and returns the need matrices A and B for the Crank-Nicolson approach.
 * 
 * @param V Matrix containing the potential of the system.
 * @param h Stepsize in the spacial dimensions x and y.
 * @param dt Stepsize in the time dimension.
 * @param M Size of the complete matrix for the system, including the boundaries.
 * @return std::vector<arma::sp_cx_dmat> 
 */
std::vector<arma::sp_cx_dmat> Matrix::create_AB(arma::mat &V, const double h, const double dt, const int M) {
    int len = M-2;
    double r = dt/(h*h);
    std::complex<double> im(0., 1.);
    arma::cx_dvec a = arma::cx_dvec(len*len);
    arma::cx_dvec b = arma::cx_dvec(len*len);

    for (int ii=0; ii < (len); ii++) {
        arma::vec temp_col = V.col(ii);

        for (int i=0; i < (len); i++) {
            std::complex<double> temp_im = im * ((dt*temp_col(i)) / 2.);

            a(pair_to_single(i, ii, len)) = 1. + 4.*r + temp_im;
            b(pair_to_single(i, ii, len)) = 1. - 4.*r - temp_im;
        }
    }

    return std::vector<arma::sp_cx_dmat>{create_mat(a, -r, len), create_mat(b, r, len)};
}

/**
 * @brief Multiplies matrices with values along diagonals -3, -1, 0, 1, 3 with a vector.
 * 
 * @param u Vector.
 * @param B Matrix.
 * @return arma::cx_dvec Result of the multiplication.
 */
arma::cx_dvec Matrix::mult_Bu(arma::cx_dvec &u, arma::sp_cx_dmat &B) {
    // trying to avoid multiplying many times by zero
    // only take the diagonals from B that matter which are 0, 1, 3, -1, and -3
    // perform Hadamard product of the diagonals with their respective parts of the u vector

    int len = u.size();

    arma::cx_dvec b = arma::cx_dvec(len);
    std::complex<double> temp_u;
    std::complex<double> temp_mat;

    // I have checked that these loops should work, by hand
    // double check that it does indeed work through an example
    // mutliplying the main diagonal of B with u
    // operator % should be multiplying elementwise
    b = B.diag() % u;

    // muliplying the lower and upper diagonals 1 and -1
    for (int i = 1; i < len; i++) {
        // diag 1
        temp_u = u(i);
        
        temp_mat = B.diag(1)(i-1);
        b(i-1) += temp_mat * temp_u;

        // diag -1
        temp_u = u(i-1);
        
        temp_mat = B.diag(-1)(i-1);
        b(i) += temp_mat * temp_u;
    }

    // muliplying the lower and upper diagonals 3 and -3
    for (int i = 3; i < len; i++) {
        // diag 3
        temp_u = u(i);
        
        temp_mat = B.diag(3)(i-3);
        b(i-3) += temp_mat * temp_u;

        // diag -3
        temp_u = u(i-3);
        
        temp_mat = B.diag(-3)(i-3);
        b(i) += temp_mat * temp_u;
    }

    return b;
}

// NOT TESTED YET
// need to figure out how to fix convergence
int Matrix::gauss_seidel(arma::sp_cx_dmat* mat, arma::cx_dvec* b, arma::cx_dvec* u_old, arma::cx_dvec* u_new, double criteria) {
    // see page 191 in Morten's lecture notes


    // proposal: change u_old to be a reference to the variable you want to update
    int lenlen = u_old->size();

    *u_new = *b;

    // diag 1 and -1
    for (int i = 1; i < lenlen; i++) {
        std::complex<double> temp_u = (u_old->at(i));
        std::complex<double> temp_mat = mat->diag(1)(i-1);
        u_new->at(i-1) -= (temp_mat * temp_u);

        temp_u = u_new->at(i-1);
        temp_mat = mat->diag(-1)(i-1);
        u_new->at(i) -= (temp_mat * temp_u);
    }

    // diag 3 and -3
    for (int i = 3; i < lenlen; i++) {
        std::complex<double> temp_u = (u_old->at(i));
        std::complex<double> temp_mat = mat->diag(3)(i-3);
        u_new->at(i-3) -= (temp_mat * temp_u);

        temp_u = u_new->at(i-3);
        temp_mat = mat->diag(-3)(i-3);
        u_new->at(i) -= (temp_mat * temp_u);
    }

    // main diag
    // operator / should be doing the operations elementwise
    *u_new /= mat->diag();

    // calc residual
    // Ax - b must be less than a tolerance
    // risky play but let's try it
    // if convergence criteria is met, then return the result
    // if not do it again, recursive style
    *u_old = *u_new;

    return (arma::approx_equal(this->mult_Bu(*u_new, *mat), *b, "reldiff", criteria)) ? 0 : gauss_seidel(mat, b, u_old, u_new, criteria); // golfing
}

arma::cx_dmat Matrix::reshape(arma::cx_dvec &u, int M) {
    arma::cx_dmat u_mat = arma::cx_dmat(M-2, M-2);

    for (int j=0; j < M-2; j++) {
        for (int i=0; i < M-2; i++) {
            u_mat(i, j) = u(Matrix::pair_to_single(i, j, M-2));
        }
    }

    return u_mat;
}
