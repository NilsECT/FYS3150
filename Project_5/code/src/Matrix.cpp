#include "Matrix.hpp"
#include <armadillo>
#include <iostream>


int Matrix::pair_to_single(const int i, const int j, const int len) {
    // len = M - 2
    return (j) * len + (i);
}

std::tuple<int, int> Matrix::single_to_pair(const int k, const int len) {
    int i = k % len;
    int j = k / len;
    
    return std::tuple<int, int>{i, j};
}

arma::sp_cx_dmat Matrix::create_tri(arma::cx_dvec &a, const double r, const int len, const int i) {
    arma::sp_cx_dmat temp = arma::sp_cx_dmat(len, len);

    for (int ii = 0; ii<len; ii++) {
        temp.at(ii, ii) = a(pair_to_single(ii, i, len));
        temp.diag(1).fill(r);
        temp.diag(-1).fill(r);
    }

    return temp;
}

arma::sp_cx_dmat Matrix::create_rdiag(const double r, const int len) {
    arma::sp_cx_dmat temp = arma::sp_cx_dmat(len, len);

    temp.diag().fill(r);

    return temp;
}

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

            a.at(pair_to_single(i, ii, len)) = 1. + 4.*r + temp_im;
            b.at(pair_to_single(i, ii, len)) = 1. - 4.*r - temp_im;
        }
    }

    return std::vector{create_mat(a, -r, len), create_mat(b, r, len)};
}

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
    std::cout << "Matrix 94" << std::endl;
    // mutliplying the main diagonal of B with u
    // operator % should be multiplying elementwise
    b = B.diag() % u;

    // muliplying the lower and upper diagonals 1 and -1
    for (int i = 1; i < len; i++) {
        // diag 1
        std::cout << "Matrix 102" << std::endl;
        temp_u = u.at(i);
        std::cout << "Matrix 104" << std::endl;
        temp_mat = B.diag(1).at(i-1);
        b.at(i-1) += temp_mat * temp_u;

        // diag -1
        std::cout << "Matrix 109" << std::endl;
        temp_u = u.at(i-1);
        std::cout << "Matrix 111" << std::endl;
        temp_mat = B.diag(-1).at(i-1);
        b.at(i) += temp_mat * temp_u;
    }

    // muliplying the lower and upper diagonals 3 and -3
    for (int i = 3; i < len; i++) {
        // diag 1
        std::cout << "Matrix 119" << std::endl;
        temp_u = u.at(i);
        std::cout << "Matrix 121" << std::endl;
        temp_mat = B.diag(3).at(i-3);
        b.at(i-3) += temp_mat * temp_u;

        // diag -1
        std::cout << "Matrix 126" << std::endl;
        temp_u = u.at(i-3);
        std::cout << "Matrix 128" << std::endl;
        temp_mat = B.diag(-3).at(i-3);
        b.at(i) += temp_mat * temp_u;
    }

    return b;
}

// NOT TESTED YET
// need to figure out how to fix convergence
// int Matrix::gauss_seidel(arma::sp_cx_dmat* mat, arma::cx_dvec* b, arma::cx_dvec* u_old, arma::cx_dvec* u_new, double criteria) {
//     // see page 191 in Morten's lecture notes


//     // proposal: change u_old to be a reference to the variable you want to update
//     int lenlen = u_old->size();

//     *u_new = *b;

//     // diag 1 and -1
//     for (int i = 1; i < lenlen; i++) {
//         std::complex<double> temp_u = (u_old->at(i));
//         std::complex<double> temp_mat = mat->diag(1).at(i-1);
//         u_new->at(i-1) -= (temp_mat * temp_u);

//         temp_u = u_new->at(i-1);
//         temp_mat = mat->diag(-1).at(i-1);
//         u_new->at(i) -= (temp_mat * temp_u);
//     }

//     // diag 3 and -3
//     for (int i = 3; i < lenlen; i++) {
//         std::complex<double> temp_u = (u_old->at(i));
//         std::complex<double> temp_mat = mat->diag(3).at(i-3);
//         u_new->at(i-3) -= (temp_mat * temp_u);

//         temp_u = u_new->at(i-3);
//         temp_mat = mat->diag(-3).at(i-3);
//         u_new->at(i) -= (temp_mat * temp_u);
//     }

//     // main diag
//     // operator / should be doing the operations elementwise
//     *u_new /= mat->diag();

//     // calc residual
//     // Ax - b must be less than a tolerance
//     // risky play but let's try it
//     // if convergence criteria is met, then return the result
//     // if not do it again, recursive style
//     *u_old = *u_new;

//     arma::cx_dmat temp = arma::cx_dmat(lenlen, 5);

//     for (int i=0; i < lenlen, i++) {
//         temp.at(0) = mat->col(i) * u_new->at(i);
//     }

//     for

//     return (arma::approx_equal(arma::sum(mat->each_row()%u_new->as_row(), 1), *b, "reldiff", criteria)) ? 0 : gauss_seidel(mat, b, u_old, u_new, criteria); // golfing
// }

arma::cx_dmat Matrix::reshape(arma::cx_dvec &u, int M) {
    arma::cx_dmat u_mat = arma::cx_dmat(M-2, M-2);

    for (int j=0; j < M-2; j++) {
        for (int i=0; i < M-2; i++) {
            u_mat.at(i, j) = u.at(Matrix::pair_to_single(i, j, M-2));
        }
    }

    return u_mat;
}
