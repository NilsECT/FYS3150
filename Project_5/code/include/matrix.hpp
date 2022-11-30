#include <armadillo>
#include <assert.h>
#include <iostream>

int pair_to_single(const int i, const int j, const int len=0) {
    // len = M - 2
    return (j) * len + (i);
}

std::tuple<int, int> single_to_pair(const int k, const int len) {
    int i = k % len;
    int j = k / len;
    
    return std::tuple<int, int>{i, j};
}

arma::cx_dmat create_tri(arma::cx_dvec &a, const double r, const int len, const int i) {
    arma::cx_dmat temp = arma::cx_dmat(len, len);

    for (int ii = 0; ii<len; ii++) {
        temp.at(ii, ii) = a(pair_to_single(ii, i, len));
        temp.diag(1).fill(r);
        temp.diag(-1).fill(r);
    }

    return temp;
}

arma::cx_dmat create_rdiag(const double r, const int len) {
    arma::cx_dmat temp = arma::cx_dmat(len, len);

    temp.diag().fill(r);

    return temp;
}

arma::cx_dmat* create_mat(arma::cx_dvec &a, const double r, const int len) {
    int lenlen = len*len;
    arma::cx_dmat *A = new arma::cx_dmat(lenlen, lenlen); // L+U+D, will be named mother matrix

    // filling and placing main diagonal matrices in mother
    for (int i = 0; i < len; i++) {
        int start = len*i;
        int end = len*i + (len-1);
        A->submat(start, start, end, end) = create_tri(a, r, len, i);
    }

    // filling and placing the off-diagonal matrices in mother
    for (int i = 1; i<(len); i++) {
        int start = len*i;
        int end = len*i + (len-1);
        A->submat(start-len, start, end-len, end) = create_rdiag(r, len);
        A->submat(start, start-len, end, end-len) = create_rdiag(r, len);
    }

    return A;
}

std::vector<arma::cx_dmat*> create_AB(arma::mat &V, const double h, const double dt, const int M) {
    int len = M-2;
    double r = dt/(h*h);
    std::complex<double> im(0, 1);
    arma::cx_dvec a = arma::cx_dvec(len*len);
    arma::cx_dvec b = arma::cx_dvec(len*len);

    for (int ii=0; ii < (len); ii++) {
        arma::vec temp_col = V.col(ii);

        for (int i=0; i < (len); i++) {
            std::complex<double> temp_im = im * ((dt*temp_col(i)) / 2);

            a.at(pair_to_single(i, ii, len)) = 1 + 4*r + temp_im;
            b.at(pair_to_single(i, ii, len)) = 1 - 4*r - temp_im;
            std::cout << i << std::endl;
        }
    }

    return std::vector{create_mat(a, -r, len), create_mat(b, r, len)};
}

arma::cx_dvec mult_Bu(arma::cx_dvec &u, arma::cx_dmat &B) {
    // trying to avoid multiplying many times by zero
    // only take the diagonals from B that matter which are 0, 1, 3, -1, and -3
    // perform Hadamard product of the diagonals with their respective parts of the u vector

    int len = u.size();

    arma::cx_dvec b = arma::cx_dvec(len);

    // I have checked that these loops should work, by hand
    // double check that it does indeed work through an example

    // mutliplying the main diagonal of B with u
    // operator % should be multiplying elementwise
    b = B.diag() % u;

    // muliplying the lower and upper diagonals 1 and -1
    for (int i = 1; i < len; i++) {
        // diag 1
        b.at(i-1) += B.diag(1).at(i-1) * u.at(i);

        // diag -1
        b.at(i) += B.diag(-1).at(i-1) * u.at(i-1);
    }

    // muliplying the lower and upper diagonals 3 and -3
    for (int i = 3; i < len; i++) {
        // diag 1
        b.at(i-3) += B.diag(3).at(i-3) * u.at(i);

        // diag -1
        b.at(i) += B.diag(-3).at(i-3) * u.at(i-3);
    }

    return b;
}

// NOT TESTED YET
arma::cx_dvec gauss_seidel(arma::cx_dmat &mat, arma::cx_dvec &b, arma::cx_dvec &u_old, double criteria) {
    // see page 191 in Morten's lecture notes
    int lenlen = u_old.size();
    int len = std::sqrt(lenlen);

    arma::cx_dvec temp_u = b;

    // diag 1 and -1
    for (int i = 1; i < lenlen; i++) {
        temp_u.at(i-1) -= (mat.diag(1).at(i-1)*u_old.at(i));
        temp_u.at(i) -= (mat.diag(-1).at(i-1)*temp_u.at(i-1));
    }

    // diag 3 and -3
    for (int i = 3; i < lenlen; i++) {
        temp_u.at(i-3) -= (mat.diag(3).at(i-3)*u_old.at(i));
        temp_u.at(i) -= (mat.diag(-3).at(i-3)*temp_u.at(i-3));
    }

    // main diag
    // operator / should be doing the operations elementwise
    temp_u /= mat.diag();

    // calc residual
    // Ax - b must be less than a tolerance
    // risky play but let's try it
    // if convergence criteria is met, then return the result
    // if not do it again, recursive style
    return (arma::approx_equal(arma::sum(mat.each_row()%temp_u.as_row(), 1), b, "reldiff", criteria)) ? temp_u : gauss_seidel(mat, b, temp_u, criteria); // golfing
}