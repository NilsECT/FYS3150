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

std::vector<arma::cx_dmat*> ALUD(arma::cx_dvec &a, const double r, const int len) {
    int lenlen = len*len;
    arma::cx_dmat *A = new arma::cx_dmat(lenlen, lenlen); // L+U+D, will be named mother matrix
    arma::cx_dmat *L = new arma::cx_dmat(lenlen, lenlen); // lower triag
    arma::cx_dmat *U = new arma::cx_dmat(lenlen, lenlen); // upper triag
    arma::cx_dmat *D = new arma::cx_dmat(lenlen, lenlen); // diagonal

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

    arma::cx_dmat lower = arma::trimatl(*A, -1); // lower triag
    arma::cx_dmat upper = arma::trimatu(*A, -1); // upper triag
    arma::cx_dmat diag = A->diag(); // diag of A

    // fill lower traingular matrix
    L->diag(-1) = A->diag(-1);
    L->diag(-3) = A->diag(-3);

    // fill upper traingular matrix
    U->diag(1) = A->diag(1);
    U->diag(3) = A->diag(3);

    // fill diagonal matrix
    D->diag() = A->diag();

    return std::vector<arma::cx_dmat*>{A, L, U, D};
}

std::vector<std::vector<arma::cx_dmat*>> create_AB(arma::mat &V, const double h, const double dt, const int M) {
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

    return std::vector{ALUD(a, -r, len), ALUD(b, r, len)};
}
