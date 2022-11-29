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

arma::mat create_tri(arma::vec &a, const int r, const int len, const int i) {
    arma::mat temp = arma::mat(len, len);

    for (int ii = 0; ii<len; ii++) {
        temp.at(ii, ii) = a(pair_to_single(ii, i, len));
        temp.diag(1).fill(r);
        temp.diag(-1).fill(r);
    }

    return temp;
}

arma::mat create_rdiag(const int r, const int len) {
    arma::mat temp = arma::mat(len, len);

    temp.diag().fill(r);

    return temp;
}

std::vector<arma::mat*> ALUD(arma::vec &a, const int r, const int len) {
    int lenlen = len*len;
    arma::mat *A = new arma::mat(lenlen, lenlen); // L+U+D, will be named mother matrix
    arma::mat *L = new arma::mat(lenlen, lenlen); // lower triag
    arma::mat *U = new arma::mat(lenlen, lenlen); // upper triag
    arma::mat *D = new arma::mat(lenlen, lenlen); // diagonal

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

    arma::mat lower = arma::trimatl(*A, -1); // lower triag
    arma::mat upper = arma::trimatu(*A, -1); // upper triag
    arma::mat diag = A->diag(); // diag of A

    // fill lower traingular matrix
    L->diag(-1) = A->diag(-1);
    L->diag(-3) = A->diag(-3);

    // fill upper traingular matrix
    U->diag(1) = A->diag(1);
    U->diag(3) = A->diag(3);

    // fill diagonal matrix
    D->diag() = A->diag();

    return std::vector<arma::mat*>{A, L, U, D};
}