#include <iostream>
#include <armadillo>

doubel get_k_l(arma::mat A, int N){

    for(int i = 0; i<N; i++){

        for(int j = i+1; j<N; j++){

            if (std::abs(A(i,j))>a_max){

                double a_max = A(i,j);
                int k = i;
                int l = j;
            }
        }
    }
    return k, l;
}


double tau(int k, int l){
    tau = (A(ll)-A(kk))/(2*A(k,l));
    return tau;
}

double t(double tau){
    if (tau > 0){
        double t = 1/(tau+std::sqrt(1+tau*tau));
    }
    if (tau < 0){
        double t = -1/(-tau+std::sqrt(1+tau*tau));
    }
    return t;
}

double c(double t){
    c = 1/std::sqrt(1+t*t)
    return c;
}

double s(double c, double t){
    s = c*t
    return s
}

