#ifndef __TDS_hpp__
#define __TDS_hpp__

class TriDiagS{
    private:
    // Declare private variables
    arma::vec g, x, v, subdiag, maindiag, supdiag;
    arma::mat A;
    int n_steps = 100;
    // n is the length of the main diagonal of A
    // N is the max range of x from 0
    int n, N;
    double v0, vn, h;
    bool fillMat;

    public:
    // Declaration of constructors and whatnot
    TriDiagS(int N, int n_steps){
        this->N = N;
        this->n_steps = n_steps;
        h = N/n_steps;
        x = arma::linspace(0, N, n_steps+1);
    }

    // CHECK WHAT V0 AND VN SHOULD BE, NOW THEY ARE THE FUNCTION EVALUATED AT THE ENDPOINTS

    // MAKE SOLVER THAT GIVES X AND V
    arma::mat solve(double (*f)(int), bool fillMat = false){
        if (fillMat) {
            // functionality not added yet
            exit(0);
        }
        else{
            symDiag();
            make_g_v(f);
            forward_elim();
            backward_elim();

            // print out x and v
        }
    }

    // forward elimination
    void forward_elim(){
        for (int i = 1; i < n; i++){

            // setting tilde(b and g)
            maindiag(i) = maindiag(i) - supdiag(i-1) * ( subdiag(i)/maindiag(i-1) );
            g(i) = g(i) - g(i-1) * ( subdiag(i)/maindiag(i-1) );
        }
    }

    // backward elimination
    void backward_elim(){

        v(n-1) = g(n-1)/maindiag(n-1);

        for (int i = n-2; i >= 0; i--){
            v(i) = (g(i) - supdiag(i) * v(i+1)) / maindiag(i);
        }
    }

    // makes the g vector
    void make_g_v(double (*f)(int)){
        set_initial(f(0), f(N));

        g = arma::vec(n);
        v = arma::vec(n);

        // filling in the values of g
        for (int i = 0; i < n; i++){

            // start and end of g are a little different
            if (i == 0){
                g(i) = v0 + h*h*f(x(i+1));
            }
            else if (i == n-1){
                g(i) = vn + h*h*f(x(i+1));
            }
            else{
                g(i) = h*h*f(x(i+1));
            }
        }
    }

    // sets initial conditions [start, end]
    void set_initial(double start, double end){
        v0 = start;
        vn = end;
    }

    void symDiag(){
        // WARNING
        // MINSTAKE IN LABELING n WILL BE HERE
        // WARNING
        n = n_steps - 1;

        // assuming it makes n number of entries if not I'm f****d, 
        // keep it PG y'kno
        maindiag = arma::vec(n).fill(2);
        supdiag, subdiag = arma::vec(n-1).fill(-1);
    }

    // make matrix A
    void set_diagonals(arma::vec subdiag, arma::vec maindiag, arma::vec supdiag){
        // setting the diagonals
        set_subdiag(subdiag);
        set_maindiag(maindiag);
        set_supdiag(supdiag);
        n = this->maindiag.size();
    }

    // filling in values of A
    void fill_A(){
        for (int i = 0; i < n; i++){
            
            //placing the maindiagonal entries
            A(i, i) = maindiag(i);
            
            // placing the supdiagonal entries
            if (i < (n-1)){
                A(i, i+1) = supdiag(i);
            }

            // placing the subdiagonal entries
            if (i > 0){
                A(i, i-1) = subdiag(i);
            }
        }
    }

    // finds the different diagonals
    void set_subdiag(arma::vec diag){
        subdiag = diag;
    }

    void set_maindiag(arma::vec diag){
        maindiag = diag;
    }

    void set_supdiag(arma::vec diag){
        supdiag = diag;
    }
};

#endif