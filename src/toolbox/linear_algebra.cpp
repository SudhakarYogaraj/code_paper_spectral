#include "toolbox/linear_algebra.hpp"
#include "templates.hpp"
#include "io/io.hpp"

using namespace std;
using namespace arma;

double ipow(double x, int e) {
    if(e == 0) return 1;
    if (e == 1) return x;
    double aux = ipow(x,e/2);
    if (e%2 == 0) return aux*aux;
    return x*aux*aux;
}

vector< vector<double> > symmetric(vector< vector<double> > A) {
    return to_std(0.5 * to_arma(A).t() + 0.5 * to_arma(A));
}

mat to_arma(const vector< vector<double> > &A) {
    mat B(A.size(), A[0].size());
    for (size_t i = 0; i < A.size(); ++i) {
        B.col(i) = conv_to<vec>::from(A[i]);
    };
    return B;
}

vector< vector<double> > to_std(const mat &A) {
    vector< vector<double> > B(A.n_rows);
    for (size_t i = 0; i < A.n_rows; ++i) {
        B[i] = conv_to< vector<double> >::from(A.col(i));
    };
    return B;
}

arma::vec to_arma_vec(const std::vector<double> &v) {
    return arma::conv_to< vec >::from(v);
}

std::vector<double> to_std_vec(const arma::vec &v) {
    typedef std::vector<double> std_vec;
    return arma::conv_to< std_vec >::from(v);
}

vector< vector<double> > cholesky(vector< vector<double> > A) {
    return to_std(chol(to_arma(A)));
}

vector< vector<double> > transpose(vector< vector<double> > A) {
    return to_std(to_arma(A).t());
}

vector<double> solve(mat A, vec b) {
    return to_std_vec(arma::solve(A, b));
}
