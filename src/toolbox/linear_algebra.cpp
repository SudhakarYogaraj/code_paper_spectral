#include "toolbox/linear_algebra.hpp"
#include "global/templates.hpp"
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

std::mat symmetric(std::mat A) {
    return to_std(0.5 * to_arma(A).t() + 0.5 * to_arma(A));
}

arma::mat to_arma(const std::mat &A) {
    arma::mat B(A.size(), A[0].size());
    for (size_t i = 0; i < A.size(); ++i) {
        B.col(i) = conv_to<arma::vec>::from(A[i]);
    };
    return B;
}

std::mat to_std(const arma::mat &A) {
    std::mat B(A.n_rows);
    for (size_t i = 0; i < A.n_rows; ++i) {
        B[i] = conv_to< vector<double> >::from(A.col(i));
    };
    return B;
}

arma::vec to_arma_vec(const std::vector<double> &v) {
    return arma::conv_to< arma::vec >::from(v);
}

std::vec to_std_vec(const arma::vec &v) {
    typedef std::vector<double> std_vec;
    return arma::conv_to< std_vec >::from(v);
}

std::mat cholesky(std::mat A) {
    return to_std(chol(to_arma(A)));
}

std::mat transpose(std::mat A) {
    return to_std(to_arma(A).t());
}

vector<double> solve(std::mat A, std::vec b) {
    return to_std_vec(arma::solve(to_arma(A), to_arma_vec(b)));
}
