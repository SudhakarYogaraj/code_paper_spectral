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

// Symmetric part of a matrix
vector< vector<double> > symmetric(vector< vector<double> > A) {
    int n = A.size();
    vector< vector <double> > result(n,vector<double>(n,0.));

    for (int i = 0.; i < n; i++) {
        for (int j = 0; j <= i; j++) {
            result[i][j] = 0.5*(A[i][j] + A[j][i]);
            result[j][i] = result[i][j];
        }
    }
    return result;
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

// Cholesky factorization of a matrix
vector< vector<double> > cholesky(vector< vector<double> > matrix) {
    mat A = to_arma(matrix);
    mat result = chol(A);
    return to_std(result);
}

vector< vector<double> > transpose(vector< vector<double> > A)
{
    int n = A.size();
    vector< vector <double> > T(n,vector<double>(n,0.));

    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            T[i][j] = A[j][i];
        }
    }
    return T;
}

vector<double> solve(vector< vector<double> > A, vector<double> b) {
    mat matrix = to_arma(A);
    vec vector = conv_to< vec >::from(b);
    vec solution = arma::solve(matrix, vector);
    typedef std::vector<double> stdvec;
    return arma::conv_to< stdvec >::from(solution);
}

void eig_qr(vector< vector<double> > a, vector< vector<double> >& v, vector<double>& l) {
    vec eigval;
    mat eigvec;
    eig_sym(eigval, eigvec, to_arma(a));
    l = conv_to< vector<double> >::from(eigval);
    v = to_std(eigvec);
}
