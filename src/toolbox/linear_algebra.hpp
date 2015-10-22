#ifndef LINEAR_ALGEBRA_H

#define LINEAR_ALGEBRA_H
#define PI 3.141592653589793238462643383279502884

#include <iostream>
#include <cmath>
#include <vector>

// Symmetric part of a matrix.
std::vector< std::vector<double> > symmetric(std::vector< std::vector<double> > A);

// Cholesky decomposition of a matrix.
std::vector< std::vector<double> > cholesky(std::vector< std::vector<double> > A);
void lu(std::vector< std::vector<double> > a, std::vector< std::vector<double> >& l, std::vector< std::vector<double> >& u);
std::vector< std::vector<double> > transpose(std::vector< std::vector<double> > A);

// Solution of linear system for SPD matrices.
std::vector<double> solve(std::vector< std::vector<double> > A, std::vector<double> b);

// Integer power.
double ipow(double x, int e);

// Eigenvalue decomposition
void qr(std::vector< std::vector<double> > a, std::vector< std::vector<double> >& q, std::vector< std::vector<double> >& r);
void eig_qr(std::vector< std::vector<double> > a, std::vector< std::vector<double> >& v, std::vector<double>& l);

#endif