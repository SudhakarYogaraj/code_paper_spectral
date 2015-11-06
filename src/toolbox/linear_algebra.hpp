#ifndef LINEAR_ALGEBRA_H

#define LINEAR_ALGEBRA_H
#define PI 3.141592653589793238462643383279502884

#include <iostream>
#include <cmath>

#include "global.hpp"

// Symmetric part of a matrix.
std::mat symmetric(std::mat A);

// Cholesky decomposition of a matrix.
std::mat cholesky(std::mat A);
std::mat transpose(std::mat A);

// Solution of linear system for SPD matrices.
std::vector<double> solve(std::mat, std::vec);

// Integer power.
double ipow(double x, int e);

arma::mat to_arma(const std::mat &A);
std::mat to_std(const arma::mat &A);

arma::vec to_arma_vec(const std::vec &v);
std::vec to_std_vec(const arma::vec &v);

#endif
