#ifndef AUX_H
#define AUX_H

#define PI 3.141592653589793238462643383279502884

#include <iostream>
#include <cmath>
#include <vector>

int bin(int n, int k);
std::vector< std::vector<double> > symmetric(std::vector< std::vector<double> > A);
std::vector< std::vector<double> > cholesky(std::vector< std::vector<double> > A);
double gaussian(std::vector<double> y, std::vector<double> sigmas);
double monomial(std::vector<int> mult, std::vector<double> x, std::vector<double> sigmas);
double monomial(std::vector<int> mult, std::vector<double> x);
std::vector<double> solve(std::vector< std::vector<double> > A, std::vector<double> b);
double ipow(double x, int e);

#endif
