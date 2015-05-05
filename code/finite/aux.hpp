#ifndef AUX_H
#define AUX_H

#include <iostream>
#include <cmath>
#include <vector>

int bin(int n, int k);
int canonicalInd(std::vector<int> alpha, int n, int degree);
double delta(int a, int b);
std::vector< std::vector<double> > symmetric(std::vector< std::vector<double> > A);
std::vector< std::vector<double> > cholesky(std::vector< std::vector<double> > A);


double monomial(std::vector<int> mult, std::vector<double> x, std::vector<double> sigmas);
std::vector<double> hcoeffs (std::vector<double> mcoeffs, int n, int d);
int mult2ind(std::vector<int> m, int d);
std::vector<int> ind2mult(int ind, int d, int n);

#endif
