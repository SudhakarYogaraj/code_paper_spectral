#ifndef AUX_H
#define AUX_H

#include <iostream>
#include <cmath>
#include <vector>

int bin(int n, int k);
std::vector< std::vector<double> > symmetric(std::vector< std::vector<double> > A);
std::vector< std::vector<double> > cholesky(std::vector< std::vector<double> > A);


double monomial(std::vector<int> mult, std::vector<double> x, std::vector<double> sigmas);
std::vector<double> mon2herm (std::vector<double> mcoeffs, int n, int d);
int mult2ind(std::vector<int> m, int d);
std::vector<int> ind2mult(int ind, int d, int n);

#endif
