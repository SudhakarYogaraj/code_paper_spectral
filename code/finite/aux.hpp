#ifndef AUX_H
#define AUX_H

#include <iostream>
#include <cmath>
#include <vector>

int bin(int n, int k);
double hermite(int n, double x, double sigma);
double hermiteM(std::vector<int> multIndex, std::vector<double> x, std::vector<double> sigmas);
int canonicalInd(std::vector<int> alpha, int n, int degree);
double delta(int a, int b);
std::vector< std::vector<double> > symmetric(std::vector< std::vector<double> > A);
std::vector< std::vector<double> > cholesky(std::vector< std::vector<double> > A);
std::vector<double> operator-(const std::vector<double>& v1, const std::vector<double>& v2);
std::vector< std::vector<double> > operator-(const std::vector< std::vector<double> >& v1, const std::vector< std::vector<double> >& v2);

#endif
