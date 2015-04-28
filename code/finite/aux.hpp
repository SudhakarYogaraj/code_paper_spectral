#ifndef AUX_H
#define AUX_H

#include <iostream>
#include <iomanip>
#include <iterator>
#include <cmath>
#include <vector>
#include <sstream>
#include <fstream>

int bin(int n, int k);
double hermite(int n, double x, double sigma);
double hermiteM(std::vector<int> multIndex, std::vector<double> x, std::vector<double> sigmas);
int canonicalInd(std::vector<int> alpha, int n, int degree);
double normVec (std::vector<double> x);
double normMat (std::vector< std::vector<double> > x);
double delta(int a, int b);
std::vector< std::vector<double> > symmetric(std::vector< std::vector<double> > A);
std::vector< std::vector<double> > cholesky(std::vector< std::vector<double> > A);
#endif
