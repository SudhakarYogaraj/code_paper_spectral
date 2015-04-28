#ifndef AUX_H
#define AUX_H

#define PI 3.141592653589793238462643383279502884

#include <iostream>
#include <iomanip>
#include <functional>
#include <cmath>
#include <vector>
#include <sstream>
#include <fstream>

int bin(int n, int k);
double hermite(int n, double x, double sigma);
double hermiteM(std::vector<int> multIndex, std::vector<double> x, std::vector<double> sigmas);
int canonicalInd(std::vector<int> alpha, int n, int degree);
double gauss_hermite_1D(std::function<double (double)> f, double sigma);
double gauss_hermite_nD(std::function<double (std::vector<double>)> f, std::vector<double> sigmas);
void print2Mats (std::vector< std::vector<double> > x, std::vector< std::vector<double> > y);
void printMat (std::vector< std::vector<double> > x);
double normVec (std::vector<double> x);
double normMat (std::vector< std::vector<double> > x);
void writeMatToFile(std::string s, std::vector< std::vector<double> > x);
void print2Vecs (std::vector<double> x, std::vector<double> y);
void printVec (std::vector<double> x);
double delta(int a, int b);
std::vector< std::vector<double> > symmetric(std::vector< std::vector<double> > A);
std::vector< std::vector<double> > cholesky(std::vector< std::vector<double> > A);
void writeToFile(std::string s, std::vector<double> x);
#endif
