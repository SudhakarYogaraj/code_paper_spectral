#ifndef COMBINATORIAL_H
#define COMBINATORIAL_H

#include <vector>
#include <iostream>
#include <algorithm>

#include "linear_algebra.hpp"

// Functions to list multi-indices
std::vector< std::vector<int> > equal_multi_indices(int n, int d);
std::vector< std::vector<int> > lower_multi_indices(int n, int d);
std::vector< std::vector<int> > interval_multi_indices(int n, int a, int b);

// Binomial coefficients
int bin(int n, int k);

// Standard gaussian
double gaussian(std::vector<double> y);

#endif
