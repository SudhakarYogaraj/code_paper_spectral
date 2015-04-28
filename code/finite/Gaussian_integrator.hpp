#ifndef GAUSSIAN_INTEGRATOR_H
#define GAUSSIAN_INTEGRATOR_H

#define PI 3.141592653589793238462643383279502884

#include <functional>
#include <vector>
#include <cmath>

double quad1d(std::function<double (double)> f, double sigma);
double quadnd(std::function<double (std::vector<double>)> f, std::vector<double> sigmas);

#endif
