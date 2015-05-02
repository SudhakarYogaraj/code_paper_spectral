#ifndef GAUSSIAN_INTEGRATOR_H
#define GAUSSIAN_INTEGRATOR_H

#define PI 3.141592653589793238462643383279502884

#include <functional>
#include <vector>
#include <cmath>
#include <iostream>

class Gaussian_integrator {
    public:
        double quadnd(std::function<double (std::vector<double>)> f, std::vector<double> sigmas);
        Gaussian_integrator(int nNodes, int nVars);

    private:
        std::vector< std::vector<double> > nodes;
        std::vector<double> weights;
};
#endif
