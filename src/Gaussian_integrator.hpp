#ifndef GAUSSIAN_INTEGRATOR_H
#define GAUSSIAN_INTEGRATOR_H

#define PI 3.141592653589793238462643383279502884

#include <functional>
#include <vector>
#include <cmath>
#include <iostream>
#include "templates.hpp"

class Gaussian_integrator {
    public:
        Gaussian_integrator(int nNodes, int nVars);

        double quadnd(std::function<double(std::vector<double>)> f);

        template<typename T, typename F> std::vector<T> quadnd(F f, const std::vector<T>& v0) {
            std::vector<T> result = v0;
            for (unsigned int i = 0; i < nodes.size(); ++i) {
                std::vector<double> x = nodes[i];
                for (int j = 0; j < nVars; ++j) {
                    x[j] = x[j] * sqrt(2);
                }
                result = result + (f(x) * weights[i]);
            }
            return result*(1./pow(sqrt(PI),nVars));
        }

        std::vector< std::vector<double> > nodes;
        std::vector<double> weights;

    private:
        int nVars;
};
#endif
