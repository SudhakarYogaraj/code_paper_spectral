#ifndef GAUSSIAN_INTEGRATOR_H
#define GAUSSIAN_INTEGRATOR_H

#define PI 3.141592653589793238462643383279502884

#include <functional>
#include <vector>
#include <cmath>
#include <iostream>

class Gaussian_integrator {
    public:
        Gaussian_integrator(int nNodes, int nVars);
        template<typename T, typename F> std::vector<T> quadnd(F f, const std::vector<double>& sigmas, const std::vector<T>& v0) {
            std::vector<T> result = v0;
            int nVars = sigmas.size();
            for (unsigned int i = 0; i < nodes.size(); ++i) {
                std::vector<double> x = nodes[i];
                for (int j = 0; j < nVars; ++j) {
                    x[j] = x[j] * (sqrt(2)*sigmas[j]);
                }
                result = result + (f(x) * weights[i]);
            }
            return result*(1./pow(sqrt(PI),nVars));
        }

        template<typename T, typename F> std::vector<T> flatquadnd(F f, const std::vector<double>& sigmas, const std::vector<T>& v0) {
            std::vector<T> result = v0;
            int nVars = sigmas.size();
            for (unsigned int i = 0; i < nodes.size(); ++i) {
                std::vector<double> x = nodes[i];
                for (int j = 0; j < nVars; ++j) {
                    x[j] = x[j] * (sqrt(2)*sigmas[j]);
                }
                double gaussian = 1.;
                for (unsigned int i = 0; i < sigmas.size(); ++i) {
                    double s = sigmas[i];
                    gaussian *= exp(-x[i]*x[i]/(2*s*s))/(sqrt(2*PI)*s);
                }
                result = result + (f(x) * weights[i] * (1/gaussian) );
            }
            return result*(1./pow(sqrt(PI),nVars));
        }

    private:
        std::vector< std::vector<double> > nodes;
        std::vector<double> weights;
};
#endif
