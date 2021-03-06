#ifndef GAUSSIAN_INTEGRATOR_H
#define GAUSSIAN_INTEGRATOR_H

#define PI 3.141592653589793238462643383279502884

#include <functional>
#include <vector>
#include <cmath>
#include <iostream>
#include "global/templates.hpp"
#include "toolbox/combinatorics.hpp"

class Gaussian_integrator {


    public:
        Gaussian_integrator(int nNodes, int nVars);

        double quadnd(std::function<double(std::vector<double>)> f);

        template<typename T, typename F> std::vector<T> quadnd(F f, const std::vector<T>& v0) {
            std::vector<T> result = v0;
            for (unsigned int i = 0; i < nodes.size(); ++i) {
                result = result + f(nodes[i]) * weights[i];
            }
            return result;
        }

        std::vector< std::vector<double> > nodes;
        std::vector<double> weights;

        static void test_integrator();

    private:

        // Number of dimensions
        int nVars;

        // Product of quadrature rules
        void quad_prod(std::vector<int> sizes, std::vector< std::vector<double> >& nodes, std::vector<double>& weights);


        // Smolyak algorithm
        void Smolyak(std::vector< std::vector<double> >& nodes, std::vector<double>& weights);
};
#endif
