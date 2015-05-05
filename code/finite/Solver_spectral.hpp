#ifndef SOLVER_SPECTRAL_H
#define SOLVER_SPECTRAL_H

#include <numeric>
#include <functional>
#include <vector>

/* TODO: Delete this after debugging (urbain, Wed 29 Apr 2015 16:31:21 BST) */
#include "io.hpp"

class Solver_spectral {
    public:

        int n_mcmc;
        int degree;
        int nvars;
        int nNodes;
        double p;
        int nBasis;

		void set(double p, int n);
        void estimator(Problem&, std::vector<double> x, std::vector<double>& fi, std::vector< std::vector<double> >& hi, double t);
};
#endif
