#ifndef SOLVER_SPECTRAL_H
#define SOLVER_SPECTRAL_H

#include <numeric>
#include <functional>
#include <vector>

class Solver_spectral {
    public:

        int n_mcmc;
        int degree;
        int nNodes;

        Solver_spectral(int degree, int nNodes);
        void estimator(Problem&, std::vector<double> x, SDE_coeffs& , double t);
};
#endif
