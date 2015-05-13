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

        Solver_spectral(int degree, int nNodes, int n_vars);
        void estimator(Problem&, std::vector<double> x, std::vector<SDE_coeffs>& , double t);

    private:
        std::vector< std::vector<double> > hermiteCoeffs;
        double basis(std::vector<int> mult, std::vector<double> x, std::vector<double> sigmas);

        /* 
         * Function that returns the vector of standard (normalized with weight of variance 1) hermite coefficients, based on the monomial expression
         */
        std::vector<double> mon2herm (std::vector<double> mcoeffs, int n, int d);
        int mult2ind(std::vector<int> m, int d);
        std::vector<int> ind2mult(int ind, int d, int n);
};
#endif
