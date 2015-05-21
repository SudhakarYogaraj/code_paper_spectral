#ifndef SOLVER_SPECTRAL_H
#define SOLVER_SPECTRAL_H

#include <numeric>
#include <string>
#include <functional>
#include <vector>
#include <assert.h>

class Solver_spectral {
    public:

        int n_mcmc;
        int degree;
        int nNodes;

        Solver_spectral(int degree, int nNodes, int n_vars, std::string poly_basis);
        void estimator(Problem&, std::vector<double> x, std::vector<SDE_coeffs>& , double t);

    private:
        // Matrix containing the coefficients of uni- and multi-
        // dimensional Hermite polynomials in terms of monomials.
        std::vector< std::vector<double> > hermiteCoeffs_1d;
        std::vector< std::vector<double> > hermiteCoeffs_nd;
        std::vector< std::vector<double> > herm_to_basis;

        // Basis used for the method
        double basis(std::vector<int> mult, std::vector<double> x, std::vector<double> sigmas);

        // Function to obtain the coefficients of the expansion of a
        // polynomial in terms of hermite polynomials, based on its
        // exansion in monomials.
        std::vector<double> basis2herm (std::vector<double> mcoeffs, int n, int d);

        // Auxiliary function to convert multi-indices to linear indices.
        int mult2ind(std::vector<int> m, int d);

        // Auxiliary function to convert linear indices to multi-indices.
        std::vector<int> ind2mult(int ind, int d, int n);

        // basis of polynomials to use
        std::string poly_basis;
};
#endif
