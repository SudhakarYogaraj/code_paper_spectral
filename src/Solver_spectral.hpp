#ifndef SOLVER_SPECTRAL_H
#define SOLVER_SPECTRAL_H

#include <numeric>
#include <string>
#include <functional>
#include <vector>

class Solver_spectral {
    public:

        int n_mcmc;
        int degree;
        int nNodes;

        Solver_spectral(int degree, int nNodes, int n_vars);
        SDE_coeffs estimator(Problem&, std::vector<double> x, double t);

    private:
        // Matrix containing the coefficients of uni- and multi-
        // dimensional Hermite polynomials in terms of monomials.
        std::vector< std::vector<double> > hermiteCoeffs_1d;
        std::vector< std::vector<double> > hermiteCoeffs_nd;
        std::vector< std::vector<double> > herm_to_basis;

        // Basis used for the method
        double monomial(std::vector<int> mult, std::vector<double> x);

        // Function to obtain the coefficients of the expansion of a
        // polynomial in terms of hermite polynomials, based on its
        // exansion in monomials.
        std::vector<double> basis2herm (std::vector<double> mcoeffs, int n, int ns);

        // Auxiliary function to convert multi-indices to linear indices.
        int mult2ind(std::vector<int> m, int ns);

        // Auxiliary function to convert linear indices to multi-indices.
        std::vector<int> ind2mult(int ind, int ns, int n);

        // basis of polynomials to use
        std::string poly_basis;

        // Calculate coefficients of Hermite polynomials.
        void hermite_coefficients (int degree, std::vector< std::vector<double> >& matrix);
};
#endif
