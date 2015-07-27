#ifndef SOLVER_SPECTRAL_H
#define SOLVER_SPECTRAL_H

#include <numeric>
#include <string>
#include <functional>
#include <vector>
#include <unordered_map>

#include "Gaussian_integrator.hpp"
#include "Problem.hpp"
#include "structures.hpp"

class Solver_spectral {
    public:

        int n_mcmc;
        int degree;
        int nNodes;

        Solver_spectral(int degree, int nNodes, int n_vars);
        SDE_coeffs estimator(Problem&, std::vector<double> x, double t);

    private:

        struct Hash_multi_index {

            // "Hashing" function for multi_indices. Associate a unique integer to
            // any multi-index with |m|_inf < base.
            size_t operator() (const std::vector<int>& multi_index) const {

                // The base, which limits the validity of the hash function.
                int base = 100;

                int tmp = 0;
                for (unsigned int i = 0; i < multi_index.size(); ++i) {
                    tmp += multi_index[i] * pow(base, i);
                }

                return std::hash<int>()(tmp);
            }
        };

        // Cached multi-indices
        std::vector< std::vector<int> > ind2mult;
        std::unordered_map<std::vector<int>, int, Hash_multi_index> mult2ind;

        // Number of dimension to solve on
        int nf;

        // Matrix containing the coefficients of uni- and multi-
        // dimensional Hermite polynomials in terms of monomials.
        std::vector< std::vector<double> > hermiteCoeffs_1d;
        std::vector< std::vector<double> > hermiteCoeffs_nd;

        double gaussian_linear_term(std::vector<double> z);
        std::vector<double> map_to_real(std::vector<double> z);

        // Basis used for the method
        double monomial(std::vector<int> mult, std::vector<double> x);

        // Calculate coefficients of Hermite polynomials.
        void hermite_coefficients (int degree, std::vector< std::vector<double> >& matrix);

        // Update variance and bias of gaussian
        void update_stats(Problem& problem, std::vector<double> var_scaling);

        std::vector<double> discretize(Problem &problem, std::vector<double> x, Gaussian_integrator& gauss, double(*f)(std::vector<double>,std::vector<double>));
        std::vector<double> project(int nf, int degree, Gaussian_integrator& gauss, std::vector<double> f_discretized, int rescale);

        // Statistics associated with the hermite functions
        std::vector<double> bias;
        std::vector<double> eig_val_cov;
        std::vector< std::vector<double> > eig_vec_cov;
        std::vector< std::vector<double> > sqrt_cov;
        double det_cov;
};
#endif
