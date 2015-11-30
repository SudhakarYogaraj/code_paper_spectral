#ifndef SOLVER_SPECTRAL_H
#define SOLVER_SPECTRAL_H

#include <numeric>
#include <string>
#include <functional>
#include <vector>
#include <unordered_map>
#include <armadillo>

#include "global/global.hpp"
#include "solvers/Solver.hpp"
#include "toolbox/Gaussian_integrator.hpp"
#include "problems/Problem.hpp"
#include "solvers/Analyser.hpp"

struct config_spectral {
    int n_nodes;
    int degree;
    std::vec scaling;
};

class Solver_spectral : public Solver {
    public:

        Solver_spectral(Problem*, Analyser*, config_spectral*);
        SDE_coeffs estimator(std::vec x, double t);
        std::vector<SDE_coeffs> full_estimator(std::vec x, double t, std::vector<int> degrees);

    private:

        config_spectral* conf;

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
        int ns;

        // Matrix containing the coefficients of uni- and multi- dimensional Hermite polynomials in terms of monomials.
        std::mat hermiteCoeffs_1d;
        std::mat hermiteCoeffs_nd;

        double gaussian_linear_term(std::vec z);
        std::vec map_to_real(std::vec z);

        // Calculate coefficients of Hermite polynomials.
        void hermite_coefficients (int degree, std::mat& matrix);

        // Update variance and bias of gaussian
        void update_stats();

        // Compute matrix of the linear system
        std::mat compute_matrix(std::vec x);

        // Compute effective coefficients from expansions in Hermite functions
        SDE_coeffs compute_averages(const std::mat& functions, const std::mat& solutions);

        // Discretize function on grid
        std::vec discretize(std::vec x, double(*f)(std::vec,std::vec));

        // Project discretized functions on monomials and hermite polynomials
        std::vec project_mon(int nf, int degree, std::vec f_discretized, int rescale);
        std::vec project_herm(int nf, int degree, std::vec f_discretized, int rescale);

        // Statistics associated with the hermite functions
        std::vec bias;
        std::vec eig_val_cov;
        std::mat eig_vec_cov;
        std::mat sqrt_cov;
        double det_cov;

        // Integrator
        Gaussian_integrator *gauss;
        Problem *problem;
        Analyser *analyser;
};
#endif
