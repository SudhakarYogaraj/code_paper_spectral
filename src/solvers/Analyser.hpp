#ifndef ANALYSER_H

#define ANALYSER_H
#define PI 3.141592653589793238462643383279502884

#include <math.h>
#include <vector>
#include <iostream>

#include "problems/Problem.hpp"

/*! Class to track the properties of the invariant measure of the fast process
 * Provides function that describe the statistics of this measure.
 */

class Analyser {
    public:

        // Constructor
        Analyser(Problem *problem);

        // Copy constructor;
        /* Analyser(const Analyser& a); */

        // Value of the slow variable to which the statistics computed correspond.
        std::vector<double> x;

        // Bias of the invariant measure at x.
        std::vector<double> bias;

        // Covariance of the invariant measure at x.
        std::vector< std::vector<double> > covariance;

        // Eigenvalues of the covariance matrix
        std::vector<double> eig_val_cov;

        // Eigenvectors of the covariant matrix
        std::vector< std::vector<double> > eig_vec_cov;

        // Square root of the covariant matrix
        std::vector< std::vector<double> > sqrt_cov;

        // Inverse of the covariant matrix.
        std::vector< std::vector<double> > inv_cov;

        // Determinant of the covariant matrix
        double det_sqrt_cov;

        // Normalization of the density
        double normalization;

        // Normalized density of the problem
        double rho(std::vector<double> x, std::vector<double> y);

        // Update stats of the invariant density
        void update_stats(std::vector<double> x);

        // Rescaling
        std::vector<double> rescale(std::vector<double> y);

    private:

        // Problem associated with Analyser
        Problem *problem;

        // Variables derived from Problem, copied for convenience.
        int nf;
};
#endif
