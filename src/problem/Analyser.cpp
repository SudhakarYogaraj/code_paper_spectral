#include "Analyser.hpp"
#include "Problem.hpp"

#include "templates.hpp"
#include "Gaussian_integrator.hpp"
#include "linear_algebra.hpp"
#include "io/io.hpp"
#include "global.hpp"

using namespace std;

// Consntructor of the analyser. The only parameter is the problem that
// the analyser will follow.
Analyser::Analyser(Problem *p) {

    // Assign variables derived from problem
    problem = p;
    nf = problem->nf;

    // Initialize other variables to reasonable values
    normalization = 1.;
    det_sqrt_cov = 1.;

    bias = vector<double>(problem->nf, 0.);
    eig_val_cov = vector<double>(problem->nf, 1.);

    covariance = vector< vector<double> > (nf, vector<double> (nf,0.));
    inv_cov = vector< vector<double> > (nf, vector<double> (nf,0.));
    sqrt_cov = vector< vector<double> > (nf, vector<double> (nf,0.));
    eig_vec_cov = vector< vector<double> > (nf, vector<double> (nf,0.));

    for (int i = 0; i < problem->nf; ++i) {
        covariance[i][i] = 1.;
        inv_cov[i][i] = 1.;
        sqrt_cov[i][i] = 1.;
        eig_vec_cov[i][i] = 1.;
    }
}

/* Analyser::Analyser(const Analyser& a) { */
/*     this->problem = a.problem; */
/*     this->nf = a.nf; */
/* } */

// Implementation of a rescaling used to approximate the invariant measure.
vector<double> Analyser::rescale(vector<double> y) {
    vector<double> result(y.size(), 0.);
    for (int i = 0; i < y.size(); ++i) {
        for (int j = 0; j < y.size(); ++j) {
            result[i] += sqrt_cov[i][j] * y[j];
        }
    }
    return (result + bias);
}

// Update statistics of the invariant measure at x.
// The outer loop serves to obtain a more accurate result.
void Analyser::update_stats(vector<double> x) {

    Gaussian_integrator gauss_plus = Gaussian_integrator(100, nf);
    Gaussian_integrator gauss = Gaussian_integrator(30, nf);

    int n_iterations = 5, i = 0;
    for (int n = 0; n < n_iterations; ++n) {

        // Normalization constant
        auto lambda = [&] (vector<double> z) -> double {
            vector<double> y = rescale(z);
            return det_sqrt_cov * problem->zrho(x,y)/gaussian(z);
        };
        normalization = gauss_plus.quadnd(lambda);

        // Calculation of the bias of 'rho'
        bias = vector<double> (nf, 0.);
        for (i = 0; i < nf; ++i) {
            auto lambda = [&] (vector<double> z) -> double {
                vector<double> y = rescale(z);
                return det_sqrt_cov * y[i] * (rho(x,y)/gaussian(z));
            };
            bias[i] = gauss.quadnd(lambda);
        }

        // Calculation of the covariance matrix
        for (i = 0; i < nf; ++i) {
            for (int j = 0; j < nf; ++j) {
                auto lambda = [&] (vector<double> z) -> double {
                    vector<double> y = rescale(z);
                    return det_sqrt_cov * (y[i] - bias[i]) * (y[j] - bias[j]) * (rho(x,y)/gaussian(z));
                };
                covariance[i][j] = gauss.quadnd(lambda);
            }
        }

        // Eigenvalue decomposition
        eig_qr(covariance, eig_vec_cov, eig_val_cov);

        sqrt_cov = eig_vec_cov;
        for (int i = 0; i < nf; i++) {
            for (int j = 0; j < nf; j++) {
                sqrt_cov[i][j] *= sqrt(eig_val_cov[j]);
            }
        }

        // Determinant of sqrt_cov
        for (det_sqrt_cov = 1., i = 0; i < nf; ++i)
            det_sqrt_cov *= sqrt(eig_val_cov[i]);

        // Inverse of covariance matrix
        for (int i = 0; i < nf; ++i) {
            vector<double> rhs_tmp(nf, 0.); rhs_tmp[i] = 1.;
            inv_cov[i] = solve(covariance, rhs_tmp);
        }
        inv_cov = transpose(inv_cov);
    }

    if(DEBUG) {
        cout << endl << "* Covariance matrix of the invariant density" << endl;
        niceMat(covariance);

        cout << endl << "* Inverse of the covariance matrix" << endl;
        niceMat(inv_cov);

        cout << endl << "* Bias of the invariant density" << endl;
        niceVec(bias);

        cout << endl << "* Eigenvectors of the covariance matrix" << endl;
        niceMat(eig_vec_cov);
    }
}

double Analyser::rho(vector<double> x, vector<double> y) {
    return problem->zrho(x,y)/normalization;
}
