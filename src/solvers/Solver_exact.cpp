#include "problem/Problem.hpp"
#include "toolbox/Gaussian_integrator.hpp"
#include "solvers/Solver_exact.hpp"

using namespace std;

Solver_exact::Solver_exact(Problem* p, Analyser* a) {
    problem = p;
    analyser = a;
}

SDE_coeffs Solver_exact::estimator(vector<double> x, double t) {
    analyser->update_stats(x);
    SDE_coeffs sde_coeffs;
    sde_coeffs.drif = soldrif(x);
    sde_coeffs.diff = soldiff(x);
    return sde_coeffs;
}

vector<double> Solver_exact::soldrif(vector<double> x) {

    int ns = problem->ns;
    int nf = problem->nf;

    vector<double> result(ns,0.);
    Gaussian_integrator gauss = Gaussian_integrator(100,nf);

    // Function to integrate to obtain exact drift.
    auto lambda = [&] (vector<double> z) -> vector<double> {
        vector<double> y = analyser->rescale(z);
        vector<double> tmp(ns, 0.);
        for (int i = 0; i < ns; ++i) {
            tmp[i] += problem->phi[i](x,y) * problem->stardiv_h(x,y);
            for (int j = 0; j < ns; ++j) {
                tmp[i] += problem->dxphi[i][j](x,y) * problem->a[j](x,y);
            }
        }
        return tmp*(analyser->rho(x,y)/gaussian(z));
    };

    result = gauss.quadnd(lambda, result) * analyser->det_sqrt_cov;
    return result;
}

vector< vector<double> > Solver_exact::soldiff(vector<double> x) {

    int ns = problem->ns;
    int nf = problem->nf;

    Gaussian_integrator gauss = Gaussian_integrator(100,nf);
    vector< vector<double> > result(ns,vector<double>(ns,0.));

    // Function to integrate to obtain exact diffusion coefficient.
    auto lambda = [&] (vector<double> z) -> vector< vector<double> > {
        vector<double> y = analyser->rescale(z);
        vector< vector<double> > tens_prod(ns, vector<double>(ns, 0.));
        for (int i = 0; i < ns; ++i) {
            for (int j = 0; j < ns; ++j) {
                tens_prod[i][j] = 2*problem->a[i](x,y)*problem->phi[j](x,y)*(analyser->rho(x,y)/gaussian(z));
            }
        }
        return tens_prod;
    };
    result = cholesky(symmetric( gauss.quadnd(lambda, result) * analyser->det_sqrt_cov ));
    return result;
}
