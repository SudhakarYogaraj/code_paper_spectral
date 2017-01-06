#include "problems/Problem.hpp"
#include "toolbox/Gaussian_integrator.hpp"
#include "solvers/Solver_exact.hpp"

using namespace std;

Solver_exact::Solver_exact(Problem* p, Analyser* a) {
    problem = p;
    analyser = a;
}

SDE_coeffs Solver_exact::estimator(vec x, double t) {
    analyser->update_stats(x);
    SDE_coeffs sde_coeffs;
    sde_coeffs.drif = soldrif(x);
    sde_coeffs.diff = soldiff(x);
    return sde_coeffs;
}

vec Solver_exact::soldrif(vec x) {

    int ns = problem->ns;
    int nf = problem->nf;

    vec result(ns,0.);
    Gaussian_integrator gauss = Gaussian_integrator(100,nf);

    // Function to integrate to obtain exact drift.
    auto lambda = [&] (vec z) -> vec {
        vec y = analyser->rescale(z);
        vec tmp(ns, 0.);
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

mat Solver_exact::soldiff(vec x) {

    int ns = problem->ns;
    int nf = problem->nf;

    Gaussian_integrator gauss = Gaussian_integrator(100,nf);
    mat result(ns,vec(ns,0.));

    // Function to integrate to obtain exact diffusion coefficient.
    auto lambda = [&] (vec z) -> mat {
        vec y = analyser->rescale(z);
        mat tens_prod(ns, vec(ns, 0.));
        for (int i = 0; i < ns; ++i) {
            for (int j = 0; j < ns; ++j) {
                tens_prod[i][j] = 2*problem->a[i](x,y)*problem->phi[j](x,y)*(analyser->rho(x,y)/gaussian(z));
            }
        }
        return tens_prod;
    };
    /* mat tmp = symmetric( gauss.quadnd(lambda, result) * analyser->det_sqrt_cov ); */
    /* cout << tmp[0][0] << "," << tmp[0][1] << endl << tmp[1][0] << "," << tmp[1][1]; */
    /* exit(0); */
    result = square_root(symmetric( gauss.quadnd(lambda, result) * analyser->det_sqrt_cov ));
    return result;
}
