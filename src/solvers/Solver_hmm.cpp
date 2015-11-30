#include "global/global.hpp"
#include "toolbox/linear_algebra.hpp"
#include "global/templates.hpp"
#include "problems/Problem.hpp"
#include "solvers/Solver_hmm.hpp"

using namespace std;

Solver_hmm::Solver_hmm(Problem *prob, config_hmm *config) {

    conf = config;
    problem = prob;
}

SDE_coeffs Solver_hmm::estimator(vec xt, double t) {

    // Vectors to store the coefficients of the sde
    SDE_coeffs sde_coeffs;

    // Initial value for the fast process at each macro time-step
    vec yInit(2*problem->nf, 0.);

    // Initialization of the coefficients of the SDE
    sde_coeffs.drif = vec(problem->ns, 0.);
    sde_coeffs.diff = mat(problem->ns, vec(problem->ns, 0.));

    default_random_engine generator;
    generator.seed(time(NULL));
    normal_distribution<double> distribution(0.0,1.0);
    int seed = (int) abs(1000*distribution(generator));
    generator.seed(seed);

    // Construction of the array that will contain the solution for the
    // fast process at each macro time-step.
    mat yAux(conf->nt + conf->n + conf->np, vec(2*problem->nf,0.));

    // Loop for ensemble average
    for (int m = 0; m < conf->M; m++) {

        // Initialization of the fast variables
        yAux[0] = yInit;

        vec drift(2*problem->nf);
        vec diffu(2*problem->nf);

        // Euler-Maruyama method for the fast processes:
        for (unsigned int j = 0; j < yAux.size() - 1; j++) {

            for (int k = 0; k < 2*problem->nf; ++k) {
                drift[k] = problem->drif[k](xt, yAux[j]);
                diffu[k] = problem->diff[k](xt, yAux[j]);
            }

            for (int k = 0; k < 2*problem->nf; k++)
            {
                yAux[j+1][k] = yAux[j][k] + drift[k]*conf->micro_dt +
                    diffu[k]*sqrt(conf->micro_dt)*distribution(generator);
            }
        }

        // Construction of auxiliary vector for efficiency.

        // sumsAux1 =approx= dax
        vector< mat > sumsAux1(conf->nt + conf->n, mat(problem->ns,vec(problem->ns,0.)));

        // sumsAux2 =approx= a
        mat sumsAux2(conf->nt + conf->n, vec(problem->ns, 0.));

        // First component of each vector
        for (int index = 0; index <= conf->np; index++) {
            for (int k = 0; k < problem->ns; ++k) {
                sumsAux2[0][k] = sumsAux2[0][k] + problem->a[k](xt, yAux[index]);
                for (int l = 0; l < problem->ns; ++l) {
                    sumsAux1[0][k][l] = sumsAux1[0][k][l] + problem->dxa[k][l](xt, yAux[index]);
                }
            }
        }

        // Recursion to obtain the other components
        for (int index = 1; index < conf->nt + conf->n; index++) {
            for (int k = 0; k < problem->ns; ++k) {
                sumsAux2[index][k] = sumsAux2[index-1][k] + problem->a[k](xt, yAux[index + conf->np]) - problem->a[k](xt, yAux[index-1]);
                for (int l = 0; l < problem->ns; ++l) {
                    sumsAux1[index][k][l] = sumsAux1[index-1][k][l] + problem->dxa[k][l](xt, yAux[index + conf->np]) - problem->dxa[k][l](xt, yAux[index-1]);
                }
            }
        }

        // Drift term by the HMM
        vec fim(problem->ns, 0.);
        vec fim1(problem->ns, 0.);
        vec fim2(problem->ns, 0.);

        for (int j = conf->nt; j < conf->nt + conf->n; j++) {

            // evaluation of data
            mat dya_j(problem->ns, vec (problem->nf));
            for (int k = 0; k < problem->ns; ++k) {
                for (int l = 0; l < problem->nf; ++l) {
                    dya_j[k][l] = problem->dya[k][l](xt, yAux[j]);
                }
            }

            // first term
            for (int i1 = 0; i1 < problem->ns; i1++) {
                for (int k = 0; k < problem->nf; k++)
                    fim1[i1] += dya_j[i1][k]*yAux[j][problem->nf+k];
            }

            // second term: improved
            for (int i1 = 0; i1 < problem->ns; i1++) {
                for (int i2 = 0; i2 < problem->ns; i2++) {
                    fim2[i1] += conf->micro_dt*problem->a[i2](xt, \
                            yAux[j])*sumsAux1[j][i1][i2];
                }
            }
        }

        for (int i1 = 0; i1 < problem->ns; i1++) {
            fim1[i1] = fim1[i1]/conf->n;
            fim2[i1] = fim2[i1]/conf->n;
            fim[i1] = fim1[i1] + fim2[i1];
        }

        // Diffusion term by the HMM
        mat him(problem->ns, vec(problem->ns,0.));

        for (int i1 = 0; i1 < problem->ns; i1++) {
            for (int i2 = 0; i2 < problem->ns; i2++) {
                for (int j = conf->nt; j < conf->nt + conf->n; j++) {
                    him[i1][i2] += conf->micro_dt*problem->a[i1](xt, \
                            yAux[j])*sumsAux2[j][i2];
                }
                him[i1][i2] = 2*him[i1][i2]/conf->n;
            }
        }

        // Incremental mean over the samples
        for (int i1 = 0; i1 < problem->ns; i1++) {
            for (int i2 = 0; i2 < problem->ns; i2++) {
                sde_coeffs.diff[i1][i2] = (sde_coeffs.diff[i1][i2]*m + him[i1][i2])/(m+1);
            }
            sde_coeffs.drif[i1] = (sde_coeffs.drif[i1]*m + fim[i1])/(m+1);
        }
    }

    // Approximate diffusion coefficient
    sde_coeffs.diff = cholesky(symmetric(sde_coeffs.diff));

    // Initial condition for next iteration and storage of y
    yInit = yAux[yAux.size()-1];

    // return coefficients of the SDE
    return sde_coeffs;
}

config_hmm Solver_hmm::sensible_conf(int p, int M) {

    config_hmm conf;

    // Precision parameter
    conf.p = p;

    // Number of replicas of the fast process
    conf.M = M;

    // Order of micro-solver (has to be 1)
    conf.l = 1;

    // Micro time step
    conf.micro_dt = pow(2.,-conf.p/conf.l);
    // conf.micro_dt = 0.05*pow(2.,-conf.p/conf.l); // KS

    // Number of micro time-steps taken into account
    // in the average to obtain the coefficients of the
    // effective equation at each macro time-step.
    conf.n = (int) 10*pow(2,conf.p*(2+1./conf.l));

    // Number of micro time-steps that are not taken
    // into account in the averages
    conf.nt = 16;

    // Number of micre time-steps used for the discretization
    // of the integrals in time.
    conf.np = (int) pow(2.,conf.p/conf.l)*conf.p;

    return conf;
}
