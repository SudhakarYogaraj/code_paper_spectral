#include "structures.hpp"
#include "toolbox.hpp"
#include "templates.hpp"
#include "Problem.hpp"
#include "Solver_hmm.hpp"

using namespace std;

Solver_hmm::Solver_hmm(double p, int M) {
    // Precision parameter
    this->p = p;

    // Order of the micro solver
    this->l = 1;

    // Macro time-step
    this->macro_dt = .01;

    // Micro time-step

    // If KS
    // this->micro_dt = 0.05*pow(2.,-this->p/this->l);

    // If Burgers
    this->micro_dt = pow(2.,-this->p/this->l);

    // Number of micro time-steps taken into account
    // in the average to obtain the coefficients of the
    // effective equation at each macro time-step.
    this->n = (int) 10*pow(2,this->p*(2+1./this->l));

    // Number of micro time-steps that are not taken
    // into account in the averages
    this->nt = 16;

    // Number of micre time-steps used for the discretization
    // of the integrals in time.
    this->np = (int) pow(2.,this->p/this->l)*this->p;

    // Number of replicas of the fast process
    this->M = M;
}

void Solver_hmm::estimator(Problem &problem, \
               vector<double> xt, \
               vector<double>& yInit, \
               SDE_coeffs& data, \
               int seed, \
               double t) {

    data.drif = vector<double>(problem.ns, 0.);
    data.diff = vector< vector<double> >(problem.ns, vector<double>(problem.ns, 0.));

    default_random_engine generator;
    normal_distribution<double> distribution(0.0,1.0);
    generator.seed(seed);

    // Construction of the array that will contain the solution for the
    // fast process at each macro time-step.
    vector< vector<double> > yAux(this->nt + this->n + this->np, \
            vector<double>(2*problem.nf,0.));

    // Loop for ensemble average
    for (int m = 0; m < this->M; m++) {

        // Initialization of the fast variables
        yAux[0] = yInit;

        vector<double> drift;
        vector<double> diffu;

        // Euler-Maruyama method for the fast processes:
        for (unsigned int j = 0; j < yAux.size() - 1; j++) {

            drift = problem.drif(xt, yAux[j]);
            diffu = problem.diff(xt, yAux[j]);

            for (int k = 0; k < 2*problem.nf; k++)
            {
                yAux[j+1][k] = yAux[j][k] + drift[k]*this->micro_dt +
                    diffu[k]*sqrt(this->micro_dt)*distribution(generator);
            }
        }

        // Construction of auxiliary vector for efficiency.

        // sumsAux1 =approx= dax
        vector< vector< vector<double> > > sumsAux1(this->nt + this->n, \
                vector< vector<double> >(problem.ns,vector<double>(problem.ns,0.)));

        // sumsAux2 =approx= a
        vector< vector<double> > sumsAux2(this->nt + this->n, \
                vector<double>(problem.ns, 0.));

        // First component of each vector
        for (int index = 0; index <= this->np; index++) {
            sumsAux2[0] = sumsAux2[0] + problem.a(xt, yAux[index]);
            sumsAux1[0] = sumsAux1[0] + problem.dax(xt, yAux[index]);
        }

        // Recursion to obtain the other components
        for (int index = 1; index < this->nt + this->n; index++) {
            sumsAux1[index] = sumsAux1[index-1] + problem.dax(xt, yAux[index + this->np]) - problem.dax(xt, yAux[index-1]);
            sumsAux2[index] = sumsAux2[index-1] + problem.a(xt, yAux[index + this->np]) - problem.a(xt, yAux[index-1]);
        }

        // Drift term by the HMM
        vector<double> fim(problem.ns, 0.);
        vector<double> fim1(problem.ns, 0.);
        vector<double> fim2(problem.ns, 0.);

        for (int j = this->nt; j < this->nt + this->n; j++) {

            // evaluation of data
            vector< vector<double> > dya_j = problem.day(xt, yAux[j]);

            // first term
            for (int i1 = 0; i1 < problem.ns; i1++) {
                for (int k = 0; k < problem.nf; k++)
                    fim1[i1] += dya_j[i1][k]*yAux[j][problem.nf+k];
            }

            // second term: improved
            for (int i1 = 0; i1 < problem.ns; i1++) {
                for (int i2 = 0; i2 < problem.ns; i2++) {
                    fim2[i1] += this->micro_dt*problem.a(xt, \
                            yAux[j])[i2]*sumsAux1[j][i1][i2];
                }
            }
        }

        for (int i1 = 0; i1 < problem.ns; i1++) {
            fim1[i1] = fim1[i1]/this->n;
            fim2[i1] = fim2[i1]/this->n;
            fim[i1] = fim1[i1] + fim2[i1];
        }

        // Diffusion term by the HMM
        vector< vector<double> > him(problem.ns, vector<double>(problem.ns,0.));

        for (int i1 = 0; i1 < problem.ns; i1++) {
            for (int i2 = 0; i2 < problem.ns; i2++) {
                for (int j = this->nt; j < this->nt + this->n; j++) {
                    him[i1][i2] += this->micro_dt*problem.a(xt, \
                            yAux[j])[i1]*sumsAux2[j][i2];
                }
                him[i1][i2] = 2*him[i1][i2]/this->n;
            }
        }

        // Incremental mean over the samples
        for (int i1 = 0; i1 < problem.ns; i1++) {
            for (int i2 = 0; i2 < problem.ns; i2++) {
                data.diff[i1][i2] = (data.diff[i1][i2]*m + him[i1][i2])/(m+1);
            }
            data.drif[i1] = (data.drif[i1]*m + fim[i1])/(m+1);
        }
    }

    // Approximate diffusion coefficient
    data.diff = cholesky(symmetric(data.diff));

    // Initial condition for next iteration and storage of y
    yInit = yAux[yAux.size()-1];
}