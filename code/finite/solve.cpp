#include "header.h"

// Function implementing the numerical method
void solve(Problem &problem, \
        Solver &solver, \
        vector<double> xt, \
        vector<double>& yInit, \
        vector<double>& fi, \
        vector< vector<double> >& hi, \
        int seed, \
        double t) { 

    default_random_engine generator; 
    normal_distribution<double> distribution(0.0,1.0);
    generator.seed(seed);

    // Construction of the array that will contain the solution for the 
    // fast process at each macro time-step. 
    vector< vector<double> > yAux(solver.nt + solver.n + solver.np, \
            vector<double>(2*problem.nf,0.)); 

    // Loop for ensemble average
    for (int m = 0; m < solver.M; m++) {

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
                yAux[j+1][k] = yAux[j][k] + drift[k]*solver.micro_dt + 
                    diffu[k]*sqrt(solver.micro_dt)*distribution(generator);
            }
        }

        // Construction of auxiliary vector for efficiency.

        // sumsAux1 =approx= dax 
        vector< vector< vector<double> > > sumsAux1(solver.nt + solver.n, \
                vector< vector<double> >(problem.d,vector<double>(problem.d,0.)));

        // sumsAux2 =approx= a
        vector< vector<double> > sumsAux2(solver.nt + solver.n, \
                vector<double>(problem.d, 0.));

        // First component of each vector
        for (int index = 0; index <= solver.np; index++) {
            for (int i1 = 0; i1 < problem.d; i1++) {
                for (int i2 = 0; i2 < problem.d; i2++) {
                    sumsAux1[0][i1][i2] += problem.dax(xt, yAux[index])[i1][i2];
                }
                sumsAux2[0][i1] += problem.a(xt, yAux[index])[i1];
            }
        }

        // Recursion to obtain the other components
        for (int index = 1; index < solver.nt + solver.n; index++) {
            for (int i1 = 0; i1 < problem.d; i1++) {
                for (int i2 = 0; i2 < problem.d; i2++) {
                    sumsAux1[index][i1][i2] = sumsAux1[index-1][i1][i2] + problem.dax(xt, \
                            yAux[index + solver.np])[i1][i2] - problem.dax(xt, yAux[index-1])[i1][i2];
                }
                sumsAux2[index][i1] = sumsAux2[index-1][i1] + problem.a(xt, \
                        yAux[index + solver.np])[i1] - problem.a(xt, yAux[index-1])[i1];
            }
        }

        // Drift term by the HMM
        vector<double> fim(problem.d, 0.);
        vector<double> fim1(problem.d, 0.);
        vector<double> fim2(problem.d, 0.);

        for (int j = solver.nt; j < solver.nt + solver.n; j++) {

            // evaluation of data
            vector< vector<double> > dya_j = problem.day(xt, yAux[j]);             

            // first term
            for (int i1 = 0; i1 < problem.d; i1++) {
                for (int k = 0; k < problem.nf; k++) 
                    fim1[i1] += dya_j[i1][k]*yAux[j][problem.nf+k];
            }

            // second term: improved 
            for (int i1 = 0; i1 < problem.d; i1++) {
                for (int i2 = 0; i2 < problem.d; i2++) {
                    fim2[i1] += solver.micro_dt*problem.a(xt, \
                            yAux[j])[i2]*sumsAux1[j][i1][i2];
                }
            }
        }

        for (int i1 = 0; i1 < problem.d; i1++) {
            fim[i1] = fim1[i1] + fim2[i1];
            fim[i1] = fim[i1]/solver.n;
        }

        // Diffusion term by the HMM
        vector< vector<double> > him(problem.d, vector<double>(problem.d,0.));

        for (int i1 = 0; i1 < problem.d; i1++) {
            for (int i2 = 0; i2 < problem.d; i2++) { 
                for (int j = solver.nt; j < solver.nt + solver.n; j++) {
                    him[i1][i2] += solver.micro_dt*problem.a(xt, \
                            yAux[j])[i1]*sumsAux2[j][i2];
                }
                him[i1][i2] = 2*him[i1][i2]/solver.n;
            }
        }

        // Incremental mean over the samples
        for (int i1 = 0; i1 < problem.d; i1++) {
            for (int i2 = 0; i2 < problem.d; i2++) {
                hi[i1][i2] = (hi[i1][i2]*m + him[i1][i2])/(m+1);
            }
            fi[i1] = (fi[i1]*m + fim[i1])/(m+1);
        }
    } 

    // Approximate diffusion coefficient
    hi = cholesky(hi);

    // Initial condition for next iteration and storage of y
    yInit = yAux[yAux.size()-1]; 
}
