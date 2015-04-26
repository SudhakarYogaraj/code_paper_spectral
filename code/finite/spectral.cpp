#include "header.h"

// Main method
void solve_spectral(Problem &problem, \
                    Solver_spectral &solver, \
                    vector<double> xt, \
                    vector<double>& fi, \
                    vector< vector<double> >& hi, \
                    int seed, \
                    double t) {

    int nBasis = solver.nBasis;
    int nf     = problem.nf;
    int degree = solver.degree;
    int ns     = problem.d;

    for (int i = 0; i < nBasis; ++i) {
        vector<int> multIndex = solver.ind2mult(i);
        /* cout << "i = " << i << endl; */
        /* cout << "multIndex: " << multIndex[0] << " " << multIndex[1] << endl; */
        /* cout << "mult2ind: " << solver.mult2ind(multIndex) << endl; */
    }

    // Eigenvalues
    vector<double> sigmas(nf, 0.);
    for (int i = 0; i < nf; ++i) {
        double aux = problem.betas[i];
        sigmas[i] = sqrt(0.5*aux*aux/problem.lambdas[i]);
    }

    // Parameters for random numbers
    default_random_engine generator; generator.seed(seed);
    normal_distribution<double> distribution(0.0,1.0);

    // Expansion of right-hand side of the Poisson equation
    vector< vector<double> > coefficients(ns, vector<double>(nBasis, 0.));
    vector< vector <vector<double> > > coefficients_dx(ns, vector< vector<double> >(ns, vector<double>(nBasis, 0.)));
    vector< vector<double> > coefficients_h(nf, vector<double>(nBasis,0.));

    for (int j = 0; j < nBasis; ++j) {
        /* cout << "j = " << j << endl; */
        vector<int> multIndex = solver.ind2mult(j);
        /* cout << "multIndex" << multIndex[0] << " " << multIndex[1] << endl; */
        double sum = 0.;

        // Monte-Carlo to compute the coefficients
        for (int k = 0; k < solver.n_mcmc; ++k) {
            vector<double> randn(nf, 0.);

            for (int l = 0; l < nf; ++l)
                randn[l] = sigmas[l]*distribution(generator);
            double h_eval = hermiteM(multIndex,randn,sigmas);
            /* cout << "h_eval" << h_eval << endl; */

            vector<double> slow_drift = problem.a(xt,randn);
            vector< vector<double> > slow_drift_dx = problem.dax(xt,randn);
            for (int l = 0; l < ns; ++l) {
                coefficients[l][j] += h_eval*slow_drift[l];
                for (int m = 0; m < ns; ++m) {
                    coefficients_dx[l][m][j] += h_eval*slow_drift_dx[l][m];
                }
            }

            vector<double> fast_drift_aux(nf, 0.);
            fast_drift_aux = problem.fast_drift_h(xt,randn);
            for (int l = 0; l < nf; ++l) {
                coefficients_h[l][j] += h_eval*fast_drift_aux[l];
            }
        }
        for (int l = 0; l < ns; ++l) {
            coefficients[l][j] /= solver.n_mcmc;
            for (int m = 0; m < ns; ++m)
                coefficients_dx[l][m][j] /= solver.n_mcmc;
        }
        for (int l = 0; l < nf; ++l) {
            coefficients_h[l][j] /= solver.n_mcmc;
        }
        /* cout << "coefficient_h " << j << ": " << coefficients_h[0][j] << endl; */
        /* cout << "coefficient_h " << j << ": " << coefficients_h[1][j] << endl << endl; */
    }

    // Solution of the Poisson equation
    vector< vector<double> > solution(ns, vector<double>(nBasis,0.));
    vector< vector < vector<double> > > solution_dx(ns, vector< vector<double> >(ns, vector<double>(nBasis, 0.)));
    vector< vector< vector<double> > > solution_dy(ns, vector< vector <double> >(nf, vector<double>(nBasis,0.)));
    for (int j = 1; j < nBasis; ++j) {
        /* cout << "multIndex: " << solver.ind2mult(j)[0] << " " << solver.ind2mult(j)[1] << endl; */
        double eig = 0.;
        for (int k = 0; k < nf; ++k) {
            eig += solver.ind2mult(j)[k]*problem.lambdas[k];
        }
        for (int l = 0; l < ns; ++l) {
            solution[l][j] = coefficients[l][j]/eig;
            for (int m = 0; m < ns; ++m)
                solution_dx[l][m][j] = coefficients_dx[l][m][j]/eig;
        }
    }

    // y-Derivatives of the solution
    for (int j = 0; j < nBasis; ++j) {
        vector<int> thisMult = solver.ind2mult(j);
        int sum = 0;
        for (int l = 0; l < nf; ++l) {
            sum += thisMult[l];
        }
        /* cout << "ThisMult: " <<  thisMult[0] << " " <<  thisMult[1] << endl; */
        if (sum < degree) {
            for (int l = 0; l < nf; ++l) {
                vector<int> newMult(nf, 0);
                for (int m = 0; m < nf; ++m) {
                    newMult[m] = thisMult[m];
                }
                newMult[l] ++;
                int newInd = solver.mult2ind(newMult);
                /* cout << "New mult: " <<  newMult[0] << " " <<  newMult[1] << endl; */
                /* cout << "New ind: " <<  newInd << "/" << nBasis << endl; */
                for (int m = 0; m < ns; ++m) {
                    solution_dy[m][l][j] = solution[m][newInd]*sqrt(newMult[l])/sigmas[l];
                }
            }
        }

        /* cout << "solution_dy 1 " << j << ": " << solution_dy[0][0][j] << endl << endl; */
        /* cout << "solution_dy 2 " << j << ": " << solution_dy[0][1][j] << endl << endl; */
    }

    // Calculation of the coefficients of the simplified equation
    vector<double> F1(ns, 0.);
    vector<double> F2(ns, 0.);
    vector< vector <double> > A0(ns, vector<double>(ns,0.));
    for (int j = 0; j < nBasis; ++j) {
        for (int k = 0; k < ns; ++k) {
            for (int l = 0; l < ns; ++l)
                F1[k] += solution_dx[k][l][j]*coefficients[k][j];
        }

        for (int k = 0; k < ns; ++k) {
            for (int l = 0; l < nf; ++l)
                F2[k] += solution_dy[k][l][j]*coefficients_h[l][j];
        }

        for (int k = 0; k < ns; ++k) {
            for (int l = 0; l < ns; ++l) {
                /* cout << 2*solution[k][j]*coefficients[l][j] << endl; */
                A0[k][l] += 2*solution[k][j]*coefficients[l][j];
            }
        }
    }

    /* cout << "First drift" << F1[0]; */
    /* cout << "Second drift" << F2[0]; */

    hi = cholesky(symmetric(A0));
    for (int i = 0; i < ns; ++i) {
        fi[i] = F1[i] + F2[i];
    }
}
