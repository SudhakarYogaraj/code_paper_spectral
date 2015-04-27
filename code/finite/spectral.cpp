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

    // Calculation of the coefficients
    for (int j = 0; j < nBasis; ++j) {
        vector<int> multIndex = solver.ind2mult(j);
        for (int k = 0; k < ns; ++k) {
            for (int l = 0; l < ns; ++l) {
                auto lambda = [&] (vector<double> y) -> double {
                    return problem.dax(xt,y)[k][l]*hermiteM(multIndex, y, sigmas);
                };
                coefficients_dx[k][l][j] = gauss_hermite_nD(lambda,sigmas);
            }
            auto lambda = [&] (vector<double> y) -> double {
                return problem.a(xt,y)[k]*hermiteM(multIndex, y, sigmas);
            };
            coefficients[k][j] = gauss_hermite_nD(lambda,sigmas);
        }
        for (int k = 0; k < nf; ++k) {
            auto lambda = [&] (vector<double> y) -> double {
                return problem.fast_drift_h(xt,y)[k]*hermiteM(multIndex, y, sigmas);
            };
            coefficients_h[k][j] = gauss_hermite_nD(lambda,sigmas);
        }
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
        int sum = accumulate(thisMult.begin(),thisMult.end(),0);
        if (sum < degree) {
            for (int l = 0; l < nf; ++l) {
                vector<int> newMult(nf, 0);
                for (int m = 0; m < nf; ++m) {
                    newMult[m] = thisMult[m];
                }
                newMult[l] ++;
                int newInd = solver.mult2ind(newMult);
                for (int m = 0; m < ns; ++m) {
                    solution_dy[m][l][j] = solution[m][newInd]*sqrt(newMult[l])/sigmas[l];
                }
            }
        }
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
