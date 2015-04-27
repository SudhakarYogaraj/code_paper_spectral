#include "header.h"

// Main method
void solve_spectral(Problem &problem, \
                    Solver_spectral &solver, \
                    vector<double> x, \
                    vector<double>& fi, \
                    vector< vector<double> >& hi, \
                    double t) {

    int nBasis = solver.nBasis;
    int nf     = problem.nf;
    int degree = solver.degree;
    int ns     = problem.d;

    // Eigenvalues
    vector<double> sigmas(nf, 0.);
    for (int i = 0; i < nf; ++i) {
        double aux = problem.betas[i];
        sigmas[i] = sqrt(0.5*aux*aux/problem.lambdas[i]);
    }

    // Expansion of right-hand side of the Poisson equation
    vector< vector<double> > coefficients(ns, vector<double>(nBasis, 0.));
    vector< vector <vector<double> > > coefficients_dx(ns, vector< vector<double> >(ns, vector<double>(nBasis, 0.)));
    vector< vector<double> > coefficients_h(nf, vector<double>(nBasis,0.));
    for (int i = 0; i < nBasis; ++i) {
        vector<int> multIndex = solver.ind2mult(i);
        for (int j = 0; j < ns; ++j) {
            for (int k = 0; k < ns; ++k) {
                auto lambda = [&] (vector<double> y) -> double {
                    return problem.dax(x,y)[j][k]*hermiteM(multIndex, y, sigmas);
                };
                coefficients_dx[j][k][i] = gauss_hermite_nD(lambda,sigmas);
            }
            auto lambda = [&] (vector<double> y) -> double {
                return problem.a(x,y)[j]*hermiteM(multIndex, y, sigmas);
            };
            coefficients[j][i] = gauss_hermite_nD(lambda,sigmas);
        }
        for (int j = 0; j < nf; ++j) {
            auto lambda = [&] (vector<double> y) -> double {
                return problem.fast_drift_h(x,y)[j]*hermiteM(multIndex, y, sigmas);
            };
            coefficients_h[j][i] = gauss_hermite_nD(lambda,sigmas);
        }
    }

    // Solution of the Poisson equation
    vector< vector<double> > solution(ns, vector<double>(nBasis,0.));
    vector< vector < vector<double> > > solution_dx(ns, vector< vector<double> >(ns, vector<double>(nBasis, 0.)));
    vector< vector< vector<double> > > solution_dy(ns, vector< vector <double> >(nf, vector<double>(nBasis,0.)));
    for (int i = 1; i < nBasis; ++i) {
        double eig = 0.;
        for (int j = 0; j < nf; ++j) {
            eig += solver.ind2mult(i)[j]*problem.lambdas[j];
        }
        for (int j = 0; j < ns; ++j) {
            solution[j][i] = coefficients[j][i]/eig;
            for (int k = 0; k < ns; ++k)
                solution_dx[j][k][i] = coefficients_dx[j][k][i]/eig;
        }
    }

    // y-Derivatives of the solution
    for (int i = 0; i < nBasis; ++i) {
        vector<int> mult = solver.ind2mult(i);
        if (accumulate (mult.begin(), mult.end(),0) == degree) continue;
        for (int j = 0; j < nf; ++j) {
            vector<int> newMult = mult;
            newMult[j] ++;
            int newInd = solver.mult2ind(newMult);
            for (int k = 0; k < ns; ++k) {
                solution_dy[k][j][i] = solution[k][newInd]*sqrt(newMult[j])/sigmas[j];
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


    hi = cholesky(symmetric(A0));
    for (int i = 0; i < ns; ++i)
        fi[i] = F1[i] + F2[i];
}
