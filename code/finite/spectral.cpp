#include "header.h"

/* TODO: Take problem as input (urbain, Sun 19 Apr 2015 13:45:33 BST) */

// Canonical integer associated with a multindex
int canonicalInd(vector<int> alpha, int n, int degree) {
    int toReturn; 
    for (int j = 0; j < n; ++j) {
        toReturn += alpha[j]*pow(degree + 1, n -j - 1);
    }
    return toReturn;
}

// Main method
void solve_spectral(Problem &problem, \
                    Solver &solver, \
                    vector<double> xt, \
                    vector<double>& fi, \
                    vector< vector<double> >& hi, \
                    int seed, \
                    double t) { 

    // Degree of polynomials approximation.
    int degree = 6;

    // Number of fast and slow variables
    int nf = problem.nf;
    int ns = problem.d;

    // Eigenvalues
    vector<double> sigmas(nf, 0.);
    for (int i = 0; i < nf; ++i) {
        double aux = problem.betas[i];
        sigmas[i] = sqrt(0.5*aux*aux/problem.lambdas[i]);
    }

    // Number of polynomials in the basis
    int nBasis = bin(degree + nf, nf);

    // Relation linear index - multiindex
    vector< vector<int> > ind2mult(nBasis, vector<int>(nf,0));
    vector<int> mult2ind(pow(degree + 1, nf), -1);
    vector<int> currentMult(nf,0);
    for (int i = 0; i < nBasis; ++i) {
        for (int j = 0; j < nf; ++j) {
            ind2mult[i][j] = currentMult[j];
        }
        mult2ind[canonicalInd(currentMult, nf, degree)] = i;
        int sum = 0;
        for (int j = 0; j < nf; ++j) {
            sum += currentMult[j];
        }
        if (sum < degree) {
            currentMult[nf-1] += 1;
        } else {
            int auxIndex = nf - 1;
            while (currentMult[auxIndex] == 0) {
                auxIndex --;
            }
            currentMult[auxIndex] = 0;
            currentMult[auxIndex-1] ++;
        } 
    }

    // Parameters for random numbers
    default_random_engine generator; generator.seed(seed);
    normal_distribution<double> distribution(0.0,1.0);

    // Expansion of right-hand side of the Poisson equation
    vector< vector<double> > coefficients(ns, vector<double>(nBasis, 0.));
    vector< vector <vector<double> > > coefficients_dx(ns, vector< vector<double> >(ns, vector<double>(nBasis, 0.)));
    vector< vector<double> > coefficients_h(nf, vector<double>(nBasis,0.));
    for (int j = 0; j < nBasis; ++j) {
        vector<int> multIndex = ind2mult[j];
        int N_mc = 100000;
        double sum = 0.;

        // Monte-Carlo to compute the coefficients
        for (int k = 0; k < N_mc; ++k) {

            vector<double> randn(nf, 0.);

            for (int l = 0; l < nf; ++l) 
                randn[l] = distribution(generator);
            double h_eval = hermiteM(multIndex,randn,sigmas);
            
            vector<double> slow_drift = problem.a(xt,randn);
            vector< vector<double> > slow_drift_dx = problem.dax(xt,randn);
            for (int l = 0; l < ns; ++l) {
                coefficients[l][j] += h_eval*slow_drift[l];
                for (int m = 0; m < ns; ++m) {
                    coefficients_dx[l][m][j] += h_eval*slow_drift_dx[l][m];
                }
            }
            vector<double> fast_drift_aux = problem.fast_drift_h(xt,randn);
            for (int l = 0; l < nf; ++l) {
                coefficients_h[l][j] += h_eval*fast_drift_aux[l];
            }
        }
        for (int l = 0; l < ns; ++l) {
            coefficients[l][j] /= N_mc;
            for (int m = 0; m < ns; ++m)
                coefficients_dx[l][m][j] /= N_mc;
        }
        for (int l = 0; l < nf; ++l) {
            coefficients_h[l][j] /= N_mc;
        }
    }

    // Solution of the Poisson equation
    vector< vector<double> > solution(ns, vector<double>(nBasis,0.));
    vector< vector < vector<double> > > solution_dx(ns, vector< vector<double> >(ns, vector<double>(nBasis, 0.)));
    vector< vector< vector<double> > > solution_dy(ns, vector< vector <double> >(nf, vector<double>(nBasis,0.)));
    for (int j = 0; j < nBasis; ++j) {
        double eig = 0.;
        for (int k = 0; k < nf; ++k) {
            eig += ind2mult[j][k]*problem.lambdas[k];
        }
        for (int l = 0; l < ns; ++l) {
            solution[l][j] = coefficients[l][j]/eig; 
            for (int m = 0; m < ns; ++m)
                solution_dx[l][m][j] = coefficients_dx[l][m][j]/eig;
        }
        vector<int> thisMult = ind2mult[j];
        int sum = 0;
        for (int l = 0; l < nf; ++l) {
            sum += thisMult[l];
        }
        if (sum < degree) {
            for (int l = 0; l < nf; ++l) {
                vector<int> newMult(nf, 0);
                for (int m = 0; m < nf; ++m) {
                    newMult[m] = thisMult[m];
                }
                newMult[l] = newMult[l] + 1;
                int newInd = mult2ind[canonicalInd(newMult, nf, degree)];
                for (int m = 0; m < ns; ++m)
                    solution_dy[m][l][j] = solution[m][newInd]*sqrt(newMult[l])/sigmas[l];
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
            for (int l = 0; l < ns; ++l) 
                A0[k][l] += 2*solution[k][j]*coefficients[l][j];
        }
    }

    hi = cholesky(A0);
    for (int i = 0; i < ns; ++i) 
        fi[i] = F1[i] + F2[i];
}
