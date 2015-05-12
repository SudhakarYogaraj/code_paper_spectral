/* TODO: Linear system solution (urbain, Sat 09 May 2015 14:32:40 BST) */
/* TODO: Integration over infinite domain: which variance for gaussian(urbain, Sat 09 May 2015 14:33:06 BST) */
/* TODO: Convergence of hermite functions? (urbain, Sat 09 May 2015 14:33:28 BST) */

#include "structures.hpp"
#include "aux.hpp"
#include "Problem.hpp"
#include "Gaussian_integrator.hpp"
#include "Solver_spectral.hpp"
#include "templates.hpp"

using namespace std;

Solver_spectral::Solver_spectral(int degree, int nNodes)
{
    this->degree = degree;
    this->nNodes = nNodes;
}

void Solver_spectral::estimator(Problem &problem, vector<double> x,  vector<SDE_coeffs>& c, double t) {

    int nf     = problem.nf;
    int ns     = problem.d;
    int nb     = bin(degree + nf, nf);

    c = vector<SDE_coeffs> (nb);
    Gaussian_integrator gauss = Gaussian_integrator(nNodes,nf);

    // Eigenvalues
    vector<double> sigmas(nf, 0.);
    for (int i = 0; i < nf; ++i) {
        sigmas[i] = problem.betas[i]*sqrt(0.5/problem.lambdas[i]);
    }

    // Variance parameter for hermite functions
    vector<double> sigmas_hf(1, 1./sqrt(2.));

    // Variance parameter for integration
    vector<double> sigmas_quad(1, 0.8);

    // Expansion of right-hand side of the Poisson equation
    vector< vector<double> > coefficients(ns, vector<double>(nb, 0.));
    vector< vector <vector<double> > > coefficients_dx(ns, vector< vector<double> >(ns, vector<double>(nb, 0.)));
    vector< vector<double> > coefficients_h(nf, vector<double>(nb,0.));
    for (int i = 0; i < nb; ++i) {

        // Graphical progression bar.
        /* cout << "["; */
        /* int bw = 103; */
        /* double progress = ( (double) i) / ( (double) nb); */
        /* int pos = bw * progress; */
        /* for (int j = 0; j < bw; ++j) { */
        /*     if (j < pos) cout << "="; */
        /*     else if (j == pos) cout << ">"; */
        /*     else cout << " "; */
        /* } */
        /* cout << "] " << int(progress * 100.0) << " %\r"; */

        vector<int> multIndex = ind2mult(i, degree, nf);

        vector<double> v0(ns,0.);
        vector<double> h0(nf,0.);
        vector< vector<double> > m0(ns, vector<double> (ns,0.));

        auto lambda = [&] (vector<double> y) -> vector<double> {
            return problem.a(x,y)*monomial(multIndex, y, sigmas_hf)*sqrt(problem.rho(x,y)*gaussian(y,sigmas_hf)); };
        auto lambda_dx = [&] (vector<double> y) -> vector< vector<double> > {
            return problem.dax(x,y)*monomial(multIndex, y, sigmas_hf)*sqrt(problem.rho(x,y)*gaussian(y,sigmas_hf)); };
        auto lambda_h = [&] (vector<double> y) -> vector<double> {
            return problem.fast_drift_h(x,y)*monomial(multIndex, y, sigmas_hf)*sqrt(problem.rho(x,y)*gaussian(y,sigmas_hf)); };

        vector<double> result = gauss.flatquadnd(lambda, sigmas_quad, v0);
        vector< vector<double> > result_dx = gauss.flatquadnd(lambda_dx, sigmas_quad, m0);
        vector<double> result_h = gauss.flatquadnd(lambda_h, sigmas_quad, h0);

        for (int j = 0; j < ns; ++j) {
            coefficients[j][i] = result[j];
            for (int k = 0; k < ns; ++k) {
                coefficients_dx[j][k][i] = result_dx[j][k];
            }
        }
        for (int j = 0; j < nf; ++j) {
            coefficients_h[j][i] = result_h[j];
        }
        cout.flush();
    }

    for (int i = 0; i < ns; ++i) {
        for (int j = 0; j < ns; ++j) {
            coefficients_dx[i][j] = mon2herm(coefficients_dx[i][j],nf,degree);
        }
        coefficients[i] = mon2herm(coefficients[i],nf,degree);
    }
    for (int i = 0; i < nf; ++i) {
        coefficients_h[i] = mon2herm(coefficients_h[i],nf,degree);
    }

    // Solution of the Poisson equation
    vector< vector<double> > solution(ns, vector<double>(nb,0.));
    vector< vector < vector<double> > > solution_dx(ns, vector< vector<double> >(ns, vector<double>(nb, 0.)));
    vector< vector< vector<double> > > solution_dy(ns, vector< vector <double> >(nf, vector<double>(nb,0.)));
    // Construction of the matrix
    vector< vector<double> > mat(nb, vector<double>(nb, 0.));
    for (int i = 0; i < nb; ++i) {
        vector<int> m1 = ind2mult(i, degree, nf);
        for (int j = 0; j < nf; ++j) {
            mat[i][i] += m1[j]/(sigmas_hf[j]*sigmas_hf[j]);
        }
        /* for (int j = 0; j < nb; ++j) { */
        /*     vector<int> m2 = ind2mult(j, degree, nf); */
        /*     auto lambda = [&] (vector<double> y) -> double { */
        /*         return monomial(m1, y, sigmas_hf) * monomial(m2, y, sigmas_hf) * gaussian(y,sigmas_hf) * problem.rho(x,y); }; */
        /*     mat[i][j] = gauss.flatquadnd(lambda, sigmas); */
        /*     /1* cout << mat[i][j] << endl; *1/ */
        /* } */
    }
    mat[0][0] = 1.;
    for (int i = 0; i < ns; ++i) {
        solution[i] = solve(mat, coefficients[i]);
        for (int j = 0; j < ns; ++j) {
            solution_dx[i][j] = solve(mat, coefficients_dx[i][j]);
        }
    }

    // y-Derivatives of the solution
    for (int i = 0; i < nb; ++i) {
        vector<int> mult = ind2mult(i,degree,nf);
        if (accumulate (mult.begin(), mult.end(),0) == degree) {
            continue;
        }
        for (int j = 0; j < nf; ++j) {
            vector<int> newMult = mult;
            newMult[j] ++;
            int newInd = mult2ind(newMult, degree);
            for (int k = 0; k < ns; ++k) {
                solution_dy[k][j][i] = solution[k][newInd]*sqrt(newMult[j])/sigmas[j];
            }
        }
    }

    // Calculation of the coefficients of the simplified equation
    vector<double> F1(ns, 0.);
    vector<double> F2(ns, 0.);
    vector< vector <double> > A0(ns, vector<double>(ns,0.));
    for (int j = 0; j < nb; ++j) {
        for (int k = 0; k < ns; ++k) {
            for (int l = 0; l < ns; ++l)
                F1[k] += solution_dx[k][l][j]*coefficients[l][j];
        }

        for (int k = 0; k < ns; ++k) {
            for (int l = 0; l < nf; ++l)
                F2[k] += solution_dy[k][l][j]*coefficients_h[l][j];
        }

        for (int k = 0; k < ns; ++k) {
            for (int l = 0; l < ns; ++l) {
                A0[k][l] += 2*solution[k][j]*coefficients[l][j];
            }
        }
        c[j].diff =  cholesky(symmetric(A0));
        c[j].drif = F1 + F2;
    }
}
