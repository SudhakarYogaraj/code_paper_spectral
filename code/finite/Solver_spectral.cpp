/* TODO: Linear system solution (urbain, Sat 09 May 2015 14:32:40 BST) */
/* TODO: Integration over infinite domain: which variance for gaussian(urbain, Sat 09 May 2015 14:33:06 BST) */
/* TODO: Convergence of hermite functions? (urbain, Sat 09 May 2015 14:33:28 BST) */

#include "structures.hpp"
#include "aux.hpp"
#include "Problem.hpp"
#include "Gaussian_integrator.hpp"
#include "Solver_spectral.hpp"
#include "templates.hpp"
#include <iomanip>

using namespace std;

Solver_spectral::Solver_spectral(int degree, int nNodes)
{
    this->degree = degree;
    this->nNodes = nNodes;
}

/* vector<double> hermite_transform(function<double(vector<double>)> f, vector<double> sigmas) */
/* } */

vector<double> aux( vector< vector<double> > mat, vector<double> coefficients, vector<double> weights) {

    /* FIXME: bad conditioning (urbain, Tue 12 May 2015 18:53:50 BST) */

    int nb = mat.size();

    for (int i = 0; i < nb; ++i) {
        cout << weights[i] << endl;
    }

    vector< vector<double> > system_mat(nb, vector<double>(nb, 0.));
    for (int i = 1; i < nb; ++i) {
        for (int j = 1; j < nb; ++j) {
            system_mat[i][j] = mat[i][j] - 2*weights[i]/weights[0]*mat[i][0] + weights[i]*weights[i]/(weights[0]*weights[0])*mat[0][0];
        }
    }
    system_mat[0][0] = 1.;


/*     // print */
/*     for (int i = 0; i < nb; ++i) { */
/*         for (int j = 0; j < nb; ++j) { */
/*             cout << setw(10) <<  mat[i][j] << " " ; */
/*         } */
/*         cout << endl; */
/*     } */
/*     cout << "*** " << endl; */

    /* // print */
    /* for (int i = 0; i < nb; ++i) { */
    /*     for (int j = 0; j < nb; ++j) { */
    /*         cout << setw(10) <<  system_mat[i][j] - mat[i][j] << " " ; */
    /*     } */
    /*     cout << endl; */
    /* } */


    vector<double> system_rhs(nb, 0.);
    for (int i = 0; i < nb; ++i) {
        system_rhs[i] = coefficients[i] - weights[i]/weights[0] * coefficients[0];
    }

    /* for (int i = 0; i < nb; ++i) { */
    /*     cout << coefficients[i] << endl; */
    /* } */
    /* cout << " --- " << endl; */

    for (int i = 0; i < nb; ++i) {
        cout << system_rhs[i] << endl;
    }

    vector<double> result = solve(system_mat, system_rhs);
    for (int i = 1; i < nb; ++i) {
        result[0] += weights[i]/weights[0] * result[i];
    }
    return result;
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
    /* vector<double> sigmas_hf(1, 1./sqrt(2.)); */
    vector<double> sigmas_hf(1, 0.7 );

    // Variance parameter for integration
    vector<double> sigmas_quad(1, 1./sqrt(2.));

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


    // weights(i,j) = int (phi_i, e^(-V) )
    vector<double> weights(nb, 0.);
    for (int i = 0; i < nb; ++i) {
        vector<int> m = ind2mult(i, degree, nf);
        auto lambda = [&] (vector<double> y) -> double {
            return monomial(m, y, sigmas_hf) * sqrt(problem.rho(x,y) * gaussian(y,sigmas_hf)); };
        weights[i] = gauss.flatquadnd(lambda, sigmas_hf);
    }

    // mat(i,j) = int ( L phi_i, phi_j)
    vector< vector<double> > mat(nb, vector<double>(nb, 0.));
    for (int i = 0; i < nb; ++i) {
        vector<int> m1 = ind2mult(i, degree, nf);
        for (int j = 0; j < nf; ++j) {
            mat[i][i] += m1[j]/(sigmas_hf[j]*sigmas_hf[j]);
        }
        for (int j = 0; j < nb; ++j) {
            vector<int> m2 = ind2mult(j, degree, nf);
            auto lambda = [&] (vector<double> y) -> double {
                double tmp = 1/(2*pow(sigmas_hf[0], 2)) - pow(y[0], 2)/(4*pow(sigmas_hf[0], 4));
                return (problem.linearTerm(x,y) - tmp) * monomial(m1, y, sigmas_hf) * monomial(m2, y, sigmas_hf) * gaussian(y,sigmas_hf); };
            mat[i][j] += gauss.flatquadnd(lambda, sigmas_hf);
        }
    }

    // print
    /* for (int i = 0; i < nb; ++i) { */
    /*     for (int j = 0; j < nb; ++j) { */
    /*         cout << setw(10) <<  mat[i][j] << " " ; */
    /*     } */
    /*     cout << endl; */
    /* } */

    for (int i = 0; i < ns; ++i) {
        /* vector<double> tmp = aux(mat, coefficients[i], weights); */
        /* mat[0][0] = 1.0; solution[i] = solve(mat, coefficients[i]); */
        /* cout << " ***" << endl; */
        /* for (int j = 0; j < nb; ++j) { */
        /*     cout << solution[0][j] - tmp[j] << endl; */
        /* } */
        /* cout << " ***" << endl; */
        /* for (int j = 0; j < nb; ++j) { */
        /*     cout << solution[0][j] << endl; */
        /* } */
        solution[i] = aux(mat, coefficients[i], weights);
        for (int j = 0; j < ns; ++j) {
            solution_dx[i][j] = aux(mat, coefficients_dx[i][j], weights);
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
