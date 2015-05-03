#include "aux.hpp"
#include "Problem.hpp"
#include "Gaussian_integrator.hpp"
#include "Solver_spectral.hpp"
#include "templates.hpp"

using namespace std;

void Solver_spectral::set(double p, int n)
{
    this->p = p;
    this->n_mcmc = 100000;
    this->degree = 10;
    this->nNodes = 20;
    this->nvars = n;

    this-> nBasis = bin(this->degree + this->nvars, this->nvars);
    this->ind2mult_aux = vector< vector<int> >(this->nBasis, vector<int>(this->nvars,0));
    this->mult2ind_aux = vector<int>(pow(degree + 1, this->nvars), -1);

    vector<int> currentMult(nvars,0);
    for (int i = 1; i < nBasis; ++i) {
        int sum = 0;
        for (int j = 0; j < nvars; ++j) {
            sum += currentMult[j];
        }
        if (sum < this->degree) {
            currentMult[nvars-1] ++;
        } else {
            int auxIndex = nvars - 1;
            while (currentMult[auxIndex] == 0) {
                auxIndex --;
            }
            currentMult[auxIndex] = 0;
            currentMult[auxIndex-1] = currentMult[auxIndex-1] + 1;
        }
        for (int j = 0; j < nvars; ++j) {
            ind2mult_aux[i][j] = currentMult[j];
        }
        mult2ind_aux[canonicalInd(currentMult, nvars, this->degree)] = i;
    }
}

int Solver_spectral::mult2ind(vector<int> alpha) {
    int canonicalInd = 0.;
    for (int j = 0; j < this->nvars; ++j) {
        canonicalInd += alpha[j]*pow(this->degree + 1, (this->nvars - 1) -j);
    }
    return this->mult2ind_aux[canonicalInd];
}

vector<int> Solver_spectral::ind2mult(int ind) {
    return ind2mult_aux[ind];
}

void Solver_spectral::estimator(Problem &problem, vector<double> x, vector<double>& fi, vector< vector<double> >& hi, double t) {

    int nBasis = this->nBasis;
    int nf     = problem.nf;
    int ns     = problem.d;
    Gaussian_integrator gauss = Gaussian_integrator(nNodes,nf);

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
        vector<int> multIndex = this->ind2mult(i);

        vector<double> v0(ns,0.);
        vector<double> h0(nf,0.);
        vector< vector<double> > m0(ns, vector<double> (ns,0.));

        auto lambda = [&] (vector<double> y) -> vector<double> {
            return problem.a(x,y)*monomial(multIndex, y, sigmas); };
        auto lambda_dx = [&] (vector<double> y) -> vector< vector<double> > {
            return problem.dax(x,y)*monomial(multIndex, y, sigmas); };
        auto lambda_h = [&] (vector<double> y) -> vector<double> {
            return problem.fast_drift_h(x,y)*monomial(multIndex, y, sigmas); };

        vector<double> result = gauss.quadnd(lambda,sigmas,v0);
        vector< vector<double> > result_dx = gauss.quadnd(lambda_dx,sigmas,m0);
        vector<double> result_h = gauss.quadnd(lambda_h,sigmas,h0);

        for (int j = 0; j < ns; ++j) {
            coefficients[j][i] = result[j];
            for (int k = 0; k < ns; ++k) {
                coefficients_dx[j][k][i] = result_dx[j][k];
            }
        }
        for (int j = 0; j < nf; ++j) {
            coefficients_h[j][i] = result_h[j];
        }
    }

    for (int i = 0; i < ns; ++i) {
        for (int j = 0; j < ns; ++j) {
            coefficients_dx[i][j] = hcoeffs(coefficients_dx[i][j],nf,degree);
        }
        coefficients[i] = hcoeffs(coefficients[i],nf,degree);
    }
    for (int i = 0; i < nf; ++i) {
        coefficients_h[i] = hcoeffs(coefficients_h[i],nf,degree);
    }

    // Solution of the Poisson equation
    vector< vector<double> > solution(ns, vector<double>(nBasis,0.));
    vector< vector < vector<double> > > solution_dx(ns, vector< vector<double> >(ns, vector<double>(nBasis, 0.)));
    vector< vector< vector<double> > > solution_dy(ns, vector< vector <double> >(nf, vector<double>(nBasis,0.)));
    for (int i = 1; i < nBasis; ++i) {
        double eig = 0.;
        for (int j = 0; j < nf; ++j) {
            eig += this->ind2mult(i)[j]*problem.lambdas[j];
        }
        for (int j = 0; j < ns; ++j) {
            solution[j][i] = coefficients[j][i]/eig;
            for (int k = 0; k < ns; ++k)
                solution_dx[j][k][i] = coefficients_dx[j][k][i]/eig;
        }
    }

    // y-Derivatives of the solution
    for (int i = 0; i < nBasis; ++i) {
        vector<int> mult = this->ind2mult(i);
        if (accumulate (mult.begin(), mult.end(),0) == degree) {
            continue;
        }
        for (int j = 0; j < nf; ++j) {
            vector<int> newMult = mult;
            newMult[j] ++;
            int newInd = this->mult2ind(newMult);
            for (int k = 0; k < ns; ++k) {
                solution_dy[k][j][i] = solution[k][newInd]*sqrt(newMult[j])/sigmas[j];
                /* cout << i << "   " <<  solution_dy[k][j][i] << endl; */
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
    }

    hi = cholesky(symmetric(A0));
    for (int i = 0; i < ns; ++i)
        fi[i] = F1[i] + F2[i];
}
