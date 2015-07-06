#include "Problem.hpp"
#include "templates.hpp"
#include "Gaussian_integrator.hpp"
#include "toolbox.hpp"
#include "io.hpp"
#include "global.hpp"

using namespace std;

void Problem::init() {

    this->init_functions();
    this->t_end = 1.;
    this->nf = 2;
    this->ns = 2;
    this->x0 = vector<double>(ns,1.2);
    this->lambdas = {2};
    this->betas   = {sqrt(2)};
    this->bias = vector<double>(nf, 0.);
    this->eig_val_cov = vector<double>(nf, 1.);
    this->sqrt_cov = vector< vector<double> > (nf, vector<double> (nf,0.));
    this->covariance = vector< vector<double> > (nf, vector<double> (nf,0.));
    this->eig_vec_cov = vector< vector<double> > (nf, vector<double> (nf,0.));
    this->inv_cov = vector< vector<double> > (nf, vector<double> (nf,0.));
    this->det_sqrt_cov = 1.;
    for (int i = 0; i < nf; ++i) {
        sqrt_cov[i][i] = 1.;
        eig_vec_cov[i][i] = 1.;
    }
}

// Calculation of the sqrt_cov and mean of the invariant measure.
// The outer loop serves to obtain a more accurate result.
void Problem::update_stats(vector<double> x) {

    Gaussian_integrator gauss = Gaussian_integrator(30, nf);

    int n_iterations = 5, i = 0;
    for (int n = 0; n < n_iterations; ++n) {

        // Calculation of the bias of 'rho'
        bias = vector<double> (nf, 0.);
        for (i = 0; i < nf; ++i) {
            auto lambda = [&] (vector<double> z) -> double {
                vector<double> y = rescale(z);
                return det_sqrt_cov * y[i] * (rho(x,y)/gaussian(z));
            };
            bias[i] = gauss.quadnd(lambda);
        }

        // Calculation of the covariance matrix
        for (i = 0; i < nf; ++i) {
            for (int j = 0; j < nf; ++j) {
                auto lambda = [&] (vector<double> z) -> double {
                    vector<double> y = rescale(z);
                    return det_sqrt_cov * (y[i] - bias[i]) * (y[j] - bias[j]) * (rho(x,y)/gaussian(z));
                };
                covariance[i][j] = gauss.quadnd(lambda);
            }
        }

        // Eigenvalue decomposition
        eig_qr(covariance, eig_vec_cov, eig_val_cov);

        sqrt_cov = eig_vec_cov;
        for (int i = 0; i < nf; i++) {
            for (int j = 0; j < nf; j++) {
                sqrt_cov[i][j] *= sqrt(eig_val_cov[j]);
            }
        }

        // Determinant of sqrt_cov
        for (det_sqrt_cov = 1., i = 0; i < nf; ++i)
            det_sqrt_cov *= sqrt(eig_val_cov[i]);

        // Inverse of covariance matrix
        for (int i = 0; i < nf; ++i) {
            vector<double> rhs_tmp(nf, 0.); rhs_tmp[i] = 1.;
            inv_cov[i] = solve(covariance, rhs_tmp);
        }
        inv_cov = transpose(inv_cov);
    }

    if(DEBUG) {
        cout << endl << "* Covariance matrix of the invariant density" << endl;
        niceMat(covariance);

        cout << endl << "* Inverse of the covariance matrix" << endl;
        niceMat(inv_cov);

        cout << endl << "* Bias of the invariant density" << endl;
        niceVec(bias);

        cout << endl << "* Eigenvectors of the covariance matrix" << endl;
        niceMat(eig_vec_cov);
    }
}

vector<double> Problem::rescale(vector<double> y) {
    vector<double> result(y.size(), 0.);
    for (int i = 0; i < y.size(); ++i) { for (int j = 0; j < y.size(); ++j) {
            result[i] += sqrt_cov[i][j] * y[j];
        }
    }
    return (result + bias);
}

vector<double> Problem::soldrif(vector<double> x) {
    vector<double> result(this->ns,0.);
    Gaussian_integrator gauss = Gaussian_integrator(100,this->nf);
    auto lambda = [&] (vector<double> z) -> vector<double> {
        vector<double> y = rescale(z);
        vector<double> tmp = phi_x(x,y) * a(x,y) + phi(x,y) * stardiv_h(x,y);
        return tmp*(rho(x,y)/gaussian(z));
    };
    result = gauss.quadnd(lambda, result) * det_sqrt_cov;
    return result;
}

vector< vector<double> > Problem::soldiff(vector<double> x) {
    Gaussian_integrator gauss = Gaussian_integrator(100,this->nf);
    vector< vector<double> > result(this->ns,vector<double>(this->ns,0.));
    auto lambda = [&] (vector<double> z) -> vector< vector<double> > {
        vector<double> y = rescale(z);
        vector< vector<double> > tens_prod(this->ns, vector<double>(this->ns, 0.));
        vector<double> slowdrift = this->a(x,y);
        vector<double> solution  = this->phi(x,y);
        for (int i = 0; i < this->ns; ++i) {
            for (int j = 0; j < this->ns; ++j) {
                tens_prod[i][j] = 2*slowdrift[i]*solution[j]*(this->rho(x,y)/gaussian(z));
            }
        }
        return tens_prod;
    };
    result = cholesky(symmetric( gauss.quadnd(lambda, result) * det_sqrt_cov ));
    return result;
}

vector< vector<double> > Problem::phi_x(vector<double> x, vector<double> y) {
    vector< vector<double> > result(this->ns,vector<double>(this->ns,0.));
    result[0][0] = -sin(y[0]*y[1])*sin(x[0]);
    result[1][1] = -sin(x[1]+y[1])*sin(y[0]+y[1]);
    return result;
}

vector< vector<double> > Problem::dax(vector<double> x, vector<double> y) {
    vector< vector<double> > result(this->ns,vector<double>(this->ns,0.));
    result[0][0] = sin(x[0])*(y[0]*cos(y[0]*y[1])*-1.1E1-y[1]*cos(y[0]*y[1])*2.1E1+(y[0]*y[0])*cos(y[0]*y[1])+(y[1]*y[1])*cos(y[0]*y[1])+(y[0]*y[0])*sin(y[0]*y[1])*1.0E1+(y[1]*y[1])*sin(y[0]*y[1])*1.0E1+y[0]*y[1]*cos(y[0]*y[1])*7.0E1-(y[0]*y[0])*y[1]*cos(y[0]*y[1])*6.0E1+(y[0]*y[0]*y[0])*y[1]*cos(y[0]*y[1])*2.0E1)*(-1.0/1.0E1);
    result[1][1] = cos(x[1]-y[0])*(-1.0/2.0)+sin(x[1]-y[0])*(2.1E1/2.0E1)+cos(x[1]+y[0]+y[1]*2.0)*(5.0/2.0)+sin(x[1]+y[0]+y[1]*2.0)*(4.3E1/2.0E1)+(y[0]*y[0])*sin(x[1]-y[0])*3.0-(y[0]*y[0]*y[0])*sin(x[1]-y[0])-y[0]*sin(x[1]+y[0]+y[1]*2.0)*(3.1E1/1.0E1)-y[1]*sin(x[1]+y[0]+y[1]*2.0)*(2.1E1/2.0E1)+(y[0]*y[0])*sin(x[1]+y[0]+y[1]*2.0)*3.0-(y[0]*y[0]*y[0])*sin(x[1]+y[0]+y[1]*2.0)-y[0]*sin(x[1]-y[0])*3.0-y[1]*sin(x[1]-y[0])*(1.0/2.0E1);
    return result;
}
double Problem::stardiv_h(vector<double> x, vector<double> y) {
    double result = 0.;
    result = cos(x[0]+x[1])*sin(y[0]+y[1])+cos(x[0])*cos(y[1])*sin(y[0])+cos(x[0]+x[1])*cos(y[0]+y[1])*(y[0]*(1.0/1.0E1)+y[1]-1.1E1/1.0E1)+cos(x[0])*cos(y[0])*cos(y[1])*(y[1]*(1.0/1.0E1)+pow(y[0]-1.0,3.0)*2.0-1.0/1.0E1);
    return result;
}

vector<double> Problem::grad(vector<double> x, vector<double> y){
    vector<double> result(this->nf);
    result[0] = y[1]*(1.0/1.0E1)+pow(y[0]-1.0,3.0)*2.0-1.0/1.0E1;
    result[1] = y[0]*(1.0/1.0E1)+y[1]-1.1E1/1.0E1;
    return result;
}

vector<double> Problem::phi(vector<double> x, vector<double> y) {
    vector<double> result(this->ns,0.);
    result[0] = sin(y[0]*y[1])*cos(x[0]);
    result[1] = cos(x[1]+y[1])*sin(y[0]+y[1]);
    return result;
}

vector<double> Problem::a(vector<double> x, vector<double> y) {
    vector<double> result(this->ns,0.);
    result[0] = cos(x[0])*(y[0]*cos(y[0]*y[1])*-1.1E1-y[1]*cos(y[0]*y[1])*2.1E1+(y[0]*y[0])*cos(y[0]*y[1])+(y[1]*y[1])*cos(y[0]*y[1])+(y[0]*y[0])*sin(y[0]*y[1])*1.0E1+(y[1]*y[1])*sin(y[0]*y[1])*1.0E1+y[0]*y[1]*cos(y[0]*y[1])*7.0E1-(y[0]*y[0])*y[1]*cos(y[0]*y[1])*6.0E1+(y[0]*y[0]*y[0])*y[1]*cos(y[0]*y[1])*2.0E1)*(1.0/1.0E1);
    result[1] = cos(x[1]-y[0])*(-2.1E1/2.0E1)-sin(x[1]-y[0])*(1.0/2.0)-cos(x[1]+y[0]+y[1]*2.0)*(4.3E1/2.0E1)+sin(x[1]+y[0]+y[1]*2.0)*(5.0/2.0)-(y[0]*y[0])*cos(x[1]-y[0])*3.0+(y[0]*y[0]*y[0])*cos(x[1]-y[0])+y[0]*cos(x[1]+y[0]+y[1]*2.0)*(3.1E1/1.0E1)+y[1]*cos(x[1]+y[0]+y[1]*2.0)*(2.1E1/2.0E1)-(y[0]*y[0])*cos(x[1]+y[0]+y[1]*2.0)*3.0+(y[0]*y[0]*y[0])*cos(x[1]+y[0]+y[1]*2.0)+y[0]*cos(x[1]-y[0])*3.0+y[1]*cos(x[1]-y[0])*(1.0/2.0E1);
    return result;
}

vector<double> Problem::fast_drift_h(vector<double> x, vector<double> y) {
    vector<double> result(this->nf);
    result[0] = cos(x[0])*cos(y[0])*cos(y[1]);
    result[1] = cos(x[0]+x[1])*cos(y[0]+y[1]);
    return result;
}

double Problem::potential(vector<double> x, vector<double> y) {
    double result;
    result = pow(y[1]-1.0,2.0)*(1.0/2.0)+pow(y[0]-1.0,4.0)*(1.0/2.0)+(y[0]*(1.0/5.0)-1.0/5.0)*(y[1]-1.0)*(1.0/2.0);
    return result;
}

double Problem::linearTerm(vector<double> x, vector<double> y){
    double result;
    result = pow(y[0]*(1.0/1.0E1)+y[1]-1.1E1/1.0E1,2.0)*(-1.0/4.0)-pow(y[1]*(1.0/1.0E1)+pow(y[0]-1.0,3.0)*2.0-1.0/1.0E1,2.0)*(1.0/4.0)+pow(y[0]-1.0,2.0)*3.0+1.0/2.0;
    return result;
}

double Problem::rho(vector<double> x, vector<double> y) {
    double result;
    result = exp(pow(y[1]-1.0,2.0)*(-1.0/2.0)-pow(y[0]-1.0,4.0)*(1.0/2.0)-(y[0]*(1.0/5.0)-1.0/5.0)*(y[1]-1.0)*(1.0/2.0))*1.846129117716609E-1;
    return result;
}

vector< vector<double> > Problem::day(vector<double> x, vector<double> y) {
    vector< vector<double> > result(this->ns,vector<double>(this->nf,0.));
    result[0][0] = cos(x[0])*(cos(y[0]*y[1])*-1.1E1+y[0]*cos(y[0]*y[1])*2.0+y[1]*cos(y[0]*y[1])*7.0E1+y[0]*sin(y[0]*y[1])*2.0E1+(y[1]*y[1]*y[1])*cos(y[0]*y[1])*1.0E1+(y[1]*y[1])*sin(y[0]*y[1])*2.1E1-(y[1]*y[1]*y[1])*sin(y[0]*y[1])-y[0]*y[1]*cos(y[0]*y[1])*1.2E2+y[0]*y[1]*sin(y[0]*y[1])*1.1E1+(y[0]*y[0])*y[1]*cos(y[0]*y[1])*7.0E1-y[0]*(y[1]*y[1])*sin(y[0]*y[1])*7.0E1-(y[0]*y[0])*y[1]*sin(y[0]*y[1])+(y[0]*y[0])*(y[1]*y[1])*sin(y[0]*y[1])*6.0E1-(y[0]*y[0]*y[0])*(y[1]*y[1])*sin(y[0]*y[1])*2.0E1)*(1.0/1.0E1);
    result[0][1] = cos(x[0])*(cos(y[0]*y[1])*-2.1E1+y[0]*cos(y[0]*y[1])*7.0E1+y[1]*cos(y[0]*y[1])*2.0+y[1]*sin(y[0]*y[1])*2.0E1-(y[0]*y[0])*cos(y[0]*y[1])*6.0E1+(y[0]*y[0]*y[0])*cos(y[0]*y[1])*3.0E1+(y[0]*y[0])*sin(y[0]*y[1])*1.1E1-(y[0]*y[0]*y[0])*sin(y[0]*y[1])+y[0]*y[1]*sin(y[0]*y[1])*2.1E1+y[0]*(y[1]*y[1])*cos(y[0]*y[1])*1.0E1-y[0]*(y[1]*y[1])*sin(y[0]*y[1])-(y[0]*y[0])*y[1]*sin(y[0]*y[1])*7.0E1+(y[0]*y[0]*y[0])*y[1]*sin(y[0]*y[1])*6.0E1-(y[0]*y[0]*y[0]*y[0])*y[1]*sin(y[0]*y[1])*2.0E1)*(1.0/1.0E1);
    result[1][0] = cos(x[1]-y[0])*(7.0/2.0)-sin(x[1]-y[0])*(2.1E1/2.0E1)+cos(x[1]+y[0]+y[1]*2.0)*(2.8E1/5.0)+sin(x[1]+y[0]+y[1]*2.0)*(4.3E1/2.0E1)+(y[0]*y[0])*cos(x[1]-y[0])*3.0-(y[0]*y[0])*sin(x[1]-y[0])*3.0+(y[0]*y[0]*y[0])*sin(x[1]-y[0])-y[0]*cos(x[1]+y[0]+y[1]*2.0)*6.0-y[0]*sin(x[1]+y[0]+y[1]*2.0)*(3.1E1/1.0E1)-y[1]*sin(x[1]+y[0]+y[1]*2.0)*(2.1E1/2.0E1)+(y[0]*y[0])*cos(x[1]+y[0]+y[1]*2.0)*3.0+(y[0]*y[0])*sin(x[1]+y[0]+y[1]*2.0)*3.0-(y[0]*y[0]*y[0])*sin(x[1]+y[0]+y[1]*2.0)-y[0]*cos(x[1]-y[0])*6.0+y[0]*sin(x[1]-y[0])*3.0+y[1]*sin(x[1]-y[0])*(1.0/2.0E1);
    result[1][1] = cos(x[1]-y[0])*(1.0/2.0E1)+cos(x[1]+y[0]+y[1]*2.0)*(1.21E2/2.0E1)+sin(x[1]+y[0]+y[1]*2.0)*(4.3E1/1.0E1)-y[0]*sin(x[1]+y[0]+y[1]*2.0)*(3.1E1/5.0)-y[1]*sin(x[1]+y[0]+y[1]*2.0)*(2.1E1/1.0E1)+(y[0]*y[0])*sin(x[1]+y[0]+y[1]*2.0)*6.0-(y[0]*y[0]*y[0])*sin(x[1]+y[0]+y[1]*2.0)*2.0;
    return result;
}

vector<double> Problem::drif(vector<double> x, vector<double> y) {
    vector<double> result(2*this->nf,0.);
    result[0] = y[1]*(-1.0/1.0E1)-pow(y[0]-1.0,3.0)*2.0+1.0/1.0E1;
    result[1] = y[0]*(-1.0/1.0E1)-y[1]+1.1E1/1.0E1;
    result[2] = y[3]*(-1.0/1.0E1)-y[2]*pow(y[0]-1.0,2.0)*6.0+cos(x[0])*cos(y[0])*cos(y[1]);
    result[3] = y[2]*(-1.0/1.0E1)-y[3]+cos(x[0]+x[1])*cos(y[0]+y[1]);
    return result;
}

vector<double> Problem::diff(vector<double> x, vector<double> y) {
    vector<double> result(2*nf,0.);
    result[0] = sqrt(2.0);
    result[1] = sqrt(2.0);
    return result;
}
double a1 (vector<double> x, vector<double> y) {
    return cos(x[0])*(y[0]*cos(y[0]*y[1])*-1.1E1-y[1]*cos(y[0]*y[1])*2.1E1+(y[0]*y[0])*cos(y[0]*y[1])+(y[1]*y[1])*cos(y[0]*y[1])+(y[0]*y[0])*sin(y[0]*y[1])*1.0E1+(y[1]*y[1])*sin(y[0]*y[1])*1.0E1+y[0]*y[1]*cos(y[0]*y[1])*7.0E1-(y[0]*y[0])*y[1]*cos(y[0]*y[1])*6.0E1+(y[0]*y[0]*y[0])*y[1]*cos(y[0]*y[1])*2.0E1)*(1.0/1.0E1);
}

double a2 (vector<double> x, vector<double> y) {
    return cos(x[1]-y[0])*(-2.1E1/2.0E1)-sin(x[1]-y[0])*(1.0/2.0)-cos(x[1]+y[0]+y[1]*2.0)*(4.3E1/2.0E1)+sin(x[1]+y[0]+y[1]*2.0)*(5.0/2.0)-(y[0]*y[0])*cos(x[1]-y[0])*3.0+(y[0]*y[0]*y[0])*cos(x[1]-y[0])+y[0]*cos(x[1]+y[0]+y[1]*2.0)*(3.1E1/1.0E1)+y[1]*cos(x[1]+y[0]+y[1]*2.0)*(2.1E1/2.0E1)-(y[0]*y[0])*cos(x[1]+y[0]+y[1]*2.0)*3.0+(y[0]*y[0]*y[0])*cos(x[1]+y[0]+y[1]*2.0)+y[0]*cos(x[1]-y[0])*3.0+y[1]*cos(x[1]-y[0])*(1.0/2.0E1);
}

double dxa11 (vector<double> x, vector<double> y) {
    return sin(x[0])*(y[0]*cos(y[0]*y[1])*-1.1E1-y[1]*cos(y[0]*y[1])*2.1E1+(y[0]*y[0])*cos(y[0]*y[1])+(y[1]*y[1])*cos(y[0]*y[1])+(y[0]*y[0])*sin(y[0]*y[1])*1.0E1+(y[1]*y[1])*sin(y[0]*y[1])*1.0E1+y[0]*y[1]*cos(y[0]*y[1])*7.0E1-(y[0]*y[0])*y[1]*cos(y[0]*y[1])*6.0E1+(y[0]*y[0]*y[0])*y[1]*cos(y[0]*y[1])*2.0E1)*(-1.0/1.0E1);
}

double dxa12 (vector<double> x, vector<double> y) {
    return 0.0;
}

double dxa21 (vector<double> x, vector<double> y) {
    return 0.0;
}

double dxa22 (vector<double> x, vector<double> y) {
    return cos(x[1]-y[0])*(-1.0/2.0)+sin(x[1]-y[0])*(2.1E1/2.0E1)+cos(x[1]+y[0]+y[1]*2.0)*(5.0/2.0)+sin(x[1]+y[0]+y[1]*2.0)*(4.3E1/2.0E1)+(y[0]*y[0])*sin(x[1]-y[0])*3.0-(y[0]*y[0]*y[0])*sin(x[1]-y[0])-y[0]*sin(x[1]+y[0]+y[1]*2.0)*(3.1E1/1.0E1)-y[1]*sin(x[1]+y[0]+y[1]*2.0)*(2.1E1/2.0E1)+(y[0]*y[0])*sin(x[1]+y[0]+y[1]*2.0)*3.0-(y[0]*y[0]*y[0])*sin(x[1]+y[0]+y[1]*2.0)-y[0]*sin(x[1]-y[0])*3.0-y[1]*sin(x[1]-y[0])*(1.0/2.0E1);
}

void Problem::init_functions() {

    fsplit[0] = a1;
    fsplit[1] = a2;

    fxsplit[0][0] = dxa11;
    fxsplit[0][1] = dxa12;
    fxsplit[1][0] = dxa21;
    fxsplit[1][1] = dxa22;
}