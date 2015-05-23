#include "Problem.hpp"
#include "templates.hpp"
#include "Gaussian_integrator.hpp"
#include "toolbox.hpp"

using namespace std;

void Problem::init() {

    this->t_end = 1.;
    this->nf = 1;
    this->d = 1;
    this->x0 = vector<double>(d,1.2);
    this->lambdas = {2};
    this->betas   = {sqrt(2)};

    this->sigmas = {0.7};
    /* this->sigmas = vector<double>(this->nf, 0.); */
    /* for (int i = 0; i < this->nf; ++i) { */
    /*     this->sigmas[i] = sqrt(betas[i]*betas[i] / (2*lambdas[i])); */
    /* } */
}

double Problem::stardiv_h(vector<double> x, vector<double> y) {
    double result = 0.;
    vector<double> fast_h = this->fast_drift_h(x,y);
    vector<double> gradient_v = this->grad(x,y);
    vector< vector<double> > gradient_h = this->grad_h(x,y);
    for (int i = 0; i < this->nf; ++i) {
        result += gradient_v[i] * fast_h[i] - gradient_h[i][i];
    }
    return result;
}

vector<double> Problem::soldrif(vector<double> x) {
    vector<double> result(this->d,0.);
    Gaussian_integrator gauss = Gaussian_integrator(100,this->nf);
    auto lambda = [&] (vector<double> y) -> vector<double> {
        double div_h = stardiv_h(x,y);
        vector<double> tmp(this->d,0.);
        vector<double> slow_drift = this->a(x,y);
        vector<double> solution = this->phi(x,y);
        vector< vector<double> > diff_phi_x = this->phi_x(x,y);
        for (int i = 0; i < this->d; ++i) {
            tmp = tmp + diff_phi_x[i]*slow_drift[i];
        }
        tmp = tmp + solution*div_h;
        return tmp*(this->rho(x,y)/gaussian(y,sigmas));
    };
    result = gauss.quadnd(lambda, sigmas, result);
    return result;
}

vector< vector<double> > Problem::soldiff(vector<double> x) {
    Gaussian_integrator gauss = Gaussian_integrator(100,this->nf);
    vector< vector<double> > result(this->d,vector<double>(this->d,0.));
    auto lambda = [&] (vector<double> y) -> vector< vector<double> > {
        vector< vector<double> > tens_prod(this->d, vector<double>(this->d, 0.));
        vector<double> slowdrift = this->a(x,y);
        vector<double> solution  = this->phi(x,y);
        for (int i = 0; i < this->d; ++i) {
            for (int j = 0; j < this->d; ++j) {
                tens_prod[i][j] = 2*slowdrift[i]*solution[j]*(this->rho(x,y)/gaussian(y,sigmas));
            }
        }
        return tens_prod;
    };
    result = cholesky(symmetric( gauss.quadnd(lambda, sigmas, result) ));
    return result;
}

vector< vector<double> > Problem::day(vector<double> x, vector<double> y) {
    vector< vector<double> > result(this->d,vector<double>(this->nf,0.));
    result[0][0] = cos(x[0])*(2*cos(y[0]) - 2*y[0]*sin(y[0]) + cos(y[0]));
    return result;
}

vector<double> Problem::drif(vector<double> x, vector<double> y) {
    vector<double> result(2*this->nf,0.);
    result[0] = -2*y[0];
    result[1] = -2*y[1] + cos(x[0])*cos(y[0]);
    return result;
}

vector<double> Problem::diff(vector<double> x, vector<double> y) {
    vector<double> result(2*nf,0.);
    result[0] = sqrt(2.);
    return result;
}

vector< vector<double> > Problem::phi_x(vector<double> x, vector<double> y) {
    vector< vector<double> > result(this->d,vector<double>(this->d,0.));
    result[0][0] = -sin(x[0])*sin(y[0]);
    return result;
}

vector< vector<double> > Problem::dax(vector<double> x, vector<double> y) {
    vector< vector<double> > result(this->d,vector<double>(this->d,0.));
    result[0][0] = -sin(x[0])*(sin(y[0])+y[0]*cos(y[0])*2.0);
    return result;
}

vector< vector<double> > Problem::grad_h(vector<double> x, vector<double> y) {
    vector< vector<double> > result(this->nf, vector<double>(this->nf, 0.));
    result[0][0] = -cos(x[0])*sin(y[0]);
    return result;
}

vector<double> Problem::grad(vector<double> x, vector<double> y){
    vector<double> result(this->nf);
    result[0] = y[0]*2.0;
    return result;
}

vector<double> Problem::phi(vector<double> x, vector<double> y) {
    vector<double> result(this->d,0.);
    result[0] = cos(x[0])*sin(y[0]);
    return result;
}

vector<double> Problem::a(vector<double> x, vector<double> y) {
    vector<double> result(this->d,0.);
    result[0] = cos(x[0])*(sin(y[0])+y[0]*cos(y[0])*2.0);
    return result;
}

vector<double> Problem::fast_drift_h(vector<double> x, vector<double> y) {
    vector<double> result(this->nf);
    result[0] = cos(x[0])*cos(y[0]);
    return result;
}

double Problem::potential(vector<double> x, vector<double> y) {
    double result = 0.;
    result = y[0]*y[0]+5.723649429247001E-1;
    return result;
}

double Problem::linearTerm(vector<double> x, vector<double> y){
    double result;
    result = -y[0]*y[0]+1.0;
    return result;
}

double Problem::rho(vector<double> x, vector<double> y) {
    double result = 0.;
    result = exp(-y[0]*y[0]-5.723649429247001E-1)*1.000000000070638;
    return result;
}
