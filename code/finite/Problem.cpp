#include "Problem.hpp"
#include "aux.hpp"

using namespace std;

void Problem::init() {

    this->t_end = 1.;
    this->nf = 1;
    this->d = 1;
    this->x0 = vector<double>(d,1.2);
    this->lambdas = {2};
    this->betas   = {sqrt(2)};
}

// Potential for gradient structure. sqrt(2) for the diff is assumed.
double Problem::potential(vector<double> x, vector<double> y) {
    return y[0]*y[0] + log(PI)/2.;
}

double Problem::rho(vector<double> x, vector<double> y) {
    return exp(-potential(x,y));
}

vector<double> Problem::grad(vector<double> x, vector<double> y){
    vector<double> result(this->nf);
    result[0] = 2.*y[0];
    return result;
}

double Problem::linearTerm(vector<double> x, vector<double> y){
    double result;
    result = 1.-y[0]*y[0];
    return result;
}

vector<double> Problem::fast_drift_h(vector<double> x, vector<double> y) {
    vector<double> result(this->nf);
    result[0] = cos(x[0])*cos(y[0]);
    return result;
}

vector< vector<double> > Problem::grad_h(vector<double> x, vector<double> y) {
    vector< vector<double> > result(this->nf, vector<double>(this->nf, 0.));
    result[0][0] = -cos(x[0])*sin(y[0]);
    return result;
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

// divide factor 2
// rho = e^(-x^2)/sqrt(pi);
// phi = 2*cos(x)*sin(y)
// dphix = -2*sin(x)*sin(y)
// dphiy = 2*cos(x)*cos(y)
// a   = cos(x)*(2*y*cos(y) + sin(y))
// dax = -sin(x)*(2*y*cos(y) + sin(y))
// day = cos(x)*(2*cos(y) - 2*y*sin(y) + cos(y))
// h   = cos(x)*cos(y)
//
// drift = -2*sin(x)*cos(x) * int (sin(y) * (2*y*cos(y) + sin(y)) * e^(-y^2)/sqrt(pi) )
//         + 2*cos(x)*cos(x) * int (cos(y)*cos(y) * e^(-y^2)/sqrt(pi) )
// diff  = SQRT 2? * cos(x)*cos(x) * int (sin(y) * (2*y*cos(y) + sin(y)) * e^(-y^2)/sqrt(pi) )
vector<double> Problem::soldrif(vector<double> x) {
    vector<double> result(this->d,0.);
    /* cout << "drift second term " <<  cos(x[0])*cos(x[0])*(1. + exp(-1.))/2. << endl; */
    result[0] = -sin(x[0])*cos(x[0])*(1 + exp(-1.0))/2. + cos(x[0])*cos(x[0])*(1. + exp(-1.))/2.;
    return result;
}

// ! Coefficient 2? Seems to have been forgotten in exact solution
vector< vector<double> > Problem::soldiff(vector<double> x) {
    vector< vector<double> > result(this->d,vector<double>(this->d,0.));
    result[0][0] = sqrt(cos(x[0])*cos(x[0])*(1 + exp(-1.0)));
    return result;
}

// Phi = 2*cos(x)*sin(y)
vector<double> Problem::a(vector<double> x, vector<double> y) {
    vector<double> result(this->d,0.);
    result[0] = cos(x[0])*(2*y[0]*cos(y[0]) + sin(y[0]));
    return result;
}

vector< vector<double> > Problem::dax(vector<double> x, vector<double> y) {
    vector< vector<double> > result(this->d,vector<double>(this->d,0.));
    result[0][0] = - sin(x[0])*(2*y[0]*cos(y[0]) + sin(y[0]));
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
