#include "Problem.hpp"
#include "aux.hpp"

using namespace std;

void Problem::init() {

    this->t_end = 1.;
    this->nf = 1;
    this->d = 1;
    this->x0 = vector<double>(d,1.2);
    this->lambdas = {1.};
    this->betas   = {1.};
}

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
    result[0] = -sin(x[0])*cos(x[0])*(1 + exp(-1.0));
    return result;
}

// ! Coefficient 2? Seems to have been forgotten in exact solution
vector< vector<double> > Problem::soldiff(vector<double> x) {
    vector< vector<double> > result(this->d,vector<double>(this->d,0.));
    result[0][0] = sqrt(2*cos(x[0])*cos(x[0])*(1 + exp(-1.0)));
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
    result[0] = -y[0];
    result[1] = -y[1];
    return result;
}

vector<double> Problem::diff(vector<double> x, vector<double> y) {
    vector<double> result(2*nf,0.);
    result[0] = 1.0;
    return result;
}

vector<double> Problem::fast_drift_h(vector<double> x, vector<double> y) {
    vector<double> result(this->nf);
    result[0] = 0.;
    return result;
}
