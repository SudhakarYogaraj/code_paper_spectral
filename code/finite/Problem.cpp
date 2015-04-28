#include "Problem.hpp"

using namespace std;

void Problem::init() {
    // Final time
    this->t_end = 1.;

    // Number of slow processes
    this->d = 1;

    // Number of fast processes
    this->nf = 2.;

    // Initial condition
    this->x0 = vector<double>(d,0.2);

    this-> lambdas = {1., 2.};
    this-> betas   = {1., 1.};
}

vector<double> Problem::soldrif(vector<double> x) {
    vector<double> result(this->d,0.);
    result[0] = x[0]*(x[0]-2)/4. + (x[0]-2)/24.;
    return result;
}

vector< vector<double> > Problem::soldiff(vector<double> x) {
    vector< vector<double> > result(this->d,vector<double>(this->d,0.));
    result[0][0] = pow(x[0]-2,2)/12.;
    result = cholesky(result);
    return result;
}


vector<double> Problem::a(vector<double> x, vector<double> y) {
    vector<double> result(this->d,0.);
    result[0] = ((x[0]-2)*y[0]*y[1]);
    return result;
}

vector< vector<double> > Problem::dax(vector<double> x, vector<double> y) {
    vector< vector<double> > result(this->d,vector<double>(this->d,0.));
    result[0][0] = (y[0]*y[1]);
    return result;
}

vector< vector<double> > Problem::day(vector<double> x, vector<double> y) {
    vector< vector<double> > result(this->d,vector<double>(this->nf,0.));
    result[0][0] = (x[0]-2)*y[1];
    result[0][1] = (x[0]-2)*y[0];
    return result;
}

vector<double> Problem::drif(vector<double> x, vector<double> y) {
    vector<double> result(2*this->nf);
    result[0] = -y[0];
    result[1] = -2*y[1];
    result[2] = -y[2] + x[0]*y[1];
    result[3] = -2*y[3] + x[0]*y[0];
    return result;
}

vector<double> Problem::diff(vector<double> x, vector<double> y) {
    vector<double> result(2*nf);
    result[0] = 1.;
    result[1] = 1.;
    result[2] = 0.;
    result[3] = 0.;
    return result;
}

vector<double> Problem::fast_drift_h(vector<double> x, vector<double> y) {
    vector<double> result(this->nf);
    result[0] = x[0]*y[1];
    result[1] = x[0]*y[0];
    return result;
}
