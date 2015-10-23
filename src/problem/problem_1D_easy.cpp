#include "problem/Problem.hpp"
#include "templates.hpp"
#include "toolbox/Gaussian_integrator.hpp"
#include "toolbox/linear_algebra.hpp"
#include "io/io.hpp"
#include "global.hpp"

using namespace std;

void Problem::init() {
    this->init_functions();
    this->x0 = vector<double> (ns, 1.2);
    this->t_end = 1.;
}
double stardiv_h_n(vector<double> x, vector<double> y){
    double result = (1.0*y[0]*(pow(y[0], 2) - 1)*cos(0.5*y[0]) + 0.5*sin(0.5*y[0]))*cos(x[0]);
    return result; 
}

double potential_n(vector<double> x, vector<double> y){
    double result = 0.25*pow(y[0], 4) - 0.5*pow(y[0], 2);
    return result; 
}

double linearTerm_n(vector<double> x, vector<double> y){
    double result = -0.25*pow(y[0], 2)*pow(pow(y[0], 2) - 1, 2) + 1.5*pow(y[0], 2) - 0.5;
    return result; 
}

double zrho_n(vector<double> x, vector<double> y){
    double result = exp(-0.25*pow(y[0], 4) + 0.5*pow(y[0], 2));
    return result; 
}

double phi0(vector<double> x, vector<double> y){
    double result = x[0]*cos(0.5*y[0]) + sin((1.0L/2.0L)*x[0]*y[0]);
    return result; 
}

double a0(vector<double> x, vector<double> y){
    double result = x[0]*(0.25*x[0]*sin((1.0L/2.0L)*x[0]*y[0]) - 0.5*pow(y[0], 3)*sin(0.5*y[0]) + 0.5*pow(y[0], 3)*cos((1.0L/2.0L)*x[0]*y[0]) + 0.5*y[0]*sin(0.5*y[0]) - 0.5*y[0]*cos((1.0L/2.0L)*x[0]*y[0]) + 0.25*cos(0.5*y[0]));
    return result; 
}

double dxphi00(vector<double> x, vector<double> y){
    double result = (1.0L/2.0L)*y[0]*cos((1.0L/2.0L)*x[0]*y[0]) + cos(0.5*y[0]);
    return result; 
}

double dxa00(vector<double> x, vector<double> y){
    double result = x[0]*(0.125*x[0]*y[0]*cos((1.0L/2.0L)*x[0]*y[0]) - 0.25*pow(y[0], 4)*sin((1.0L/2.0L)*x[0]*y[0]) + 0.25*pow(y[0], 2)*sin((1.0L/2.0L)*x[0]*y[0]) + 0.25*sin((1.0L/2.0L)*x[0]*y[0])) + 0.25*x[0]*sin((1.0L/2.0L)*x[0]*y[0]) - 0.5*pow(y[0], 3)*sin(0.5*y[0]) + 0.5*pow(y[0], 3)*cos((1.0L/2.0L)*x[0]*y[0]) + 0.5*y[0]*sin(0.5*y[0]) - 0.5*y[0]*cos((1.0L/2.0L)*x[0]*y[0]) + 0.25*cos(0.5*y[0]);
    return result; 
}

double dya00(vector<double> x, vector<double> y){
    double result = x[0]*(0.125*pow(x[0], 2)*cos((1.0L/2.0L)*x[0]*y[0]) - 0.25*x[0]*pow(y[0], 3)*sin((1.0L/2.0L)*x[0]*y[0]) + 0.25*x[0]*y[0]*sin((1.0L/2.0L)*x[0]*y[0]) - 0.25*pow(y[0], 3)*cos(0.5*y[0]) - 1.5*pow(y[0], 2)*sin(0.5*y[0]) + 1.5*pow(y[0], 2)*cos((1.0L/2.0L)*x[0]*y[0]) + 0.25*y[0]*cos(0.5*y[0]) + 0.375*sin(0.5*y[0]) - 0.5*cos((1.0L/2.0L)*x[0]*y[0]));
    return result; 
}

double dyv0(vector<double> x, vector<double> y){
    double result = 1.0*pow(y[0], 3) - 1.0*y[0];
    return result; 
}

double h0(vector<double> x, vector<double> y){
    double result = cos(x[0])*cos(0.5*y[0]);
    return result; 
}

double drif0(vector<double> x, vector<double> y){
    double result = -1.0*pow(y[0], 3) + 1.0*y[0];
    return result; 
}

double diff0(vector<double> x, vector<double> y){
    double result = sqrt(2);
    return result; 
}

double drif1(vector<double> x, vector<double> y){
    double result = -y[1]*(3.0*pow(y[0], 2) - 1.0) + cos(x[0])*cos(0.5*y[0]);
    return result; 
}

double diff1(vector<double> x, vector<double> y){
    double result = sqrt(2);
    return result; 
}

void Problem::init_functions() {

    ns = 1;
    nf = 1;

    dyv = vector<double (*) (vector<double> x,             vector<double> y)> (nf);
    h = vector<double (*) (vector<double> x,             vector<double> y)> (nf);
    a = vector<double (*) (vector<double> x,             vector<double> y)> (ns);
    dxa = vector< vector<double (*) (vector<double> x,             vector<double> y)> >(ns, vector<double (*) (vector<double> x,             vector<double> y)> (ns));
    dya = vector< vector<double (*) (vector<double> x,             vector<double> y)> >(ns, vector<double (*) (vector<double> x,             vector<double> y)> (nf));
    phi = vector<double (*) (vector<double> x,             vector<double> y)> (ns);
    dxphi = vector< vector<double (*) (vector<double> x,             vector<double> y)> >(ns, vector<double (*) (vector<double> x,             vector<double> y)> (ns));
    drif = vector<double (*) (vector<double> x,             vector<double> y)> (2*nf);
    diff = vector<double (*) (vector<double> x,             vector<double> y)> (2*nf);

    stardiv_h = stardiv_h_n;

    zrho = zrho_n;

    linearTerm = linearTerm_n;

    potential = potential_n;

    dyv[0] = dyv0;

    h[0] = h0;

    a[0] = a0;

    dxa[0][0] = dxa00;

    dya[0][0] = dya00;

    phi[0] = phi0;

    dxphi[0][0] = dxphi00;

    drif[0] = drif0;
    drif[1] = drif1;

    diff[0] = diff0;
    diff[1] = diff1;

}