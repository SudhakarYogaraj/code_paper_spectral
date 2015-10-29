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
    double result = (1.0*y[0]*(pow(y[0], 2) - 1)*cos(y[0]) + sin(y[0]))*cos(x[0]);
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
    double result = x[0]*cos(5*y[0]) + sin((1.0L/2.0L)*x[0]*y[0]);
    return result; 
}

double a0(vector<double> x, vector<double> y){
    double result = x[0]*(0.25*x[0]*sin((1.0L/2.0L)*x[0]*y[0]) - 5.0*pow(y[0], 3)*sin(5*y[0]) + 0.5*pow(y[0], 3)*cos((1.0L/2.0L)*x[0]*y[0]) + 5.0*y[0]*sin(5*y[0]) - 0.5*y[0]*cos((1.0L/2.0L)*x[0]*y[0]) + 25.0*cos(5*y[0]));
    return result; 
}

double dxphi00(vector<double> x, vector<double> y){
    double result = (1.0L/2.0L)*y[0]*cos((1.0L/2.0L)*x[0]*y[0]) + cos(5*y[0]);
    return result; 
}

double dxa00(vector<double> x, vector<double> y){
    double result = x[0]*(0.125*x[0]*y[0]*cos((1.0L/2.0L)*x[0]*y[0]) - 0.25*pow(y[0], 4)*sin((1.0L/2.0L)*x[0]*y[0]) + 0.25*pow(y[0], 2)*sin((1.0L/2.0L)*x[0]*y[0]) + 0.25*sin((1.0L/2.0L)*x[0]*y[0])) + 0.25*x[0]*sin((1.0L/2.0L)*x[0]*y[0]) - 5.0*pow(y[0], 3)*sin(5*y[0]) + 0.5*pow(y[0], 3)*cos((1.0L/2.0L)*x[0]*y[0]) + 5.0*y[0]*sin(5*y[0]) - 0.5*y[0]*cos((1.0L/2.0L)*x[0]*y[0]) + 25.0*cos(5*y[0]);
    return result; 
}

double dya00(vector<double> x, vector<double> y){
    double result = x[0]*(0.125*pow(x[0], 2)*cos((1.0L/2.0L)*x[0]*y[0]) - 0.25*x[0]*pow(y[0], 3)*sin((1.0L/2.0L)*x[0]*y[0]) + 0.25*x[0]*y[0]*sin((1.0L/2.0L)*x[0]*y[0]) - 25.0*pow(y[0], 3)*cos(5*y[0]) - 15.0*pow(y[0], 2)*sin(5*y[0]) + 1.5*pow(y[0], 2)*cos((1.0L/2.0L)*x[0]*y[0]) + 25.0*y[0]*cos(5*y[0]) - 120.0*sin(5*y[0]) - 0.5*cos((1.0L/2.0L)*x[0]*y[0]));
    return result; 
}

double dyv0(vector<double> x, vector<double> y){
    double result = 1.0*pow(y[0], 3) - 1.0*y[0];
    return result; 
}

double h0(vector<double> x, vector<double> y){
    double result = cos(x[0])*cos(y[0]);
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
    double result = -y[1]*(3.0*pow(y[0], 2) - 1.0) + cos(x[0])*cos(y[0]);
    return result; 
}

double diff1(vector<double> x, vector<double> y){
    double result = sqrt(2);
    return result; 
}

void Problem::init_functions() {

    ns = 1;
    nf = 1;

    stardiv_h = stardiv_h_n;
    zrho = zrho_n;
    linearTerm = linearTerm_n;
    potential = potential_n;
    dyv = {dyv0};
    h = {h0};
    a = {a0};
    dxa = {{dxa00}};
    dya = {{dya00}};
    phi = {phi0};
    dxphi = {{dxphi00}};
    drif = {drif0, drif1};
    diff = {diff0, diff1};
}