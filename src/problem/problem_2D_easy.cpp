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
    double result = (y[0]*(0.4*pow(y[0], 2) + 0.4*pow(y[1], 2) + 2)*cos(y[0])*cos(y[1]) + y[1]*(0.4*pow(y[0], 2) + 0.4*pow(y[1], 2) + 2)*cos(y[0] + y[1]) + sin(y[0])*cos(y[1]) + sin(y[0] + y[1]))*cos(x[0]);
    return result; 
}

double potential_n(vector<double> x, vector<double> y){
    double result = pow(y[0], 2) + pow(y[1], 2) + 0.1*pow(pow(y[0], 2) + pow(y[1], 2), 2);
    return result; 
}

double linearTerm_n(vector<double> x, vector<double> y){
    double result = -0.125*pow(4, 2)*pow(y[0], 2)*pow(0.4*pow(y[0], 2) + 0.4*pow(y[1], 2) + 2, 2) - 0.125*pow(4, 2)*pow(y[1], 2)*pow(0.4*pow(y[0], 2) + 0.4*pow(y[1], 2) + 2, 2) + pow(4, 2)*(0.1*pow(y[0], 2) + 0.3*pow(y[1], 2) + 0.5) + pow(4, 2)*(0.3*pow(y[0], 2) + 0.1*pow(y[1], 2) + 0.5);
    return result; 
}

double zrho_n(vector<double> x, vector<double> y){
    double result = exp(-pow(y[0], 2) - pow(y[1], 2) - 0.1*pow(pow(y[0], 2) + pow(y[1], 2), 2));
    return result; 
}

double phi0(vector<double> x, vector<double> y){
    double result = cos(x[0] + y[0] + y[1]);
    return result; 
}

double a0(vector<double> x, vector<double> y){
    double result = pow(4, 2)*(-0.2*pow(y[0], 3)*sin(x[0] + y[0] + y[1]) - 0.2*pow(y[0], 2)*y[1]*sin(x[0] + y[0] + y[1]) - 0.2*y[0]*pow(y[1], 2)*sin(x[0] + y[0] + y[1]) - 1.0*y[0]*sin(x[0] + y[0] + y[1]) - 0.2*pow(y[1], 3)*sin(x[0] + y[0] + y[1]) - 1.0*y[1]*sin(x[0] + y[0] + y[1]) + 1.0*cos(x[0] + y[0] + y[1]));
    return result; 
}

double dxphi00(vector<double> x, vector<double> y){
    double result = -sin(x[0] + y[0] + y[1]);
    return result; 
}

double dxa00(vector<double> x, vector<double> y){
    double result = pow(4, 2)*(-0.2*pow(y[0], 3)*cos(x[0] + y[0] + y[1]) - 0.2*pow(y[0], 2)*y[1]*cos(x[0] + y[0] + y[1]) - 0.2*y[0]*pow(y[1], 2)*cos(x[0] + y[0] + y[1]) - 1.0*y[0]*cos(x[0] + y[0] + y[1]) - 0.2*pow(y[1], 3)*cos(x[0] + y[0] + y[1]) - 1.0*y[1]*cos(x[0] + y[0] + y[1]) - 1.0*sin(x[0] + y[0] + y[1]));
    return result; 
}

double dxphi01(vector<double> x, vector<double> y){
    double result = 0;
    return result; 
}

double dxa01(vector<double> x, vector<double> y){
    double result = 0;
    return result; 
}

double dya00(vector<double> x, vector<double> y){
    double result = pow(4, 2)*(-0.2*pow(y[0], 3)*cos(x[0] + y[0] + y[1]) - 0.2*pow(y[0], 2)*y[1]*cos(x[0] + y[0] + y[1]) - 0.6*pow(y[0], 2)*sin(x[0] + y[0] + y[1]) - 0.2*y[0]*pow(y[1], 2)*cos(x[0] + y[0] + y[1]) - 0.4*y[0]*y[1]*sin(x[0] + y[0] + y[1]) - 1.0*y[0]*cos(x[0] + y[0] + y[1]) - 0.2*pow(y[1], 3)*cos(x[0] + y[0] + y[1]) - 0.2*pow(y[1], 2)*sin(x[0] + y[0] + y[1]) - 1.0*y[1]*cos(x[0] + y[0] + y[1]) - 2.0*sin(x[0] + y[0] + y[1]));
    return result; 
}

double dya01(vector<double> x, vector<double> y){
    double result = pow(4, 2)*(-0.2*pow(y[0], 3)*cos(x[0] + y[0] + y[1]) - 0.2*pow(y[0], 2)*y[1]*cos(x[0] + y[0] + y[1]) - 0.2*pow(y[0], 2)*sin(x[0] + y[0] + y[1]) - 0.2*y[0]*pow(y[1], 2)*cos(x[0] + y[0] + y[1]) - 0.4*y[0]*y[1]*sin(x[0] + y[0] + y[1]) - 1.0*y[0]*cos(x[0] + y[0] + y[1]) - 0.2*pow(y[1], 3)*cos(x[0] + y[0] + y[1]) - 0.6*pow(y[1], 2)*sin(x[0] + y[0] + y[1]) - 1.0*y[1]*cos(x[0] + y[0] + y[1]) - 2.0*sin(x[0] + y[0] + y[1]));
    return result; 
}

double phi1(vector<double> x, vector<double> y){
    double result = sin(x[1])*sin(y[0] + y[1]);
    return result; 
}

double a1(vector<double> x, vector<double> y){
    double result = pow(4, 2)*(0.2*pow(y[0], 3)*cos(y[0] + y[1]) + 0.2*pow(y[0], 2)*y[1]*cos(y[0] + y[1]) + 0.2*y[0]*pow(y[1], 2)*cos(y[0] + y[1]) + 1.0*y[0]*cos(y[0] + y[1]) + 0.2*pow(y[1], 3)*cos(y[0] + y[1]) + 1.0*y[1]*cos(y[0] + y[1]) + 1.0*sin(y[0] + y[1]))*sin(x[1]);
    return result; 
}

double dxphi10(vector<double> x, vector<double> y){
    double result = 0;
    return result; 
}

double dxa10(vector<double> x, vector<double> y){
    double result = 0;
    return result; 
}

double dxphi11(vector<double> x, vector<double> y){
    double result = sin(y[0] + y[1])*cos(x[1]);
    return result; 
}

double dxa11(vector<double> x, vector<double> y){
    double result = pow(4, 2)*(0.2*pow(y[0], 3)*cos(y[0] + y[1]) + 0.2*pow(y[0], 2)*y[1]*cos(y[0] + y[1]) + 0.2*y[0]*pow(y[1], 2)*cos(y[0] + y[1]) + 1.0*y[0]*cos(y[0] + y[1]) + 0.2*pow(y[1], 3)*cos(y[0] + y[1]) + 1.0*y[1]*cos(y[0] + y[1]) + 1.0*sin(y[0] + y[1]))*cos(x[1]);
    return result; 
}

double dya10(vector<double> x, vector<double> y){
    double result = pow(4, 2)*(-0.2*pow(y[0], 3)*sin(y[0] + y[1]) - 0.2*pow(y[0], 2)*y[1]*sin(y[0] + y[1]) + 0.6*pow(y[0], 2)*cos(y[0] + y[1]) - 0.2*y[0]*pow(y[1], 2)*sin(y[0] + y[1]) + 0.4*y[0]*y[1]*cos(y[0] + y[1]) - 1.0*y[0]*sin(y[0] + y[1]) - 0.2*pow(y[1], 3)*sin(y[0] + y[1]) + 0.2*pow(y[1], 2)*cos(y[0] + y[1]) - 1.0*y[1]*sin(y[0] + y[1]) + 2.0*cos(y[0] + y[1]))*sin(x[1]);
    return result; 
}

double dya11(vector<double> x, vector<double> y){
    double result = pow(4, 2)*(-0.2*pow(y[0], 3)*sin(y[0] + y[1]) - 0.2*pow(y[0], 2)*y[1]*sin(y[0] + y[1]) + 0.2*pow(y[0], 2)*cos(y[0] + y[1]) - 0.2*y[0]*pow(y[1], 2)*sin(y[0] + y[1]) + 0.4*y[0]*y[1]*cos(y[0] + y[1]) - 1.0*y[0]*sin(y[0] + y[1]) - 0.2*pow(y[1], 3)*sin(y[0] + y[1]) + 0.6*pow(y[1], 2)*cos(y[0] + y[1]) - 1.0*y[1]*sin(y[0] + y[1]) + 2.0*cos(y[0] + y[1]))*sin(x[1]);
    return result; 
}

double dyv0(vector<double> x, vector<double> y){
    double result = 0.4*y[0]*(pow(y[0], 2) + pow(y[1], 2)) + 2*y[0];
    return result; 
}

double h0(vector<double> x, vector<double> y){
    double result = cos(x[0])*cos(y[0])*cos(y[1]);
    return result; 
}

double dyv1(vector<double> x, vector<double> y){
    double result = 0.4*y[1]*(pow(y[0], 2) + pow(y[1], 2)) + 2*y[1];
    return result; 
}

double h1(vector<double> x, vector<double> y){
    double result = cos(x[0])*cos(y[0] + y[1]);
    return result; 
}

double drif0(vector<double> x, vector<double> y){
    double result = -0.4*y[0]*(pow(y[0], 2) + pow(y[1], 2)) - 2*y[0];
    return result; 
}

double diff0(vector<double> x, vector<double> y){
    double result = 4;
    return result; 
}

double drif1(vector<double> x, vector<double> y){
    double result = -0.4*y[1]*(pow(y[0], 2) + pow(y[1], 2)) - 2*y[1];
    return result; 
}

double diff1(vector<double> x, vector<double> y){
    double result = 4;
    return result; 
}

double drif2(vector<double> x, vector<double> y){
    double result = -0.8*y[0]*y[1]*y[3] - y[2]*(1.2*pow(y[0], 2) + 0.4*pow(y[1], 2) + 2) + cos(x[0])*cos(y[0])*cos(y[1]);
    return result; 
}

double diff2(vector<double> x, vector<double> y){
    double result = 4;
    return result; 
}

double drif3(vector<double> x, vector<double> y){
    double result = -0.8*y[0]*y[1]*y[2] - y[3]*(0.4*pow(y[0], 2) + 1.2*pow(y[1], 2) + 2) + cos(x[0])*cos(y[0] + y[1]);
    return result; 
}

double diff3(vector<double> x, vector<double> y){
    double result = 4;
    return result; 
}

void Problem::init_functions() {

    ns = 2;
    nf = 2;
    s = 4;

    stardiv_h = stardiv_h_n;
    zrho = zrho_n;
    linearTerm = linearTerm_n;
    potential = potential_n;
    dyv = {dyv0, dyv1};
    h = {h0, h1};
    a = {a0, a1};
    dxa = {{dxa00, dxa01}, {dxa10, dxa11}};
    dya = {{dya00, dya01}, {dya10, dya11}};
    phi = {phi0, phi1};
    dxphi = {{dxphi00, dxphi01}, {dxphi10, dxphi11}};
    drif = {drif0, drif1, drif2, drif3};
    diff = {diff0, diff1, diff2, diff3};
}
