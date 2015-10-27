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
    double result = (1.0L/2.0L)*(6*pow(y[0], 5)*cos(y[0] - y[1]) + 6*pow(y[0], 5)*cos(y[0] + y[1]) + 12*pow(y[0], 4)*y[1]*cos(y[0] + y[1]) + 12*pow(y[0], 3)*pow(y[1], 2)*cos(y[0] - y[1]) + 12*pow(y[0], 3)*pow(y[1], 2)*cos(y[0] + y[1]) + 24*pow(y[0], 2)*pow(y[1], 3)*cos(y[0] + y[1]) - 6*pow(y[0], 2)*cos(y[0] - y[1]) - 6*pow(y[0], 2)*cos(y[0] + y[1]) + 6*y[0]*pow(y[1], 4)*cos(y[0] - y[1]) + 6*y[0]*pow(y[1], 4)*cos(y[0] + y[1]) + 24*y[0]*y[1]*cos(y[0] + y[1]) + 12*pow(y[1], 5)*cos(y[0] + y[1]) + 6*pow(y[1], 2)*cos(y[0] - y[1]) + 6*pow(y[1], 2)*cos(y[0] + y[1]) + sin(y[0] - y[1]) + 3*sin(y[0] + y[1]))*cos(x[0]);
    return result; 
}

double potential_n(vector<double> x, vector<double> y){
    double result = (pow(y[1], 2) + pow(y[0] - 1, 2))*(pow(y[0] + 1.0L/2.0L, 2) + pow(y[1] - 1.0L/2.0L*sqrt(3), 2))*(pow(y[0] + 1.0L/2.0L, 2) + pow(y[1] + (1.0L/2.0L)*sqrt(3), 2));
    return result; 
}

double linearTerm_n(vector<double> x, vector<double> y){
    double result = -9.0*pow(y[0], 10) - 36.0*pow(y[0], 8)*pow(y[1], 2) + 18.0*pow(y[0], 7) - 54.0*pow(y[0], 6)*pow(y[1], 4) + 18.0*pow(y[0], 5)*pow(y[1], 2) - 36.0*pow(y[0], 4)*pow(y[1], 6) + 9.0*pow(y[0], 4) - 18.0*pow(y[0], 3)*pow(y[1], 4) - 9.0*pow(y[0], 2)*pow(y[1], 8) + 54.0*pow(y[0], 2)*pow(y[1], 2) - 18.0*y[0]*pow(y[1], 6) + 9.0*pow(y[1], 4) - pow(y[1], 2)*(9.0*pow(y[0], 8) + 36.0*pow(y[0], 6)*pow(y[1], 2) + 36.0*pow(y[0], 5) + 54.0*pow(y[0], 4)*pow(y[1], 4) + 72.0*pow(y[0], 3)*pow(y[1], 2) + 36.0*pow(y[0], 2)*pow(y[1], 6) + 36.0*pow(y[0], 2) + 36.0*y[0]*pow(y[1], 4) + 9.0*pow(y[1], 8));
    return result; 
}

double zrho_n(vector<double> x, vector<double> y){
    double result = exp(-(pow(y[1], 2) + pow(y[0] - 1, 2))*(pow(y[0] + 1.0L/2.0L, 2) + pow(y[1] - 1.0L/2.0L*sqrt(3), 2))*(pow(y[0] + 1.0L/2.0L, 2) + pow(y[1] + (1.0L/2.0L)*sqrt(3), 2)));
    return result; 
}

double phi0(vector<double> x, vector<double> y){
    double result = cos(x[0] + y[0] + y[1]);
    return result; 
}

double a0(vector<double> x, vector<double> y){
    double result = -6.0*pow(y[0], 5)*sin(x[0] + y[0] + y[1]) - 6.0*pow(y[0], 4)*y[1]*sin(x[0] + y[0] + y[1]) - 12.0*pow(y[0], 3)*pow(y[1], 2)*sin(x[0] + y[0] + y[1]) - 12.0*pow(y[0], 2)*pow(y[1], 3)*sin(x[0] + y[0] + y[1]) + 6.0*pow(y[0], 2)*sin(x[0] + y[0] + y[1]) - 6.0*y[0]*pow(y[1], 4)*sin(x[0] + y[0] + y[1]) - 12.0*y[0]*y[1]*sin(x[0] + y[0] + y[1]) - 6.0*pow(y[1], 5)*sin(x[0] + y[0] + y[1]) - 6.0*pow(y[1], 2)*sin(x[0] + y[0] + y[1]) + 2.0*cos(x[0] + y[0] + y[1]);
    return result; 
}

double dxphi00(vector<double> x, vector<double> y){
    double result = -sin(x[0] + y[0] + y[1]);
    return result; 
}

double dxa00(vector<double> x, vector<double> y){
    double result = -6.0*pow(y[0], 5)*cos(x[0] + y[0] + y[1]) - 6.0*pow(y[0], 4)*y[1]*cos(x[0] + y[0] + y[1]) - 12.0*pow(y[0], 3)*pow(y[1], 2)*cos(x[0] + y[0] + y[1]) - 12.0*pow(y[0], 2)*pow(y[1], 3)*cos(x[0] + y[0] + y[1]) + 6.0*pow(y[0], 2)*cos(x[0] + y[0] + y[1]) - 6.0*y[0]*pow(y[1], 4)*cos(x[0] + y[0] + y[1]) - 12.0*y[0]*y[1]*cos(x[0] + y[0] + y[1]) - 6.0*pow(y[1], 5)*cos(x[0] + y[0] + y[1]) - 6.0*pow(y[1], 2)*cos(x[0] + y[0] + y[1]) - 2.0*sin(x[0] + y[0] + y[1]);
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
    double result = -6.0*pow(y[0], 5)*cos(x[0] + y[0] + y[1]) - 6.0*pow(y[0], 4)*y[1]*cos(x[0] + y[0] + y[1]) - 30.0*pow(y[0], 4)*sin(x[0] + y[0] + y[1]) - 12.0*pow(y[0], 3)*pow(y[1], 2)*cos(x[0] + y[0] + y[1]) - 24.0*pow(y[0], 3)*y[1]*sin(x[0] + y[0] + y[1]) - 12.0*pow(y[0], 2)*pow(y[1], 3)*cos(x[0] + y[0] + y[1]) - 36.0*pow(y[0], 2)*pow(y[1], 2)*sin(x[0] + y[0] + y[1]) + 6.0*pow(y[0], 2)*cos(x[0] + y[0] + y[1]) - 6.0*y[0]*pow(y[1], 4)*cos(x[0] + y[0] + y[1]) - 24.0*y[0]*pow(y[1], 3)*sin(x[0] + y[0] + y[1]) - 12.0*y[0]*y[1]*cos(x[0] + y[0] + y[1]) + 12.0*y[0]*sin(x[0] + y[0] + y[1]) - 6.0*pow(y[1], 5)*cos(x[0] + y[0] + y[1]) - 6.0*pow(y[1], 4)*sin(x[0] + y[0] + y[1]) - 6.0*pow(y[1], 2)*cos(x[0] + y[0] + y[1]) - 12.0*y[1]*sin(x[0] + y[0] + y[1]) - 2.0*sin(x[0] + y[0] + y[1]);
    return result; 
}

double dya01(vector<double> x, vector<double> y){
    double result = -6.0*pow(y[0], 5)*cos(x[0] + y[0] + y[1]) - 6.0*pow(y[0], 4)*y[1]*cos(x[0] + y[0] + y[1]) - 6.0*pow(y[0], 4)*sin(x[0] + y[0] + y[1]) - 12.0*pow(y[0], 3)*pow(y[1], 2)*cos(x[0] + y[0] + y[1]) - 24.0*pow(y[0], 3)*y[1]*sin(x[0] + y[0] + y[1]) - 12.0*pow(y[0], 2)*pow(y[1], 3)*cos(x[0] + y[0] + y[1]) - 36.0*pow(y[0], 2)*pow(y[1], 2)*sin(x[0] + y[0] + y[1]) + 6.0*pow(y[0], 2)*cos(x[0] + y[0] + y[1]) - 6.0*y[0]*pow(y[1], 4)*cos(x[0] + y[0] + y[1]) - 24.0*y[0]*pow(y[1], 3)*sin(x[0] + y[0] + y[1]) - 12.0*y[0]*y[1]*cos(x[0] + y[0] + y[1]) - 12.0*y[0]*sin(x[0] + y[0] + y[1]) - 6.0*pow(y[1], 5)*cos(x[0] + y[0] + y[1]) - 30.0*pow(y[1], 4)*sin(x[0] + y[0] + y[1]) - 6.0*pow(y[1], 2)*cos(x[0] + y[0] + y[1]) - 12.0*y[1]*sin(x[0] + y[0] + y[1]) - 2.0*sin(x[0] + y[0] + y[1]);
    return result; 
}

double phi1(vector<double> x, vector<double> y){
    double result = sin(x[1])*sin(y[0] + y[1]);
    return result; 
}

double a1(vector<double> x, vector<double> y){
    double result = (6.0*pow(y[0], 5)*cos(y[0] + y[1]) + 6.0*pow(y[0], 4)*y[1]*cos(y[0] + y[1]) + 12.0*pow(y[0], 3)*pow(y[1], 2)*cos(y[0] + y[1]) + 12.0*pow(y[0], 2)*pow(y[1], 3)*cos(y[0] + y[1]) - 6.0*pow(y[0], 2)*cos(y[0] + y[1]) + 6.0*y[0]*pow(y[1], 4)*cos(y[0] + y[1]) + 12.0*y[0]*y[1]*cos(y[0] + y[1]) + 6.0*pow(y[1], 5)*cos(y[0] + y[1]) + 6.0*pow(y[1], 2)*cos(y[0] + y[1]) + 2.0*sin(y[0] + y[1]))*sin(x[1]);
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
    double result = (6.0*pow(y[0], 5)*cos(y[0] + y[1]) + 6.0*pow(y[0], 4)*y[1]*cos(y[0] + y[1]) + 12.0*pow(y[0], 3)*pow(y[1], 2)*cos(y[0] + y[1]) + 12.0*pow(y[0], 2)*pow(y[1], 3)*cos(y[0] + y[1]) - 6.0*pow(y[0], 2)*cos(y[0] + y[1]) + 6.0*y[0]*pow(y[1], 4)*cos(y[0] + y[1]) + 12.0*y[0]*y[1]*cos(y[0] + y[1]) + 6.0*pow(y[1], 5)*cos(y[0] + y[1]) + 6.0*pow(y[1], 2)*cos(y[0] + y[1]) + 2.0*sin(y[0] + y[1]))*cos(x[1]);
    return result; 
}

double dya10(vector<double> x, vector<double> y){
    double result = (-6.0*pow(y[0], 5)*sin(y[0] + y[1]) - 6.0*pow(y[0], 4)*y[1]*sin(y[0] + y[1]) + 30.0*pow(y[0], 4)*cos(y[0] + y[1]) - 12.0*pow(y[0], 3)*pow(y[1], 2)*sin(y[0] + y[1]) + 24.0*pow(y[0], 3)*y[1]*cos(y[0] + y[1]) - 12.0*pow(y[0], 2)*pow(y[1], 3)*sin(y[0] + y[1]) + 36.0*pow(y[0], 2)*pow(y[1], 2)*cos(y[0] + y[1]) + 6.0*pow(y[0], 2)*sin(y[0] + y[1]) - 6.0*y[0]*pow(y[1], 4)*sin(y[0] + y[1]) + 24.0*y[0]*pow(y[1], 3)*cos(y[0] + y[1]) - 12.0*y[0]*y[1]*sin(y[0] + y[1]) - 12.0*y[0]*cos(y[0] + y[1]) - 6.0*pow(y[1], 5)*sin(y[0] + y[1]) + 6.0*pow(y[1], 4)*cos(y[0] + y[1]) - 6.0*pow(y[1], 2)*sin(y[0] + y[1]) + 12.0*y[1]*cos(y[0] + y[1]) + 2.0*cos(y[0] + y[1]))*sin(x[1]);
    return result; 
}

double dya11(vector<double> x, vector<double> y){
    double result = (-6.0*pow(y[0], 5)*sin(y[0] + y[1]) - 6.0*pow(y[0], 4)*y[1]*sin(y[0] + y[1]) + 6.0*pow(y[0], 4)*cos(y[0] + y[1]) - 12.0*pow(y[0], 3)*pow(y[1], 2)*sin(y[0] + y[1]) + 24.0*pow(y[0], 3)*y[1]*cos(y[0] + y[1]) - 12.0*pow(y[0], 2)*pow(y[1], 3)*sin(y[0] + y[1]) + 36.0*pow(y[0], 2)*pow(y[1], 2)*cos(y[0] + y[1]) + 6.0*pow(y[0], 2)*sin(y[0] + y[1]) - 6.0*y[0]*pow(y[1], 4)*sin(y[0] + y[1]) + 24.0*y[0]*pow(y[1], 3)*cos(y[0] + y[1]) - 12.0*y[0]*y[1]*sin(y[0] + y[1]) + 12.0*y[0]*cos(y[0] + y[1]) - 6.0*pow(y[1], 5)*sin(y[0] + y[1]) + 30.0*pow(y[1], 4)*cos(y[0] + y[1]) - 6.0*pow(y[1], 2)*sin(y[0] + y[1]) + 12.0*y[1]*cos(y[0] + y[1]) + 2.0*cos(y[0] + y[1]))*sin(x[1]);
    return result; 
}

double dyv0(vector<double> x, vector<double> y){
    double result = (2*y[0] - 2)*(pow(y[0] + 1.0L/2.0L, 2) + pow(y[1] - 1.0L/2.0L*sqrt(3), 2))*(pow(y[0] + 1.0L/2.0L, 2) + pow(y[1] + (1.0L/2.0L)*sqrt(3), 2)) + (2*y[0] + 1)*(pow(y[1], 2) + pow(y[0] - 1, 2))*(pow(y[0] + 1.0L/2.0L, 2) + pow(y[1] - 1.0L/2.0L*sqrt(3), 2)) + (2*y[0] + 1)*(pow(y[1], 2) + pow(y[0] - 1, 2))*(pow(y[0] + 1.0L/2.0L, 2) + pow(y[1] + (1.0L/2.0L)*sqrt(3), 2));
    return result; 
}

double h0(vector<double> x, vector<double> y){
    double result = cos(x[0])*cos(y[0])*cos(y[1]);
    return result; 
}

double dyv1(vector<double> x, vector<double> y){
    double result = 2*y[1]*(pow(y[0] + 1.0L/2.0L, 2) + pow(y[1] - 1.0L/2.0L*sqrt(3), 2))*(pow(y[0] + 1.0L/2.0L, 2) + pow(y[1] + (1.0L/2.0L)*sqrt(3), 2)) + (2*y[1] - sqrt(3))*(pow(y[1], 2) + pow(y[0] - 1, 2))*(pow(y[0] + 1.0L/2.0L, 2) + pow(y[1] + (1.0L/2.0L)*sqrt(3), 2)) + (2*y[1] + sqrt(3))*(pow(y[1], 2) + pow(y[0] - 1, 2))*(pow(y[0] + 1.0L/2.0L, 2) + pow(y[1] - 1.0L/2.0L*sqrt(3), 2));
    return result; 
}

double h1(vector<double> x, vector<double> y){
    double result = cos(x[0])*cos(y[0] + y[1]);
    return result; 
}

double drif0(vector<double> x, vector<double> y){
    double result = -(2*y[0] - 2)*(pow(y[0] + 1.0L/2.0L, 2) + pow(y[1] - 1.0L/2.0L*sqrt(3), 2))*(pow(y[0] + 1.0L/2.0L, 2) + pow(y[1] + (1.0L/2.0L)*sqrt(3), 2)) - (2*y[0] + 1)*(pow(y[1], 2) + pow(y[0] - 1, 2))*(pow(y[0] + 1.0L/2.0L, 2) + pow(y[1] - 1.0L/2.0L*sqrt(3), 2)) - (2*y[0] + 1)*(pow(y[1], 2) + pow(y[0] - 1, 2))*(pow(y[0] + 1.0L/2.0L, 2) + pow(y[1] + (1.0L/2.0L)*sqrt(3), 2));
    return result; 
}

double diff0(vector<double> x, vector<double> y){
    double result = sqrt(2);
    return result; 
}

double drif1(vector<double> x, vector<double> y){
    double result = -2*y[1]*(pow(y[0] + 1.0L/2.0L, 2) + pow(y[1] - 1.0L/2.0L*sqrt(3), 2))*(pow(y[0] + 1.0L/2.0L, 2) + pow(y[1] + (1.0L/2.0L)*sqrt(3), 2)) - (2*y[1] - sqrt(3))*(pow(y[1], 2) + pow(y[0] - 1, 2))*(pow(y[0] + 1.0L/2.0L, 2) + pow(y[1] + (1.0L/2.0L)*sqrt(3), 2)) - (2*y[1] + sqrt(3))*(pow(y[1], 2) + pow(y[0] - 1, 2))*(pow(y[0] + 1.0L/2.0L, 2) + pow(y[1] - 1.0L/2.0L*sqrt(3), 2));
    return result; 
}

double diff1(vector<double> x, vector<double> y){
    double result = sqrt(2);
    return result; 
}

double drif2(vector<double> x, vector<double> y){
    double result = -y[2]*(2*(2*y[0] - 2)*(2*y[0] + 1)*(pow(y[0] + 1.0L/2.0L, 2) + pow(y[1] - 1.0L/2.0L*sqrt(3), 2)) + 2*(2*y[0] - 2)*(2*y[0] + 1)*(pow(y[0] + 1.0L/2.0L, 2) + pow(y[1] + (1.0L/2.0L)*sqrt(3), 2)) + 2*pow(2*y[0] + 1, 2)*(pow(y[1], 2) + pow(y[0] - 1, 2)) + 2*(pow(y[1], 2) + pow(y[0] - 1, 2))*(pow(y[0] + 1.0L/2.0L, 2) + pow(y[1] - 1.0L/2.0L*sqrt(3), 2)) + 2*(pow(y[1], 2) + pow(y[0] - 1, 2))*(pow(y[0] + 1.0L/2.0L, 2) + pow(y[1] + (1.0L/2.0L)*sqrt(3), 2)) + 2*(pow(y[0] + 1.0L/2.0L, 2) + pow(y[1] - 1.0L/2.0L*sqrt(3), 2))*(pow(y[0] + 1.0L/2.0L, 2) + pow(y[1] + (1.0L/2.0L)*sqrt(3), 2))) - y[3]*(2*y[1]*(2*y[0] + 1)*(pow(y[0] + 1.0L/2.0L, 2) + pow(y[1] - 1.0L/2.0L*sqrt(3), 2)) + 2*y[1]*(2*y[0] + 1)*(pow(y[0] + 1.0L/2.0L, 2) + pow(y[1] + (1.0L/2.0L)*sqrt(3), 2)) + (2*y[0] - 2)*(2*y[1] - sqrt(3))*(pow(y[0] + 1.0L/2.0L, 2) + pow(y[1] + (1.0L/2.0L)*sqrt(3), 2)) + (2*y[0] - 2)*(2*y[1] + sqrt(3))*(pow(y[0] + 1.0L/2.0L, 2) + pow(y[1] - 1.0L/2.0L*sqrt(3), 2)) + (2*y[0] + 1)*(2*y[1] - sqrt(3))*(pow(y[1], 2) + pow(y[0] - 1, 2)) + (2*y[0] + 1)*(2*y[1] + sqrt(3))*(pow(y[1], 2) + pow(y[0] - 1, 2))) + cos(x[0])*cos(y[0])*cos(y[1]);
    return result; 
}

double diff2(vector<double> x, vector<double> y){
    double result = sqrt(2);
    return result; 
}

double drif3(vector<double> x, vector<double> y){
    double result = -y[2]*(2*y[1]*(2*y[0] + 1)*(pow(y[0] + 1.0L/2.0L, 2) + pow(y[1] - 1.0L/2.0L*sqrt(3), 2)) + 2*y[1]*(2*y[0] + 1)*(pow(y[0] + 1.0L/2.0L, 2) + pow(y[1] + (1.0L/2.0L)*sqrt(3), 2)) + (2*y[0] - 2)*(2*y[1] - sqrt(3))*(pow(y[0] + 1.0L/2.0L, 2) + pow(y[1] + (1.0L/2.0L)*sqrt(3), 2)) + (2*y[0] - 2)*(2*y[1] + sqrt(3))*(pow(y[0] + 1.0L/2.0L, 2) + pow(y[1] - 1.0L/2.0L*sqrt(3), 2)) + (2*y[0] + 1)*(2*y[1] - sqrt(3))*(pow(y[1], 2) + pow(y[0] - 1, 2)) + (2*y[0] + 1)*(2*y[1] + sqrt(3))*(pow(y[1], 2) + pow(y[0] - 1, 2))) - y[3]*(4*y[1]*(2*y[1] - sqrt(3))*(pow(y[0] + 1.0L/2.0L, 2) + pow(y[1] + (1.0L/2.0L)*sqrt(3), 2)) + 4*y[1]*(2*y[1] + sqrt(3))*(pow(y[0] + 1.0L/2.0L, 2) + pow(y[1] - 1.0L/2.0L*sqrt(3), 2)) + 2*(2*y[1] - sqrt(3))*(2*y[1] + sqrt(3))*(pow(y[1], 2) + pow(y[0] - 1, 2)) + 2*(pow(y[1], 2) + pow(y[0] - 1, 2))*(pow(y[0] + 1.0L/2.0L, 2) + pow(y[1] - 1.0L/2.0L*sqrt(3), 2)) + 2*(pow(y[1], 2) + pow(y[0] - 1, 2))*(pow(y[0] + 1.0L/2.0L, 2) + pow(y[1] + (1.0L/2.0L)*sqrt(3), 2)) + 2*(pow(y[0] + 1.0L/2.0L, 2) + pow(y[1] - 1.0L/2.0L*sqrt(3), 2))*(pow(y[0] + 1.0L/2.0L, 2) + pow(y[1] + (1.0L/2.0L)*sqrt(3), 2))) + cos(x[0])*cos(y[0] + y[1]);
    return result; 
}

double diff3(vector<double> x, vector<double> y){
    double result = sqrt(2);
    return result; 
}

void Problem::init_functions() {

    ns = 2;
    nf = 2;

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
    dyv[1] = dyv1;

    h[0] = h0;
    h[1] = h1;

    a[0] = a0;
    a[1] = a1;

    dxa[0][0] = dxa00;
    dxa[0][1] = dxa01;
    dxa[1][0] = dxa10;
    dxa[1][1] = dxa11;

    dya[0][0] = dya00;
    dya[0][1] = dya01;
    dya[1][0] = dya10;
    dya[1][1] = dya11;

    phi[0] = phi0;
    phi[1] = phi1;

    dxphi[0][0] = dxphi00;
    dxphi[0][1] = dxphi01;
    dxphi[1][0] = dxphi10;
    dxphi[1][1] = dxphi11;

    drif[0] = drif0;
    drif[1] = drif1;
    drif[2] = drif2;
    drif[3] = drif3;

    diff[0] = diff0;
    diff[1] = diff1;
    diff[2] = diff2;
    diff[3] = diff3;

}