#include "Gaussian_integrator.hpp"

using namespace std;

double Gaussian_integrator::quad1d(function<double (double)> f, double sigma) {
    double result = 0.;
    double x_aux;
    for (unsigned int i = 0; i < nodes.size(); ++i) {
        x_aux = sqrt(2)*sigma*this->nodes[i];
        result += this->weights[i]*f(x_aux);
    }
    return result/sqrt(PI);
     
}

double Gaussian_integrator::quadnd(function<double (vector<double>)> f, vector<double> sigmas) {
    double result = 0.;
    int size = sigmas.size();
    if (size == 1) {
        auto lambda = [f,this](double x) -> double {
            vector<double> aux_vec(1,x);
            return f(aux_vec);
        };
        return this->quad1d(lambda, sigmas[0]);
    }
    else {
        auto lambda = [f,size,sigmas,this](vector<double> x) -> double {
            auto nested_lambda = [f,&x](double y) -> double {
                x.push_back(y);
                double result = f(x);
                x.pop_back();
                return result;
            };
            return this->quad1d(nested_lambda, sigmas[size-1]);
        };
        sigmas.pop_back();
        return this->quadnd(lambda, sigmas);
    }
    return result;
}
