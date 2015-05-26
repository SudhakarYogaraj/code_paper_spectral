#ifndef PROBLEM_H
#define PROBLEM_H

#define PI 3.141592653589793238462643383279502884

#include <math.h>
#include <vector>
#include <iostream>

class Problem {
    public:

        // Final time
        double t_end;

        // Number of fast processes
        int nf;

        // Number of slow processes
        int ns;

        // Initial condition for the slow process
        std::vector<double> x0;

        // Potential for gradient structure
        std::vector<double> phi(std::vector<double> x, std::vector<double> y);
        std::vector< std::vector<double> > phi_x(std::vector<double> x, std::vector<double> y);
        double potential(std::vector<double> x, std::vector<double> y);
        std::vector<double> grad(std::vector<double> x, std::vector<double> y);
        double linearTerm(std::vector<double> x, std::vector<double> y);
        double rho(std::vector<double> x, std::vector<double> y);

        // Drift coefficient of the slow process
        std::vector<double> a(std::vector<double> x, std::vector<double> y); // f0 in Pavliotis-Stuart

        // Derivatives of a(,)
        std::vector< std::vector<double> > dax(std::vector<double> x, std::vector<double> y);
        std::vector< std::vector<double> > day(std::vector<double> x, std::vector<double> y);

        // FOR SPECTRAL
        // Non-leading order part of drift in the fast process
        std::vector<double>  fast_drift_h(std::vector<double> x, std::vector<double> y);
        std::vector< std::vector<double> > grad_h(std::vector<double> x, std::vector<double> y);
        double stardiv_h(std::vector<double> x, std::vector<double> y);
        std::vector<double> lambdas;
        std::vector<double> betas;
        std::vector<double> sigmas;

        // FOR HMM
        // Drift and diffusion terms of the fast system, in its transformed
        // version. The first nf components correspond to the initial
        // variables, whereas the last nf components correspend to the
        // auxiliary variables.
        std::vector<double> drif(std::vector<double> x, std::vector<double> y);
        std::vector<double> diff(std::vector<double> x, std::vector<double> y);

        // Drift and diffusion coefficients of the solution
        std::vector<double> soldrif(std::vector<double> x);
        std::vector< std::vector<double> > soldiff(std::vector<double> x);

        // Initialization of the problem
        void init();
};
#endif
