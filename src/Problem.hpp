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

        // FOR HMM
        // Drift and diffusion terms of the fast system, in its transformed
        // version. The first nf components correspond to the initial
        // variables, whereas the last nf components correspend to the
        // auxiliary variables.

        // Initialization of the problem
        void init();
        void init_functions();

        double (*potential) (std::vector<double> x, std::vector<double> y);
        double (*linearTerm) (std::vector<double> x, std::vector<double> y);
        double (*zrho) (std::vector<double> x, std::vector<double> y);
        double (*stardiv_h) (std::vector<double> x, std::vector<double> y);
        std::vector<double (*) (std::vector<double> x, std::vector<double> y)> dyv;
        std::vector<double (*) (std::vector<double> x, std::vector<double> y)> h;
        std::vector<double (*) (std::vector<double> x, std::vector<double> y)> a;
        std::vector< std::vector<double (*) (std::vector<double> x, std::vector<double> y)> > dxa;
        std::vector< std::vector<double (*) (std::vector<double> x, std::vector<double> y)> > dya;
        std::vector<double (*) (std::vector<double> x, std::vector<double> y)> phi;
        std::vector< std::vector<double (*) (std::vector<double> x, std::vector<double> y)> > dxphi;
        std::vector<double (*) (std::vector<double> x, std::vector<double> y)> drif;
        std::vector<double (*) (std::vector<double> x, std::vector<double> y)> diff;
};
#endif
