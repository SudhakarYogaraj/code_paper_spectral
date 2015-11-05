#ifndef PROBLEM_H

#define PROBLEM_H
#define PI 3.141592653589793238462643383279502884

#include <math.h>
#include <vector>
#include <iostream>

typedef std::vector<double> std_vec;

class Problem {
    public:

        // Final time
        double t_end;

        // Number of fast processes
        int nf;

        // Number of slow processes
        int ns;

        // Initial condition for the slow process
        std_vec x0;

        // FOR HMM
        // Drift and diffusion terms of the fast system, in its transformed
        // version. The first nf components correspond to the initial
        // variables, whereas the last nf components correspend to the
        // auxiliary variables.

        // Initialization of the problem
        void init();
        void init_functions();


        double (*potential) (std_vec x, std_vec y);
        double (*linearTerm) (std_vec x, std_vec y);
        double (*zrho) (std_vec x, std_vec y);
        double (*stardiv_h) (std_vec x, std_vec y);
        std::vector<double (*) (std_vec x, std_vec y)> dyv;
        std::vector<double (*) (std_vec x, std_vec y)> h;
        std::vector<double (*) (std_vec x, std_vec y)> a;
        std::vector< std::vector<double (*) (std_vec x, std_vec y)> > dxa;
        std::vector< std::vector<double (*) (std_vec x, std_vec y)> > dya;
        std::vector<double (*) (std_vec x, std_vec y)> phi;
        std::vector< std::vector<double (*) (std_vec x, std_vec y)> > dxphi;
        std::vector<double (*) (std_vec x, std_vec y)> drif;
        std::vector<double (*) (std_vec x, std_vec y)> diff;
};
#endif
