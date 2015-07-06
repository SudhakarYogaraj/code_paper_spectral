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
        std::vector<double> sigmas;

        // Potential for gradient structure

        // FOR SPECTRAL
        // Non-leading order part of drift in the fast process
        std::vector<double> bias;
        std::vector<double> eig_val_cov;
        std::vector< std::vector<double> > eig_vec_cov;
        std::vector< std::vector<double> > sqrt_cov;
        std::vector< std::vector<double> > covariance;
        std::vector< std::vector<double> > inv_cov;
        double det_sqrt_cov;

        double rho(std::vector<double> x, std::vector<double> y);
        double normalization;

        void update_stats(std::vector<double> x);
        std::vector<double> rescale(std::vector<double> y);

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

        std::vector<double> soldrif(std::vector<double> x);
        std::vector< std::vector<double> > soldiff(std::vector<double> x);
};
#endif
