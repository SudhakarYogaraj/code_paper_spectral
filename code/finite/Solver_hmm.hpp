#ifndef SOLVER_HMM_H
#define SOLVER_HMM_H

#include <random>
#include <vector>
#include "Problem.hpp"
#include "aux.hpp"

class Solver_hmm {
    public:

        // Macro and micro time-steps
        double macro_dt;
        double micro_dt;

        // Parameters of the estimator
        int n;
        int nt;
        int np;
        int M;

        // Order of accuracy of the micro-solver
        double l;

        // Precision parameter of the solver
        double p;

        void estimator(Problem&, std::vector<double> x, std::vector<double>& yInit, std::vector<double>& fi, std::vector< std::vector<double> >& hi, int seed, double t);
        void set(double, int);
};
#endif
