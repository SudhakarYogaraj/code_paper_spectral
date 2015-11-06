#ifndef SOLVER_HMM_H
#define SOLVER_HMM_H

#include <random>
#include <vector>

#include "problem/Problem.hpp"
#include "solvers/Solver.hpp"

struct config_hmm {

    // Macro and micro time-steps
    double macro_dt;
    double micro_dt;

    // Order of the micro solver
    double l;

    // Parameters of the estimator
    int n;
    int nt;
    int np;
    int M;

    // Precision parameter of the solver
    double p;
};


class Solver_hmm : public Solver {

    public:
        Solver_hmm(Problem*, config_hmm*);
        SDE_coeffs estimator(std::vector<double> x, double t);
        static config_hmm sensible_conf(int p, int M);

    private:
        Problem *problem;
        config_hmm *conf;
};
#endif
