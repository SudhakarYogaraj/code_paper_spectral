#ifndef SOLVER_HMM_H
#define SOLVER_HMM_H

#include <random>
#include <vector>

#include "problem/Problem.hpp"
#include "solvers/Solver.hpp"

class Solver_hmm : Solver {
    public:

        Solver_hmm(Problem *prob, double p, int M);

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

        SDE_coeffs estimator(std::vector<double> x, double t);

    private:
        Problem* problem;
};
#endif
