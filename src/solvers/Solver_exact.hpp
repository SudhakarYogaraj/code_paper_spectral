#ifndef SOLVER_EXACT_H
#define SOLVER_EXACT_H

#include "global.hpp"
#include "problem/Problem.hpp"
#include "problem/Analyser.hpp"
#include "solvers/Solver.hpp"

class Solver_exact : public Solver {

    public:
        Solver_exact(Problem *p, Analyser *a);
        std::vector<double> soldrif(std::vector<double> x);
        std::vector< std::vector<double> > soldiff(std::vector<double> x);
        SDE_coeffs estimator(std::vector<double> x, double t);

    private:
        Problem *problem;
        Analyser *analyser;
};
#endif
