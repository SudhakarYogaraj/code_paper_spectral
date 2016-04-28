#ifndef SOLVER_EXACT_H
#define SOLVER_EXACT_H

#include "global/global.hpp"
#include "problems/Problem.hpp"
#include "solvers/Analyser.hpp"
#include "solvers/Solver.hpp"

class Solver_exact : public Solver {

    public:
        Solver_exact(Problem *p, Analyser *a);
        SDE_coeffs estimator(std::vec x, double t);

    private:
        Problem *problem;
        Analyser *analyser;
        std::vec soldrif(std::vec x);
        std::mat soldiff(std::vec x);
};
#endif
