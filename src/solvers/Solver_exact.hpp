#ifndef SOLVER_EXACT_H
#define SOLVER_EXACT_H

#include <vector>
#include "problem/Problem.hpp"
#include "problem/Analyser.hpp"

class Solver_exact {

    public:
        Solver_exact(Problem *p, Analyser *a);
        std::vector<double> soldrif(std::vector<double> x);
        std::vector< std::vector<double> > soldiff(std::vector<double> x);

    private:
        Problem *problem;
        Analyser *analyser;
};
#endif
