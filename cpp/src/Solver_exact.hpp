#ifndef SOLVER_EXACT_H
#define SOLVER_EXACT_H

#include <vector>
#include "Problem.hpp"

class Solver_exact {

    public:
        std::vector<double> soldrif(Problem problem, std::vector<double> x);
        std::vector< std::vector<double> > soldiff(Problem problem, std::vector<double> x);
};
#endif
