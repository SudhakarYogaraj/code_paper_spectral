#ifndef SOLVER_EXACT_H
#define SOLVER_EXACT_H

#include <vector>
#include "Problem.hpp"
#include "Analyser.hpp"

class Solver_exact {

    public:
        std::vector<double> soldrif(Problem& problem, Analyser& analyser, std::vector<double> x);
        std::vector< std::vector<double> > soldiff(Problem& problem, Analyser& analyser, std::vector<double> x);
};
#endif
