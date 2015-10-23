#ifndef SOLVER_H
#define SOLVER_H

#include "structures.hpp"

class Solver {
    public:
        virtual SDE_coeffs estimator(std::vector<double> x, double t) = 0;
};

#endif
