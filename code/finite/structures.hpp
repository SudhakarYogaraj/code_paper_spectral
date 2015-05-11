#ifndef STRUCTURES_H
#define STRUCTURES_H

#include <vector>

struct SDE_coeffs {
    std::vector<double> drif;
    std::vector< std::vector<double> > diff;
};

#endif
