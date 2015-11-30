#ifndef ERROR_HMM_H
#define ERROR_HMM_H

#include <iostream>
#include <iomanip>
#include <fstream>

#include "global/global.hpp"

namespace tests {
    void error_hmm(std::vector<double> x, Problem*, Analyser*);
}

#endif
