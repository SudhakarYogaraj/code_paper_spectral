#ifndef ERROR_SPECTRAL_H
#define ERROR_SPECTRAL_H

#include <iostream>
#include <iomanip>
#include <vector>
#include <fstream>

#include "global/global.hpp"
#include "problem/Problem.hpp"
#include "problem/Analyser.hpp"

namespace tests {
    void error_spectral(std::vector<double> x, Problem*, Analyser*);
}

#endif
