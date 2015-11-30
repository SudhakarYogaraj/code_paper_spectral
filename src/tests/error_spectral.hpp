#ifndef ERROR_SPECTRAL_H
#define ERROR_SPECTRAL_H

#include <iostream>
#include <iomanip>
#include <vector>
#include <fstream>

#include "global/global.hpp"
#include "problems/Problem.hpp"
#include "problems/Analyser.hpp"

namespace tests {
    void error_spectral(std::vector<double> x, Problem*, Analyser*);
}

#endif
