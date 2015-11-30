#ifndef TIME_INTEGRATION_H
#define TIME_INTEGRATION_H

#include <iostream>
#include <iomanip>
#include <vector>
#include <fstream>

#include "global/global.hpp"

namespace tests {
    void integrate(Problem*, Solver*, int seed, std::string s);
}

#endif
