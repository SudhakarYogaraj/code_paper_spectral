#ifndef TIME_INTEGRATION_H
#define TIME_INTEGRATION_H

#include <iostream>
#include <iomanip>
#include <vector>
#include <fstream>

namespace tests {
    void integrate(Problem *problem, Solver_exact *se, Solver_spectral *ss, Solver_hmm *sh);
}

#endif
