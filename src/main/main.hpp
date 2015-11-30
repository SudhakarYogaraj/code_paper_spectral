#ifndef MAIN_H
#define MAIN_H

#include "global/global.hpp"
#include "io/tictoc.hpp"
#include "global/templates.hpp"
#include "io/io.hpp"
#include "toolbox/combinatorics.hpp"

// Problem and analyser
#include "problems/Problem.hpp"
#include "solvers/Analyser.hpp"

// Solvers
#include "solvers/Solver_spectral.hpp"
#include "solvers/Solver_hmm.hpp"
#include "solvers/Solver_exact.hpp"

// Tests
#include "tests/error_spectral.hpp"
#include "tests/error_hmm.hpp"
#include "tests/time_integration.hpp"

#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>
#include <random>

#endif
