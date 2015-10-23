#ifndef MAIN_H
#define MAIN_H

#include "structures.hpp"
#include "io/tictoc.hpp"
#include "templates.hpp"
#include "io/io.hpp"
#include "global.hpp"
#include "toolbox/combinatorics.hpp"

// Problem and analyser
#include "problem/Problem.hpp"
#include "problem/Analyser.hpp"

// Solvers
#include "solvers/Solver_spectral.hpp"
#include "solvers/Solver_hmm.hpp"
#include "solvers/Solver_exact.hpp"

// Tests
#include "tests/error_spectral.hpp"
#include "tests/error_hmm.hpp"

#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>
#include <random>

#endif
