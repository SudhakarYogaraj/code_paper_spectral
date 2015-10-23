#include "main.hpp"

using namespace std;

int main(int argc, char* argv[]) {

    // Initialization of the problem and helper analyser
    Problem problem;
    problem.init();
    Analyser analyser(&problem);


    // Initialization of the default solvers
    Solver_exact solver_exact(&problem, &analyser);
    Solver_hmm solver_hmm = Solver_hmm(&problem, 4, 1);
    Solver_spectral solver_spectral = Solver_spectral(&problem, &analyser, 15, 100);

    // Test of the solvers
    tests::error_spectral(problem.x0, &problem, &analyser);
    tests::error_hmm(problem.x0, &problem, &analyser);
    tests::integrate(&problem, &solver_exact, &solver_spectral, &solver_hmm);
}
