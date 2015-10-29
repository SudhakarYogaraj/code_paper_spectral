#include "main.hpp"

using namespace std;

int main(int argc, char* argv[]) {

    // Initialization of the problem and helper analyser
    Problem problem;
    problem.init();
    Analyser analyser(&problem);

    // Configuration for the HMM solver
    config_hmm conf_hmm = Solver_hmm::sensible_conf(4,1);

    // Initialization of the default solvers
    Solver_exact solver_exact(&problem, &analyser);
    Solver_spectral solver_spectral = Solver_spectral(&problem, &analyser, 15, 100);
    Solver_hmm solver_hmm = Solver_hmm(&problem, &conf_hmm);

    // Test of the solvers
    tests::error_spectral(problem.x0, &problem, &analyser);
    tests::error_hmm(problem.x0, &problem, &analyser);
    tests::integrate(&problem, &solver_exact, 0, "exact");
    tests::integrate(&problem, &solver_spectral, 0, "spectral");
    tests::integrate(&problem, &solver_hmm, 0, "hmm");
}
