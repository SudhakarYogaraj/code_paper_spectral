#ifndef SOLVER_SPECTRAL_H
#define SOLVER_SPECTRAL_H

#include <numeric>
#include <functional>
#include <vector>
#include "aux.hpp"
#include "Problem.hpp"

class Solver_spectral {
    public:

        int n_mcmc;
        int degree;
        int nvars;
        double p;
        int nBasis;
        std::vector< std::vector<int> > ind2mult_aux;
        std::vector<int> mult2ind_aux;

		void set(double p, int n);
		int mult2ind(std::vector<int>);
		std::vector<int> ind2mult(int);
        void estimator(Problem&, std::vector<double> x, std::vector<double>& fi, std::vector< std::vector<double> >& hi, double t);
};
#endif
