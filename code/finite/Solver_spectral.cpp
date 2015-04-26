#include "header.h"

void Solver_spectral::set(double p, int n)
{
    this->p = p;
    this->n_mcmc = 100000;
    this->degree = 4;
    this->nvars = n;

    this-> nBasis = bin(this->degree + this->nvars, this->nvars);
    this->ind2mult_aux = vector< vector<int> >(this->nBasis, vector<int>(this->nvars,0));
    this->mult2ind_aux = vector<int>(pow(degree + 1, this->nvars), -1);

    vector<int> currentMult(nvars,0);
    for (int i = 1; i < nBasis; ++i) {
        int sum = 0;
        for (int j = 0; j < nvars; ++j) {
            sum += currentMult[j];
        }
        if (sum < this->degree) {
            currentMult[nvars-1] ++;
        } else {
            int auxIndex = nvars - 1;
            while (currentMult[auxIndex] == 0) {
                auxIndex --;
            }
            currentMult[auxIndex] = 0;
            currentMult[auxIndex-1] = currentMult[auxIndex-1] + 1;
        }
        for (int j = 0; j < nvars; ++j) {
            ind2mult_aux[i][j] = currentMult[j];
        }
        mult2ind_aux[canonicalInd(currentMult, nvars, this->degree)] = i;
    }
}

int Solver_spectral::mult2ind(vector<int> alpha) {
    int canonicalInd = 0.;
    for (int j = 0; j < this->nvars; ++j) {
        canonicalInd += alpha[j]*pow(this->degree + 1, (this->nvars - 1) -j);
    }
    return this->mult2ind_aux[canonicalInd];
}

vector<int> Solver_spectral::ind2mult(int ind) {
    return ind2mult_aux[ind];
}
