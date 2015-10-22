#include "toolbox/combinatorics.hpp"

using namespace std;

// Binomial coefficients
int bin(int n, int k) {
    int res = 1;
    if ( k > n - k )
        k = n - k;
    for (int i = 0; i < k; ++i) {
        res *= (n - i);
        res /= (i + 1);
    }
    return res;
}

// Standard normal gaussian
double gaussian(vector<double> y) {
    double result = 1.;
    for (unsigned int i = 0; i < y.size(); ++i) {
        result *= exp(-y[i]*y[i]/2)/(sqrt(2*PI));
    }
    return result;
}

// Enumeration of the n-dimensional multi-indices i such that |i| = d
vector< vector<int> > equal_multi_indices(int n, int d) {

    // Univariate case
    if (n == 1)
        return vector< vector<int> > (1, vector<int> (1, d));

    // Case where |i| = 0
    if (d == 0)
        return vector< vector<int> > (1, vector<int> (n, 0));

    // Initialization of result
    vector< vector<int> > result(0);

    // Iteration on the first variable
    for (int i = 0; i <= d; ++i) {
        vector< vector<int> > aux = equal_multi_indices(n-1, d-i);
        for (unsigned int j = 0; j < aux.size(); ++j) {
            aux[j].push_back(i);
            result.push_back(aux[j]);
        }
    }

    return result;
}

// Enumeration of the n-dimensional multi-indices i such that a <= |i| <= b
vector< vector<int> > interval_multi_indices(int n, int a, int b) {

    if (b == a)
        return equal_multi_indices(n, a);

    // Auxiliary results
    vector< vector<int> > r_lower = interval_multi_indices(n, a, b-1);
    vector< vector<int> > r_equal = equal_multi_indices(n, b);

    // Concatenation of the vectors
    r_lower.insert(r_lower.end(), r_equal.begin(), r_equal.end());

    return r_lower;
}

// Enumeration of the n-dimensional multi-indices i such that |i| <= d
vector< vector<int> > lower_multi_indices(int n, int d) {
    return interval_multi_indices(n, 0, d);
}
