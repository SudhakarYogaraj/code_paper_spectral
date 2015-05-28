#include "templates.hpp"

using namespace std;

vector< vector<double> > operator*(const vector< vector<double> >& mat1, const vector< vector<double> >& mat2) {
    vector< std::vector<double> > result (mat1.size(), vector<double> (mat1.size(), 0.));
    for (int i = 0; i < mat1.size(); i++) {
        for (int j = 0; j < mat1.size(); j++) {
            for (int k = 0; k < mat1.size(); ++k) {
                result[i][j] += mat1[i][k] * mat2[k][j];
            }
        }
    }
    return result;
}

vector<double> operator*(const vector< vector<double> >& mat, const vector<double>& vec) {
    vector<double> result(mat.size(),0.);
    for (int i = 0; i < mat.size(); i++) {
        for (int j = 0; j < mat.size(); j++) {
            result[i] += mat[i][j] * vec[j];
        }
    }
    return result;
}

double operator*(const vector<double>& vec1, const vector<double>& vec2) {
    double result = 0.;
    for (int i = 0; i < vec1.size(); ++i) {
        result += vec1[i] * vec2[i];
    }
    return result;
}
