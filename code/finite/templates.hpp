#ifndef TEMPLATES_H
#define TEMPLATES_H

#include <vector>

template<class type> std::vector<type> operator-(const std::vector<type>& v1, const std::vector<type>& v2) {
    std::vector<type> result(v1.size());
    for (unsigned int i = 0; i < v1.size(); ++i) {
        result[i] = v1[i] - v2[i];
    }
    return result;
}

template<class type> std::vector<type> operator+(const std::vector<type>& v1, const std::vector<type>& v2) {
    std::vector<type> result(v1.size());
    for (unsigned int i = 0; i < v1.size(); ++i) {
        result[i] = v1[i] + v2[i];
    }
    return result;
}

template<class type> double fabs(std::vector<type> vec) {
    double result = 0.;
    for (unsigned int i = 0; i < vec.size(); ++i) {
        result += fabs(vec[i]);
    }
    return result;
}

template<class type> std::vector<type> operator*(const std::vector<type>& vec, const double& x) {
    std::vector<type> result = vec;
    for (unsigned int i = 0; i < vec.size(); ++i) {
        result[i] = result[i]*x;
    }
    return result;
}

std::vector< std::vector<double> > operator*(const std::vector< std::vector<double> >& mat1, const std::vector< std::vector<double> >& mat2);
std::vector<double> operator*(const std::vector< std::vector<double> >& mat, const std::vector<double>& vec);
double operator*(const std::vector<double>& vec1, const std::vector<double>& vec2);

#endif
