// Includes
#include <iostream>
#include <iomanip>
#include <sstream>
#include <vector>
#include <vector>
#include <vector>
#include <cmath>
#include <random>
#include <algorithm>
#include <iterator>
#include <fstream>
#include <time.h>
#include <functional>
#include <stack>
#include <typeinfo>

#define PI 3.141592653589793238463

// Namespace
using namespace std;
stack<clock_t> tictoc_stack;

double ipow(double x, int e) {
    if(e == 0) return 1;
    if (e == 1) return x;
    double aux = ipow(x,e/2);
    if (e%2 == 0) return aux*aux;
    return x*aux*aux;
}
void tic() {
    tictoc_stack.push(clock());
}

int bin(int n, int k) {
    int res = 1;

    if ( k > n - k )
        k = n - k;

    for (int i = 0; i < k; ++i)
    {
        res *= (n - i);
        res /= (i + 1);
    }
    return res;
}

int mult2ind(vector<int> m, int d) {
    int l = m.size() - 1; int i;
    for(i = l; i > 0 & m[i] == 0; i--);

    if (i == 0 & m[0] == 0)
        return 0;

    int s = 0;
    for (int j = 0; j < m.size(); ++j)
        s += m[j];

    int dr = d - s; if(dr < 0) return -1;
    int vr = l - i;
    m[i] = m[i] - 1;
    /* cout << "dr " << dr << " vr " << vr << endl; */
    /* cout << "bin(dr + vr, vr)" <<  bin(dr + vr, vr) << endl; */
    return bin(dr + vr + 1, vr) + mult2ind(m, d);
}

vector<int> ind2mult(int ind, int d, int n) {
    vector<int> m(n,0); int s = 0;
    for (int i = 0; i < ind; ++i) {
        if (s < d) {
            m[n-1] ++; s++;
        } else {
            int j; for(j = n-1; m[j] == 0; j--);
            s -= m[j]; m[j] = 0; m[j-1] ++; s ++;
        }
    }
    return m;
}


void toc() {
    std::cout << "Time elapsed: "
        << ((double)(clock() - tictoc_stack.top())) / CLOCKS_PER_SEC
        << std::endl;
    tictoc_stack.pop();
}

static vector<double> nodes_30 = { 0.2011285765488714855458, 0.6039210586255523077782, 1.008338271046723461805, 1.415527800198188511941, 1.826741143603688038836, 2.243391467761504072473, 2.667132124535617200571, 3.099970529586441748689, 3.544443873155349886925, 4.003908603861228815228, 4.483055357092518341887, 4.988918968589943944486, 5.533147151567495725118, 6.138279220123934620395, 6.863345293529891581061};
static vector<double> weights_30 = { 0.3863948895418138625556, 0.2801309308392126674135, 0.1467358475408900997517, 0.05514417687023425116808, 0.01470382970482668351528, 0.00273792247306765846299, 3.48310124318685523421E-4, 2.9387252289229876415E-5, 1.57909488732471028835E-6, 5.10852245077594627739E-8, 9.178580424378528209E-10, 8.10618629746304420399E-12, 2.87860708054870606219E-14, 2.8103336027509037088E-17, 2.90825470013122622941E-21};

template<class type> std::vector<type> operator-(const std::vector<type>& v1, const std::vector<type>& v2) {
    std::vector<type> result(v1.size());
    for (unsigned int i = 0; i < v1.size(); ++i) {
        result[i] = v1[i] - v2[i];
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

template<class type> std::vector<type> operator+(const std::vector<type>& v1, const std::vector<type>& v2) {
    std::vector<type> result(v1.size());
    for (unsigned int i = 0; i < v1.size(); ++i) {
        result[i] = v1[i] + v2[i];
    }
    return result;
}

template<typename T, typename F> vector<T> quad1d(F f, double s, vector<T> v0) {
    std::vector<T> result = v0;
    double x_aux;
    for (unsigned int i = 0; i < nodes_30.size(); ++i) {
        x_aux = sqrt(2)*s*nodes_30[i];
        result = result + f(x_aux) * weights_30[i];
        result = result + f(-x_aux) * weights_30[i];
    }
    return result * (1./sqrt(PI));
}

template<typename T, typename F> vector<T> quadnd(F f, vector<double> s, vector<T> v0) {
    std::vector<T> result = v0;
    int size = s.size();
    if (size == 1) {
        auto lambda = [f](double x) -> vector<T> {
            vector<double> aux_vec(1,x);
            return f(aux_vec);
        };
        return quad1d(lambda, s[0], v0);
    }
    else {
        auto lambda = [f,size,s,v0](vector<double> x) -> vector<T> {
            auto nested_lambda = [f,&x](double y) -> vector<T> {
                x.push_back(y);
                vector<T> result = f(x);
                x.pop_back();
                return result;
            };
            return quad1d(nested_lambda, s[size-1], v0);
        };
        s.pop_back();
        return quadnd(lambda, s);
    }
    return result;
}

vector<double> func(vector<double> x) {
    vector<double> result = {x[0]*x[0], x[1]*x[1]*x[1]*x[1]};
    return result;
}

int main(int argc, char *argv[])
{
    double a = 2.;
    /* auto autofunc =  [a](double x) -> vector< vector<double> > { */
    /*     vector< vector<double> > result = {{1,x*x*x*x},{x,x*x*x*x*x*x}}; */
    /*     return result; */
    /* }; */

    vector<double> s = {1.,2.};
    vector<double> v0 = {0.,0.};
    vector<double> result = quadnd(func,s,v0);
    /* vector<double> v2 = times(v,2.); */

    /* for (int i = 0; i < 3; ++i) { */
    /*     for (int j = 0; j < 3; ++j) { */
    /*         for (int k = 0; k < 3; ++k) { */
    /*             vector<int> m = {i,j,k}; */
    /*             int ind = mult2ind(m,4); */
    /*             cout << i << j << k << " = " << ind << endl; */
    /*             if(ind >= 0) { */
    /*                 m = ind2mult(ind,4,3); */
    /*                 cout << m[0] << m[1] << m[2] << endl; */
    /*             } */
    /*         } */
    /*     } */
    /* } */
    return 0;
}
