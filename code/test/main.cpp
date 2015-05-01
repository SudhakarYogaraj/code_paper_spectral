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

int main(int argc, char *argv[])
{
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            for (int k = 0; k < 3; ++k) {
                vector<int> m = {i,j,k};
                int ind = mult2ind(m,4);
                cout << i << j << k << " = " << ind << endl;

                if(ind >= 0) {
                    m = ind2mult(ind,4,3);
                    cout << m[0] << m[1] << m[2] << endl;
                }
            }
        }
    }
    return 0;
}
