#include "io/tictoc.hpp"
using namespace std;

stack<clock_t> tictoc_stack;

void tic() {
    tictoc_stack.push(clock());
}

double toc() {
    double toReturn = ((double)(clock() - tictoc_stack.top())) / CLOCKS_PER_SEC;
    tictoc_stack.pop();
    return toReturn;
}
