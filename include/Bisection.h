#ifndef BISECTION
#define BISECTION

#include <CP_Root.h>

struct interval {
    mpz c;
    int k;
    poly p;
};

struct rootInterval {
    mpz c;
    int k, h;
    rootInterval(mpz c_ = 0, int k_ = 0, int h_ = 0) : c(c_), k(k_), h(h_) {}
};

int pErrorBound = 100;

const int iterationLimit = 100;

#endif