#pragma once

#include <CP_Root.h>

#include <gmpxx.h>

#include <functional>
#include <fstream>
#include <utility>
#include <string>
#include <vector>
#include <chrono>

namespace testlib {

    constexpr double defaultErrorBound = 1e-15;

    struct poly {
        std::vector<double> coefficients;
        std::vector<double> roots;
        int degree;
        double e;
    } ;

    typedef std::function<void(const std::vector<double>&, double, int&, std::vector<CP_Root>&)> pSolveFunc;

    typedef std::chrono::high_resolution_clock clockType;
    typedef std::chrono::time_point<clockType> timePoint;

    class testMainApp {
    private:
        std::string path, base;
        int flag, mode;
        poly p;
        std::vector<CP_Root> res;
        timePoint timer;
        std::vector<std::pair<std::string, timePoint> > stopList;
        std::ofstream fout, ferr, fraw, fans;

        void getPoly();
        void readPoly(const std::string& fn);
        void prettyPoly(std::ostream& os);
        void rawPoly(std::ostream& os);
        void prettyAns(std::ostream& os);
        void logPolynomial();
        void compareRoots(const std::vector<CP_Root>& v1, const std::vector<CP_Root>& v2, std::vector<int>& r);
        void compareResult(std::vector<int>& r1, std::vector<int>& r2);

    public:
        ~testMainApp();
        void openFile(std::string path = "../data/", std::string base = "poly");
        void setMode(int mode);
        void logCaseNumber(int i, int nTests);
        void restartTimer();
        void lapTimer(const std::string& name);
        void displayTimer(std::ostream& os);
        bool test(std::ostream& os, pSolveFunc gb_solve);
        bool testMultiple(std::ostream& os, std::vector<pSolveFunc> gb_solves);
    };

    extern testMainApp app;
}

// void writePolynomial(std::ostream& os, const testlib::poly& p);
// testlib::poly generatePolynomial(int mode);
