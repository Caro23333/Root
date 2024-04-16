// 通用定义和声明 by Caro23333 1108

#ifndef CP_ROOT_H
#define CP_ROOT_H

#include <vector>
#include <string>
#include <unordered_map>
#include <utility>
#include <algorithm>
#include <iostream>
#include <cmath>
#include <queue>
#include <cstdlib>
#include <tuple>
#include <cassert>
#include <ctime>

#include "gmp.h"
#include "gmpxx.h"

typedef mpq_class mpq;
typedef mpz_class mpz;
typedef mpf_class mpf;
typedef std::vector<mpq_class> poly;
typedef unsigned long long ULL;
typedef long double LD;

class CP_Root
{
    public:
        int m_repetitions; // 当前根的重数: 要求不小于1。
        double m_errorBound; // 当前根的实际误差上界: 要求不小于0。
        double m_root; // 计算出来的实根，要求|m_root-精确实根|不大于m_errorBound。
        CP_Root() : m_repetitions(1), m_errorBound(0), m_root(0){ }
        CP_Root(int repetitions, double errorBound, double root)
        {
            m_repetitions = repetitions;
            m_errorBound = errorBound;
            m_root = root;
        } 
        virtual ~CP_Root() { }
        bool operator < (const CP_Root t) const
            { return m_root < t.m_root; }
}; 

extern double gb_getNext(double x);
extern double gb_getPrevious(double x);
extern double gb_toDouble(ULL x);
extern ULL gb_toULL(double x);

// 求解函数
extern void gb_solve2(const std::vector<double> &a, double e, int &flag, std::vector<CP_Root> &r);
extern void gb_solve3(const std::vector<double> &a, double e, int &flag, std::vector<CP_Root> &r);
extern void gb_solve4(const std::vector<double> &a, double e, int &flag, std::vector<CP_Root> &r);

// 近似约分，给定精度要求、近似方向
extern mpq gb_approxReduction(mpq x, int eps, bool greater);

// 获取多项式 a 对多项式 b 取模并取反的结果
extern poly gb_getMod(const poly &a, const poly &b);

// 求多项式的导数
extern poly gb_getDiff(const poly &a);

extern std::vector<poly> squareDecomposition(poly f);

extern double mpq_to_double(const mpq &x);

extern std::ostream & operator << (std::ostream & out, poly x);

extern poly reduce(poly f);

extern poly operator + (poly f, poly g);

extern poly operator - (poly f, poly g);

extern poly derivative(poly f);

extern poly negative(poly f);

extern poly opposite(poly f);

extern poly mod(poly f, poly g, poly * d = NULL);

extern poly div(poly f, poly g);

extern poly gcd(poly f, poly g);

extern int variate(poly f);

extern mpq calculate(poly f, mpq x);

extern poly reverse(poly f);

extern poly divbyx(poly f);

extern mpq simpleUpperBound(poly f);

extern inline double to_double(const mpq_class& x);

#endif
