#include <iostream>
using namespace std;
#include <math.h>
#include "CP_RootQuadratic.h"

// gb_solveQuadratic计算一元二次方程的实根:
// 输入多项式的系数a: a[2]*x^2+...+a[1]*x+a[0]。
// 输入求根的容差e: 要求e不小于0，求出的实根的误差应当不大于e。
// 输出flag: -1(表示输入有误)、
//           0(该多项式不存在实根)、
//           1(计算出所有实根，而且全部在容差范围之内)、
//           3(只计算出部分实根，还有部分计算不出来。计算出来的都在容差范围之内)、
//           5(该多项式存在实根，但都没有计算出来)、
// 输出r: r.size()是计算出来的实根的个数，r的各个元素是实根及其计算误差界，
//        要求r的各个元素对应的实根各不相等。

void gb_solveQuadratic(const vector<double> &a, double e, int &flag, vector<CP_Root> &r)
{
    r.clear();
    if (e < 0)
    {
        flag = -1;
        return;
    } // if结束
    int n = a.size();
    if (n != 3)
    {
        flag = -1;
        return;
    } // if结束
    if (a[2] == 0)
    {
        flag = -1;
        return;
    } // if结束
    double d = a[1] * a[1] - 4 * a[0] * a[2];
    if (d < 0)
    {
        flag = 0;
        return;
    } // if结束
    double r1, r2, t;
    t = 2 * a[2];
    if (d == 0)
    {
        r1 = -a[1] / t;
        if ((!isinf(r1)) && (!isnan(r1)))
        {
            r.resize(1);
            r[0].m_repetitions = 2;
            r[0].m_root = r1;
            r[0].m_errorBound = 0; // 待计算: 实际的误差界是多少?
            flag = 1;
            return;
        } // if结束
        flag = 5;
        return;
    } // 外部if结束
    d = sqrt(d);
    r1 = (-a[1] + d) / t;
    r2 = (-a[1] - d) / t;
    if ((!isinf(r1)) && (!isnan(r1)))
    {
        if ((!isinf(r2)) && (!isnan(r2)))
        {
            if (r1 == r2)
            {
                r.resize(1);
                r[0].m_repetitions = 2;
                r[0].m_root = r1;
                r[0].m_errorBound = 0; // 待计算: 实际的误差界是多少?
                flag = 1;
                return;
            } // if结束
            r.resize(2);
            r[0].m_repetitions = 1;
            r[0].m_root = r1;
            r[0].m_errorBound = 0; // 待计算: 实际的误差界是多少?
            r[1].m_repetitions = 1;
            r[1].m_root = r2;
            r[1].m_errorBound = 0; // 待计算: 实际的误差界是多少?
            flag = 1;
            return;
        } // if结束
        r.resize(1);
        r[0].m_repetitions = 1;
        r[0].m_root = r1;
        r[0].m_errorBound = 0; // 待计算: 实际的误差界是多少?
        flag = 3;
        return;
    } // if结束
    if ((!isinf(r2)) && (!isnan(r2)))
    {
        r.resize(1);
        r[0].m_repetitions = 1;
        r[0].m_root = r2;
        r[0].m_errorBound = 0; // 待计算: 实际的误差界是多少?
        flag = 3;
        return;
    } // if结束
    flag = 5;
} // 函数gb_solveQuadratic结束
