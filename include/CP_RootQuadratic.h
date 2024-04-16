#ifndef CP_ROOTQUADRATIC_H
#define CP_ROOTQUADRATIC_H
#include "CP_Root.h"

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
extern void gb_solveQuadratic(const vector<double> &a, double e, int &flag, vector<CP_Root> &r);
#endif
