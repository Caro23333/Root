#ifndef STURM_NEWTON
#define STURM_NEWTON

#include <CP_Root.h>

#define ITER_LIM 20
#define PREC 256

class CP_Interval 
{
    public:
        double l;
        double r;
        CP_Interval(double rl = 0, double rr = 0): l(rl), r(rr) {}
        virtual ~CP_Interval() {}
}; 

// 用 mpf 类型求多项式的值
extern mpq gb_getValue(const poly &a, mpq x);

// 判断多项式值 aVal 是否满足容差 e 的要求，其导数值为 dVal
extern bool gb_judge(mpq aVal, mpq dVal, mpq e);

// 以 eps 的容差判断 x 与零的大小关系
extern int gb_compareWithZero(mpq x, mpq eps);

// 找到与 x 最近的使得多项式 a 非零的数。步长为 e，容差为 eps，方向为 dir
extern mpq gb_findNotZero(const poly &a, mpq x, mpq e, mpq eps, int dir);

// 计算 x 在 s 中的变号数
extern int gb_calcVar(const std::vector<poly> &s, mpq x, int eps);

// 援引 Sturm 定理计算 (a, b) 上的单根数。
extern int gb_getRootNumber(const std::vector<poly> &s, double l, double r, int eps);

// 牛顿法，初值为 x，容差为 e
extern mpq gb_newtonMethod(const poly &a, const poly &d, mpq x0, mpq e, int eps);

#endif