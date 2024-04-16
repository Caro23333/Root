#ifndef CP_ROOT_H
#define CP_ROOT_H
#include <vector>
#include <string>
#include <unordered_map>
#include "gmpxx.h"
#define MAXN 1005
#define ITER_LIM 20
#define LEN 1e-5

using namespace std;
typedef mpq_class mpq;
typedef vector<mpq_class> poly;

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
    {
        return m_root < t.m_root;
    }
}; // 类CP_Root定义结束

class CP_Interval
{
    public:
        double l;
        double r;
        CP_Interval(double rl = 0, double rr = 0): l(rl), r(rr) {}
        virtual ~CP_Interval() {}
}; // 求解时的二分区间

typedef unsigned long long ULL;


// 检查根的精确性
// 输入多项式的系数a: a[n-1]*x^(n-1)+...+a[1]*x+a[0], 其中n=a.size();
// 输入待检测的实根系列r;
// 输入nSpace: 是每行句首的空格数。
extern double gb_checkAccuracy(ostream& os, const vector<double> &a, const vector<double> &r,
    int nSpace);

// 比较浮点数a和b的绝对值
// 返回值: -1(|a|<|b|), 0(|a|==|b|), 1(|a|>|b|)
extern int gb_compareAbs(double a, double b);

// 统计在向量a中总共有多少个元素的值等于v。
extern int gb_coutVectorValue(const vector<int> &a, int v);

// 多项式展开
// 将多项式: r[0]*(x-r[1])*...*(x-r[n-1]),
// 展开为: a[n-1]*x^(n-1)+...+a[1]*x+a[0], 其中n=a.size()=r.size()。
// 要求n>=2。
extern void gb_expand(const vector<double> &r, vector<double> &a);

// 展开在r中的重根，从而将所有的根保存在a中。
extern void gb_expandRoot(const vector<CP_Root> &r, vector<double> &a);

// 从文件file当中读取一行字符并保存到字符串s当中。
extern void gb_getLine(string &s, ifstream& file);

// 计算下一个浮点数
extern double gb_getNext(double x);

// 计算前一个浮点数
extern double gb_getPrevious(double x);

// 计算多项式的值: a[n-1]*x^(n-1)+...+a[1]*x+a[0], 其中n=a.size()。
extern double gb_getValue(const vector<double> &a, double x);

// 返回: | |a|-|b| |
extern double gb_getValueCompareAbs(double a, double b);

// 如果c是回车符，则返回true；否则，返回false。
extern bool gb_isCharReturn(char c);

// 匹配两个浮点数序列a和b，输入e是容差。
// 输出va: 如果va[i]==-1，则表明在序列b中没有与a[i]相匹配的元素；否则，
//         va[i]记录与a[i]相匹配的序列b的下标值，其中va.size()==a.size();
// 输出vb: 如果vb[i]==-1，则表明在序列a中没有与b[i]相匹配的元素；否则，
//         vb[i]记录与b[i]相匹配的序列a的下标值，其中vb.size()==b.size();
// 输入ve: 如果va[i]==-1，则ve[i]=-1; 
//         否则，ve[i]=1: 如果||a[i]-b[va[i]]||<=e, 
//               ve[i]=0: 如果||a[i]-b[va[i]]||>e, 
// 输出maxd: 如果maxd=-1，表示没有找到最大误差处；否则，maxd为最大误差。
// 输出maxa: 如果maxa=-1，表示没有找到最大误差处；否则，最大误差所对应的序列a的下标。
// 输出maxb: 如果maxb=-1，表示没有找到最大误差处；否则，最大误差所对应的序列b的下标。
extern void gb_match(const vector<double> a, const vector<double> b, double e,
    vector<int> &va, vector<int> &vb, vector<int> &ve,
    double &maxd, int &maxa, int &maxb);

// 取最大值
extern double gb_max(double a, double b);

// 取最小值
extern double gb_min(double a, double b);

// 产生随机浮点数: 有效位不超过dataBitNum（要求dataBitNum<=1000）
// 指数在区间[exponentLow, exponentUp]。
extern double gb_randDouble(int dataBitNum, int exponentLow, int exponentUp);

// 产生在区间[lowerbound, upperbound]内的随机整数
extern int gb_randInt(int lowerbound, int upperbound);

// 输出匹配情况
// 输入a和b: 用来匹配的两个浮点数序列a和b，
// 输入e: 容差。
// 输入va: 如果va[i]==-1，则表明在序列b中没有与a[i]相匹配的元素；否则，
//         va[i]记录与a[i]相匹配的序列b的下标值，其中va.size()==a.size();
// 输入vb: 如果vb[i]==-1，则表明在序列a中没有与b[i]相匹配的元素；否则，
//         vb[i]记录与b[i]相匹配的序列a的下标值，其中vb.size()==b.size();
// 输入ve: 如果va[i]==-1，则ve[i]=-1; 
//         否则，ve[i]=1: 如果||a[i]-b[va[i]]||<=e, 
//               ve[i]=0: 如果||a[i]-b[va[i]]||>e, 
// 输入maxd: 如果maxd=-1，表示没有找到最大误差处；否则，maxd为最大误差。
// 输入maxa: 如果maxa=-1，表示没有找到最大误差处；否则，最大误差所对应的序列a的下标。
// 输入maxb: 如果maxb=-1，表示没有找到最大误差处；否则，最大误差所对应的序列b的下标。
// 输入nSpace: 是每行句首的空格数。
extern void gb_showMatch(ostream& os,
    const vector<double> a, const vector<double> b, double e,
    const vector<int> &va, const vector<int> &vb, const vector<int> &ve,
    double maxd, int maxa, int maxb, int nSpace);

// 输出多项式: a[n-1]*x^(n-1)+...+a[1]*x+a[0], 其中n=a.size()。
// 输入nSpace: 是每行句首的空格数。
extern void gb_showPolynomial(ostream& os, const vector<double> &a, int nSpace);

// 输出多项式: r[0]*(x-r[1])*...*(x-r[n-1]), 其中n=r.size()。
// 输入nSpace: 是每行句首的空格数。
extern void gb_showPolynomialByFactor(ostream& os, const vector<double> &r, int nSpace);

// 输出根的信息r。
// 输入nSpace: 是每行句首的空格数。
extern void gb_showRoots(ostream& os, const vector<CP_Root> &r, int nSpace);

// 输出n个空格。
extern void gb_showSpace(ostream& os, int n);

// 计算多项式所有实根的函数gb_solve的说明如下:
// 输入多项式的系数a: a[n-1]*x^(n-1)+...+a[1]*x+a[0], 其中n=a.size();
//                    要求a[n-1]不等于0。
// 输入求根的容差e: 要求e不小于0，求出的实根的误差应当不大于e。
// 输出flag: -1(表示输入有误)、
//           0(该多项式不存在实根)、
//           1(计算出所有实根，而且全部在容差范围之内)、
//           2(计算出所有实根，但部分或全部不在容差范围之内)、
//           3(只计算出部分实根，还有部分计算不出来。计算出来的都在容差范围之内)、
//           4(只计算出部分实根，还有部分计算不出来。但部分实根不在容差范围之内)、
//           5(该多项式存在实根，但都没有计算出来)、
//           6(没能判断出该多项式是否存在实根)。
// 输出r: r.size()是计算出来的实根的个数，r的各个元素是实根及其计算误差界，
//        要求r的各个元素对应的实根各不相等。
extern void gb_solve(const vector<double> &a, double e, int &flag, vector<CP_Root> &r);

// 将字符串s转换为浮点数
extern double gb_toDouble(const string &s);

// 将ULL数的内存解析为浮点数
extern double gb_toDouble(ULL x);

// 将浮点数d转换为字符串s
extern void gb_toString(string &s, double d);

// 将浮点数的内存解析为ULL数
extern ULL gb_toULL(double x);

extern void gb_solve2(const vector<double> &a, double e, int &flag, vector<CP_Root> &r);
extern void gb_solve3(const vector<double> &a, double e, int &flag, vector<CP_Root> &r);
extern void gb_solve4(const vector<double> &a, double e, int &flag, vector<CP_Root> &r);

// 近似约分，给定精度要求、近似方向
extern mpq gb_approxReduction(mpq x, int eps, int dir);

// 获取多项式在某一点的值
extern mpq gb_getValue(const poly &a, mpq x, int eps = 0);

// 获取多项式 a 对多项式 b 取模并取反的结果
extern poly gb_getMod(const poly &a, const poly &b);

// 求多项式的导数
extern poly gb_getDiff(const poly &a);

// 判断多项式值 aVal 是否满足容差 e 的要求，其导数值为 dVal
extern bool gb_judge(mpq aVal, mpq dVal, mpq e);

// 计算多项式在 x 处的重根数
extern int gb_getRepetition(const vector<poly> &da, mpq x, mpq e);

// 以 eps 的容差判断 x 与零的大小关系
extern int gb_compareWithZero(mpq x, mpq eps);

// 找到与 x 最近的使得多项式 a 非零的数。步长为 e，容差为 eps，方向为 dir
extern mpq gb_findNotZero(const poly &a, mpq x, mpq e, mpq eps, int dir);

// 计算 x 在 s 中的变号数
extern int gb_calcVar(const vector<poly> &s, mpq x, int eps);

// 援引 Sturm 定理计算 (a, b) 上的单根数。
extern int gb_getRootNumber(const vector<poly> &s, double l, double r, int eps);

// 牛顿法，初值为 x，容差为 e
extern mpq gb_newtonMethod(const poly &a, const poly &d, mpq x0, mpq e, int eps);

// 利用一阶导数估计当前根的误差
extern mpq gb_getError(const poly &a, const poly &d, mpq x);

extern std :: vector<poly> squareDecomposition(poly f);

extern inline double to_double(const mpq_class& x);

#endif
