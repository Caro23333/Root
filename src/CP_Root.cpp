#include <iostream>
#include <fstream>
#include <queue>
#include <cstdlib>
#include <algorithm>
#include <unordered_map>
#include <assert.h>
#include "CP_Root.h"
#include "testlib.h"

using namespace std;

// 检查根的精确性
// 输入多项式的系数a: a[n-1]*x^(n-1)+...+a[1]*x+a[0], 其中n=a.size();
// 输入待检测的实根系列r;
// 输入nSpace: 是每行句首的空格数。
double gb_checkAccuracy(ostream& os, const vector<double> &a, const vector<double> &r,
    int nSpace)
{
    int nr = r.size();
    int nn = 0; // nn记录误差大于0的实根个数
    int i, k, idm;
    double d, dp, dn, dtp, dtn, ds;
    double dm = -1;
    string s;
    for (i = 0; i < nr; i++)
    {
        d = gb_getValue(a, r[i]);
        if (d != 0)
        {
            nn++;
            ds = (d < 0 ? -d : d);
            if (dm < ds)
            {
                dm = ds;
                idm = i;
            } // if结束
            gb_showSpace(os, nSpace+4);
            os << "[" << i << "]root: ";
            gb_toString(s, r[i]);
            os << s << endl;
            gb_showSpace(os, nSpace + 8);
            os << "The result of substitution is not 0 but ";
            gb_toString(s, d);
            os << s << endl;

            dp = gb_getPrevious(r[i]);
            dtp = gb_getValue(a, r[i]);
            if (dtp == 0)
            {
                gb_showSpace(os, nSpace + 8);
                os << "The value of polynomial with the previous floating point number p is 0, where p = ";
                gb_toString(s, dp);
                os << s << endl;
            }
            else
            {
                k = gb_compareAbs(dtp, d);
                if (k == -1)
                {
                    gb_showSpace(os, nSpace + 8);
                    os << "The value of polynomial with the previous floating point number p is less, where p = ";
                    gb_toString(s, dp);
                    os << s << endl;
                    gb_showSpace(os, nSpace + 8);
                    os << "The value of polynomial with p is ";
                    gb_toString(s, dtp);
                    os << s << endl;
                    gb_showSpace(os, nSpace + 8);
                    dp = gb_getValueCompareAbs(dtp, d);
                    os << "The absolute value of the difference between two values is: ";
                    gb_toString(s, dp);
                    os << s << endl;
                } // if (k == -1)结束
            } // if/else结束

            dn = gb_getNext(r[i]);
            dtn = gb_getValue(a, r[i]);
            if (dtn == 0)
            {
                gb_showSpace(os, nSpace + 8);
                os << "The value of polynomial with the next floating point number q is 0, where q = ";
                gb_toString(s, dn);
                os << s << endl;
            }
            else
            {
                k = gb_compareAbs(dtn, d);
                if (k == -1)
                {
                    gb_showSpace(os, nSpace + 8);
                    os << "The value of polynomial with the next floating point number q is less, where q = ";
                    gb_toString(s, dn);
                    os << s << endl;
                    gb_showSpace(os, nSpace + 8);
                    os << "The value of polynomial with q is ";
                    gb_toString(s, dtn);
                    os << s << endl;
                    gb_showSpace(os, nSpace + 8);
                    dn = gb_getValueCompareAbs(dtn, d);
                    os << "The absolute value of the difference between two values is: ";
                    gb_toString(s, dn);
                    os << s << endl;
                } // if (k == -1)结束
            } // if/else结束
        } // if (d != 0)结束
    } // for结束
    gb_showSpace(os, nSpace + 4);
    os << "The polynomial equals to 0 is not satisfied with " << nn << " roots in total." << endl;
    if (nn > 1)
    {
        gb_showSpace(os, nSpace + 4);
        os << "The maximum absolute value of the polynomial is: ";
        gb_toString(s, dm);
        os << s << endl;
        gb_showSpace(os, nSpace + 4);
        os << "With the root[" << idm << "]." << endl;
    } // if结束
    return dm;
} // 函数gb_checkAccuracy结束

// 比较浮点数a和b的绝对值
// 返回值: -1(|a|<|b|), 0(|a|==|b|), 1(|a|>|b|)
int gb_compareAbs(double a, double b)
{
    if (a < 0)
        a = -a;
    if (b < 0)
        b = -b;
    if (a < b)
        return -1;
    if (a == b)
        return 0;
    return 1;
} // 函数gb_compareAbs结束

// 统计在向量a中总共有多少个元素的值等于v。
int gb_coutVectorValue(const vector<int> &a, int v)
{
    int n = a.size();
    int i;
    int r = 0;
    for (i = 0; i < n; i++)
        if (a[i] == v)
            r++;
    return r;
} // 函数gb_coutVectorValue结束

// 多项式展开
// 将多项式: r[0]*(x-r[1])*...*(x-r[n-1]),
// 展开为: a[n-1]*x^(n-1)+...+a[1]*x+a[0], 其中n=a.size()=r.size()。
// 要求n>=2。
void gb_expand(const vector<double> &r, vector<double> &a)
{
    int n = r.size();
    if (n <= 1)
    {
        a.clear();
        return;
    } // if结束
    int m = a.size();
    if (m != n)
        a.resize(n);
    int i, j;
    a[0] = - r[0] * r[1];
    a[1] = r[0];
    double t;
    for (i = 2; i < n; i++)
    {
        a[i] = a[i - 1];
        t = -r[i];
        for (j = i - 1; j>0; j--)
        {
            a[j] = a[j - 1] + a[j] * t;
        } // 内部for结束
        a[0] *= t;
    } // 外部for结束
} // 函数gb_expand结束

// 展开在r中的重根，从而将所有的根保存在a中。
void gb_expandRoot(const vector<CP_Root> &r, vector<double> &a)
{
    int n = r.size();
    if (n <= 0)
    {
        a.clear();
        return;
    } // if结束
    int na = n;
    int i, j, k;
    for (i = 0; i < n; i++)
        if (r[i].m_repetitions > 1)
            na += (r[i].m_repetitions - 1);
    i = a.size();
    if (i != na)
        a.resize(na);
    for (i = 0, k = 0; i < n; i++)
        for (j = 0; j < r[i].m_repetitions; j++)
            a[k++] = r[i].m_root;
} // 函数gb_expandRoot结束

// 从文件file当中读取一行字符并保存到字符串s当中。
void gb_getLine(string &s, ifstream& file)
{
    char c;
    s.clear();
    do
    { // 首先跳过回车与换行符
        c = (char)(file.peek());
        if (!gb_isCharReturn(c))
            break;
        file.get();
        if (file.eof())
            return;
    } while (true);
    do
    { // 接着读取数据
        file.get(c);
        if (!file.good())
            break;
        if (!gb_isCharReturn(c))
            s += c;
        else break;
    } while (true);
} // 函数gb_getLine结束

// 计算下一个浮点数
double gb_getNext(double x)
{
    ULL a = gb_toULL(x);
    ULL sign = a >> 63;
    ULL exponent = (a >> 52) & 0x7ff;
    ULL fraction = a & 0xfffffffffffff;

    if (exponent != 0 && exponent != 0x7ff) // 规格化的
    {

        if (sign == 0 && exponent == 0x7fe && fraction == 0xfffffffffffff) // 最大的正规格化数
            return gb_toDouble(0x7ff0000000000000);
        if (sign == 1 && exponent == 0x001 && fraction == 0) // 最大的负规格化数
            return gb_toDouble(0x800fffffffffffff);
        if (sign) // 负数
        {
            if (fraction == 0)
                exponent--, fraction = 0xfffffffffffff;
            else
                fraction--;
        }
        else // 正数
        {
            fraction++;
            if (fraction == 0x10000000000000)
                exponent++, fraction = 0;
        }
    }
    else if (exponent == 0) // 非规格化的
    {
        if (sign == 0 && fraction == 0xfffffffffffff) // 最大的正非规格化数
            return gb_toDouble(0x0010000000000000);
        if (sign == 1 && fraction == 0) // 最大的负非规格化数
            return gb_toDouble(0x0000000000000000);
        if (sign) // 负数
            fraction--;
        else // 正数
            fraction++;
    }
    else if (exponent == 0x7ff && fraction == 0) // 无穷大
    {
        if (sign) // 负无穷大
            return gb_toDouble(0xffefffffffffffff);
        else // 正无穷大
            return gb_toDouble(0xfff0000000000000);
    }
    else // NaN
        return x;
    return gb_toDouble((sign << 63) | (exponent << 52) | fraction);
} // 函数gb_getNext结束

// 计算前一个浮点数
double gb_getPrevious(double x)
{
    ULL a = gb_toULL(x);
    ULL sign = a >> 63;
    ULL exponent = (a >> 52) & 0x7ff;
    ULL fraction = a & 0xfffffffffffff;

    if (exponent != 0 && exponent != 0x7ff) // 规格化的
    {

        if (sign == 0 && exponent == 0x001 && fraction == 0) // 最小的正规格化数
            return gb_toDouble(0x000fffffffffffff);
        if (sign == 1 && exponent == 0x7fe && fraction == 0xfffffffffffff) // 最小的负规格化数
            return gb_toDouble(0xfff0000000000000);
        if (sign) // 负数
        {
            fraction++;
            if (fraction == 0x10000000000000)
                exponent++, fraction = 0;
        }
        else // 正数
        {
            if (fraction == 0)
                exponent--, fraction = 0xfffffffffffff;
            else
                fraction--;
        }
    }
    else if (exponent == 0) // 非规格化的
    {
        if (sign == 0 && fraction == 0) // 最小的正非规格化数
            return gb_toDouble(0x8000000000000000);
        if (sign == 1 && fraction == 0xfffffffffffff) // 最小的负非规格化数
            return gb_toDouble(0x8010000000000000);
        if (sign) // 负数
            fraction++;
        else // 正数
            fraction--;
    }
    else if (exponent == 0x7ff && fraction == 0) // 无穷大
    {
        if (sign) // 负无穷大
            return gb_toDouble(0x7ff0000000000000);
        else // 正无穷大
            return gb_toDouble(0x7fefffffffffffff);
    }
    else // NaN
        return x;
    return gb_toDouble((sign << 63) | (exponent << 52) | fraction);
} // 函数gb_getPrevious结束

// 计算多项式的值: a[n-1]*x^(n-1)+...+a[1]*x+a[0], 其中n=a.size()。
double gb_getValue(const vector<double> &a, double x)
{
    int n = a.size();
    int i;
    double t = 1;
    double r = 0;
    for (i = 0; i < n; i++)
    {
        r += (t*a[i]);
        t *= x;
    } // for结束
    return r;
} // 函数gb_getValue结束

// 如果c是回车符，则返回true；否则，返回false。
bool gb_isCharReturn(char c)
{
    if ((c == '\r') || (c == '\n') || (c == '\f'))
        return true;
    return false;
} // gb_isCharReturn函数结束

// 返回: | |a|-|b| |
double gb_getValueCompareAbs(double a, double b)
{
    if (a < 0)
        a = -a;
    if (b < 0)
        b = -b;
    double d = a - b;
    if (d < 0)
        d = -d;
    return d;
} // 函数gb_getValueCompareAbs结束

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
void gb_match(const vector<double> a, const vector<double> b, double e,
        vector<int> &va, vector<int> &vb, vector<int> &ve,
        double &maxd, int &maxa, int &maxb)
{
    maxd = -1;
    maxa = -1;
    maxb = -1;
    int na = a.size();
    int nb = b.size();
    int i, j, k, n;
    if (na > 0)
    {
        n = va.size();
        if (n != na)
            va.resize(na);
        n = ve.size();
        if (n != na)
            ve.resize(na);
        for (i = 0; i < na; i++)
        {
            va[i] = -1;
            ve[i] = -1;
        } // for结束
    }
    else
    {
        va.clear();
        ve.clear();
    } // if/else结束

    if (nb > 0)
    {
        n = vb.size();
        if (n != nb)
            vb.resize(nb);
        for (i = 0; i < nb; i++)
            vb[i] = -1;
    }
    else vb.clear();

    if ((na <= 0) || (nb <= 0))
        return;

    double d, dk;
    for (j = 0; j < nb; j++)
    {
        k = -1;
        for (i = 0; i < na; i++)
        {
            if (va[i] == -1)
            {
                if (a[i] >= b[j])
                    d = a[i] - b[j];
                else d = b[j] - a[i];
                if (k == -1)
                {
                    k = i;
                    dk = d;
                }
                else if (dk>d)
                {
                    k = i;
                    dk = d;
                } // if/else if结束
            } //if结束
        } // 内部
        if (k >= 0)
        {
            va[k] = j;
            vb[j] = k;
            if (dk <= e)
                ve[k] = 1;
            else ve[k] = 0;
            if (maxd < dk)
            {
                maxd = dk;
                maxa = k;
                maxb = j;
            } // if结束
        } // if结束
    } // 外部for结束
} // 函数gb_match结束

// 取最大值
double gb_max(double a, double b)
{
    if (a >= b)
        return a;
    return b;
} // 函数gb_max结束

// 取最小值
double gb_min(double a, double b)
{
    if (a <= b)
        return a;
    return b;
} // 函数gb_max结束


// 产生随机浮点数: 有效位不超过dataBitNum（要求dataBitNum<=1000）
// 指数在区间[exponentLow, exponentUp]。
double gb_randDouble(int dataBitNum, int exponentLow, int exponentUp)
{
    string s;
    char sa[1010];
    int i, k;
    k = gb_randInt(0, 1);
    if (k == 0)
        s += "-";
    k = gb_randInt(1, 9);
    sa[0] = '0' + (char)k;
    sa[1] = '.';
    if (dataBitNum > 1000)
        dataBitNum = 1000;
    dataBitNum += 1;
    for (i = 2; i < dataBitNum; i++)
    {
        k = gb_randInt(0, 9);
        sa[i] = '0' + (char)k;
    } // for结束
    sa[i] = 0;
    s += sa;
    k = gb_randInt(exponentLow, exponentUp);
    if (k != 0)
    {
        snprintf(sa, 1010, "e%d", k);
        s += sa;
    } // if结束
    double r = gb_toDouble(s);
    return r;
} // 函数gb_randDouble结束

// 产生在区间[lowerbound, upperbound]内的随机整数
int gb_randInt(int lowerbound, int upperbound)
{
    if (lowerbound == upperbound)
        return lowerbound;
    int result;
    if (lowerbound < upperbound)
        result = (rand() % (upperbound - lowerbound + 1)) + lowerbound;
    else result = (rand() % (lowerbound - upperbound + 1)) + upperbound;
    return result;
} // 函数gb_randInt结束

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
void gb_showMatch(ostream& os,
    const vector<double> a, const vector<double> b, double e,
    const vector<int> &va, const vector<int> &vb, const vector<int> &ve,
    double maxd, int maxa, int maxb, int nSpace)
{
    gb_showSpace(os, nSpace);
    string s;
    os << "The maximum error of the roots is: ";
    gb_toString(s, maxd);
    os << s << ", " << endl;
    if ((maxa != -1) && (maxb != -1))
    {
        gb_showSpace(os, nSpace+4);
        os << "The root from input which makes maximum error is: ";
        gb_toString(s, a[maxa]);
        os << s << "," << endl;
        gb_showSpace(os, nSpace + 4);
        os << "While the calculated root is: ";
        gb_toString(s, b[maxb]);
        os << s << "." << endl;
    } // if结束
    int nva = va.size();
    int i;
    i = gb_coutVectorValue(va, -1);
    if (i > 0)
    {
        gb_showSpace(os, nSpace);
        os << "There are totally " << i << " roots from input have(has) not been figured out. They are(is), respectively: " << endl;
        for (i = 0; i < nva; i++)
        {
            if (va[i] == -1)
            {
                gb_showSpace(os, nSpace+4);
                os << "[" << i << "]: ";
                gb_toString(s, a[i]);
                os << s << endl;
            } // if结束
        } // for结束
    } // if结束
    int nvb = vb.size();
    i = gb_coutVectorValue(vb, -1);
    if (i > 0)
    {
        gb_showSpace(os, nSpace);
        os << "There are(is) totally " << i << " unexpected roots have(has) been figured out. They are(is), respectively: " << endl;
        for (i = 0; i < nvb; i++)
        {
            if (vb[i] == -1)
            {
                gb_showSpace(os, nSpace + 4);
                os << "[" << i << "]: ";
                gb_toString(s, b[i]);
                os << s << endl;
            } // if结束
        } // for结束
    } // if结束
    int ne = ve.size();
    i = gb_coutVectorValue(ve, 0);
    double d;
    if (i > 0)
    {
        gb_showSpace(os, nSpace);
        os << "There are(is) " << i << " pair(s) of roots which exceeded the error bound. They are, respectively: " << endl;
        for (i = 0; i < ne; i++)
        {
            if (ve[i] == 0)
            {
                gb_showSpace(os, nSpace + 4);
                os << "[" << i << "]the root from input is: ";
                gb_toString(s, a[i]);
                os << s << endl;
                gb_showSpace(os, nSpace + 4);
                os << "[" << i << "]the calculated root is ";
                gb_toString(s, b[va[i]]);
                os << s << endl;
                d = a[i] - b[va[i]];
                if (d < 0)
                    d = -d;
                gb_showSpace(os, nSpace + 4);
                os << "[" << i << "]with an error of: ";
                gb_toString(s, d);
                os << s << endl;
            } // if结束
        } // for结束
    } // if结束
} // 函数gb_showMatch结束

// 输出多项式: a[n-1]*x^(n-1)+...+a[1]*x+a[0], 其中n=a.size()。
// 输入nSpace: 是每行句首的空格数。
void gb_showPolynomial(ostream& os, const vector<double> &a, int nSpace)
{
    int n = a.size();
    if (n <= 0)
    {
        gb_showSpace(os, nSpace);
        os << "The polynomial is empty." << endl;
    } // if 结束
    string s;
    int i, flag = 0;
    for (i = n-1; i >= 0; i--)
    {
        if (a[i] != 0)
        {
            gb_showSpace(os, nSpace);
            if ((i != (n - 1)) && (flag==1))
                os << "+";
            if (i>1)
                os << "x^" << i << "*";
            else if (i==1)
                os << "x" << "*";
            gb_toString(s, a[i]);
            os << '(' << s << ')' << endl;
            flag = 1;
        } // if结束
    } // for结束
    if (flag == 0)
        os << "0" << endl;
} // 函数gb_showPolynomial结束

// 输出多项式: r[0]*(x-r[1])*...*(x-r[n-1]), 其中n=r.size()。
// 输入nSpace: 是每行句首的空格数。
void gb_showPolynomialByFactor(ostream& os, const vector<double> &r, int nSpace)
{
    int n = r.size();
    if (n <= 0)
    {
        gb_showSpace(os, nSpace);
        os << "The polynomial is empty." << endl;
    } // if 结束
    gb_showSpace(os, nSpace);
    string s;
    gb_toString(s, r[0]);
    os << "(" << s << ")" << endl;
    int i;
    for (i = 1; i < n; i++)
    {
        gb_showSpace(os, nSpace);
        if (r[i] == 0)
        {
            os << "*x" << endl;
            continue;
        } // if结束
        os << "*(x-(";
        gb_toString(s, r[i]);
        os << s;
        os << "))" << endl;
    } // for结束
} // 函数gb_showPolynomialByFactor结束

// 输出根的信息r。
// 输入nSpace: 是每行句首的空格数。
void gb_showRoots(ostream& os, const vector<CP_Root> &r, int nSpace)
{
    int n = r.size();
    gb_showSpace(os, nSpace);
    os << "There are(is) " << n << " different roots in total." << endl;
    int i;
    string s;
    for (i = 0; i < n; i++)
    {
        gb_showSpace(os, nSpace);
        os << "The value of " << (i + 1) << "th root is: ";
        gb_toString(s, r[i].m_root);
        os << s << "." << endl;
        gb_showSpace(os, nSpace+4);
        os << "which repeats for " << r[i].m_repetitions << " times." << endl;
        gb_showSpace(os, nSpace + 4);
        os << "The calculated error bound is: ";
        gb_toString(s, r[i].m_errorBound);
        os << s << "." << endl;
    } // for结束
} // 函数gb_showRoots结束

// 输出n个空格。
void gb_showSpace(ostream& os, int n)
{
    int i;
    for (i = 0; i < n; i++)
        os << ' ';
} // 函数gb_showSpace结束

// 将字符串s转换为浮点数
double gb_toDouble(const string &s)
{
    double d = strtod(s.c_str(), nullptr);
    return d;
} // 函数gb_toDouble结束

// 将ULL数的内存解析为浮点数
double gb_toDouble(ULL x)
{
    void *a = &x;
    return *((double *)a);
} // 函数gb_toDouble结束

// 将浮点数d转换为字符串s
void gb_toString(string &s, double d)
{
    char str[1100];
    string s1, s2;
    snprintf(str, 1100, "%.1000e", d);
    s = str;
    int i, j, k;
    i = s.rfind('e');
    if (i == -1)
    { // 没有找到字母'e'
        j = s.rfind('.');
        if (j == -1)
            return; // 表明没有小数点
        j = s.find_last_not_of('0');
        if (j == -1)
            return; // 表明末尾不是0
        s = s.substr(0, j + 1);
        return;
    } // if结束
    s1 = s.substr(0, i); // 获取到小数部分
    s2 = s.substr(i); // 获取到指数部分
    j = s1.rfind('.');
    if (j >= 0)
    { // 表明有小数点
        k = s1.find_last_not_of('0');
        if (k >= 0)
        {
            s1 = s1.substr(0, k + 1);
            if (j == k)
                s1 += "0";
            s = s1 + s2;
        } // 内部if结束
    } // 外部if结束
} // 函数gb_toString结束

// 将浮点数的内存解析为ULL数
ULL gb_toULL(double x)
{
    void *a = &x;
    return *((ULL *)a);
} // 函数gb_toULL结束

double gb_ABS(double x) // 获取浮点数的绝对值
{
    return x >= 0 ? x : -x;
}

int gb_compareWithZero(double x, double eps) // 以 eps 的容差范围判定 x 与 0 的大小关系. < = >  -->  -1 0 1
{
   // double lim = gb_min(EPS, gb_ABS(hi * 1e-6));
    if(gb_ABS(x) <= eps)
        return 0;
    else if(x > eps)
        return 1;
    else 
        return -1;
}

int gb_getRepetition(const vector<double> *da, double x, double e) // 求出多项式在 x 处重根的数量(至少为 1)
{
    int n = da[0].size();
    double nowVal = gb_getValue(da[1], x), nxtVal;
    for(int i = 1; i <= n - 1; i++)
    {
        if(i == n - 1)
        {
            if(da[i][0] != 0)
                return i;
        }
        else
        {
            nxtVal = gb_getValue(da[i + 1], x);
            if(nxtVal != 0 && gb_compareWithZero(nowVal / nxtVal, e) != 0)
                return i;
        }
    }
    return n;
}

double gb_findNotZero(const vector<double> &a, double x, double e, double eps, int dir) // 找到相邻的使多项式非零的浮点数. dir = 1, -1  --> next, previous
{
    double res = x;
    double val = gb_getValue(a, x);
    int flag = gb_compareWithZero(val, eps);
    while(flag == 0)
    {
        //res = (dir == 1 ? gb_getNext(res) : gb_getPrevious(res));
        res = (dir == 1 ? res + e : res - e); // 在已知为根的数附近找下一个数时，以容差 e 为步长，否则以最近浮点数为步长
        val = gb_getValue(a, res);
        flag = gb_compareWithZero(val, eps);
    }
    return res;
}

int gb_calcSignVar(const vector<double> *t, unordered_map<double, int> &v, double x) // 求出 x 处 Sturm 序列的变号次数
{
    if(v.find(x) != v.end())
        return v[x];
    int cnt = 0, lstFlag = gb_compareWithZero(gb_getValue(t[0], x), 0);
    for(int i = 1; t[i].size(); i++)
    {
        int nowFlag = gb_compareWithZero(gb_getValue(t[i], x), 0);
        if(nowFlag == 0)
            continue;
        if(nowFlag * lstFlag < 0)
            cnt++;
        lstFlag = nowFlag;
    }
    v[x] = cnt;
    return cnt;
}

int gb_getRootNumber(const vector<double> *t, unordered_map<double, int> &v, double l, double r) // 利用 Sturm Theorem 计算区间内的实根数
    { return gb_ABS(gb_calcSignVar(t, v, l) - gb_calcSignVar(t, v, r)); }

void gb_negPolyMod(const vector<double> &a, const vector<double> &b, vector<double> &res) // 求出多项式 a mod b 再取负的结果，存储在 res 
{
    vector<double> tmp = a;
    int n = a.size(), m = b.size();
    assert(a[n - 1] != 0 && b[m - 1] != 0);
    res.clear();
    for(int i = n - 1; i >= m - 1; i--)
        for(int j = 0; j <= m - 1; j++)
            tmp[i - m + 1 + j] -= tmp[i] * b[j] / b[m - 1]; 
    int deg = 0;
    for(int i = 0; i <= m - 2; i++)
        if(gb_compareWithZero(tmp[i], 1e-40) != 0)
            deg = i;
    for(int i = 0; i <= deg; i++)
        res.push_back(-tmp[i]);
}

CP_Root gb_NewtonMethod(const vector<double> &a, const vector<double> &d, double x0, double e) // 以 x0 为初值在一个单根区间内进行牛顿迭代法
{
    double x = x0;
    double aVal = gb_getValue(a, x), dVal = gb_getValue(d, x);
    if(dVal == 0)
    {
        x = gb_findNotZero(d, x, e, 0, 1); 
        aVal = gb_getValue(a, x);
        dVal = gb_getValue(d, x);
    }
    int cnt = 0;
    while(aVal != 0 && cnt < ITER_LIM)
    {
        cnt++;
        x = x - aVal / dVal;
        aVal = gb_getValue(a, x);
        dVal = gb_getValue(d, x);
        if(aVal == 0)
            return CP_Root(1, 0, x);
        if(dVal == 0)
        {
            x = gb_findNotZero(d, x, e, 0, 1);
            aVal = gb_getValue(a, x);
            dVal = gb_getValue(d, x);
        }
    }
    return CP_Root(1, gb_ABS(aVal / dVal), x);
}
