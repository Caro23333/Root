#include <iostream>
#include <fstream>
#include <queue>
#include <cstdlib>
#include <algorithm>
#include <unordered_map>
#include <vector>
#include "CP_Root.h"
#include "testlib.h"
#include <gmpxx.h>
#include <ctime>
#include <cmath>

// Sturm - Newton Method with High-Precision Rational Number
// by 吕敬一 on 20221009
using namespace std;
typedef mpq_class mpq;
typedef mpz_class mpz;
typedef vector<mpq_class> poly;

mpq gb_approxReduction(mpq x, int eps, int dir) 
{
    int flag = gb_compareWithZero(x, mpq_class(0));
    if(!flag)
        return x;
    x = abs(x);
	mpz den = x.get_den();
	mpz td = den / (x.get_num() / x.get_den() + 1);
	int dig = mpz_sizeinbase(td.get_mpz_t(), 2);
	int mov = dig - eps - 3;
	if(mov > 0)
    {
		mpz num = x.get_num();
		num >>= mov;
		den >>= mov;
		if(dir * flag == 1) 
            ++num;
		else 
            ++den;
		return flag * mpq_class(num, den);
	} 
    else 
        return flag * x;
}

mpq gb_getValue(const poly &a, mpq x, int eps)
{
    mpq res = mpq_class(0);
    int n = a.size() - 1;
    for(int i = n; i >= 0; i--)
    {
        res = res * x + a[i];
        if(eps != 0)
            res = gb_approxReduction(res, eps, gb_compareWithZero(res, mpq_class(0)));
    }
  //  if(eps == 0 && x >= -73 && x <= -71)
    //    cout << "getValue " << x << " " << res << endl;
    return res;
}

poly gb_getMod(const poly &a, const poly &b)
{
    poly res, tmp;
    res.clear();
    int n = a.size() - 1, m = b.size() - 1;
    for(int i = 0; i <= n; i++)
        tmp.push_back(a[i]);
    for(int i = n; i >= m; i--)
    {
        mpq rt = tmp[i] / b[m];
        for(int j = i - m; j <= i; j++)
            tmp[j] -= rt * b[j - i + m];
    }
    for(int i = 0; i <= m - 1; i++)
        res.push_back(0 - tmp[i]);
    for(int i = m - 1; i >= 0; i--)
    {
        if(res[i] == 0)
            res.pop_back();
        else
            break;
    }
    return res;
}

poly gb_getDiff(const poly &a)
{
    poly res;
    res.clear();
    int n = a.size() - 1;
    for(int i = 0; i <= n - 1; i++)
        res.push_back(a[i + 1] * (i + 1));
    return res;
}

int gb_compareWithZero(mpq x, mpq eps)
{
    if(abs(x) <= eps)
        return 0;
    else if(x > eps)
        return 1;
    else
        return -1;
}

bool gb_judge(mpq aVal, mpq dVal, mpq e)
    { return (dVal == 0 && aVal == 0) || (dVal != 0 && gb_compareWithZero(aVal / dVal, e) == 0); }

int gb_getRepetition(const vector<poly> &da, mpq x, mpq e)
{
    int n = da[0].size();
    mpq nowVal = gb_getValue(da[1], x), nxtVal;
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
            if(!gb_judge(nowVal, nxtVal, e))
                return i;
        }
    }
    return n;
}

mpq gb_findNotZero(const poly &a, mpq x, mpq e, mpq eps, int dir)
{
    mpq val = gb_getValue(a, x), res = x;
    int flag = gb_compareWithZero(val, eps);
    while(flag != 0)
    {
        res += dir * e;
        val = gb_getValue(a, res);
        flag = gb_compareWithZero(val, eps);
    }
    return res;
}

int gb_calcVar(const vector<poly> &s, mpq x, int eps)
{
    int lstFlag = gb_compareWithZero(gb_getValue(s[0], x), mpq_class(0));
    int m = s.size(), res = 0;
    for(int i = 1; i <= m - 1; i++)
    {
        int nowFlag = gb_compareWithZero(gb_getValue(s[i], x), mpq_class(0));
        if(lstFlag != 0 && nowFlag * lstFlag < 0)
            res++;
        if(nowFlag != 0)
            lstFlag = nowFlag;
    }
    return res;
}

int gb_getRootNumber(const vector<poly> &s, mpq l, mpq r, int eps)
    { return gb_calcVar(s, l, eps) - gb_calcVar(s, r, eps); }

mpq gb_newtonMethod(const poly &a, const poly &d, mpq x0, mpq e, int eps)
{
    mpq x = x0;
    mpq aVal = gb_getValue(a, x0), dVal = gb_getValue(d, x0);
    if(dVal == 0)
    {
        x = gb_findNotZero(d, x, e, mpq_class(0), 1);
        aVal = gb_getValue(a, x0), dVal = gb_getValue(d, x0);
    }
    mpq del = aVal / dVal;
    int cnt = 0;
    while(cnt < ITER_LIM)
    {
        cnt++;
        x0 = x;
        x = x - del;
        aVal = gb_getValue(a, x), dVal = gb_getValue(d, x);
        if(dVal == 0)
        {
            x = gb_findNotZero(d, x, e, mpq_class(0), 1);
            aVal = gb_getValue(a, x), dVal = gb_getValue(d, x);
        }
        if(aVal == 0)
            return x;
        del = aVal / dVal;
        if(gb_compareWithZero(abs(del), e / 10) == 0)
            return x;
        if(del < 0)
            x = gb_approxReduction(x, eps * 2, 1), del = gb_approxReduction(del, eps * 2, -1);
        else
            x = gb_approxReduction(x, eps * 2, -1), del = gb_approxReduction(del, eps * 2, 1);
    }
    return x;
}

mpq gb_getError(const poly &a, const poly &d, mpq x)
{
    mpq aVal = gb_getValue(a, x), dVal = gb_getValue(d, x);
    if(dVal == 0)
        return aVal;
    else
        return abs(aVal / dVal);
}

void gb_solveSqrFree(poly &p, double e, vector<CP_Root> &r, int k) // k - 重根数
{
    cerr << "solve " << k << endl;
    int n = p.size() - 1;
    const double maxLen = 0.1; // 最大迭代初始区间长度
    mpq re = 1;
    int epsk = 0;
    while(re >= e / 4)
        re /= 2, epsk++;
    for(int i = 0; i <= n; i++)
        p[i] /= p[n];
    mpq bound = abs(p[0]);
    for(int i = 1; i <= n; i++) 
    {
        mpq val = abs(p[i]);
        if (val > bound)   
            bound = val;
    }
    bound += 1;
    double dBound = 1;
    while(dBound <= bound)
        dBound *= 2;
    //cout << bound.get_d() << endl;
    // 求出多项式的若干阶导数
    vector<poly> dp;
    dp.push_back(p);
    for(int i = 1; i <= n; i++)
        dp.push_back(gb_getDiff(dp[i - 1]));
    // 求出多项式的 Sturm 序列
    vector<poly> st;
    st.push_back(p);
    st.push_back(dp[1]);
    while(st[st.size() - 1].size() > 1)
        st.push_back(gb_getMod(st[st.size() - 2], st[st.size() - 1]));
    // 其他准备工作
    queue<CP_Interval> q;
    q.push(CP_Interval(-dBound, dBound));
    vector<mpq_class> res;
    // 开始求解
    while(!q.empty())
    {
        CP_Interval now = q.front();
        q.pop();
        mpq lVal = gb_getValue(p, mpq_class(now.l)), dlVal = gb_getValue(dp[1], mpq_class(now.l));
        mpq rVal = gb_getValue(p, mpq_class(now.r)), drVal = gb_getValue(dp[1], mpq_class(now.r));
        if(now.l >= now.r)
            continue;
        if(gb_judge(lVal, dlVal, 0))
        {
            res.push_back(mpq_class(now.l));
            now.l += max(re.get_d(), 1e-14);
         /////   cerr << "testLeft " << now.l << " " << mpq_class(now.l) << endl;
        }
        if(gb_judge(rVal, drVal, 0))
        {
            res.push_back(mpq_class(now.r));
            now.r -= max(re.get_d(), 1e-14);
          //  cerr << "testRight " << now.r << " " << endl;
        }
        if(now.l >= now.r)
            continue;
        int rootCnt = gb_getRootNumber(st, mpq_class(now.l), mpq_class(now.r), epsk);
       // cerr << "test " << now.l << " " << now.r << " " << rootCnt << endl;
        double mid = (now.l + now.r) / 2.0;
        if(rootCnt == 0)
            continue;
        else if(rootCnt >= 2 || (rootCnt == 1 && now.r - now.l > maxLen))
        {
            q.push(CP_Interval(now.l, mid));
            q.push(CP_Interval(mid, now.r));
        }
        else
        {
            res.push_back(gb_newtonMethod(p, dp[1], mpq_class(mid), re, epsk));
            cerr << "newton done" << endl;
            //system("pause");
        }
    }
    
    sort(res.begin(), res.end());
    vector<mpq> :: iterator ed = unique(res.begin(), res.end());
    vector<mpq> :: iterator it;
    // 将根加入答案列表
    for(it = res.begin(); it != ed; it++)
        r.push_back(CP_Root(k, 0, to_double(*it)));
    res.clear();
}

void gb_solve3(const vector<double> &a, double e, int &flag, vector<CP_Root> &r)
{
    // 准备工作，转换为 mpq_class 类型，进行平方因子分解
    int n = a.size() - 1;
    if(a[n] == 0 || e <= 0)
    {
        flag = -1;
        return;
    }
    poly f(n + 1);
    for(int i = 0; i <= n; i++)
        f[i] = mpq_class(a[i]);
    vector<poly> dec = squareDecomposition(f);
    cout << "square " << dec.size() << std::endl;
	for (auto t : dec) {
		for (int i = 0; i < (int) t.size(); ++i)
		    cerr << to_double(t[i]) << ',';
	}
    //dec.push_back(f);
    r.clear();

    // 依次求根
    for(int i = 0; i < dec.size(); i++)
        gb_solveSqrFree(dec[i], e, r, i + 1);

    //  后处理
    sort(r.begin(), r.end());
    if(r.size() == 0)
        flag = 0;
    else   
        flag = 1;
}
