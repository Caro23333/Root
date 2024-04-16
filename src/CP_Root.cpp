#include <CP_Root.h>

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

// 将ULL数的内存解析为浮点数
double gb_toDouble(ULL x)
{
    void *a = &x;
    return *((double *)a);
} // 函数gb_toDouble结束

// 将浮点数的内存解析为ULL数
ULL gb_toULL(double x)
{
    void *a = &x;
    return *((ULL *)a);
} // 函数gb_toULL结束

double mpq_to_double(const mpq &x) {
    if (sgn(x) == 0) return 0;
    double b = x.get_d(), a = sgn(x) > 0 ? gb_getNext(b) : gb_getPrevious(b);
    return abs(x - b) < abs(x - a) ? b : a;
}

mpq gb_approxReduction(mpq x, int eps, bool greater) 
{
    if(x == 0)
        return mpq(0);
    bool flag = x < 0;
    mpq y = abs(x);
    mpz den = y.get_den();
	mpz td = den / (y.get_num() / y.get_den() + 1);
	int dig = mpz_sizeinbase(td.get_mpz_t(), 2);
	int mov = dig - eps - 3;
	if (mov > 0) {
		mpz num = y.get_num();
		num >>= mov;
		den >>= mov;
		if (greater ^ flag) ++num;
		else ++den;
		return (flag ? -1 : 1) * mpq(num, den);
	} else return x;
}

std::ostream & operator << (std::ostream & out, poly x) {
	for (int i = 0; i < (int) x.size(); ++i)
		out << x[i].get_d() << ',';
	return out;
}

poly reduce(poly f) {
	while (f.size() > 1 && f.back() == 0) f.pop_back();
	return f;
}

poly operator + (poly f, poly g) {
	if (f.size() < g.size()) std::swap(f, g);
	for (int i = 0; i < (int) g.size(); ++i) f[i] += g[i];
	return reduce(f);
}

poly operator - (poly f, poly g) {
	for (int i = 0; i < (int) g.size(); ++i) g[i] = -g[i];
	return reduce(f + g);
}
poly derivative(poly f) {
	if (f.size() == 1) {
		f[0] = 0;
		return f;
	}
	const int sz = f.size();
	for (int i = 0; i + 1 < sz; ++i)
		f[i] = f[i + 1] * (i + 1);
	f.pop_back();
	return f;
}

poly negative(poly f) {
	const int sz = f.size();
	for (int i = sz - 2; i >= 0; i -= 2)
		f[i] = -f[i];
	return f;
}

poly opposite(poly f) {
    const int sz = f.size();
    for (int i = sz - 1; i >= 0; i--)
        f[i] = -f[i];
    return f;
}

poly mod(poly f, poly g, poly * d) {
	const int sf = f.size(), sg = g.size();
	if (d != NULL) {
		if (sf >= sg)
			d->resize(sf + 1 - sg);
		else {
			d->resize(1);
			d->operator[](0) = 0;
		}
	}
	for (int i = sf - 1; i >= sg - 1; --i) {
		if (f[i] != 0) {
			mpq v = f[i] / g[sg - 1];
			if (d != NULL)
				d->operator[](i - sg + 1) = v;
			for (int j = sg - 1, k = i; ~j; --j, --k)
				f[k] -= v * g[j];
		}
		f.pop_back();
	}
	if (sg == 1) {
		return poly({mpq(0)});
	}
	return reduce(f);
}

// may be very low-efficienct
poly div(poly f, poly g) {
	poly ret;
	mod(f, g, &ret);
	return reduce(ret);
}

poly gcd(poly f, poly g) {
	while (!(g.size() == 1 && g[0] == 0)) {
		poly t = mod(f, g);
		f = g; g = t;
	}
	return f;
}

int variate(poly f) {
	const int sz = f.size();
	int sgn = f.back() > 0, cnt = 0;
	for (int i = sz - 1; ~i; --i) {
		if (f[i] == 0) continue;
		int sg = f[i] > 0;
		if (sg != sgn) {
			sgn = sg;
			++cnt;
		}
	}
	return cnt;
}

mpq calculate(poly f, mpq x) {
	mpq ret = 0;
	if (x == 0) return f[0];
	const int sz = f.size();
	if (x == 1) {
		for (int i = 0; i < sz; ++i)
			ret += f[i];
		return ret;
	}
	for (int i = sz - 1; ~i; --i)
		ret = ret * x + f[i];
	return ret;
}

poly reverse(poly f) {
	std::reverse(f.begin(), f.end());
	return reduce(f);
}

poly divbyx(poly f) {
	const int sz = f.size();
	for (int i = 0; i + 1 < sz; ++i)
		f[i] = f[i + 1];
	f.pop_back();
	return f;
}

// using Yun's Algorithm to decompose poly
// result[i] means polynomial constructed by (i+1)-multiple roots
std::vector<poly> squareDecomposition(poly f) {
	std::vector<poly> result;
	poly a, b, c, d;
	b = f, d = derivative(f);
	bool fir = true;
	while (b.size() > 1) {
		a = gcd(b, d);
		if (fir) fir = false;
		else result.push_back(a);
		b = div(b, a);
		c = div(d, a);
		d = c - derivative(b);
	}
	return result;
}

mpq simpleUpperBound(poly f) {
	const int sz = f.size();
	mpq ret(0);
	for (int i = 0; i + 1 < sz; ++i)
    {
        mpq val = abs(f[i] / f[sz - 1]);
		ret = std::max(val, ret);
    }
	return ret + 1;
}