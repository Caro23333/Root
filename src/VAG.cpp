/*******************************
* VAG bisection method by panjq22(panjq22@mails.tsinghua.edu.cn)
* Created on 2022/10/02
*******************************/
#include "VAG.h"
#include <queue>
#include <cmath>
#include <utility>
#include <tuple>
#include <algorithm>
#include <iostream>
#include <cassert>

typedef std::vector<mpq_class> poly;

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

poly mod(poly f, poly g, poly * d = NULL) {
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
			mpq_class v = f[i] / g[sg - 1];
			if (d != NULL)
				d->operator[](i - sg + 1) = v;
			for (int j = sg - 1, k = i; ~j; --j, --k)
				f[k] -= v * g[j];
		}
		f.pop_back();
	}
	if (sg == 1) {
		return poly({mpq_class(0)});
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

mpq_class calculate(poly f, mpq_class x) {
	mpq_class ret = 0;
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

poly shift(poly f, mpq_class x) {
	const int sz = f.size();
	for (int T = sz - 1; ~T; --T)
		for (int i = T; i + 1 < sz; ++i)
			f[i] += f[i + 1] * x;
	return f;
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
		// std::cout << "?? " << a << std::endl;
		if (fir) fir = false;
		else result.push_back(a);
		b = div(b, a);
		c = div(d, a);
		d = c - derivative(b);
		// std::cout << "=> " << b << ';' << c << ';' << d << std::endl;
	}
	return result;
}

struct mobius {
	static mpq_class rootlimit;
	mpz_class a, b, c, d;
	mobius() {}
	mobius(mpz_class _a, mpz_class _b, mpz_class _c, mpz_class _d) {
		a = _a, b = _b, c = _c, d = _d;
	}
	mobius operator() (mobius t) {
		return mobius(
				a * t.a + b * t.c, a * t.b + b * t.d,
				c * t.a + d * t.c, c * t.b + d * t.d);
	}
	mpq_class zero() {
		if (d == 0) return rootlimit;
		return mpq_class(b, d);
	}
	mpq_class infinity() {
		if (c == 0) return rootlimit;
		return mpq_class(a, c);
	}
	mpq_class operator() (mpz_class x) {
		return mpq_class(a * x + b, c * x + d);
	}
} ;
mpq_class mobius::rootlimit;

std::ostream & operator << (std::ostream & out, mobius x) {
	out << "(" << x.a.get_d() << "x+" << x.b.get_d() << ")";
	std::cout << "/";
	out << "(" << x.c.get_d() << "x+" << x.d.get_d() << ")";
	return out;
}

mpq_class simpleUpperBound(poly f) {
	const int sz = f.size();
	mpq_class ret(0);
	for (int i = 0; i + 1 < sz; ++i)
		ret += abs(f[i]);
	ret /= abs(f[sz - 1]);
	if (ret < 1) ret = 1;
	return ret + 1;
}
mpz_class polyLowerBound(poly f, int eps) {
	f = reverse(f);
	const int sz = f.size();
	if (f.back() < 0) {
		for (int i = 0; i < sz; ++i)
			f[i] = -f[i];
	}
	// local max quadratic
	// low-precision may cause low-efficiency
	double ma = -1;
	for (int i = 0; i < sz; ++i)
		if (f[i] < 0) {
			for (int j = i + 1; j < sz; ++j)
				if (f[j] > 0) {
					mpq_class tv = -f[i] / f[j];
					f[j] /= 2;
					double v = pow(tv.get_d(), 1. / (j - i));
					ma = std::max(ma, v);
				}
		}
	// -1 to prevent corner case
	return (int) (std::floor(1 / ma) - 1);
}

typedef std::pair<mpq_class, mpq_class> interval;

interval operator + (interval a, interval b) {
	return interval(a.first + b.first, a.second + b.second);
}

interval operator * (interval a, interval b) {
	mpq_class l, r, t;
	l = r = a.first * b.first;
	t = a.second * b.first;
	l = std::min(l, t), r = std::max(r, t);
	t = a.first * b.second;
	l = std::min(l, t), r = std::max(r, t);
	t = a.second * b.second;
	l = std::min(l, t), r = std::max(r, t);
	return interval(l, r);
}

interval intervalValue(poly f, mpq_class L, mpq_class R) {
//	if (L < 0 && R > 0) {
//		interval a = intervalValue(f, L, 0);
//		interval b = intervalValue(f, 0, R);
//		if (a.first > b.first) std::swap(a, b);
//		if (a.second < b.first) {
//			std::cout << "???" << std::endl;
//			abort();
//		}
//		return interval(a.first, b.second);
//	}
	const int sz = f.size();
	interval ret(0, 0);
	for (int i = sz - 1; ~i; --i) {
	//	std::cout << "???: " << L.get_d() << ' ' << R.get_d() << '\n';
	//	std::cout << "???: " << mpq_class(L * ret.first).get_d() << ' ' << mpq_class(R * ret.first).get_d() << '\n';
	//	std::cout << "???: " << (ret * interval(L, R)).second.get_d() << '\n';
		ret = ret * interval(L, R) + interval(f[i], f[i]);
	//	std::cout << f[i].get_d() << " : " << ret.first.get_d() << ' ' << ret.second.get_d() << '\n';
	}
	return ret;
}

mpq_class round(mpq_class x, int eps, int greater) {
	mpz_class den = x.get_den();
	int dig = mpz_sizeinbase(den.get_mpz_t(), 2);
	int mov = dig - eps - 2;
	if (mov > 0) {
		mpz_class num = x.get_num();
		num >>= mov;
		den >>= mov;
		if (greater) ++num;
		else ++den;
		return mpq_class(num, den);
	} else return x;
}

mpq_class intervalNewton(poly f, mpq_class L, mpq_class R, int eps) {
	mpq_class limit(1);
	limit >>= eps;
	poly df = derivative(f);
//	std::cout << '\n' << f << '\n';
	// std::cout << "----------NEWTON " << L.get_d() << ' ' << R.get_d() << '\n';
	while ((R - L) > limit) {
		L = round(L, eps, 0);
		R = round(R, eps, 1);
		mpq_class mid = (L + R) / 2;
		mpq_class val = calculate(f, mid);
		if (val == 0) return mid;
		// std::cout << "NOW " << L.get_d() << ' ' << R.get_d() << '\n';
		interval it = intervalValue(df, L, R);
		// std::cerr << "GETFIR " << it.first.get_d() << ' ' << it.second.get_d() << std::endl;
		if (it.first == 0) it.first -= 1;
		if (it.second == 0) it.second += 1;
		if (it.first < 0 && it.second > 0) {
			// for lhs
			// - val / -0
			if (val > 0) std::swap(it.first, it.second);
			mpq_class tr = std::min(R, mpq_class(mid - val / it.first));
			mpq_class tl = std::max(L, mpq_class(mid - val / it.second));
			// std::cout << "GENE " << "(" << L.get_d() << "," << tr.get_d() << ") (" << tl.get_d() << "," << R.get_d() << ")" << '\n';
			if (L <= tr) {
				mpq_class t = intervalNewton(f, L, tr, eps);
				if (t != mobius::rootlimit) return t;
			}
			if (tl <= R) {
				mpq_class t = intervalNewton(f, tl, R, eps);
				if (t != mobius::rootlimit) return t;
			}
			// std::cout << "NOTFOUND" << std::endl;
			return mobius::rootlimit;
		}
		it.first = mid - val / it.first;
		it.second = mid - val / it.second;
		if (it.first > it.second)
			std::swap(it.first, it.second);
		// std::cerr << "GET " << it.first.get_d() << ' ' << it.second.get_d() << std::endl;
		// no solution
		if (it.second < L || R < it.first)
			return mobius::rootlimit;
		L = std::max(L, it.first);
		R = std::min(R, it.second);
	}
	return L;
}

interval VAG_bisect(poly f, mobius M, int eps, bool left = false, bool right = false) {
	int var = variate(f);
	
	if (var == 0) return interval(1, 0);

	if (left && right) {
		mpq_class a = M.zero();
		mpq_class b = M.infinity();
		if (a > b) std::swap(a, b);
		return interval(a, b);
	}

	mpz_class lb = polyLowerBound(f, eps);

	if (lb > 1) {
		f = shift(f, lb);
		M = M(mobius(1, lb, 0, 1));
		return VAG_bisect(f, M, eps, true, right);
	}

	if (calculate(f, 1) == 0) {
		return interval(M(1), M(1));
	}

	interval t = VAG_bisect(shift(reverse(f), 1), M(mobius(0, 1, 1, 1)), eps, true, left);
	if (t.first <= t.second) return t;
	t = VAG_bisect(shift(f, 1), M(mobius(1, 1, 0, 1)), eps, true, right);
	return t;
}

poly VAG_algorithm(poly init_f, int eps) {
	// std::cout << "VAG ALGORITHM" << std::endl;
	// std::cout << init_f << " : " << eps << std::endl;
	// isolation
	typedef std::pair<poly, mobius> qt;
	std::queue<qt> que;
	std::vector<interval> list;
	if (init_f[0] == 0) {
		init_f = divbyx(init_f);
		list.emplace_back(0, 0);
	}
	mobius::rootlimit = simpleUpperBound(init_f);
	que.emplace(init_f, mobius(1, 0, 0, 1));
	while (!que.empty()) {
		poly f; mobius M;
		std::tie(f, M) = que.front(); que.pop();
		// std::cout << "\nGETF " << M << " (" << M.zero().get_d() << ", " << M.infinity().get_d() << ") \n" << f << std::endl;
		// std::cin.get();

		/*
		// divided by (x - 0)
		if (f[0] == 0) {
			f = divbyx(f);
			list.emplace_back(M.zero(), M.zero());
			std::cout << "FOUND " << M.zero().get_d() << std::endl;
		}
		*/
		
		int var = variate(f);
		// std::cout << "GETVAR " << var << std::endl;

		if (var == 0) continue;
		if (var == 1) {
			// std::cout << "BISECT " << M.zero() << ' ' << M.infinity() << '\n';
			interval it = VAG_bisect(f, M, eps);
			list.push_back(it);
			// std::cout << "FOUND " << it.first.get_d() << ' ' << it.second.get_d() << std::endl;
			continue;
		}

		mpz_class lb = polyLowerBound(f, eps);

		// std::cout << "GETLB " << lb << std::endl;
		// move root to [0, 1)
		if (lb > 1) {
			f = shift(f, lb);
			M = M(mobius(1, lb, 0, 1));
		}
		
		// check if 1 is root
		if (calculate(f, 1) == 0) {
			list.emplace_back(M(1), M(1));
			// std::cout << "FOUND " << M(1).get_d() << std::endl;
		}

		// get isolation in (0, 1)
		que.emplace(shift(reverse(f), 1), M(mobius(0, 1, 1, 1)));
		// get isolation in (1, +inf)
		que.emplace(shift(f, 1), M(mobius(1, 1, 0, 1)));
	}
	std::cout << "ISOLATED" << std::endl;
//	for (int i = 0; i < (int) list.size(); ++i) {
//		std::cout << list[i].first.get_d() << ' ' << list[i].second.get_d() << "\n";
//	}
	
	poly ret; ret.resize(list.size());
	for (int i = 0; i < (int) list.size(); ++i) {
		// outside or inside
		// which is faster?
		mpq_class l = list[i].first;
		mpq_class r = list[i].second;
		if (l > r) std::swap(l, r);
		if (l == r) ret[i] = l;
		else {
			ret[i] = intervalNewton(init_f, l, r, eps);
		}
		std::cout << ret[i].get_d() << std::endl;
	}
	return ret;
}

std::vector<root_t> solve_VAG(poly f, int eps) {
	std::vector<poly> dec = squareDecomposition(f);
	// std::cout << "square" << std::endl;
	// for (auto t : dec) {
	// 	std::cout << t << '\n';
	// }
	std::vector<root_t> ret;
	for (int T = 0; T < (int) dec.size(); ++T) {
		poly res = VAG_algorithm(dec[T], eps);
		// found (T+1)-multiple root t
		for (auto t : res) {
			ret.emplace_back(t.get_d(), T + 1);
		}
		// std::cin.get();
		res = VAG_algorithm(negative(dec[T]), eps);
		for (auto t : res) {
			if (t == 0) continue;
			ret.emplace_back(-t.get_d(), T + 1);
		}
	}
	std::sort(ret.begin(), ret.end());
	std::cout << "FINAL:" << std::endl;
	std::cout << "RET: " << ret.size() << '\n';
	for (auto t : ret)
		std::cout << t.first.get_d() << ' ' << t.second << '\n';
	return ret;
}
