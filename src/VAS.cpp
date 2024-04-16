/*******************************
* VAS bisection method by panjq22(panjq22@mails.tsinghua.edu.cn)
* Created on 2022/10/02
*******************************/
#include <CP_Root.h>
#include <VAS.h>

mpq mobius::rootlimit;
int rootlimitlog2;

std::ostream & operator << (std::ostream & out, mobius x) {
	out << "(" << x.a.get_d() << "x+" << x.b.get_d() << ")";
	out << "/";
	out << "(" << x.c.get_d() << "x+" << x.d.get_d() << ")";
	return out;
}

mpz polyLowerBound(poly f, int eps) {
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
					f[j] /= 2;
					mpq tv = -f[i] / f[j];
					double v = pow(tv.get_d(), 1. / (j - i));
					ma = std::max(ma, v);
				}
		}
	// -1 to prevent corner case
	return (int) (std::floor(1 / ma) - 1);
}

typedef std::pair<mpq, mpq> interval;

interval operator + (interval a, interval b) {
	return interval(a.first + b.first, a.second + b.second);
}

interval operator * (interval a, interval b) {
	mpq l, r, t;
	l = r = a.first * b.first;
	t = a.second * b.first;
	l = std::min(l, t), r = std::max(r, t);
	t = a.first * b.second;
	l = std::min(l, t), r = std::max(r, t);
	t = a.second * b.second;
	l = std::min(l, t), r = std::max(r, t);
	return interval(l, r);
}

interval round(interval x, int eps) {
	if (x.first < 0) x.first = -gb_approxReduction(-x.first, eps, 1);
	else x.first = gb_approxReduction(x.first, eps, 0);
	if (x.second < 0) x.second = -gb_approxReduction(x.second, eps, 0);
	else x.second = gb_approxReduction(x.second, eps, 1);
	return x;
}
interval intervalValue(poly f, mpq L, mpq R, int eps) {
	const int sz = f.size();
	interval ret(0, 0);
	for (int i = sz - 1; ~i; --i) {
		ret = ret * interval(L, R) + interval(f[i], f[i]);
	}
	return ret;
}

clock_t shift_clock;
unsigned shift_count;
unsigned sign_max;
poly shift(poly f, mpq x, int approx = 0) {
	shift_clock -= clock();
	const int sz = f.size();
	bool can = false;
	if (approx) {
		can = true;
		std::vector<mpf> g;
		g.resize(sz);
		unsigned prec = approx + rootlimitlog2 + sz * sz;
		mpf xt = mpf(x, prec * 2);
		for (int i = 0; i < sz; ++i)
			g[i] = mpf(f[i], prec * 2);
		for (int T = sz - 1; ~T; --T)
			for (int i = T; i + 1 < sz; ++i)
				g[i] += g[i + 1] * xt;
		if (can) {
			for (int i = 0; i < sz; ++i)
				f[i] = g[i];
		}
	}
	if (!can) {
		++shift_count;
		for (int T = sz - 1; ~T; --T)
			for (int i = T; i + 1 < sz; ++i)
				f[i] += f[i + 1] * x;
		if (!approx) {
			for (int i = 0; i < sz; ++i) {
				mpz t = f[i].get_den();
				unsigned sz = mpz_sizeinbase(t.get_mpz_t(), 2);
				sign_max = std::max(sign_max, sz);
			}
		}
	}
	shift_clock += clock();
	return f;
}

static int totalNewton;
mpq intervalNewton(poly f, mpq L, mpq R, int eps) {
	mpq limit(1);
	limit >>= eps;
	poly df = derivative(f);
	L = gb_approxReduction(L, eps, 0);
	R = gb_approxReduction(R, eps, 1);
	int var = variate(shift(f, L, eps)) - variate(shift(f, R, eps));
	if (!var) return mobius::rootlimit;
	while ((R - L) > limit) {
		++totalNewton;
		L = gb_approxReduction(L, eps, 0);
		R = gb_approxReduction(R, eps, 1);
		mpq mid = (L + R) / 2;
		mpq val = calculate(f, mid);
		if (val == 0) return mid;
		interval it = intervalValue(df, L, R, eps);
		if (it.first == 0) it.first -= 1;
		if (it.second == 0) it.second += 1;
		if (it.first < 0 && it.second > 0) {
			// for lhs
			// - val / -0
			if (val > 0) std::swap(it.first, it.second);
			mpq tr = std::min(R, mpq(mid - val / it.first));
			mpq tl = std::max(L, mpq(mid - val / it.second));
			if (L <= tr) {
				mpq t = intervalNewton(f, L, tr, eps);
				if (t != mobius::rootlimit) return t;
			}
			if (tl <= R) {
				mpq t = intervalNewton(f, tl, R, eps);
				if (t != mobius::rootlimit) return t;
			}
			return mobius::rootlimit;
		}
		it.first = mid - val / it.first;
		it.second = mid - val / it.second;
		if (it.first > it.second)
			std::swap(it.first, it.second);
		// no solution
		if (it.second < L || R < it.first) {
			return mobius::rootlimit;
		}
		L = std::max(L, it.first);
		R = std::min(R, it.second);
	}
	return L;
}

interval VAS_bisect(poly f, mobius M, int eps, bool left = false, bool right = false) {
	int var = variate(f);
	
	if (var == 0) return interval(1, 0);

	if (1 || (left && right)) {
		mpq a = M.zero();
		mpq b = M.infinity();
		if (a > b) std::swap(a, b);
		return interval(a, b);
	}

	mpz lb = polyLowerBound(f, eps);

	if (lb > 1) {
		f = shift(f, lb);
		M = M(mobius(1, lb, 0, 1));
		return VAS_bisect(f, M, eps, true, right);
	}

	if (calculate(f, 1) == 0) {
		return interval(M(1), M(1));
	}

	interval t = VAS_bisect(shift(reverse(f), 1), M(mobius(0, 1, 1, 1)), eps, true, left);
	if (t.first <= t.second) return t;
	t = VAS_bisect(shift(f, 1), M(mobius(1, 1, 0, 1)), eps, true, right);
	return t;
}

poly VAS_algorithm(poly init_f, int eps) {
	// isolation
	typedef std::pair<poly, mobius> qt;
	std::queue<qt> que;
	std::vector<interval> list;
	if (init_f[0] == 0) {
		init_f = divbyx(init_f);
		list.emplace_back(0, 0);
	}
	mobius::rootlimit = simpleUpperBound(init_f);
	mpz tdiv = mobius::rootlimit.get_num()
				/ mobius::rootlimit.get_den();
	rootlimitlog2 = mpz_sizeinbase(tdiv.get_mpz_t(), 2) + 1;
	que.emplace(init_f, mobius(1, 0, 0, 1));
	while (!que.empty()) {
		poly f; mobius M;
		std::tie(f, M) = que.front(); que.pop();
		
		int var = variate(f);

		if (var == 0) continue;
		if (var == 1) {
			interval it = VAS_bisect(f, M, eps);
			list.push_back(it);
			continue;
		}

		mpz lb = polyLowerBound(f, eps);
		// move root to [0, 1)
		if (lb > 1) {
			f = shift(f, lb);
			M = M(mobius(1, lb, 0, 1));
		}
		if (calculate(f, 1) == 0) {
			list.emplace_back(M(1), M(1));
		}
		que.emplace(shift(reverse(f), 1), M(mobius(0, 1, 1, 1)));
		que.emplace(shift(f, 1), M(mobius(1, 1, 0, 1)));
	}
	
	poly ret; ret.resize(list.size());
	for (int i = 0; i < (int) list.size(); ++i) {
		// outside or inside
		// which is faster?
		mpq l = list[i].first;
		mpq r = list[i].second;
		if (l > r) std::swap(l, r);
		if (l == r) ret[i] = l;
		else {
			ret[i] = intervalNewton(init_f, l, r, eps);
		}
	}
	return ret;
}

std::vector<root_t> VAS_Algorithm(poly f, int eps) {
	std::vector<poly> dec = squareDecomposition(f);
	std::vector<root_t> ret;
	for (int T = 0; T < (int) dec.size(); ++T) {
		poly res = VAS_algorithm(dec[T], eps);
		// found (T+1)-multiple root t
		for (auto t : res) {
			ret.emplace_back(mpq_to_double(t), T + 1);
		}
		// std::cin.get();
		res = VAS_algorithm(negative(dec[T]), eps);
		for (auto t : res) {
			if (t == 0) continue;
			ret.emplace_back(mpq_to_double(-t), T + 1);
		}
	}
	std::sort(ret.begin(), ret.end());
	return ret;
}

void gb_solve4(const std::vector<double> & f, double e, int & flag, std::vector<CP_Root> & r) {
	totalNewton = 0;
	flag = 1;
	poly nf;
	nf.resize(f.size());
	for (int i = 0; i < (int) f.size(); ++i)
		nf[i] = f[i];
	auto ret = VAS_Algorithm(nf);
	r.resize(ret.size());
	for (int i = 0; i < (int) ret.size(); ++i) {
		r[i].m_repetitions = ret[i].second;
		r[i].m_root = ret[i].first.get_d();
		r[i].m_errorBound = e;
	}
}
