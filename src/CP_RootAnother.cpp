#include <CP_Root.h>
#include <Bisection.h>
// Finding the roots of n-degree polynomials
// by 叶隽希 20220925

poly gb_scale(const poly &a) {
    poly b(a);
    int n = b.size() - 1;
    mpz d = 1;
    for (int i = n; ~i; i--) {
        b[i] *= d;
        d *= 2;
    }
    return b;
}
poly gb_inversion(const poly &a) {
    poly b(a);
    std::reverse(b.begin(), b.end());
    return b;
}
poly gb_translation(const poly &a, mpq r = 1) {
    int n = a.size() - 1;
    poly b; b.resize(n + 1);
    for (int i = n; ~i; i--) {
        b[i] = a[i];
        for (int j = i; j < n; j++)
            b[j] += b[j + 1] * r;
    }
    return b;
}
poly gb_shift(const poly &a) {
    int n = a.size() - 1;
    poly b; b.reserve(n);
    for (int i = 1; i <= n; i++)
        b.emplace_back(a[i]);
    return b;
}
int gb_variation(const poly& a) {
    int n = a.size() - 1;
    int ans = 0, s1 = sgn(a[0]);
    for (int i = 1; i <= n; i++) {
        int s2 = sgn(a[i]);
        if (s1 && s2 && s1 != s2) s1 = s2, ans++;
        else if (!s1) s1 = s2;
    }
    return ans;
}
void print(const poly &a) {
    for (int n = a.size() - 1, i = 0; i <= n; i++)
        std::cerr << a[i] << "\n "[i != n];
}
poly gb_derivative(const poly &a) {
    poly d;
    int n = a.size() - 1;
    for (int i = 1; i <= n; i++) {
        d.push_back(a[i] * i);
    }
    return d;
}

mpq gb_newton2(const poly &a, const poly &d, mpq x0) {
    for (int i = 0; i < iterationLimit; i++) {
        x0 = gb_approxReduction(x0, pErrorBound, 0);
        x0 -= calculate(a, x0) / calculate(d, x0);
    }
    return x0;
}
void gb_bisection(const poly &a, std::vector<rootInterval>& ans) {
    std::queue<interval> q;
    q.push({0, 0, a});
    ans.clear();
    while (!q.empty()) {
        interval u = q.front(); q.pop();
        while (!u.p.empty() && u.p[0] == 0) {
            u.p = gb_shift(u.p);
            ans.emplace_back(rootInterval(u.c, u.k, 0));
        }
        if (u.p.empty()) continue;
        int v = gb_variation(gb_translation(gb_inversion(u.p)));
        if (v == 0) continue;
        if (v == 1 && u.k >= pErrorBound) ans.emplace_back(rootInterval(u.c, u.k, 1));
        else {
            q.push({u.c * 2, u.k + 1, gb_scale(u.p)});
            q.push({u.c * 2 + 1, u.k + 1, gb_translation(gb_scale(u.p))});
        }
    }
}

inline mpq norm(const mpq& x, const mpq& bound) {
    return bound * (x * 2 - 1);
}

void _gb_solve2(const poly & a, double e, int& flag, std::vector<CP_Root> &r) {
    r.clear();
    int n = a.size() - 1;
    mpq an = mpq(a[n]), bound = abs(a[0] / an);
    for (int i = 1; i <= n; i++) {
        mpq val = abs(a[i] / an);
        bound = std::max(val, bound);
    }
    bound += 1;
    pErrorBound = (log(bound.get_d()) - log(e)) / log(2) * 2;
    
    poly b, der;
    mpq tmp = 1;
    for (int i = 0; i <= n; i++) {
        b.push_back(a[i] * tmp);
        tmp *= bound * 2;
    }
    b = gb_translation(b, -0.5);
    der = gb_derivative(b);

    std::vector<rootInterval> res;
    gb_bisection(b, res);

    int m = res.size();
    mpz d = 1; int k = 0;
    for (int i = 0; i < m; i++) {
        for (; k < res[i].k; k++)
            d *= 2;
        if (res[i].h == 0) {
            r.emplace_back(CP_Root(1, 0, mpq_to_double(2 * bound * (mpq(res[i].c, d)) - bound)));
        }
        else {
            mpq L = mpq(res[i].c, d), R = mpq(res[i].c + 1, d);
            mpq x0 = gb_newton2(b, der, (L + R) / 2);
            r.emplace_back(CP_Root(1, 0, mpq_to_double(norm(x0, bound))));
        }
    }
}


void gb_solve2(const std::vector<double>& a, double e, int& flag, std::vector<CP_Root> &r) {
	poly nf;
	nf.resize(a.size());
	for (int i = 0; i < (int) a.size(); ++i)
		nf[i] = a[i];
	r.clear();
	std::vector<poly> dec = squareDecomposition(nf);
	if (dec.size() == 1) dec[0] = nf;
	for (int T = 0; T < (int) dec.size(); ++T) {
		std::vector<CP_Root> d;
		_gb_solve2(dec[T], e, flag, d);
		for (auto t : d) {
			t.m_repetitions = T + 1;
			r.push_back(t);
		}
	}
	std::sort(r.begin(), r.end());
}
