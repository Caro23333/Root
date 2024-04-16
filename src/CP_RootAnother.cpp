#include <CP_Root.h>
#include <gmpxx.h>

#include <algorithm>
#include <iostream>
#include <vector>
#include <cmath>
#include <queue>

struct interval {
    mpz_class c;
    int k;
    std::vector<mpq_class> p;
};
struct rootInterval {
    mpz_class c;
    int k, h;
    rootInterval(mpz_class c_ = 0, int k_ = 0, int h_ = 0) : c(c_), k(k_), h(h_) {}
};
// Finding the roots of n-degree polynomials
// by 叶隽希 20220925

inline double to_double(const mpq_class& x) {
    if (sgn(x) == 0) return 0;
    double b = x.get_d(), a = sgn(x) > 0 ? gb_getNext(b) : gb_getPrevious(b);
    std::cerr << "x - b = " << mpf_class(abs(x - b)) << '\n';
    std::cerr << "x - a = " << mpf_class(abs(x - a)) << '\n';
    std::cerr << "Choose " << (abs(x - b) < abs(x - a) ? b : a) << '\n';
    return abs(x - b) < abs(x - a) ? b : a;
}

std::vector<mpq_class> gb_scale(const std::vector<mpq_class> &a) {
    std::vector<mpq_class> b(a);
    int n = b.size() - 1;
    mpz_class d = 1;
    for (int i = n; ~i; i--) {
        b[i] *= d;
        d *= 2;
    }
    return b;
}
std::vector<mpq_class> gb_inversion(const std::vector<mpq_class> &a) {
    std::vector<mpq_class> b(a);
    std::reverse(b.begin(), b.end());
    return b;
}
std::vector<mpq_class> gb_translation(const std::vector<mpq_class> &a, mpq_class r = 1) {
    int n = a.size() - 1;
    std::vector<mpq_class> b; b.resize(n + 1);
    for (int i = n; ~i; i--) {
        b[i] = a[i];
        for (int j = i; j < n; j++)
            b[j] += b[j + 1] * r;
    }
    return b;
}
std::vector<mpq_class> gb_shift(const std::vector<mpq_class> &a) {
    int n = a.size() - 1;
    std::vector<mpq_class> b; b.reserve(n);
    for (int i = 1; i <= n; i++)
        b.emplace_back(a[i]);
    return b;
}
// inline int sgn(const mpq_class& x) {
//     if (std::abs(x) <= 1e-6) return 0;
//     if (x > 0) return 1;
//     return -1;
// }
int gb_variation(const std::vector<mpq_class>& a) {
    int n = a.size() - 1;
    int ans = 0, s1 = sgn(a[0]);
    for (int i = 1; i <= n; i++) {
        int s2 = sgn(a[i]);
//         std::cerr << s1 << ' ' << s2 << std::endl;
        if (s1 && s2 && s1 != s2) s1 = s2, ans++;
        else if (!s1) s1 = s2;
    }
    return ans;
}
mpq_class gb_eval(const std::vector<mpq_class> &a, mpq_class x0) {
    int n = a.size() - 1;
    mpq_class res = 0;
    for (int i = n; ~i; i--) {
        res = res * x0 + a[i];
    }
    return res;
}
void print(const std::vector<mpq_class> &a) {
    for (int n = a.size() - 1, i = 0; i <= n; i++)
        std::cerr << a[i] << "\n "[i != n];
}
std::vector<mpq_class> gb_derivative(const std::vector<mpq_class> &a) {
    std::vector<mpq_class> d;
    int n = a.size() - 1;
    for (int i = 1; i <= n; i++) {
        d.push_back(a[i] * i);
    }
    return d;
}
mpq_class roundAnother(const mpq_class& x, int eps, int greater) {
	mpz_class den = x.get_den();
	mpz_class td = den / (x.get_num() / x.get_den() + 1);
	int dig = mpz_sizeinbase(td.get_mpz_t(), 2);
	int mov = dig - eps - 3;
	if (mov > 0) {
		mpz_class num = x.get_num();
		num >>= mov;
		den >>= mov;
		if (greater) ++num;
		else ++den;
		return mpq_class(num, den);
	} else return x;
}

int pErrorBound = 100;

int iterationLimit = 100;

mpq_class gb_newton2(const std::vector<mpq_class> &a, const std::vector<mpq_class> &d, mpq_class x0) {
    for (int i = 0; i < iterationLimit; i++) {
        x0 = roundAnother(x0, pErrorBound, 0);
        x0 -= gb_eval(a, x0) / gb_eval(d, x0);
        // std::cerr << "x0(" << i << ")=" << x0 << std::endl;
    }
    return x0;
}
void gb_bisection(const std::vector<mpq_class> &a, std::vector<rootInterval>& ans) {
    std::cerr << "Starting bisection.\n";
    std::queue<interval> q;
    q.push({0, 0, a});
    ans.clear();
    while (!q.empty()) {
        interval u = q.front(); q.pop();
        // std::cerr << "Operating on interval " << '(' << u.c << ", " << u.k << ")\n";
//         std::cerr << "u.p = ";
//         for (int n = u.p.size() - 1, i = 0; i <= n; i++)
//             std::cerr << u.p[i].get_d() << "\n "[i != n];
        while (!u.p.empty() && u.p[0] == 0) {
            u.p = gb_shift(u.p);
            ans.emplace_back(rootInterval(u.c, u.k, 0));
        }
        if (u.p.empty()) continue;
        int v = gb_variation(gb_translation(gb_inversion(u.p)));
//         std::cerr << "res = "; print(gb_translation(gb_inversion(u.p)));
//         std::cerr << "v = " << v << "\n";
        if (v == 0) continue;
        if (v == 1 && u.k >= pErrorBound) ans.emplace_back(rootInterval(u.c, u.k, 1));
        else {
            q.push({u.c * 2, u.k + 1, gb_scale(u.p)});
            q.push({u.c * 2 + 1, u.k + 1, gb_translation(gb_scale(u.p))});
        }
    }
//     int m = ans.size() - 1;
//     for (int i = 0; i <= m; i++) {
//         std::cerr << "(" << ans[i].c << " + " << ans[i].h << ") / 2^" << ans[i].k << "\n";
//     }
}

inline mpq_class norm(const mpq_class& x, const mpq_class& bound) {
    return bound * (x * 2 - 1);
}

void gb_solve2(const std::vector<double>& a, double e, int& flag, std::vector<CP_Root> &r) {
    std::cerr << "Starting gb_solve2()\n";
    r.clear();
    int n = a.size() - 1;
    mpq_class an = mpq_class(a[n]), bound = abs(a[0] / an);
    for (int i = 1; i <= n; i++) {
        mpq_class val = abs(a[i] / an);
        bound = std::max(val, bound);
    }
    bound += 1;
    std::cerr << "bound = " << bound.get_d() << '\n';
    pErrorBound = (log(bound.get_d()) - log(e)) / log(2) * 2;
    std::cerr << "pErrorBound = " << pErrorBound << '\n';
    // std::cerr << "(raw) bound = " << bound << '\n';
    
    std::cerr << "Standarizing polynomial...\n";
    std::vector<mpq_class> b, der;
    mpq_class tmp = 1;
    for (int i = 0; i <= n; i++) {
        b.push_back(a[i] * tmp);
        tmp *= bound * 2;
    }
    b = gb_translation(b, -0.5);
    der = gb_derivative(b);
    // print(b); print(der);

    std::vector<rootInterval> res;
    gb_bisection(b, res);

    std::cerr << "Starting Newton...\n";

    int m = res.size();
    mpz_class d = 1; int k = 0;
    for (int i = 0; i < m; i++) {
        for (; k < res[i].k; k++)
            d *= 2;
        if (res[i].h == 0) {
            r.emplace_back(CP_Root(1, 0, to_double(2 * bound * (mpq_class(res[i].c, d)) - bound)));
            std::cerr << "Found border root " << to_double(2 * bound * (mpq_class(res[i].c, d)) - bound) << '\n';
        }
        else {
            mpq_class L = mpq_class(res[i].c, d), R = mpq_class(res[i].c + 1, d);
            std::cerr << "Found isolated root in (raw) [" << mpf_class(L) << ", " << mpf_class(R) << "]\n";
            std::cerr << "Found isolated root in       [" << norm(L, bound).get_d() << ", " << norm(R, bound).get_d() << "]\n";
            std::cerr << "Starting Newton at (raw) " << mpf_class(norm((L + R) / 2, bound)) << '\n';
            std::cerr << "Starting Newton at       " << norm((L + R) / 2, bound).get_d() << '\n';
            std::cerr << "Starting Newton at       " << gb_getNext(norm((L + R) / 2, bound).get_d()) << '\n';
            std::cerr << "Starting Newton at       " << gb_getPrevious(norm((L + R) / 2, bound).get_d()) << '\n';
            mpq_class x0 = gb_newton2(b, der, (L + R) / 2);
            // std::cerr << "x0 = " << (L + R) / 2 << '\n';
            // std::cerr << "gb = " << gb_newton2(b, der, (L + R) / 2) << '\n';
            std::cerr << "Ending Newton at (raw) " << mpf_class(norm(x0, bound)) << '\n';
            std::cerr << "Ending Newton at       " << norm(x0, bound).get_d() << '\n';
            std::cerr << "Ending Newton at       " << gb_getNext(norm(x0, bound).get_d()) << '\n';
            std::cerr << "Ending Newton at       " << gb_getPrevious(norm(x0, bound).get_d()) << '\n';
            // mpq_class x1 = gb_newton2(b, der, L);
            // mpq_class x2 = gb_newton2(b, der, R);
            // std::cerr << "Newton Method solves the root (raw) " << x1.get_d() << ' ' << x0.get_d() << ' ' << x2.get_d() << std::endl;
            // std::cerr << "Newton Method solves the root       " << norm(x1, bound).get_d() << ' ' << norm(x0, bound).get_d() << ' ' << norm(x2, bound).get_d() << std::endl;

            r.emplace_back(CP_Root(1, 0, to_double(norm(x0, bound))));
        }
    }
}
/* int main() {
    std::cerr.setf(std::ios_base::fixed); std::cerr.precision(10);
    std::cout.setf(std::ios_base::fixed); std::cout.precision(10);
    // std::vector<double> a({200000000, -302000200, 103000302, -1000103, 1});
    // std::vector<double> r({1, 2, 100, 1000000});
    std::vector<double> a({3, -4, 1});
    std::fprintf(stderr, "Work on polynomial:");
    for (int i = 0, n = a.size(); i < n; i++) {
        std::fprintf(stderr, " %.6lfx^%d",  a[i], i);
    }
    std::fprintf(stderr, "\n");
    std::vector<CP_Root> res;
    int tmp;
    gb_solve2(a, 0, tmp, res);
    for (int i = 0, n = res.size(); i < n; i++) {
        std::cerr << (res[i].m_root == 1.) << std::endl;
        std::string st; gb_toString(st, res[i].m_root);
        // std::cout << "Root #" << i << ": " << res[i].m_root << std::endl;
        std::cout << "Root #" << i << ": " << st << std::endl;
    }
} */
