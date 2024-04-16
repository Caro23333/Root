/*******************************
* VAS bisection method by panjq22(panjq22@mails.tsinghua.edu.cn)
* Created on 2022/10/02
*******************************/
#ifndef VAS_H
#define VAS_H
#include "gmp.h"
#include "gmpxx.h"
#include "CP_Root.h"
#include <utility>
#include <vector>

// (root, k-multiple)
typedef std::pair<mpq_class, int> root_t;

// root-finding by VAS real-root isolating
// poly: polynomial(x) = 0, poly[i] means [x^i] poly(x)
// eps: absolute error 2^{-eps}
// guarantees all real roots are found
extern std::vector<root_t> VAS_Algorithm(
		std::vector<mpq_class> poly,
		int eps = 52);

struct mobius {
	static mpq rootlimit;
	mpz a, b, c, d;
	mobius() {}
	mobius(mpz _a, mpz _b, mpz _c, mpz _d) {
		a = _a, b = _b, c = _c, d = _d;
	}
	mobius operator() (mobius t) {
		return mobius(
				a * t.a + b * t.c, a * t.b + b * t.d,
				c * t.a + d * t.c, c * t.b + d * t.d);
	}
	mpq zero() {
		if (d == 0) return rootlimit;
		return mpq(b, d);
	}
	mpq infinity() {
		if (c == 0) return rootlimit;
		return mpq(a, c);
	}
	mpq operator() (mpz x) {
		return mpq(a * x + b, c * x + d);
	}
} ;

#endif
