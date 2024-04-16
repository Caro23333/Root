/*******************************
* VAG bisection method by panjq22(panjq22@mails.tsinghua.edu.cn)
* Created on 2022/10/02
*******************************/
#ifndef VAG_H
#define VAG_H
#include "gmp.h"
#include "gmpxx.h"
#include <utility>
#include <vector>

// (root, k-multiple)
typedef std::pair<mpq_class, int> root_t;

// root-finding by VAG real-root isolating
// poly: polynomial(x) = 0, poly[i] means [x^i] poly(x)
// eps: absolute error 2^{-eps}
// guarantees all real roots are found
extern std::vector<root_t> solve_VAG(
		std::vector<mpq_class> poly,
		int eps = 52);
#endif
