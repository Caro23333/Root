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
#endif
