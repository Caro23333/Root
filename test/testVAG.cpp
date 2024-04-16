#include <bits/stdc++.h>
#include "VAG.h"
#include "testlib.h"

int main() {
	std::ifstream fin(getDataFilePath("poly_deg_1_20.in"));
	std::vector<polynomial> polys = readDataSet(fin);
    fin.close();
	std::cout << std::setprecision(20);
    for (auto p : polys) {
		std::vector<mpq_class> f; f.resize(p.coefficients.size());
		for (int i = 0; i < (int) p.coefficients.size(); ++i)
			f[i] = mpq_class(p.coefficients[i]);
		std::cout << "TEST STARTED" << std::endl;
		std::vector<root_t> ret = solve_VAG(f);
		std::cout << "TEST SUCCEED" << std::endl;
		std::cout << "-------- Result --------" << std::endl;
		std::cout << "Roots:" << std::endl;
		int ct = 0;
        for(auto root : ret)
            for(int j = 1; j <= root.second; j++) {
				std::cout << root.first << " (" << p.roots[ct] << ")" << std::endl;
				++ct;
			}
		// std::cin.get();
    }
    return 0;
}
