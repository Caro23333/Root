#include <bits/stdc++.h>
#include "VAS.h"
#include "testlib.h"

int main() {
	std::ifstream fin("data/poly_single.in");
	std::vector<double> p;
	std::vector<CP_Root> res;
	std::string buf;
	int n;
	fin >> buf; fin >> n;
	fin >> buf;
	p.resize(n + 1);
	for (int i = 0; i <= n; ++i)
		fin >> p[i];
    fin.close();
	std::cout << std::setprecision(20);
	int flag;
	gb_solve4(p, 1, flag, res);
	std::cout << "GET ROOT:\n";
	for (auto t : res)
		std::cout << t.m_root << "\n";
    return 0;
}
