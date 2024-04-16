#include <fstream>
#include "testlib.h"
#include <iostream>

int main() {
	std::ofstream fout(getDataFilePath("poly_deg_1_20.in"));
	generateDataSet(fout, 50);
	fout.close();
	std::ifstream fin(getDataFilePath("poly_deg_1_20.in"));
	auto t = readDataSet(fin);
	fin.close();
	for (auto poly : t) {
		for (auto x : poly.roots)
			std::cout << x << ' ';
		std::cout << '\n';
		for (auto x : poly.coefficients)
			std::cout << x << ' ';
		std::cout << '\n';
		std::cout << '\n';
	}
	return 0;
}
