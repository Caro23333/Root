#include "testlib.h"
#include "CP_ROOT.h"

#include <iostream>
#include <fstream>
#include <string>
#include <chrono>
#include <cmath>

constexpr int nTests = 1;

using testlib::app;

// mode 读入 / 硬编码 自行修改
// 我知道 mode 很怪但我也没办法
// 相信你会用这个 timer，我知道这个也很怪但是能用就行（
// data reuse 用的是 .rawin!!! 数据记得备份，因为即使 reuse 也会重新写一份 raw
// reuse 暂时只支持一组数据 请把需要 delete 的那一行删掉
// errorBound 现在是 1e-15!
int main() {
    app.openFile("../data/", "poly");

    std::cout.setf(std::ios_base::fixed); std::cout.precision(63);
    std::cerr.setf(std::ios_base::fixed); std::cerr.precision(63);

    int success = 0;

    int mode; std::cin >> mode;
    // int mode = 0;

    app.setMode(mode);

    std::vector<testlib::pSolveFunc> gb_solves = {gb_solve2, gb_solve3, gb_solve4};

    for (int i = 1; i <= (mode ? nTests : 1); i++) {
        std::cout << "### Running Test Case #" << i << "\n";

        app.logCaseNumber(i, nTests);
        bool res = app.testMultiple(std::cout, gb_solves);
        // bool res = app.test(std::cout, gb_solve4);
        if (res) success++;
        system("pause");
    }
    std::cout << "Test success " << success << "/" << nTests << "\n";
    
}
