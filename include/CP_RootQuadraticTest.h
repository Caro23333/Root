#ifndef CP_ROOTQUADRATICTEST_H
#define CP_ROOTQUADRATICTEST_H
#include <fstream>
#include "CP_RootQuadratic.h"

class CP_RootQuadraticTest
{
public:
    int m_caseTotal;
    vector<int> m_caseError; // 不满足容差的案例
    double m_maxError; // 计算出来的实根与对应给定的实根之间的最大误差
    int m_maxErrorCaseID;
    double m_maxValueInput; // 在给定的实根所对应的多项式的值中的最大值
    int m_maxValueInputCaseID;
    double m_maxValueSolve; // 在计算得到的实根所对应的多项式的值中的最大值
    int m_maxValueSolveCaseID;
    string m_fileInput;
    string m_fileOutput;
    CP_RootQuadraticTest();
    virtual ~CP_RootQuadraticTest() {}
    void mb_buildCases(const char *fileOut);
    void mb_reportBegin(ostream& os);
    void mb_reportEnd(ostream& os);
    void mb_test(const char *fileIn, const char *fileOut);
}; // 类CP_RootQuadraticTest定义结束

// 手工测试gb_solveQuadratic计算一元二次方程的实根:
extern void gb_testSolveQuadraticByIo();

#endif
