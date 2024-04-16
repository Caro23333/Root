#include <iostream>
using namespace std;
#include "CP_Root.h"
#include "CP_RootTest.h"
#include "CP_RootQuadraticTest.h"
#include "testlib.h"
#include <string>
#include <vector>

vector<double> a;
vector<CP_Root> r;

void output(std :: ofstream &out, double x)
{
    string res;
    gb_toString(res, x);
    out << res << " ";
}

int main(int argc, char* args[])
{
    int n;
    ifstream fin(getDataFilePath("poly_deg_1_20.in"));
    ofstream fout(getDataFilePath("poly_deg_1_20.out"));
    fout << "******* AUTOMATIC TEST *******" << endl;
    vector<polynomial> polys = readDataSet(fin);
    int cnt = 0;
    for(auto p : polys)
    {
        cnt++;
        int flag = 0;
        gb_solve(p.coefficients, p.e, flag, r);
        fout << "-------- Result #" << cnt << " --------" << endl;
        fout << "Calculated roots:" << endl;
        for(auto root : r)
            for(int j = 1; j <= root.m_repetitions; j++)
                output(fout, root.m_root);
        fout << endl;
        fout << "Exact roots:" << endl;
        for(auto root : p.roots)
            output(fout, root);
        fout << endl << endl;
    }
    fout << "******* DONE *******" << endl;
    fin.close();
    fout.close();
    return 0; // 返回0表明程序运行成功
} // main函数结束
