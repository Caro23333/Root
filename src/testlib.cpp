#include <testlib.h>

#include <gmpxx.h>

#include <algorithm>
#include <iomanip>
#include <random>
#include <chrono>
#include <limits>
#include <vector>
#include <set>

namespace testlib {

	std::mt19937 rand;
    inline unsigned _rand(unsigned n) {
        return rand() % n;
    }

    double randomDouble() {
        const static double umax = std::numeric_limits<unsigned>::max(); // may lose precision?
        return testlib::rand() / umax;
    }

    int randomInteger(int L, int R) { // [L, R)
        unsigned len = R - L;
        return testlib::rand() % len + L;
    }

    void expandPoly(testlib::poly& p, const std::vector<mpq_class>& roots) {
        p.roots.clear(); p.coefficients.clear();
        std::vector<mpq_class> f; f.push_back(1);
        for (unsigned i = 0; i < roots.size(); i++) {
            p.roots.push_back(roots[i].get_d());
            f.push_back(1);
            for (int j = i; j > 0; j--)
                f[j] = -roots[i] * f[j] + f[j - 1];
            f[0] *= -roots[i];
        }
        for (unsigned i = 0; i < f.size(); i++)
            p.coefficients.push_back(f[i].get_d());
    }

    testlib::poly genPoly0(unsigned deg) {
        testlib::poly ret; ret.degree = deg;
        std::vector<int> vec;
        for (int i = -10; i <= 10; i++)
            vec.push_back(i);
        if (deg > vec.size())
            throw "[UserError] genPoly0(deg): deg > 21";
        std::random_shuffle(vec.begin(), vec.end(), testlib::_rand); // to prevent multiple roots
        std::vector<mpq_class> roots;
        for (unsigned i = 0; i < deg; i++) {
            roots.push_back(vec[i]);
        }
        expandPoly(ret, roots);
        ret.e = defaultErrorBound;
        return ret;
    }

    testlib::poly genPolyInt(int deg, int L, int R) { // [L, R)
        testlib::poly ret; ret.degree = deg;
        std::set<int> vis;
        if (R - L < deg)
            throw "[UserError] genPolyInt(deg, L, R): deg > R - L";
        std::vector<mpq_class> roots;
        for (int i = 0; i < deg; i++) {
            int x = randomInteger(L, R);
            for (; vis.find(x) != vis.end(); x = randomInteger(L, R));
            vis.insert(x); roots.push_back(x);
        }
        expandPoly(ret, roots); ret.e = defaultErrorBound;
        return ret;
    }

    testlib::poly genPolyDouble(int deg, int R) {
        testlib::poly ret; ret.degree = deg;
        std::vector<mpq_class> roots;
        for (int i = 0; i < deg; i++) {
            mpq_class x = (mpq_class(testlib::rand(), std::numeric_limits<unsigned>::max()) * 2 - 1) * R;
            roots.push_back(x); // 相信你的人品不会抽中两个近似实根
        }
        expandPoly(ret, roots); ret.e = defaultErrorBound;
        return ret;
    }

    testlib::poly genPolySpec() {
        std::vector<mpq_class> roots = {41, 40, 91, -60, -92};
        testlib::poly ret; ret.degree = roots.size();
        expandPoly(ret, roots); ret.e = defaultErrorBound;
        return ret;
    }

    testMainApp::~testMainApp() {
        fout.close();
        ferr.close();
        fans.close();
        fraw.close();
    }

    void testMainApp::getPoly() {
        if (mode == 0) {
            // return genPoly0(randomInteger(2, 10));
            // return genPolySpec();
            // return genPolyDouble(20, 10);
            readPoly("../data/poly.rawin");
        }
        else if (mode == 1) {
            p = genPolyInt(randomInteger(30, 31), -100, 100);
        }
        else if (mode == 2) {
            p = genPolyDouble(randomInteger(2, 20), 10);
        }
        else if (mode == 3) {
            p = genPolyInt(2, -100, 100);
        }
        else p = genPoly0(randomInteger(2, 4));
        std::sort(p.roots.begin(), p.roots.end());
    }

    void testMainApp::openFile(std::string path, std::string base) {
        this->path = path; this->base = base;
        fout.open(path + base + ".in");  fout << std::fixed << std::setprecision(53);
        ferr.open(path + base + ".out"); ferr << std::fixed << std::setprecision(53);
        fans.open(path + base + ".ans"); fans << std::fixed << std::setprecision(53);
        fraw.open(path + base + ".raw"); fraw << std::fixed << std::setprecision(53);
    }

    void testMainApp::prettyPoly(std::ostream& os) {
        os << "#DEGREE " << p.degree << "\n";
        os << "#COEFFICIENTS ";
        for (int i = 0; i <= p.degree; i++)
            os << p.coefficients[i] << "\n "[i != p.degree];
        int rt = p.roots.size();
        os << "#REAL_ROOTS " << rt << "\n";
        for (int i = 0; i < rt; i++)
            os << "    ROOTS#" << (i + 1) << " " << p.roots[i] << "\n";
        os << "#ERRORBOUND " << p.e << std::endl; // refresh?
        // todo: multiple roots
    }

    void testMainApp::rawPoly(std::ostream& os) {
        os << p.degree << "\n";
        for (int i = 0; i <= p.degree; i++) {
            double x = p.coefficients[i];
            os << *(unsigned long long*)&x << " ";
        }
        os << "\n";
        os << p.roots.size() << "\n";
        for (auto r : p.roots) {
            double x = r;
            os << *(unsigned long long*)&x << " ";
        }
        os << "\n";
        os << p.e << "\n\n";
    }

    void testMainApp::prettyAns(std::ostream& os) {
        os << "#ROOT_CALCULATED " << res.size() << "\n";
        for (unsigned i = 0; i < res.size(); i++) {
            os << "    ROOT#" << (i + 1) << " " << res[i].m_root << "\n";
        }
    }

    void testMainApp::readPoly(const std::string& fn) { // 暂只支持一组，请删除编号
        std::ifstream is(fn);
        testlib::poly ret; is >> ret.degree;
        unsigned long long x;
        for (int i = 0; i <= ret.degree; i++) {
            is >> x;
            ret.coefficients.push_back(*(double*)&x);
        }
        int n; is >> n;
        for (int i = 0; i < n; i++) {
            is >> x;
            ret.roots.push_back(*(double*)&x);
        }
        // is >> ret.e;
        ret.e = defaultErrorBound;
        p = ret;
    }

    void testMainApp::logPolynomial() {
        prettyPoly(fout);
        prettyPoly(fans);
        rawPoly(fraw);
    }

    void testMainApp::compareRoots(const std::vector<CP_Root>& v1, const std::vector<CP_Root>& v2, std::vector<int>& r) {
        r.clear();
        for (unsigned i = 0; i < v1.size(); i++) {
            // Assuming that there's no multiple roots
            bool flag = 0;
            for (unsigned j = 0; j < v2.size(); j++) {
                if (std::abs(v1[i].m_root - v2[j].m_root) <= p.e) {
                    flag = 1;
                    break;
                }
            }
            if (!flag) r.push_back(i);
        }
    }

    void testMainApp::compareResult(std::vector<int>& r1, std::vector<int>& r2) {
        std::vector<CP_Root> v;
        for (unsigned i = 0; i < p.roots.size(); i++)
            v.push_back(CP_Root(1, 0, p.roots[i]));
        compareRoots(res, v, r1); compareRoots(v, res, r2);
    }

    void testMainApp::logCaseNumber(int i, int nTests) {
        fraw << "*** Test Case " << i << " (PLEASE DELETE)\n";
        fout << "### Test Case " << i << "/" << nTests << std::endl;
        ferr << "### Test Case " << i << "/" << nTests << std::endl;
        fans << "### Test Case " << i << "/" << nTests << std::endl;
    }

    const static std::string ret_msg[] = {
        "Invalid input",
        "No real roots",
        "Found every root within error bound",
        "Found every root, but some exceed the error bound",
        "Found some roots within error bound",
        "Fount some roots, but some exceed the error bound",
        "Real root exists but not found",
        "Cannot decide whether there is real root or not"
    };

    bool flagMessage = 0;

    void testMainApp::setMode(int mode) {
        this->mode = mode;
    }

    bool testMainApp::test(std::ostream& os, pSolveFunc gb_solve) {
        getPoly(); // 这个用 mode 的想法不太好?

        logPolynomial();

        restartTimer();
        gb_solve(p.coefficients, p.e, flag, res);
        lapTimer("Algorithm ends");
        displayTimer(os);

        prettyAns(fans);

        if (flagMessage) {
            ferr << "gb_solve() returns with value " << flag << "\n";
            if (-1 <= flag && flag <= 6)
                os << "gb_solve() returns with value \"" << ret_msg[flag + 1] << "\"\n";
            else os << "gb_solve() returns invalid flag\n";
        }

        std::vector<int> r1, r2;
        compareResult(r1, r2);

        os << "There are " << r1.size() << " calculated root(s) that exceed the error bound.\n";
        ferr << "There are " << r1.size() << " calculated root(s) that exceed the error bound.";
        if (r1.size()) ferr << " They are:";
        ferr << "\n";
        for (unsigned i = 0; i < r1.size(); i++)
            ferr << "    Root #" << (r1[i] + 1) << ": " << res[r1[i]].m_root << "\n";

        os << "There are " << r2.size() << " real root(s) that are not found.\n";
        ferr << "There are " << r2.size() << " real root(s) that are not found.";
        if (r2.size()) ferr << " They are:";
        ferr << "\n";
        for (unsigned i = 0; i < r2.size(); i++)
            ferr << "    Root #" << (r2[i] + 1) << ": " << p.roots[r2[i]] << "\n";

        return r1.empty() && r2.empty();
        // TODO: evaluate f(x0)
        // TODO: multiple root support
        // TODO: distinguish fake roots with not-so-precise ones
    }

    bool testMainApp::testMultiple(std::ostream& os, std::vector<pSolveFunc> gb_solves) {
        getPoly(); // 这个用 mode 的想法不太好?

        logPolynomial();

        std::vector<std::vector<CP_Root> > ress;

        std::vector<CP_Root> vvv;
        for (unsigned i = 0; i < p.roots.size(); i++)
            vvv.push_back(CP_Root(1, 0, p.roots[i]));
        ress.push_back(vvv);

        int szSolve = gb_solves.size();
        for (int i = 0; i < szSolve; i++) {

            const pSolveFunc& gb_solve = gb_solves[i];
            restartTimer();
            gb_solve(p.coefficients, p.e, this->flag, res);
            lapTimer("Algorithm #" + std::to_string(i + 1) + " ends");
            displayTimer(os);

            ress.push_back(res);

            os << "Test #" << i + 1 << std::endl;
            for(int j = 0; j < (int) res.size(); j++)
                os << res[j].m_root << std::endl;
            //os << endl;  
        }

        bool flag = 1;
        std::vector<int> r;
        for (int i = 0; i <= szSolve; i++) {
            for (int j = 0; j <= szSolve; j++) {
                if (i != j) {
                    compareRoots(ress[i], ress[j], r);
                    if (!r.empty())
                        os << "#" << i << " has " << r.size() << " roots not in #" << j << "\n";
                    if (i && j) flag &= r.empty();
                }
            }
        }
        
        return flag;

    }

    void testMainApp::restartTimer() {
        timer = std::chrono::high_resolution_clock::now();
        stopList.clear();
    }

    void testMainApp::lapTimer(const std::string& name) {
        stopList.push_back(std::make_pair(name, std::chrono::high_resolution_clock::now()));
    }

    void testMainApp::displayTimer(std::ostream& os) {
        for (auto it : stopList) {
            // os << it.first << ": " << std::chrono::duration_cast<std::chrono::microseconds>(it.second - timer).count() << '\n';
            os << it.first << ": " << std::chrono::duration_cast<std::chrono::milliseconds>(it.second - timer).count() << "ms\n";
        }
    }
    
    testMainApp app;
}

namespace {
	class initializer {
		public:
			initializer() {
				initRNGSeed();
			}
		private:
			void initRNGSeed() {
                testlib::rand.seed(std::chrono::high_resolution_clock::now().time_since_epoch().count() + std::random_device{}());
			}
	} init;
}
