#include <CP_Root.h>
#include <functional>
#include <iomanip>
#include <ctime>

typedef std::function<void(const std::vector<double>&, double, int&, std::vector<CP_Root>&)> pSolveFunc;

int main()
{
    std::cout << "****************************************************" << std::endl;
    std::cout << "*                                                  *" << std::endl;
    std::cout << "*      Polynomial Equation Real Roots Solving      *" << std::endl;
    std::cout << "*                                                  *" << std::endl;
    std::cout << "****************************************************" << std::endl;
    std::cout << std::endl;

    int n;
    double e;
    std::vector<double> p;
    std::cout << "Please input the degree of polynomial: ";
    std::cin >> n;
    if(n <= 0)
    {
        std::cout << "Exception: Non-positive degree." << std::endl;
        system("pause");
        return 0;
    }
    if(n == 1)
    {
        std::cout << "Assume the polynomial is a[0] + a[1]x, ";
        std::cout << "Please input the coefficients of polynomial, in the order of a[0], a[1]:" << std::endl;
    } 
    else
    {
        std::cout << "Assume the polynomial is a[0] + ... + a[" << n << "]x^" << n << ", ";
        std::cout << "Please input the coefficients of polynomial, in the order of a[0], ... a[" << n << "]:" << std::endl;
    }
    for(int i = 0; i <= n; i++)
    {
        double x;
        std::cin >> x;
        p.push_back(x);
    }

    std::vector<pSolveFunc> gb_solves = {gb_solve2, gb_solve3, gb_solve4};
    std::cout << "Please decide the type of algorithm(s) that you want to apply." << std::endl;
    std::cout << "     1 - Descartes-Bisection Method" << std::endl;
    std::cout << "     2 - Sturm-Newton Method" << std::endl;
    std::cout << "     3 - VAS Method" << std::endl;
    std::cout << "     0 - All" << std::endl;
    int type;
    std::cin >> type;
    if(type > 3 || type < 0)
    {
        std::cout << "Exception: Invalid algorithm type" << std::endl;
        system("pause");
        return 0;
    }

    std::cout << "Please input the maximum accepcted error(e.g. 1e-6): ";
    std::cin >> e;

    int prec;
    std::cout << "Please determine the output precision(by decimal digits, should not exceed 18): ";
    std::cin >> prec;
    if(prec < 0 || prec > 18)
    {
        std::cout << "Exception: Invalid output precision" << std::endl;
        return 0;
    }

    std::vector<CP_Root> ans;
    int flag = 0;
    if(type > 0)
    {
        gb_solves[type - 1](p, e, flag, ans);
        std::sort(ans.begin(), ans.end());
        std::cout << "Calculation succeeded. The roots are:" << std::endl;
        for(int i = 0; i < (int) ans.size(); i++)
            std::cout << std::fixed << std::setprecision(prec) << ans[i].m_root << ", count: " << ans[i].m_repetitions << std::endl;
    }
    else
    {
        for(int k = 0; k <= 2; k++)
        {
            int startc = clock();
            gb_solves[k](p, e, flag, ans);
            int endc = clock();
            std::cout << "Calculation succeeded. Algorithm #" << k + 1 << " costed " << endc - startc << " ms. The roots are:" << std::endl;
            for(int i = 0; i < (int) ans.size(); i++)
                std::cout << std::fixed << std::setprecision(prec) << ans[i].m_root << ", count: " << ans[i].m_repetitions << std::endl;
        }
    }
    system("pause");
    return 0;
}
