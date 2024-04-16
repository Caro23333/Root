#include <iostream>
#include <iomanip>
#include <bitset>
unsigned long long to_ull(double x) {
    return *reinterpret_cast<unsigned long long*>(&x);
}
std::bitset<64> to_bitset(double x) {
    return std::bitset<64>(to_ull(x));
}
double to_double(unsigned long long x) {
    return *reinterpret_cast<double*>(&x);
}
void print(double x) {
    std::cout << to_bitset(x) << std::endl;
}
int main() {
    std::cout << std::fixed << std::setprecision(50);
    double x = 9;
    double y = to_double(to_ull(x) - 1);
    std::cout << x << '\n'; print(9);
    std::cout << y << '\n'; print(y);
}
