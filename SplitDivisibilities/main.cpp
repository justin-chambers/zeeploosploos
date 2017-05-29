/**
 * Consider the number 48.
 * There are five pairs of integers a and b (a≤b) such that a×b=48: (1,48), (2,24), (3,16), (4,12) and (6,8).
 * It can be seen that both 6 and 8 have 4 divisors.
 * So of those five pairs one consists of two integers with the same number of divisors.
 *
 * In general:
 * Let C(n) be the number of pairs of positive integers a×b=n, (a≤b) such that a and b have the same number of divisors;
 * so C(48)=1.
 *
 * You are given C(10!)=3: (1680, 2160), (1800, 2016) and (1890,1920).
 *
 * Find C(100!)
 */

#include <iostream>
#include <tuple>
#include <vector>
#include <boost/multiprecision/cpp_int.hpp>
namespace mp = boost::multiprecision;

mp::cpp_int factorial(mp::cpp_int n);
void factor(mp::cpp_int x, std::vector<std::tuple<mp::cpp_int, mp::cpp_int>>& tpls);
void display_tuples(std::vector<std::tuple<mp::cpp_int, mp::cpp_int>> tpls);

int main() {
    mp::cpp_int n = factorial(100);
    std::cout << "n = " << n << std::endl;
    std::vector<std::tuple<mp::cpp_int, mp::cpp_int>>factors;
    std::vector<std::tuple<mp::cpp_int, mp::cpp_int>>cfactors;
    factor(n,factors);
    int i = 0, bigC = 0;
    for(i; i<(int)factors.size();++i) {
        std::vector<std::tuple<mp::cpp_int, mp::cpp_int>>factorsA;
        std::vector<std::tuple<mp::cpp_int, mp::cpp_int>>factorsB;
        factor(std::get<0>(factors[i]),factorsA);
        factor(std::get<1>(factors[i]),factorsB);
        if(factorsA.size() == factorsB.size()) {
            ++bigC;
            cfactors.push_back(factors[i]);
        }
    }
    std::cout << "C(n) = " << bigC <<std::endl;
    display_tuples(cfactors);
    return 0;
}

mp::cpp_int factorial(mp::cpp_int n) {
    mp::cpp_int result = 1;
    for(mp::cpp_int i=n;i>0;--i) {
        result *= i;
    }
    return result;
}

void factor(mp::cpp_int x, std::vector<std::tuple<mp::cpp_int, mp::cpp_int>>& tpls) {
    mp::cpp_int i = 0, quotient = x;
    for(i = 1; i<quotient; ++i) {
        if((x%i)==0) {
            tpls.push_back(std::make_tuple(i,x/i));
            quotient = x/i;
        }
    }
}

void display_tuples(std::vector<std::tuple<mp::cpp_int, mp::cpp_int>> tpls) {
    int i = 0;
    for(i; i < (int)tpls.size(); ++i) {
        std::cout << "(" << std::get<0>(tpls[i]) << ", " << std::get<1>(tpls[i]) << ") ";
        if((i+1) % 10 == 0) std::cout << std::endl;
    }
    std::cout << std::endl;
}