#include <iostream>
#include "McInt.hpp"

int main()
{
    double analytic_result = 5.0/6.0;
    double result, standartDeviation;
    std::cout << "\n";
    for(int N=1000; N<=10000000; N *= 10)
    {
        print("# of MC-steps", N);
        monteCarloIntegration(result, f1, N, standartDeviation);
        print("result", result);
        print("numerical - analytic", std::abs(result - analytic_result));
        print("standart_deviation", standartDeviation);
        std::cout << "\n";
    }
    return 0;
}
