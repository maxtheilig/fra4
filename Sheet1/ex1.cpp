#include <iostream>
#include "McInt.hpp"

int main()
{
    double analytic_result = 5.0/6.0;
    double result, standartDeviation;
    std::cout << "\n";
    for(unsigned int N=1000; N<=10000000; N *= 10)
    {
        std::cout << "# of MC-steps = " << N << std::endl;
        monteCarloIntegration(result, f1, N, standartDeviation);
        std::cout << "result = " << result << std::endl;
        std::cout << "numerical - analytic = " << std::abs(result - analytic_result) << std::endl;
        std::cout << "standart_deviation = " << standartDeviation << std::endl;
        std::cout << "\n";
    }
    return 0;
}
