#include <iostream>
#include "ranlxd.h"
#include "McInt.hpp"
#include <cmath>

int main()
{
    double analytic_result[3] = {3.81256653333958, 1.80862037922022, 1.115564383765492};
    double c[3] = {1.0, 0.1, 0.01};
    double result, standartDeviation;
    int stepsNeeded;
    double accuracy = 0.001;
    double a = 0.0; double b = 100.0;
    for(unsigned int i=0; i<3; ++i) {
        std::cout << "\n";
        std::cout << "a = " << c[i] << "\n";
        monteCarloIntegration(result, f2, accuracy, standartDeviation, c[i], stepsNeeded, a, b);
        std::cout << "# steps = " << stepsNeeded << "\n";
        std::cout << "result=" << result << std::endl;
        std::cout << "numerical - analytic = " << std::abs(result - analytic_result[i]) << std::endl;
        std::cout << "standart_deviation = " << standartDeviation << std::endl;
        std::cout << "\n";
    }
    std::cout << "time = " << time(NULL) << "\n";
    return 0;
}
