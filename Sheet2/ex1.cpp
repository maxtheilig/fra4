#include <iostream>
#include "metro.hpp"

int main()
{
    double analytic_result = 0.886226925452758; //sqrt(pi)/2
    double result, acceptanceRate;
    uint steps = 10000000;
    double epsilon[] = {0.01, 0.02, 0.04, 0.08, 0.16, 0.32, 0.64, 1.28, 2.56, 5.12, 10.24, 40.96};
    double startingPoint = 0.0;
    double volumeP = sqrt(3.141592653589793238463);
    for(uint i=0; i<11; ++i)
    {
        std::cout << "\n";
        print("epsilon", epsilon[i]);
        MetropolisHastings(&result, xSquared, expMinusXSquared, volumeP, steps, startingPoint, epsilon[i], &acceptanceRate);
        print("result", result);
        print("numerical - analytic", std::abs(result - analytic_result));
        print("acceptanceRate", acceptanceRate);
        std::cout << "\n";
    }
    return 0;
}
