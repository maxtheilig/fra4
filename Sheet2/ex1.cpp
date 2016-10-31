#include <iostream>
#include "metro.hpp"

int main()
{
    double analytic_result = 0.886226925452758; //sqrt(pi)/2
    double result, acceptanceRate, uncertainty, standartDeviation;
    uint steps = 10000000;
    //double epsilon[] = {0.01, 0.02, 0.04, 0.08, 0.16, 0.32, 0.64, 1.25};
    double epsilon[] = {0.0005, 0.55, 20000.0};
    double startingPoint = 0.0;
    double volumeP = sqrt(3.141592653589793238463);
    /*for(uint i=0; i<(sizeof(epsilon)/sizeof(*epsilon)); ++i)
    {
        std::cout << "\n";
        std::cout << "################################" << "\n";
        std::cout << "####### epsilon = " << epsilon[i] << " ########" << "\n";
        std::cout << "################################" << "\n";
        MetropolisHastings(&result, xSquared, expMinusXSquared, volumeP, steps, startingPoint, epsilon[i], &acceptanceRate, & uncertainty, & standartDeviation);
        print("result", result);
        print("numerical - analytic", std::abs(result - analytic_result));
        print("acceptanceRate", acceptanceRate);
        print("uncertainty", uncertainty);
        print("standartDeviation", standartDeviation);
        std::cout << "\n";
    }*/
    MetropolisHastings(&result, xSquared, expMinusXSquared, volumeP, startingPoint, epsilon[1], analytic_result);
    print("result", result);
    
    return 0;
}
