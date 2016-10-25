#include <iostream>
#include "McInt.hpp"

int main()
{
    std::cout << "\n";
    std::cout << "###############################################################\n";
    std::cout << "##### Monte-Carlo-Integration without importance sampling #####\n";
    std::cout << "###############################################################\n";
    double analytic_result[3] = {3.81256653333958, 1.80862037922022, 1.115564383765492};//b=100
    double c[3] = {1.0, 0.1, 0.01};
    double result, standartDeviation;
    int stepsNeeded;
    double accuracy = 0.001;
    double a = 0.0; double b = 100;
    for(unsigned int i=0; i<3; ++i) {
        print("a", c[i]);
        monteCarloIntegration(result, f2, accuracy, standartDeviation, c[i], stepsNeeded, a, b);
        print("# steps needed", stepsNeeded);
        print("result", result);
        print("numerical - analytic", std::abs(result - analytic_result[i]));
        print("standart_deviation", standartDeviation);
        std::cout << "\n";
    }
    
    std::cout << "############################################################\n";
    std::cout << "##### Monte-Carlo-Integration with importance sampling #####\n";
    std::cout << "############################################################\n";
    for(unsigned int i=0; i<3; ++i) {
        print("a", c[i]);
        monteCarloIntegration(result, f3, accuracy, standartDeviation, c[i], stepsNeeded);
        print("# steps needed", stepsNeeded);
        print("result", result);
        print("numerical - analytic", std::abs(result - analytic_result[i]));
        print("standart_deviation", standartDeviation);
        std::cout << "\n";
    }
    return 0;
}
