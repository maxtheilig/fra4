#include <iostream>
#include <fstream>
#include "ex2.hpp"

int main()
{
    double analytic_result = 0.886226925452758; //sqrt(pi)/2
    double result, uncertainty;
    std::cout << "\n";
    std::ofstream myfile;
    myfile.open ("ex2c.dat");
    for(uint N=1000; N<=1000000; N *= 10)
    {
        print("# of MC-steps", N);
        MonteCarloIntegrationGaussian(&result, xSquared, N, &uncertainty);
        print("result", result);
        print("numerical - analytic", std::abs(result - analytic_result));
        print("uncertainty", uncertainty);
        myfile << N << "\t" << uncertainty << "\n";
        std::cout << "\n";
    }
    myfile.close();
    return 0;
}
