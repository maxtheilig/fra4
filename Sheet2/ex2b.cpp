#include <iostream>
#include <fstream>
#include "ex2.hpp"

int main()
{
    double x;
    rlxd_init(1, time(NULL));
    std::ofstream myfile;
    myfile.open ("sigma_1.dat");
    for(uint i=0; i<1000000; ++i)
    {
        x = gaussianPRNG(1.0);
        myfile << x << "\n";
    }
    myfile.close();
    
    std::ofstream myfile1;
    myfile1.open ("sigma_01.dat");
    for(uint i=0; i<1000000; ++i)
    {
        x = gaussianPRNG(0.1);
        myfile1 << x << "\n";
    }
    myfile1.close();
    
    return 0;
}
