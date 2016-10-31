#ifndef EX2_HPP
#define EX2_HPP

#include <iostream>
#include <cmath>
#include "../ranlux/ranlxd.hpp"
#include "../ranlux/ranlxd.cpp"
#include <time.h>
#include "../tools.hpp"

double prng_double()
{
    double tmp;
    ranlxd(&tmp, 1);
    return tmp;
}

double gaussianPRNG(double sigma = 1.0)
{
    double a = sigma * sqrt(-2.0 * log(prng_double()));
    double b = cos(2 * 3.141592653589793238463 * prng_double());
    return a * b;
}

double xSquared(double x)
{
    return x * x;
}

void MonteCarloIntegrationGaussian(double * res, double (*f)(double), uint steps, double * uncertainty)
{
    rlxd_init(1, time(NULL));
    double volumeP = sqrt(3.141592653589793238463);
    double E_f = 0.0; //Expectation value of f
    double E_f_squared = 0.0; //Expectation value of f^2
    double x, standartDeviation;
    for(uint i=0; i<steps; ++i)
    {
        x = gaussianPRNG(sqrt(0.5));
        E_f += f(x);
        E_f_squared += f(x) * f(x);
    }
    E_f = E_f * volumeP / steps;
    E_f_squared = E_f_squared * volumeP * volumeP /steps;
    standartDeviation = sqrt((E_f_squared - E_f * E_f) / (steps - 1));
    *uncertainty = standartDeviation / E_f;
    *res = E_f;
}

#endif
