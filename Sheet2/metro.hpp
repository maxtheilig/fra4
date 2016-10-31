#ifndef METRO_H
#define METRO_H

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

double xSquared(double x)
{
    return x * x;
}

double expMinusXSquared(double x)
{
    double tmp = exp(-x * x);
    return tmp;
}

double relativePrecision(double exactResult, double numericalResult)
{
    double tmp = 1.0 - exactResult / numericalResult;
    return tmp;
}

void MetropolisHastings(double * res, double (*f)(double), double (*p)(double), double volumeP, uint steps, double startingPoint, double epsilon, double * acceptanceRate, double * uncertainty, double * standartDeviation)
{
    rlxd_init(1, time(NULL));
    int acceptedUpdates = 0;
    double x, y, delta, deltaP, res_tmp, res_tmp_squared, r, E_f, E_f_squared;
    res_tmp = 0.0; res_tmp_squared = 0.0;
    x = startingPoint;
    for(uint i=0; i<steps; ++i)
    {
        delta = 2.0 * prng_double() - 1.0;
        y = x + delta * epsilon;
        //print("newPoint", y);
        deltaP = p(y) / p(x);
        //print("deltaP", deltaP);
        if(deltaP >=1) {
            x = y;
            ++acceptedUpdates;
        }
        else {
            r = prng_double();
            if(r < deltaP){x = y; ++acceptedUpdates;}
        }
        res_tmp += f(x);
        res_tmp_squared += f(x) * f(x);
    }
    E_f = volumeP * res_tmp / steps; E_f_squared = volumeP * volumeP * res_tmp_squared / steps;
    *standartDeviation = sqrt((E_f_squared - E_f * E_f) / (double)(steps - 1));
    *uncertainty = *standartDeviation / E_f;
    *acceptanceRate = (double) acceptedUpdates / (double) steps;
    *res = E_f;
}

void MetropolisHastings(double * res, double (*f)(double), double (*p)(double), double volumeP, double startingPoint, double epsilon, double exactResult)
{
    rlxd_init(1, time(NULL));
    double res_tmp = 0.0;
    double deltaP, x, y, delta, r;
    uint N = 10;
    x = startingPoint;
    for(uint i=0; i<N; ++i)
    {
        delta = 2.0 * prng_double() - 1.0;
        y = x + delta * epsilon;
        //print("newPoint", y);
        deltaP = p(y) / p(x);
        //print("deltaP", deltaP);
        if(deltaP >=1) {
            x = y;
        }
        else {
            r = prng_double();
            if(r < deltaP){
                x = y;
            }
        }
        res_tmp += f(x);
    }
    uint stepcounter = N;
    *res = volumeP * res_tmp / stepcounter;
    while(std::abs(relativePrecision(exactResult, *res)) >= 0.01)
    {
        ++stepcounter;
        delta = 2.0 * prng_double() - 1.0;
        y = x + delta * epsilon;
        //print("newPoint", y);
        deltaP = p(y) / p(x);
        //print("deltaP", deltaP);
        if(deltaP >=1) {
            x = y;
        }
        else {
            r = prng_double();
            if(r < deltaP){
                x = y;
            }
        }
        res_tmp += f(x);
        *res = volumeP * res_tmp / stepcounter;
    }
    std::cout << "For N=" << stepcounter << " and epsilon=0.55 the relative precision is smaller then 0.1%" << "\n";
}

#endif
