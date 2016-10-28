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

void MetropolisHastings(double * res, double (*f)(double), double (*p)(double), double volumeP, uint steps, double startingPoint, double epsilon, double * acceptanceRate)
{
    rlxd_init(1, time(NULL));
    int acceptedUpdates = 0;
    double x, y, delta, deltaP, res_tmp, r;
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
    }
    *acceptanceRate = (double) acceptedUpdates / (double) steps;
    *res = volumeP * res_tmp / steps;
}


#endif
