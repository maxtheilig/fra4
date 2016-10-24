#ifndef MCINT_H
#define MCINT_H

#include <cmath>
#include "ranlxd.h"
#include <time.h>

double prng_double();
double f1(double);
double f2(double x, double a = 1., double A = 1.5, double B = 0.5);
//void monteCarloIntegration(double & res, double (*f)(double), int steps, double & standartDeviation, bool zeroToOne = true, double a = 0., double b = 1.);

double prng_double()
{
    double tmp;
    ranlxd(&tmp, 1);
    return tmp;
}

double f1(double x)
{
    return x*x+3.0*x-1;
}

double f2(double x, double a, double A, double B)
{
    double tmp = (A * sqrt(x) + B) * a / (x * x + a * a);
    return tmp;
}

void monteCarloIntegration(double & res, double (*f)(double), int steps, double & standartDeviation, bool zeroToOne = true, double a = 0.0, double b = 1.0)
{
    double E_f = 0.0; //Expectation value of f
    double E_f_squared = 0.0; //Expectation value of f^2
    rlxd_init(1, time(NULL));
    for(unsigned int i=0; i<steps; ++i)
    {
        double x = prng_double();
        if(!zeroToOne) {
            x = (b - a) * x + a;}
        E_f += f(x);
        E_f_squared += f(x) * f(x);
    }
    E_f /= steps;
    E_f_squared /= steps;
    standartDeviation = sqrt((E_f_squared - E_f * E_f) / (steps - 1));
    if(zeroToOne) {
        res = E_f;}
    else {
        res = (b - a) * E_f;}
}

void monteCarloIntegration(double & res, double (*f)(double,double,double,double), double accuracy, double & standartDeviation, double c, int & stepsNeeded, double a = 0.0, double b = 1.0)
{
    double A = 1.5; double B = 0.5;
    double tmp1 = 0.0; //sum f(x_i)
    double tmp2 = 0.0; //sum f(x_i)^2
    double E_f = 0.0; //Expectation value of f
    double E_f_squared = 0.0; //Expectation value of f^2
    standartDeviation = 1.0;
    int stepcounter = 0;
    rlxd_init(1, time(NULL));
    while(standartDeviation / res >= accuracy)
    {
        ++stepcounter;
        //std::cout << "# of steps=" << stepcounter << std::endl;
        double x = prng_double();
        x = (b - a) * x + a;
        //std::cout << "randomNumber=" << x << "\n";
        tmp1 += f(x,c,A,B);
        //std::cout << "f(x)=" << f(x,c,A,B) << "\n";
        tmp2 += f(x,c,A,B) * f(x,c,A,B);
        //std::cout << "sum f(x_i)=" << tmp1 << "\n";
        //std::cout << "sum f^2(x_i)=" << tmp2 << "\n";
        E_f = tmp1 / stepcounter;
        //std::cout << "E[f]=" << E_f << std::endl;
        E_f_squared = tmp2 / stepcounter;
        //std::cout << "E[f^2]=" << E_f_squared << std::endl;
        //std::cout << "E^2[f]=" << E_f * E_f << std::endl;
        if(stepcounter < 2) {
            standartDeviation = 10000.;}
        else {
            standartDeviation = sqrt((E_f_squared - E_f * E_f) / (stepcounter - 1));}
        //std::cout << "standartDeviaton=" << standartDeviation << std::endl;
        res = (b - a) * E_f;
        //std::cout << "result_tmp=" << res << std::endl;
    }
    stepsNeeded = stepcounter;
}

#endif
