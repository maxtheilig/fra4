#ifndef EX2_HPP
#define EX2_HPP

#include <cmath>
#include "../ranlux/ranlxd.hpp"
#include "../ranlux/ranlxd.cpp"
#include <time.h>
#include "../tools.hpp"

void Metropolis(double * result, double(*f)(double), uint steps, double omega, double epsilon, double *x, double * acceptanceRate, uint Nt, bool improvedDerivative);

double prng_double();

double LocalAction(double x0, double x1, double omega);

double LocalAction(double x0, double x1, double x2, double omega);

double xSquared(double x);

double analyticResult(double omega, uint Nt);

double gaussianPRNG(double s = 1.0);

#endif
