#include "random.hpp"

double prng_double()
{
    double tmp;
    ranlxd(&tmp, 1);
    return tmp;
}

double gaussianPRNG(double s)
{
    double a = s * sqrt(-2.0 * log(prng_double()));
    double b = cos(2 * 3.141592653589793238463 * prng_double());
    return a * b;
}
