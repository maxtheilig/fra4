#ifndef RANDOM_HPP
#define RANDOM_HPP

#include <cmath>
#include "ranlux/ranlxd.hpp"
#include "ranlux/ranlxd.cpp"

double prng_double();

double gaussianPRNG(double s = 1.0);

#endif
