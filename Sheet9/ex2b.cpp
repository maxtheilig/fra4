#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include "../ranlux/ranlxd.cpp"

#define Nt 10
#define lambda1 1.0
#define lambda2 1.0

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

double expectationValue(double (*func)(double), std::vector<double> in)
{
    double result = 0.0;
    for(unsigned int i=0; i<in.size(); ++i)
    {
        result += func(in.at(i));
    }
    result /= in.size();
    return result;
}

double xSquared(double x)
{
    return x * x;
}

void overrelax_update(std::vector<double> *x, int idx)
{
    double xip1, xim1;
    
    if(idx == (Nt-1))
        xip1 = x->at(0);
    else
        xip1 = x->at(idx+1);
    
    if(idx == 0)
        xim1 = x->at(Nt-1);
    else
        xim1 = x->at(idx-1);
    
    x->at(idx) = 1.0/(2.0+lambda1+lambda2) * (xip1 + xim1) - x->at(idx);
}

void overrelax_sweep(std::vector<double> *x)
{
    int idx;
    for(unsigned int i=0; i<Nt; ++i)
    {
        idx = round(prng_double() * (Nt - 1));
        overrelax_update(x, idx);
    }
}

int main()
{
    double result;
    int Nstep, Ntherm, Nskip;
    Nstep = 1000000; Ntherm = 0; Nskip = 0;
    rlxd_init(1,time(NULL));
    std::vector<double> x;
    for(unsigned int i=0; i<Nt; ++i) x.push_back(gaussianPRNG(1.0));
    
    for(unsigned int i=0; i<Ntherm; ++i)
    {
        overrelax_sweep(&x);
    }
    for(unsigned int i=0; i<Nstep; ++i)
    {
        for(unsigned int k=0; k<Nskip; ++k)
        {
            overrelax_sweep(&x);
        }
        overrelax_sweep(&x);
        result += expectationValue(xSquared, x);
        //std::cout << "result = " << result << "\n";
    }
    result /= Nstep;
    std::cout << "<x^2> = " << result << "\n";
    
    return 0;
}
