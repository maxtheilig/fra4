#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <cmath>
#include "../ranlux/ranlxd.cpp"

#define Nt 10
#define omega 5.0

double analyticResult()
{
    double R = 1.0 + (omega * omega) / 2.0 - omega * sqrt(1.0 + (omega * omega) / 4.0);
    double tmp = 1.0 / (4.0 * omega * sqrt(1.0 + (omega * omega) / 4.0)) * (1.0 + pow(R, Nt)) / (1.0 - pow(R, Nt));
    return tmp;
}

void save(std::vector<double> in, std::string filename)
{
    std::ofstream out;
    out.open (filename);
    for(unsigned int i=0; i<in.size(); ++i)
        out << in.at(i) << "\n";
    out.close();
}

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
    
    x->at(idx) = 1.0/(2.0+pow(omega,2)) * (xip1 + xim1) - x->at(idx);
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

void update_heatbath(std::vector<double> *x, int idx)
{
    double y, xip1, xim1;
    y = gaussianPRNG(1/sqrt(2.0));
    /*if(idx == (Nt-1))
        xip1 = x->at(0);
    else
        xip1 = x->at(idx+1);
    if(idx == 0)
        xim1 = x->at(Nt-1);
    else
        xim1 = x->at(idx-1);*/
    int idx_neigh, idx_lower_neigh;
    if(idx == Nt - 1) idx_neigh = 0;
    else idx_neigh = idx + 1;
    if(idx == 0) idx_lower_neigh = Nt - 1;
    else idx_lower_neigh = idx - 1;

    x->at(idx) = 1.0/sqrt(2.0+pow(omega,2)) * y + 1.0/(2.0+pow(omega,2)) * (x->at(idx_neigh) + x->at(idx_lower_neigh));
}

void heatbath_sweep(std::vector<double> *x)
{
    int idx;
    for(unsigned int i=0; i<x->size(); ++i)
    {
        idx = round(prng_double() * (x->size() - 1));
        update_heatbath(x, idx);
    }
}

int main()
{
    std::ostringstream oss; std::string filename;
    double result;
    int Nstep, Ntherm, Nskip, Nor;
    Nstep = 100000; Ntherm = 10000; Nskip = 0; Nor = 0;
    rlxd_init(1,time(NULL));
    for(Nor=0; Nor<1; ++Nor)
    {
    std::vector<double> x, resultVector;
    for(unsigned int i=0; i<Nt; ++i) x.push_back(gaussianPRNG(1.0));
    
    for(unsigned int i=0; i<Ntherm; ++i)
        heatbath_sweep(&x);
    for(unsigned int i=0; i<Nstep; ++i)
    {
        for(unsigned int k=0; k<Nskip; ++k)
            heatbath_sweep(&x);
        heatbath_sweep(&x);
        for(unsigned int m=0; m<Nor; ++m)
            overrelax_sweep(&x);
        resultVector.push_back(expectationValue(xSquared, x));
        result += expectationValue(xSquared, x);
    }
    oss.str("");
    oss << "Nor_" << Nor << ".txt";
    filename = oss.str();
    result /= Nstep;
    std::cout << "<x^2> = " << result << "\n";
    std::cout << "analytic result = " << analyticResult() << "\n";
    std::cout << "delta = " << std::abs(result-analyticResult()) << "\n";
    //save(resultVector, filename);
    }
    
    return 0;
}

    


