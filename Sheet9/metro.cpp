#include <iostream>
#include <fstream>
#include <vector>
#include "../random.cpp"

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
    {
        out << in.at(i) << "\n";
    }
    out.close();
}

void save(std::vector<double> x, std::vector<double> y, std::string filename)
{
    std::ofstream out;
    out.open (filename);
    for(unsigned int i=0; i<x.size(); ++i)
    {
        out << x.at(i) << "\t" << y.at(i) << "\n";
    }
    out.close();
}

void save(std::vector<double> x, std::vector<double> y, std::vector<double> z, std::string filename)
{
    std::ofstream out;
    out.open (filename);
    for(unsigned int i=0; i<x.size(); ++i)
    {
        out << x.at(i) << "\t" << y.at(i) << "\t" << z.at(i) << "\n";
    }
    out.close();
}

double deltaS(std::vector<double> x, int idx, double xnew)
{
    int idx_neigh, idx_lower_neigh;
    if(idx == Nt - 1) idx_neigh = 0;
    else idx_neigh = idx + 1;
    if(idx == 0) idx_lower_neigh = Nt - 1;
    else idx_lower_neigh = idx - 1;
    double tmp = (2 + pow(omega,2))*(pow(x[idx],2) - pow(xnew,2)) + 2.0 * (x.at(idx_neigh) + x.at(idx_lower_neigh))*(x.at(idx) - xnew);
    return tmp;
}

double expectationValue(double (*func)(double), std::vector<double> in)
{
    double result = 0.0;
    for(unsigned int i=0; i<Nt; ++i)
    {
        result += func(in.at(i));
    }
    result /= Nt;
    return result;
}
double xSquared(double x)
{
    return x * x;
}

int metropolis_step(std::vector<double> *x, int idx, double epsilon)
{
    double delta, r, h, xnew;
    delta = 2.0 * prng_double() - 1.0;
    xnew = x->at(idx) + delta * epsilon;
    r = exp(deltaS(*x, idx, xnew));
    h = prng_double();
    if(r >=1 || h <= r) {
        x->at(idx) = xnew;
        return 1;
    }
    else
        return 0;
}

double metropolis_sweep(std::vector<double> *x, double epsilon)
{
    int idx;
    double acceptancedUpdates = 0;
    for(uint j=0; j<Nt; ++j)
    {
        idx = round(prng_double() * (Nt - 1));
        acceptancedUpdates += metropolis_step(x, idx, epsilon);
    }
    return (double)acceptancedUpdates/Nt;
}

int main()
{
    double result, acceptanceRate, nominator, denominator;//, W0, fW0, f, WW, g, g2;
    
    int steps = 100000;
    int Ntherm = 1000;
    int Nskip = 0;
    double epsilon = 0.285;
    
    rlxd_init(1,time(NULL));
    std::vector<double> x;
    for(uint i=0; i<Nt; ++i) x.push_back(gaussianPRNG(1.0));
    std::vector<double> resultVector;

    result = 0.0; acceptanceRate = 0.0;
    
    for(unsigned int n=0; n<Ntherm; ++n)
    {
        metropolis_sweep(&x, epsilon);
    }
    for(unsigned int n=0; n<steps; ++n)
    {
        for(unsigned int j=0; j<Nskip; ++j)
        {
            metropolis_sweep(&x, epsilon);
        }
        acceptanceRate += metropolis_sweep(&x, epsilon);
        resultVector.push_back(expectationValue(xSquared, x));
        result += expectationValue(xSquared, x); //normal Monte Carlo
        
    }
    result /= steps; //normal Monte Carlo
    acceptanceRate /= steps;
    std::cout << "acceptanceRate = " << acceptanceRate << std::endl;
    std::cout << "<x^2> = " << result << "\n";
    std::cout << "analytic result = " << analyticResult() << "\n";
    std::cout << "delta = " << std::abs(result-analyticResult()) << "\n";
    //save(resultVector, "metro4.txt");
    
    return 0;
}



