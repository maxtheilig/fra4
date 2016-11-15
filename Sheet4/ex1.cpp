#include <iostream>
#include <fstream>
#include <vector>
#include "../random.cpp"

#define Nt 10

void save(std::vector<double> in, std::string filename)
{
    std::ofstream out;
    out.open (filename);
    for(unsigned int i=0; i<100; ++i)
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

double deltaS(std::vector<double> x, int idx, double xnew, double lambda1, double lambda2)
{
    int idx_neigh, idx_lower_neigh;
    if(idx == Nt - 1) idx_neigh = 0;
    else idx_neigh = idx + 1;
    if(idx == 0) idx_lower_neigh = Nt - 1;
    else idx_lower_neigh = idx - 1;
    double tmp = (2 + lambda1)*(pow(x[idx],2) - pow(xnew,2)) + lambda2*(pow(x[idx],4) - pow(xnew,4)) + 2.0 * (x.at(idx_neigh) + x.at(idx_lower_neigh))*(x.at(idx) - xnew);
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

double xFunc(double x)
{
    return x;
}

double xAbsolute(double x)
{
    double result;
    if(x<0) result = -x;
    else result = x;
    return result;
}

double xSquared(double x)
{
    return x * x;
}

double tunneling(std::vector<double> x)
{
    int tmp = 0;
    for(unsigned int i=0; i<(x.size()-1); ++i)
    {
        if((x.at(i)<0.0 && x.at(i+1)>0.0) || (x.at(i)>0.0 && x.at(i+1)<0.0))
        {
            ++tmp;
            //std::cout << "++tunneling events at i=" << i << "\t#tunneling=" << tmp << "\n";
        }
    }
    if((x.at(x.size()-1)<0 && x.at(0)>0) || (x.at(x.size()-1)>0 && x.at(0)<0))
    {
        ++tmp;
        //std::cout << "++tunneling events at i=" << x.size() << "\t#tunneling=" << tmp << "\n";
    }
    double result = (double) tmp / (double) x.size();
    return result;
}


void metropolis(double * result, double * acceptanceRate, std::vector<double> x, double(*func)(double), int steps, int Ntherm, double lambda1, double lambda2, double epsilon)
{
    double res_tmp, x0, x1, delta, xnew, h, r;
    res_tmp = 0.0;
    int idx, acceptedUpdates;
    acceptedUpdates = 0.0;
    for(uint n=0; n<Ntherm; ++n) //thermalization
    {
        for(uint j=0; j<Nt; ++j)
        {
            delta = 2.0 * prng_double() - 1.0;
            idx = round(prng_double() * (Nt - 1));
            xnew = x[idx] + delta * epsilon;
            r = exp(deltaS(x, idx, xnew, lambda1, lambda2));
            h = prng_double();
            if(r >=1 || h <= r) {
                x[idx] = xnew;
            }
        }
    }
    for(uint n=0; n<steps; ++n)
    {
        for(uint j=0; j<Nt; ++j)
        {
            delta = 2.0 * prng_double() - 1.0;
            idx = round(prng_double() * (Nt - 1));
            xnew = x[idx] + delta * epsilon;
            r = exp(deltaS(x, idx, xnew, lambda1, lambda2));
            h = prng_double();
            if(r >=1 || h <= r) {
                x[idx] = xnew; ++acceptedUpdates;
            }
            
        }
        //res_tmp += expectationValue(func, x);
        res_tmp += tunneling(x) * tunneling(x);
        //std::cout << "res_tmp = " << res_tmp << "\n";
    }
    *acceptanceRate = (double) acceptedUpdates / (double) (steps * Nt);
    *result = res_tmp / (double) (steps);
}

int main()
{
    double result, acceptanceRate;
    int steps = 100000;
    int Ntherm = 10000;
    //int Nt = 10;
    std::vector<double> epsilon = {0.13, 0.22, 0.27, 0.285, 0.285, 0.285};
    std::vector<double> lambda1 = {-10.0, -3.0, -2.5, -2.0, -1.0, 1.0};//double lambda1 = 1.0;
    double lambda2 = 1.0;
    
    rlxd_init(1,time(NULL));
    std::vector<double> x;
    for(uint i=0; i<Nt; ++i) x.push_back(gaussianPRNG(1.0));
    std::vector<double> resultVector;
    
    for(unsigned int i=0; i<lambda1.size(); ++i)
    {
        std::cout << "lamda1 = " << lambda1.at(i) << std::endl;
        metropolis(&result, &acceptanceRate, x, xSquared, steps, Ntherm, lambda1.at(i), lambda2, epsilon.at(i));
        resultVector.push_back(result);
        std::cout << "result = " << result << std::endl;
        std::cout << "resultvector = " << resultVector.at(i) << std::endl;
        std::cout << "acceptanceRate = " << acceptanceRate << std::endl;
        std::cout << "\n";
    }
    save(lambda1, resultVector, "test.txt");
    
    return 0;
}


