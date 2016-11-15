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

double W(std::vector<double> x, double lambda1, double lambda2, double lambda10, double lambda20)
{
    double result; double x2 = 0.0; double x4 = 0.0;
    for(unsigned int i=0; i<x.size(); ++i)
    {
        x2 += pow(x.at(i),2); x4 += pow(x.at(i),4);
    }
    result = (lambda10 - lambda1)*x2 + (lambda20 - lambda2)*x4;
    return result;
}

void metropolisReweighting(double * result, double * acceptanceRate, std::vector<double> x, double(*func)(double), int steps, int Ntherm, double lambda1, double lambda2, double lambda10, double lambda20, double epsilon)
{
    double res_tmp, x0, x1, delta, xnew, h, r, nominator, denominator;
    res_tmp = 0.0; nominator = 0.0; denominator = 0.0;
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
        if( (n+1)%100 == 0 )
        {
            nominator += expectationValue(func, x) * W(x,lambda1,lambda2,lambda10,lambda20);
            //nominator += pow(W(x,lambda1,lambda2,lambda10,lambda20),3);
            denominator += W(x,lambda1,lambda2,lambda10,lambda20);
        }
        
    }
    *acceptanceRate = (double) acceptedUpdates / (double) (steps * Nt);
    *result = nominator / denominator;
}

void metropolisReweighting2(double * result, double * acceptanceRate, std::vector<double> x, double(*func)(double), int steps, int Ntherm, double lambda1, double lambda2, double lambda10, double lambda20, double epsilon)
{
    double res_tmp, x0, x1, delta, xnew, h, r, nominator, denominator;
    res_tmp = 0.0; nominator = 0.0; denominator = 0.0;
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
        //nominator += expectationValue(func, x) * W(x,lambda1,lambda2,lambda10,lambda20);
        nominator += W(x,lambda1,lambda2,lambda10,lambda20) * W(x,lambda1,lambda2,lambda10,lambda20);
        denominator += W(x,lambda1,lambda2,lambda10,lambda20);
        
    }
    *acceptanceRate = (double) acceptedUpdates / (double) (steps * Nt);
    *result = nominator / denominator;
    *result *= *result;
}


int main()
{
    double result1, result2, acceptanceRate1, acceptanceRate2;
    int steps = 100000;
    int Ntherm = 10000;
    //int Nt = 10;
    std::vector<double> epsilon = {0.13, 0.22, 0.27, 0.285, 0.285, 0.285};
    std::vector<double> lambda1 = {-10.0, -3.0, -2.5, -2.0, -1.0, 1.0}; //double lambda1 = 1.0;
    double lambda2 = 1.0;
    double lambda10 = 0.0;
    double lambda20 = 1.0;
    
    rlxd_init(1,time(NULL));
    std::vector<double> x;
    for(uint i=0; i<Nt; ++i) x.push_back(gaussianPRNG(1.0));
    std::vector<double> resultVector;
    
    for(unsigned int i=0; i<lambda1.size(); ++i)
    {
        std::cout << "lambda1 = " << lambda1.at(i) << std::endl;
        metropolisReweighting(&result1, &acceptanceRate1, x, xSquared, steps, Ntherm, lambda1.at(i), lambda2, lambda10, lambda20, epsilon.at(i));
        //metropolisReweighting2(&result2, &acceptanceRate1, x, xSquared, steps, Ntherm, lambda1.at(i), lambda2, lambda10, lambda20, epsilon.at(i));
        resultVector.push_back(result1);
        std::cout << "result = " << result1 << std::endl;
        std::cout << "resultvector = " << resultVector.at(i) << std::endl;
        std::cout << "acceptanceRate1 = " << acceptanceRate1 << std::endl;
        std::cout << "acceptanceRate2 = " << acceptanceRate2 << std::endl;
        std::cout << "\n";
    }
    save(lambda1, resultVector, "xsquared_rew.txt");
    
    return 0;
}


