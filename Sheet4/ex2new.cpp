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

int metropolis_step(std::vector<double> *x, int idx, double epsilon, double lambda1, double lambda2)
{
    double delta, r, h, xnew;
    delta = 2.0 * prng_double() - 1.0;
    xnew = x->at(idx) + delta * epsilon;
    r = exp(deltaS(*x, idx, xnew, lambda1, lambda2));
    h = prng_double();
    if(r >=1 || h <= r) {
        x->at(idx) = xnew;
        return 1;
    }
    else
        return 0;
}

double metropolis_sweep(std::vector<double> *x, double epsilon, double lambda1, double lambda2)
{
    int idx;
    double acceptancedUpdates = 0;
    for(uint j=0; j<Nt; ++j)
    {
        idx = round(prng_double() * (Nt - 1));
        acceptancedUpdates += metropolis_step(x, idx, epsilon, lambda1, lambda2);
    }
    return (double)acceptancedUpdates/Nt;
}

double errorReweighting(std::vector<std::vector<double>> configs, double (*func)(double), double lambda1, double lambda2, double lambda10, double lambda20)
{
    int steps = configs.size();
    std::cout << "steps in error function = " << steps << "\n";
    double W0, fW0, f, WW, g, g2;
    W0 = 0.0; fW0 = 0.0;
    for(unsigned int n=0; n<steps; ++n)
    {
        W0 += W(configs.at(n), lambda1, lambda2, lambda10, lambda20);
        fW0 += expectationValue(func, configs.at(n)) * W(configs.at(n), lambda1, lambda2, lambda10, lambda20);
    }
    // calculate g = f * W / W0 - W * fW0 / pow(W0, 2) for error estimation
    W0 /= steps;
    fW0 /= steps;
    for(unsigned int n=0; n<steps; ++n)
    {
        f = expectationValue(func, configs.at(n));
        WW = W(configs.at(n), lambda1, lambda2, lambda10, lambda20);
        g += f * WW / W0 - WW * fW0 / pow(W0, 2);
        g2 += pow(g, 2);
    }
    g /= steps; //expectation value of g
    g2 /= steps; //expectation value of g^2
    return sqrt( (1.0 / (steps - 1)) * (g2 - pow(g, 2)));
}

int main()
{
    double result, acceptanceRate, nominator, denominator;//, W0, fW0, f, WW, g, g2;
    
    int steps = 100000;
    int Ntherm = 10000;
    int Nskip = 10;
    
    std::vector<double> epsilon = {0.13, 0.22, 0.27, 0.285, 0.285, 0.285};
    std::vector<double> lambda1 = {-10.0, -3.0, -2.5, -2.0, -1.0, 1.0}; //double lambda1 = 1.0;
    double lambda2 = 1.0;
    double lambda10 = 0.0;
    double lambda20 = 1.0;
    
    rlxd_init(1,time(NULL));
    std::vector<double> x;
    for(uint i=0; i<Nt; ++i) x.push_back(gaussianPRNG(1.0));
    std::vector<double> resultVector, error;
    std::vector<std::vector<double>> configurations(steps, std::vector<double>(Nt));
    
    for(unsigned int i=0; i<lambda1.size(); ++i)
    {
        result = 0.0; acceptanceRate = 0.0; nominator = 0.0; denominator = 0.0; //W0 = 0.0; fW0 = 0.0;
        std::cout << "lambda1 = " << lambda1.at(i) << std::endl;
        for(unsigned int n=0; n<Ntherm; ++n)
        {
            metropolis_sweep(&x, epsilon.at(i), lambda1.at(i), lambda2);
        }
        for(unsigned int n=0; n<steps; ++n)
        {
            for(unsigned int j=0; j<Nskip; ++j){
                metropolis_sweep(&x, epsilon.at(i), lambda1.at(i), lambda2);}
            acceptanceRate += metropolis_sweep(&x, epsilon.at(i), lambda1.at(i), lambda2);
            //store configurations
            configurations.at(n) = x;
            //measurement of observable
            //result += expectationValue(xSquared, x); //normal Monte Carlo
            nominator += expectationValue(xSquared, x) * W(x, lambda1.at(i), lambda2, lambda10, lambda20); //reweighting
            denominator += W(x, lambda1.at(i), lambda2, lambda10, lambda20); //reweighting
            //W0 += W(x, lambda1.at(i), lambda2, lambda10, lambda20);
            //fW0 += expectationValue(xSquared, x) * W(x, lambda1.at(i), lambda2, lambda10, lambda20);
        }
        //result /= steps; //normal Monte Carlo
        result = nominator / denominator; //reweighting
        resultVector.push_back(result);
        acceptanceRate /= steps;
        std::cout << "result = " << result << std::endl;
        std::cout << "acceptanceRate = " << acceptanceRate << std::endl;
        
        // calculate g = f * W / W0 - W * fW0 / pow(W0, 2) for error estimation
        /*W0 /= steps; fW0 /= steps;
        for(unsigned int n=0; n<steps; ++n)
        {
            f = expectationValue(xSquared, configurations.at(n));
            WW = W(configurations.at(n), lambda1.at(i), lambda2, lambda10, lambda20);
            g += f * WW / W0 - WW * fW0 / pow(W0, 2);
            g2 += pow(g, 2);
        }
        g /= steps; //expectation value of g
        g2 /= steps; //expectation value of g^2
        error.push_back(sqrt( (1.0 / (steps - 1)) * (g2 - pow(g, 2))));*/
        error.push_back(errorReweighting(configurations, xSquared, lambda1.at(i), lambda2, lambda10, lambda20));
        std::cout << "error = " << error.at(i) << "\n";

        std::cout << "\n";
    }
    //std::vector<double> z;
    //for(uint i=0; i<Nt; ++i) z.push_back(prng_double());
    save(lambda1, resultVector, error, "testerror.txt");

    return 0;
}












