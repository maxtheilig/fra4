#include <iostream>
#include <fstream>
#include "ex2.hpp"

//#define Nt 10

int main()
{
    double result, acceptanceRate;
    uint steps = 100000;
    double epsilon[4] = {0.1, 0.2, 0.3, 0.4};
    double omega[4] = {5.0, 2.5, 1.0, 0.5};
    const uint Nt[4] = {10, 20, 50, 100};
    for(uint j=0; j<1; ++j)
    {
        std::cout << "\n";
        print("Nt", Nt[j]);
        print("omega", omega[j]);
        double *startingPoint = new double[Nt[j]];
        double result_ana = analyticResult(omega[j], Nt[j]);
        for(uint i=0; i<Nt[j]; ++i) startingPoint[i] = gaussianPRNG(result_ana);
        Metropolis(&result, xSquared, steps, omega[j], epsilon[j], startingPoint, &acceptanceRate, Nt[j], false);
        delete[] startingPoint;
        //myfile << i << "\t" << result << "\n";
        print("result", result);
        print("analytic-result", result_ana);
        print("numerical - analytic", std::abs(result - result_ana));
        print("acceptanceRate", acceptanceRate);
        std::cout << "\n";
    }
    //myfile.close();
    return 0;
}

void Metropolis(double * result, double(*f)(double), uint steps, double omega, double epsilon, double *x, double * acceptanceRate, uint Nt, bool improvedDerivative)
{
    rlxd_init(1, time(NULL));
    double res_tmp, tmp, x0, x1, delta, x_new, h, r, S_old, S_new;
    res_tmp = 0.0;
    uint idx, idx_neigh, idx_lower_neigh, acceptedUpdates;
    acceptedUpdates = 0.0;
    std::ofstream myfile;
    myfile.open ("metro_autocor2.txt");
    for(uint n=0; n<steps; ++n)
    {
        for(uint j=0; j<Nt; ++j)
        {
            delta = 2.0 * prng_double() - 1.0;
            idx = round(prng_double() * (Nt - 1));
            if(improvedDerivative)
            {
                if(idx == Nt - 1) idx_neigh = 0;
                else idx_neigh = idx + 1;
                if(idx == 0) idx_lower_neigh = Nt - 1;
                else idx_lower_neigh = idx - 1;
                x_new = x[idx] + delta * epsilon;
                S_old = LocalAction(x[idx], x[idx_neigh], x[idx_lower_neigh], omega);
                S_new = LocalAction(x_new, x[idx_neigh], x[idx_lower_neigh], omega);
            }
            else
            {
                if(idx == Nt - 1) idx_neigh = 0;
                else idx_neigh = idx + 1;
                x_new = x[idx] + delta * epsilon;
                S_old = LocalAction(x[idx], x[idx_neigh], omega);
                S_new = LocalAction(x_new, x[idx_neigh], omega);
            }
            r = exp(S_old - S_new);
            if(r >=1) {
                x[idx] = x_new;
                ++acceptedUpdates;
                myfile << x[idx] << "\n";
            }
            else {
                h = prng_double();
                if(h < r){x[idx] = x_new; ++acceptedUpdates; myfile << x[idx] << "\n";}
            }
            tmp = 0.0;
            for(uint i=0; i<Nt; ++i)
            {
                tmp = tmp + f(x[i]);
            }
            res_tmp = res_tmp + tmp / Nt;
        }
    }
    *acceptanceRate = (double) acceptedUpdates / (double) (steps * Nt);
    *result = res_tmp / (double) (steps * Nt);
    myfile.close();
}

double LocalAction(double x0, double x1, double omega)
{
    double result, kineticTerm, potential;
    kineticTerm = (x1 - x0) * (x1 - x0);
    potential = omega * omega * x0 * x0;
    result = kineticTerm + potential;
    return result;
}

double LocalAction(double x0, double x1, double x2, double omega)
{
    double result, kineticTerm, potential;
    kineticTerm = ((x1 - x2) / 2.0) * ((x1 - x2) / 2.0);
    potential = omega * omega * x0 * x0;
    result = kineticTerm + potential;
    return result;
}

double prng_double()
{
    double tmp;
    ranlxd(&tmp, 1);
    return tmp;
}

double xSquared(double x)
{
    return x * x;
}

double analyticResult(double omega, uint Nt)
{
    double R = 1.0 + (omega * omega) / 2.0 - omega * sqrt(1.0 + (omega * omega) / 4.0);
    double tmp = 1.0 / (4.0 * omega * sqrt(1.0 + (omega * omega) / 4.0)) * (1.0 + pow(R, Nt)) / (1.0 - pow(R, Nt));
    return tmp;
}

double gaussianPRNG(double s)
{
    double a = s * sqrt(-2.0 * log(prng_double()));
    double b = cos(2 * 3.141592653589793238463 * prng_double());
    return a * b;
}


