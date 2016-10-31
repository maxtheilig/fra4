#ifndef TOOLS_HPP
#define TOOLS_HPP

#include <iostream>
#include <string>

typedef unsigned int uint;

void print(std::string name, double x)
{
    std::cout << name << " = " << x << "\n";
}

void print(std::string name, int x)
{
    std::cout << name << " = " << x << "\n";
}

void print(std::string name, uint x)
{
    std::cout << name << " = " << x << "\n";
}

#endif
