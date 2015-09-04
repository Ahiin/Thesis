#pragma once

#include <functional>

class PotentialTrace
{
public:
    virtual double operator()(double x) = 0;
};