#pragma once

class EnergySpace
{
public:
    virtual double operator()(double a1, double b1, double a2, double b2) = 0;
};