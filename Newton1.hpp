#ifndef NEWTON1_HPP
#define NEWTON1_HPP

#include "BaseUnaryNonLinearSolver.hpp"
#include <complex>

using namespace std;

class Newton1:public BaseUnaryNonLinearSolver
{
public:
    Newton1(double toler = 1e-06,unsigned int iterlimit = 50):BaseUnaryNonLinearSolver(toler,iterlimit){}    
    double solve(double);
    virtual double df(double) = 0;
};

class Complex_Newton1:public BaseUnaryNonLinearSolver
{
public:
    Complex_Newton1(double toler = 1e-08,unsigned int iterlimit = 50):BaseUnaryNonLinearSolver(toler,iterlimit){}        
    complex<double> solve(complex<double>);
    double f(double){return 0.0;}
    virtual complex<double> df(complex<double>) = 0;
    virtual complex<double> f(complex<double>) = 0;
};

#endif
