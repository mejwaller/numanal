#ifndef MULLER1_HPP
#define MULLER1_HPP

#include "BaseUnaryNonLinearSolver.hpp"
#include <complex>

using namespace std;

class Muller1:public BaseUnaryNonLinearSolver
{
public:
    Muller1(double toler = 1e-06,unsigned int iterlimit = 50):BaseUnaryNonLinearSolver(toler,iterlimit){}    
    double solve(double,double,double);
};

class Complex_Muller1:public BaseUnaryNonLinearSolver
{
public:
    Complex_Muller1(double toler = 1e-08,unsigned int iterlimit = 50):BaseUnaryNonLinearSolver(toler,iterlimit){}    
    complex<double> solve(complex<double>,complex<double>,complex<double>);
    double f(double){return 0;}
    virtual complex<double> f(complex<double>) = 0;
};


#endif
