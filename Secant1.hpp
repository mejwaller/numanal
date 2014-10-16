#ifndef SECANT1_HPP
#define SECANT1_HPP

#include "BaseUnaryNonLinearSolver.hpp"

class Secant1:public BaseUnaryNonLinearSolver
{
public:
    Secant1(double toler = 1e-08,unsigned int iterlimit = 50):BaseUnaryNonLinearSolver(toler,iterlimit){}
    double solve(double,double);

};

#endif
