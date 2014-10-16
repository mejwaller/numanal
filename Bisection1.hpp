#ifndef BISECTION1_HPP
#define BISECTION1_HPP

#include "BaseUnaryNonLinearSolver.hpp"

class Bisection1:public BaseUnaryNonLinearSolver
{
public:
    Bisection1(double toler = 1e-08,unsigned int iterlimit = 50):BaseUnaryNonLinearSolver(toler,iterlimit){}
    double solve(double,double);
};

#endif



