#ifndef FALSEPOS1_HPP
#define FALSEPOS1_HPP

#include "BaseUnaryNonLinearSolver.hpp"

class FalsePos1:public BaseUnaryNonLinearSolver
{
public:
    FalsePos1(double toler = 1e-08,unsigned int iterlimit = 50):BaseUnaryNonLinearSolver(toler,iterlimit){}
    double solve(double,double);

};

#endif
