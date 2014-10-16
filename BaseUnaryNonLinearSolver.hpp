#ifndef BASEUNARYNONLINEARSOLVER_HPP
#define BASEUNARYNONLINEARSOLVER_HPP

#include <float.h>
#include <math.h>

#ifdef FCMP
#   error "FCMP already defined"
#else
#   define FCMP(aa, bb) ((fabs((double)(aa) - (double)(bb))/\
           (FLT_EPSILON+fabs((double)(aa))+fabs((double)(bb))) < FLT_EPSILON))
#endif 

class BaseUnaryNonLinearSolver
{
public:
    BaseUnaryNonLinearSolver(double toler = 1e-08,unsigned int iterlimit = 50):tol(toler),iterlim(iterlimit){}
    virtual double solve(double,double){ return 0.0; };
    virtual double solve(double){ return 0.0; };
protected:
    virtual double f(double) = 0;    
    double tol;
    unsigned int iterlim;

};
#endif





















