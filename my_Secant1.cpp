#include "Secant1.hpp"
#include <math.h>

class my_Secant1:public Secant1
{
public:
    my_Secant1(double tol = 1e-08,unsigned int it = 100):Secant1(tol,it){}
protected:
    virtual double f(double);
};

double my_Secant1::f(double x)
{    
    //double retval = x*x*x + x*x - 3.0*x - 3.0;
    double retval = 3.0*x + sin(x) - exp(x);
    return retval;
}


int main( int argc, char** argv )
{
    my_Secant1 secant;
    secant.solve(0.0,1.0);
    return 0;
}
