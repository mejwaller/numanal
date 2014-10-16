#include "Bisection1.hpp"
#include <math.h>

class my_func:public Bisection1
{
public:
    my_func(double tol = 1e-07,unsigned int it = 100):Bisection1(tol,it){}
protected:
    virtual double f(double);
};

double my_func::f(double x)
{
    /*double a = 123*3.141592654/180.0;
    double retval = 9.0*cos(a+x)/pow(sin(a+x),2.0) + 7.0*cos(x)/pow(sin(x),2.0);
    */
    double retval = x*x*x + x*x - 3.0*x - 3.0;
    return retval;
}


int main( int argc, char** argv )
{
    my_func fn;
    fn.solve(1.0,2.0);
    return 0;
}

