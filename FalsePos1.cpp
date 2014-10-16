#include "FalsePos1.hpp"
#include <iostream>
#include <math.h>

using namespace std;

double FalsePos1::solve(double xlower, double xupper)
{

    double x0 = xlower,x1 = xupper, x2;
    unsigned int iteration = 1;
    if( fabs(f(x0)) < fabs(f(x1)) )
    {
        x2 = x0;
        x0 = x1;
        x1 = x2;        
    }
    cout << "iteration\tx0\tx1\tx2\tf(x2)" << endl;
    for(;;)
    {        
        x2 = x0 - f(x0)*(x0-x1)/(f(x0)-f(x1));
        if(f(x2)*f(x0) < 0.0)
        {
            x1 = x2;
        }
        else
        {
            x0 = x2;
        }        
        cout << iteration << "\t\t" << x0 << "\t" << x1 << "\t" << x2 << "\t" << f(x2) << endl;
        if(fabs(f(x2)) < tol)
        {
            cout << "FalsePos1: Root found to specified tolerance of " << tol << endl;
            cout << "Root is " << x2 << endl;
            return x2;
        }
        if(iteration > iterlim)
        {
            cout << "Iteration limit exceeded...root not found to required tolerance." << endl;            
            cout << "Root is " << x2 << endl;
            return  x2;
        }

        iteration+=1;
    }
}
