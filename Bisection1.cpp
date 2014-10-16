#include "Bisection1.hpp"
#include <iostream>
#include <math.h>

using namespace std;

double Bisection1::solve(double xlower, double xupper)
{
    double f1,f2,x1 = xlower,x2 = xupper;
    double error;
    unsigned int iteration = 1;    
    cout << "iteration\tx1\tx2\tx3\tf(x3)" << endl;
    for(;;)    
    {
        f1 = f(x1);
        f2 = f(x2);
        if(f1*f2 > 0.0)
        {
            cout << "Values do not bracket a root! Aborting" << endl;
            return 0.0;
        }
        else
        {
            double x3 = (x1+x2)/2.0;
            double f3 = f(x3);
            cout << iteration << "\t\t" << x1 << "\t" << x2 << "\t" << x3 << "\t" << f3 << endl;
            if(f3*f1 < 0.0)
            {
                x2 = x3;
            }
            else
            {
                x1 = x3;
            }
            
            if( fabs(x1-x2)/2.0 < tol)
            {
                cout << "Root found to specified tolerance of " << tol << endl;
                error = (1.0/pow(2.0,(double)iteration))*fabs(xupper-xlower);
                cout << "Root is " << x3 << " +/- " << error << endl;

                return x3;
            }
            if(FCMP(f3,0.0))
            {
                cout << "Root found as f(root) not distinguishable from zero" << endl;
                error = (1.0/pow(2.0,(double)iteration))*fabs(xupper-xlower);
                cout << "Root is " << x3 << " +/- " << error << endl;
                return x3;
            }
            if(iteration > iterlim)
            {
                cout << "Iteration limit exceeded...root not found to required tolerance." << endl;
                error = (1.0/pow(2.0,(double)iteration))*fabs(xupper-xlower);
                cout << "Root is " << x3 << " +/- " << error << endl;
                return  x3;
            }

            iteration+=1;
        }
    }
}
    


