#include "Newton1.hpp"
#include <iostream>
#include <math.h>

double Newton1::solve(double xguess)
{

    double x0 = xguess, x1 = x0;
    double fx = f(x0);
    double dfx = df(x0);    
    unsigned int iteration = 1;
    //cout << "iteration\tx0\tx1\tf(x1)" << endl;
    for(;;)
    {
        //cout << iteration-1 << "\t\t" << x0 << "\t" << x1 << "\t" << fx << endl;
        if(FCMP(dfx,0.0))
        {
            cout << "Newton1: Can't find root - f'(" << x1 << ") == 0 at iteration " << iteration-1 << "!" << endl;
            return 0.0;
        }
        if(FCMP(fx,0.0))
        {
            cout << "Newton1: Root found after " << iteration-1 << "iterations (f(" << x1 << ") == 0):" << x0 << endl;
            return x1;
        }        
        x0 = x1;
        fx = f(x0);
        dfx = df(x0);
        x1 = x0 - fx/dfx;
        if(fabs(x0-x1) < tol)
        {
            cout << "Newton1: Root found at iteration " << iteration-1 << ": " << x1 << endl;
            return x1;
        }
        iteration+=1;
    }    
    
}

complex<double> Complex_Newton1::solve(complex<double> xguess)
{
    complex<double> x0 = xguess, x1 = x0;
    complex<double> fx = f(x0);
    complex<double> dfx = df(x0);
    unsigned int iteration = 1;
    cout << "iteration\tx0\tx1\tf(x1)" << endl;
    for(;;)
    {
        //cout << iteration-1 << "\t\t" << x0 << "\t" << x1 << "\t" << fx << endl;
        if(FCMP(abs(dfx),0.0))
        {
            cout << "Newton1: Can't find root - f'(" << x1 << ") == 0 at iteration " << iteration-1 << "!" << endl;
            return 0.0;
        }
        if(FCMP(abs(fx),0.0))
        {
            cout << "Newton1: Root found after " << iteration-1 << "iterations (f(" << x1 << ") == 0):" << x0 << endl;
            return x1;
        }        
        x0 = x1;
        fx = f(x0);
        dfx = df(x0);
        x1 = x0 - fx/dfx;
        if( fabs(x0.real()-x1.real()) < tol && fabs(x0.imag() - x1.imag()) < tol)
        {
            cout << "Newton1: Root found at iteration " << iteration-1 << ": " << x1 << endl;
            return x1;
        }
        iteration+=1;
    }    
    
}
