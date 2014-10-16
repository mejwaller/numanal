#include "Muller1.hpp"
#include <iostream>
#include <math.h>

double Muller1::solve(double x,double y,double z)
{
    double x0 = x;
    double x1 = y;
    double x2 = z;
    double f0 = f(x0);
    double f1 = f(x1);
    double f2 = f(x2);
    double h1 = x1-x0;
    double h2 = x0-x2;
    double h2overh1 = h2/h1;
    double a = (h2overh1*f1 - f0*(1.0+h2overh1)+f2)/ (h2overh1*h1*h1*(1.0+h2overh1));
    double b = (f1-f0-a*h1*h1)/h1;
    double c = f0;
    double xr;
    double toppart = b*b - 4.0*a*c;
    double retval;
    unsigned int i = 0;

    //cout << "f(x0)\tf(x1)\tf(x2)\th1\th2\th2overh1\ta\tb" << endl;

    for(;;)
    {
        //cout << f0 << "\t" << f1 << "\t" << f2 << "\t" << h1 << "\t" << h2 << "\t" << h2overh1 << "\t" << a << "\t" << b << endl;
        i++;
        if(toppart < 0.0)
        {
            cout << "Muller1 - error - complex root encountered. Aborting." << endl;
            return 0.0;
        }
        else
        {
            if(b > 0.0)
            {
                xr = x0 - (2.0*c/ (b+sqrt(toppart)));
            }
            else
            {
                xr = x0 - (2.0*c/(b-sqrt(toppart)));
            }

            retval = f(xr);

            //cout << "Current root estimate is " << xr << " at iteration " << i << endl;

            if(fabs(retval) < tol)
            {
                cout << "Muller1: Root found to required tolerance = " << xr << endl;
                return xr;
            }
            else
            {
                if(i > iterlim)
                {
                    cout << "Max iteration limit exceeded - aborting with root estimate at " << xr <<  endl;
                    return xr;
                }
                if(xr > x0)
                {
                    x2 = x0;
                    x0 = xr;
                }
                else
                {
                    x1 = x0;
                    x0 = xr;
                }
            }

            f0 = f(x0);
            f1 = f(x1);
            f2 = f(x2);
            h1 = x1-x0;
            h2 = x0-x2;
            h2overh1 = h2/h1;
            a = (h2overh1*f1 - f0*(1.0+h2overh1)+f2)/ (h2overh1*h1*h1*(1.0+h2overh1));
            b = (f1-f0-a*h1*h1)/h1;
            c = f0;            
            toppart = b*b - 4.0*a*c;
        }
    }
}

complex<double> Complex_Muller1::solve(complex<double> x,complex<double> y,complex<double> z)
{
    complex<double> x0 = x;
    complex<double> x1 = y;
    complex<double> x2 = z;
    complex<double> f0 = f(x0);
    complex<double> f1 = f(x1);
    complex<double> f2 = f(x2);
    complex<double> h1 = x1-x0;
    complex<double> h2 = x0-x2;
    complex<double> h2overh1 = h2/h1;
    complex<double> a = (h2overh1*f1 - f0*(1.0+h2overh1)+f2)/ (h2overh1*h1*h1*(1.0+h2overh1));
    complex<double> b = (f1-f0-a*h1*h1)/h1;
    complex<double> c = f0;
    complex<double> xr;
    complex<double> toppart = b*b - 4.0*a*c;
    complex<double> retval;
    unsigned int i = 0;

    //cout << "f(x0)\tf(x1)\tf(x2)\th1\th2\th2overh1\ta\tb" << endl;

    for(;;)
    {
        //cout << f0 << "\t" << f1 << "\t" << f2 << "\t" << h1 << "\t" << h2 << "\t" << h2overh1 << "\t" << a << "\t" << b << endl;
        i++;        

        if(abs(b) > 0.0)
        {
            xr = x0 - (2.0*c/ (b+sqrt(toppart)));
        }
        else
        {
            xr = x0 - (2.0*c/(b-sqrt(toppart)));
        }

        retval = f(xr);

        //cout << "Current root estimate is " << xr << " at iteration " << i << endl;

        if(fabs(retval.real()) < tol && fabs(retval.imag()) < tol)
        {
            cout << "Complex_Muller1: Root found to required tolerance = " << xr << endl;
            return xr;
        }
        else
        {
            if(i > iterlim)
            {
                cout << "Max iteration limit exceeded - aborting with root estimate at " << xr <<  endl;
                return xr;
            }
            if(abs(xr) > abs(x0))
            {
                x2 = x0;
                x0 = xr;
            }
            else
            {
                x1 = x0;
                x0 = xr;
            }
        }

        f0 = f(x0);
        f1 = f(x1);
        f2 = f(x2);
        h1 = x1-x0;
        h2 = x0-x2;
        h2overh1 = h2/h1;
        a = (h2overh1*f1 - f0*(1.0+h2overh1)+f2)/ (h2overh1*h1*h1*(1.0+h2overh1));
        b = (f1-f0-a*h1*h1)/h1;
        c = f0;            
        toppart = b*b - 4.0*a*c;
    }    

}
