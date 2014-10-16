#include "Horner1.hpp"
#include <iostream>
#include <math.h>

/** Divides the polynomial with coefficients in a,  by x-x1, where x1 is an input value.  
    Coefficents of reduced polynomial returned in vector b. */
vector<double> Horner1::dividePoly(vector<double> a, double x1)
{
    vector<double> b(a.size());
    
    b[0] = a[0];
    for(unsigned long int i = 1; i<a.size();i++)
    {
        b[i] = a[i]+b[i-1]*x1;
    }

    cout << "Original polynomial of: ";

    for(unsigned long int j = 0;j<a.size()-1;j++)
    {
        cout << a[j] << "x^" << a.size()-j-1 << " ";
    }

    cout << a[a.size()-1] << endl;

    cout << endl;

    cout << "divided by x - " << x1 << " is" << endl;

    for(unsigned long int k = 0; k < b.size()-2; k++)
    {
        cout << b[k] << "x^" << b.size()-k-2 << " ";
    }
    cout << b[b.size()-2] << endl;
    cout << "with remainder " << b[b.size()-1] << endl;

    double res = 0.0;    

    //vector<double> xpows(a.size());
    

    for(unsigned long int l = 0; l < a.size();l++)
    {
        res += a[l]*pow(x1,(double)a.size()-l-1);
        cout << a[l] << " * " << pow(x1,(double)a.size()-l-1) << endl;
    }

    cout << "P(x1) is: " << res << endl;

    return b;
}


