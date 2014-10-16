#include "Bairstow1.hpp"
#include "BaseUnaryNonLinearSolver.hpp"
#include <iostream>

vector< vector<float> > Bairstow1::factorPoly(const vector<float> &acoeffs,float rest,float sest,float toler,int iterlimit)
{    
    r = rest;
    s = sest;
    tol = toler;
    iterlim = iterlimit;    
    a=acoeffs;
    n = a.size()-1;    
    vector< vector<float> > nil(0);
    if(n<3)
    {
        cout << "Polynomial of degree 3 or higher required! Aborting." << endl;        
        return nil;
    }
    b.resize(a.size());
    c.resize(a.size());
    for(unsigned long int j= 0; j<c.size(); j++)
    {
        b[j] = 0.0;
        c[j]= 0.0;
    }

    cout << "Original polynomial is:" << endl;
    for(int i = 0; i <= n; i++)
    {
        cout << a[i] << "x^" << n-i << " + ";
    }
    cout << endl;

    for(;;)
    {
        computeBandCarrays();      
        denom = c[n-2]*c[n-2] - c[n-1]*c[n-3];
        if(!FCMP(denom,0.0))
        {
            findRandS();
            checkTolMet();
            if(flag == true)
            {
                return resvec;
            }
        }       
        else
        {
            changeRandSandK();
        }
        if(k>iterlim)
        {
            cout << "Number of iterations exceeded...roots not found." << endl;
            return nil;
        }
    }    

    return resvec;
}

void Bairstow1::computeBandCarrays()
{
    b[0] = c[0] = a[0];
    b[1] = a[1] + r*b[0];
    c[1] = b[1] + r*c[0];
    for(int i = 2; i <= n; i++)
    {
        b[i] = a[i]+r*b[i-1]+s*b[i-2];
        c[i] = b[i]+r*c[i-1] + s*c[i-2];
    }
}

void Bairstow1::findRandS()
{
    delR = (-b[n-1]*c[n-2] - -b[n]*c[n-3])/denom;
    delS = (c[n-2]*-b[n] - c[n-1]*-b[n-1])/denom;    
    r+=delR;
    s+=delS;
}

void Bairstow1::checkTolMet()
{    
    if(fabs(delR) + fabs(delS) < tol)
    {        
        //cout << endl << "Factors are: " << endl;
        //cout << "x^2 - " << r << "x - " << s << endl;
        vector<float> fac1(3);
        vector<float> fac3(2);
        vector<float> fac2(3);
        fac1.push_back(1.0);
        fac1.push_back(r);
        fac1.push_back(s);
        resvec.push_back(fac1);
        n-=2;        
        switch(n)
        {
        case -1:
            flag = true;
            break;
        case 0:
            flag = true;
            break;
        case 1:
            //cout << b[n-1] << "x + " << b[n] << endl;                        
            fac3.push_back(b[n-1]);
            fac3.push_back(-b[n]);
            resvec.push_back(fac3);
            flag = true;
            break;
        case 2:
            //cout << b[n-2] << "x^2 + " << b[n-1] << "x + " << b[n] << endl;            
            fac2.push_back(b[n-2]);
            fac2.push_back(-b[n-1]);
            fac2.push_back(-b[n]);
            resvec.push_back(fac2);
            flag = true;
            break;
        default:
            reducePolynomial();
            break;
        }

    }
    else
    {
        k+=1;
    }    
}

void Bairstow1::reducePolynomial()
{    
    k = 1;
    for(int i = 2;i<=n+2;i++)
    {
        a[i-2] = b[i-2];
    }
}

void Bairstow1::changeRandSandK()
{
    r+=1.0;
    s+=1.0;
    k=1;
}




        

