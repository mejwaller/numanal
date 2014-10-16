#include "Interpolate.hpp"
#include "Matrix1.hpp"
#include <iostream>

#include <stdio.h>

float LagrangeInterpolate(float val,const vector<float> &x, const vector<float> &fx)
{
    if(x.size() != fx.size())
    {
        cout << "LagrangeInterpolate: Error - x and fx vectors different sizes!" << endl;
        return 0.0;
    }
    float sum = 0.0,p;
    for(unsigned long int i =0; i<x.size();i++)
    {
        p=1.0;
        for(unsigned long int j = 0; j<x.size(); j++)
        {
            if(j!=i)
            {
                p*=((val-x[j])/(x[i]-x[j]));
            }            
        }
        sum+=(p*fx[i]);
    }
    return sum;
}

DividedDiffInterpolater::DividedDiffInterpolater(const vector<float> &xs, const vector<float> &fx)
{
    
    a = fx;
    x = xs;
    float temp1,temp2;

    for(unsigned long int j = 1; j < a.size(); j++)
    {
        temp1 = a[j-1];
        for(unsigned long int k = j; k<a.size();k++)
        {         
            temp2 = a[k];         
            a[k] = (a[k] - temp1)/(xs[k]-xs[k-j]);
            temp1 = temp2;            
        }
    }    
    
}

float DividedDiffInterpolater::interpolate(const float &val)
{
    float sum = 0.0;
    
    for(unsigned long int i = a.size()-1; i>=1; i--)
    {
        if(i <= 0)
        {
            break;
        }
        //printf("\na.size() is %d, i = %d, a[%d] = %f, x[%d] = %f\n",a.size(),i,i,a[i],i-1,x[i-1]);
        sum = (sum+a[i])*(val - x[i-1]);
    }
    sum+=a[0];
    return sum;
}

CubicSplineInterpolater::CubicSplineInterpolater(const vector<float> &xs, const vector <float> &ys,unsigned int condition)
{
    if(xs.size() != ys.size())
    {
        cout << "Error - number of data points != number of function points" << endl;
    }
    else
    {
        x = xs;
        fx = ys;
        h = vector<float>(xs.size());
        S = vector<float>(xs.size(),0.0);

        //Remember the end conditions dictate the values of S[0] and S[n] so the trida system
        //to solve is just for S[1] - S[n-1] and so the vectors are of size n-2.
        //EXCEPT for condition 1 where S[1] and S[n-1] are given as functiosn of estimates of the derivatives
        //and so in the case need to be solved for to.
        if(condition == 1)
        {            
            f = vector<float>(xs.size());        
            lowerdiag = vector<float>(xs.size());
            middlediag = vector<float>(xs.size());
            upperdiag = vector<float>(xs.size());

            for(unsigned long int i = 0; i< h.size(); i++)
            {
               h[i] = xs[i+1] - xs[i];            
            }
            for(unsigned long int j = 1; j< f.size()-1; j++)
            {
               f[j] = 6.0*( ((ys[j+1]-ys[j])/h[j]) - ((ys[j]-ys[j-1])/h[j-1]) );
            }
        }
        else
        {
            f = vector<float>(xs.size()-2);        
            lowerdiag = vector<float>(xs.size()-2);
            middlediag = vector<float>(xs.size()-2);
            upperdiag = vector<float>(xs.size()-2);

            for(unsigned long int i = 0; i< h.size(); i++)
            {
               h[i] = xs[i+1] - xs[i];            
            }
            for(unsigned long int j = 0; j< f.size(); j++)
            {
               f[j] = 6.0*( ((ys[j+2]-ys[j+1])/h[j+1]) - ((ys[j+1]-ys[j])/h[j]) );
            }
        }
        
        lowerdiag[0] = 0.0;
        upperdiag[f.size()-1] = 0.0;        
        switch(condition)
        {
        case 3://S0 and Sn are linear extrapolations
        {   
            middlediag[0] = (h[0]+h[1])*(h[0]+2.0*h[1]);            
            upperdiag[0] = (h[1]*h[1] - h[0]*h[0])/h[1];
            for(unsigned long int i = 1; i<f.size()-1; i++)
            {                
                upperdiag[i] = h[i+1];
                middlediag[i] = 2.0*(h[i]+h[i+1]);
                lowerdiag[i] = h[i];
            }
            unsigned long int k = f.size()-1;
            middlediag[k] = (h[k+1]+h[k])*(h[k+1]+2.0*h[k])/h[k];
            lowerdiag[k] = (h[k]*h[k]-h[k+1]*h[k+1])/h[k];
            vector<float> res = tridag(lowerdiag,middlediag,upperdiag,f);

            for(unsigned long int j = 0; j<res.size(); j++)
            {
                S[j+1] = res[j];                
            }
            S[0] = ((h[0]+h[1])*S[1] - h[0]*S[2])/h[1];
            unsigned long int n = xs.size()-1;
            S[n] = ((h[n-2]+h[n-1])*S[n-1] - h[n-1]*S[n-2])/h[n-2];
        }
        break;    
        case 2://S[0] = S[1], S[n] = S[n-1]
        {
            middlediag[0] = 3.0*h[0]+2.0*h[1];
            upperdiag[0] = h[1];
            for(unsigned int i = 1; i< f.size()-1; i++)
            {
                lowerdiag[i] = h[i];
                middlediag[i] = 2.0*(h[i]+h[i+1]);
                upperdiag[i] = h[i+1];
            }
            unsigned long int k = f.size()-1;
            middlediag[k] = 2.0*h[k]+3.0*h[k+1];
            lowerdiag[k] = h[k];
            vector<float> res = tridag(lowerdiag,middlediag,upperdiag,f);

            for(unsigned long int j = 0; j<res.size(); j++)
            {
                S[j+1] = res[j];                
            }
            S[0] = S[1];
            S[xs.size()-1] = S[xs.size()-2];
        }
        break;
        case 1://f'(x0) = A, f'(xn) = B
        {
            //estimate A and B numercailly by forawrd and backward differences respectively
            float A = (ys[1]-ys[0])/(xs[1]-xs[0]);
            float B = (ys[ys.size()-1]-ys[ys.size()-2])/(xs[xs.size()-1]-xs[xs.size()-2]);            

            middlediag[0] = 2.0*h[0];
            upperdiag[0] = h[1];
            f[0] = 6.0*( (ys[1]-ys[0])/(xs[1]-xs[0]) - A );
            f[f.size()-1] = 6.0*(B - (ys[ys.size()-1]-ys[ys.size()-2])/(x[xs.size()-1]-xs[xs.size()-2]) );                        

            for(unsigned int i = 1; i< f.size()-1;i++)
            {
                lowerdiag[i] = h[i-1];
                middlediag[i] = 2.0*(h[i-1]+h[i]);
                upperdiag[i] = h[i];
            }
            unsigned long int k = f.size()-1;
            lowerdiag[k] = h[k-1];
            middlediag[k] = 2.0*h[k];

            S = tridag(lowerdiag,middlediag,upperdiag,f);            
        }
        break;

        case 0://default - natural spline
        default:
        {            
            for(unsigned long int i = 0; i< f.size(); i++)
            {
                middlediag[i] = 2.0*(h[i]+h[i+1]);
                if(i!=f.size()-1)
                {
                    upperdiag[i] = h[i+1];
                }
                if(i!=0)
                {
                    lowerdiag[i] = h[i];
                }
            }
            vector<float> res = tridag(lowerdiag,middlediag,upperdiag,f);            

            for(unsigned long int j = 0; j<res.size(); j++)
            {
                S[j+1] = res[j];                
            }

        }
        break;
        }
    }
}

float CubicSplineInterpolater::Interpolate(float xval)
{    
    unsigned long int interval = 0;
    for(unsigned long int i = 0; i< x.size()-1; i++)
    {
        if(xval >= x[i] && xval <= x[i+1])
        {
            interval = i;
            break;
        }
    }

    float a = (S[interval+1]-S[interval])/(6.0*h[interval]);
    float b = S[interval]/2.0;
    float c = (fx[interval+1]-fx[interval])/h[interval] - (2.0*h[interval]*S[interval] + h[interval]*S[interval+1])/6.0;
    float d = fx[interval];

    float retval = a*(xval-x[interval])*(xval-x[interval])*(xval-x[interval]) + b*(xval-x[interval])*(xval-x[interval]) + c*(xval-x[interval]) + d;
            
    return retval;
}
