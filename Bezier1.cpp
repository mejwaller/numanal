#include "Bezier1.hpp"
#include <iostream>
#include <math.h>

Bezier1::Bezier1(const vector<float> &thexs, const vector<float> &theys,unsigned int theorder,float thestep)
    :x(thexs),y(theys),order(theorder),step(thestep),morepoints(false)
{    

    if(order == 0)
    {
        order = (unsigned long)thexs.size()-1;
    }
    
    if( (thexs.size()-1)%order == 0)
    {
        generateCurveSeries();
    }
    else
    {
        cout << "Bezier 1 error - number of points (" << thexs.size() << ") and order of bezier curve (" << order << ") are incompatabale." << endl;
    }    
}

void Bezier1::generateCurveSeries()
{
    vector<float> xseries;
    vector<float> yseries;

    unsigned long int loops = (x.size()-1)/order;

    if(x.size()!=y.size())
    {
        cout << "Bezier1 error - x and y series not the same!" << endl;
        return;
    }    
    float nbang = factorial(order);
    for(unsigned long int i = 0; i<loops;i++)
    {
        for(double u = 0.0; u < 1.0; u+=step)
        {
            float xval = 0.0,yval = 0.0;
            for(unsigned long int j = 0; j <= order; j++)
            {
                float jbang = factorial(j);
                float n_jbang = factorial(order-j);
                float comb = nbang/(jbang*n_jbang);
                xval+=comb *pow(1.0-(double)u,(double)order-(double)j)*pow((double)u,(double)j)*x[(i*order)+j];
                yval+=comb *pow(1.0-(double)u,(double)order-(double)j)*pow((double)u,(double)j)*y[(i*order)+j];
            }
            curveXSeries.push_back(xval);
            curveYSeries.push_back(yval);
        }
    }
}






