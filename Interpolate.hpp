#ifndef INTERPOLATE_HPP
#define INTERPOLATE_HPP

#include <vector>

using namespace std;

float LagrangeInterpolate(float,const vector<float> &, const vector<float> &);

class DividedDiffInterpolater
{
public:
    DividedDiffInterpolater(const vector<float> &, const vector<float> &);
    float interpolate(const float &);
private:
    vector<float> a;
    vector<float> x;
};

class CubicSplineInterpolater
{
public:
    CubicSplineInterpolater(const vector<float> &, const vector<float> &,unsigned int condition = 0);
    float Interpolate(float);
private:
    vector<float> h,f,S,x,fx;
    vector <float> lowerdiag,middlediag,upperdiag;    
};

#endif