#ifndef BEZIER1_HPP
#define BEZIER1_HPP

#include "MathUtils.hpp"
#include <vector>

using namespace std;

class Bezier1
{
public:
    Bezier1(const vector<float> &thexs, const vector<float> &theys, unsigned int theorder = 0,float thestep = 0.01);
    vector<float> getCurveXSeries() const { return curveXSeries; }
    vector<float> getCurveYSeries() const { return curveYSeries; }
private:
    void generateCurveSeries();
    vector<float> curveXSeries,curveYSeries,x,y;
    unsigned int order;
    float step;
    bool morepoints;
};


#endif //BEZIER1_HPP
