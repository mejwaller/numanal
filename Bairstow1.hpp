#ifndef BAIRSTOW1_HPP
#define BAIRSTOW1_HPP

#include <vector>

using namespace std;

class Bairstow1{
public:
    Bairstow1():k(0),flag(false){};
    vector< vector<float> > factorPoly(const vector<float> &, float rest = 0.0, float sest = 0.0,float tol = 0.001,int iterlimit = 20);
private:
    void computeBandCarrays();
    void findRandS();
    void checkTolMet();
    void reducePolynomial();
    void changeRandSandK();

    vector<float> a;//storage for original array of coefficients of polynomial
    vector<float> b;
    vector<float> c;
    float tol,r,s,denom,delR,delS;
    int n;//degree of polynomial    
    int k;//??
    int iterlim;
    bool flag;
    vector< vector<float> > resvec;
};

#endif
