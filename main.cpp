#include "Bisection1.hpp"
#include "Secant1.hpp"
#include "FalsePos1.hpp"
#include "Newton1.hpp"
#include "Muller1.hpp"
#include "Horner1.hpp"
#include "Bairstow1.hpp"
#include "Matrix1.hpp"
#include "Interpolate.hpp"
#include "ChartWidget1.hpp"
#include "Bezier1.hpp"
#include <iostream>
#include <math.h>
#include <vector>
#include <qapplication.h>

using namespace std;
/*
 * Bisection single var non-linear eqn. solver
 */
class my_Bisection1:public Bisection1
{
public:
    my_Bisection1(double tol = 1e-07,unsigned int it = 100):Bisection1(tol,it){}
protected:
    virtual double f(double);
};

double my_Bisection1::f(double x)
{    
    //double retval = x*x*x + x*x - 3.0*x - 3.0;
    double retval = 3.0*x + sin(x) - exp(x);
    return retval;
    return retval;
}

/*
 * Secant single var non-linear eqn. solver
 */
class my_Secant1:public Secant1
{
public:
    my_Secant1(double tol = 1e-08,unsigned int it = 100):Secant1(tol,it){}
protected:
    virtual double f(double);
};

double my_Secant1::f(double x)
{    
    //double retval = x*x*x + x*x - 3.0*x - 3.0;
    //double retval = 3.0*x + sin(x) - exp(x);
    //double retval = (x-1.0)*(exp(x-1)-1);
    //return retval;
    return (x+1.0)*(x+1.0)*(x+1.0);
}

/*
 * False Position single var non-linear eqn. solver
 */
class my_FalsePos1:public FalsePos1
{
public:
    my_FalsePos1(double tol = 1e-08,unsigned int it = 100):FalsePos1(tol,it){}
protected:
    virtual double f(double);
};

double my_FalsePos1::f(double x)
{    
    //double retval = x*x*x + x*x - 3.0*x - 3.0;
    double retval = 3.0*x + sin(x) - exp(x);
    return retval;
}

/*
 * Newton single var non-linear eqn. solver
 */
class my_Newton1:public Newton1
{
    public:
    my_Newton1(double tol = 1e-08,unsigned int it = 50):Newton1(tol,it){}
protected:
    virtual double f(double);
    virtual double df(double);
};

double my_Newton1::f(double x)
{    
    //double retval = x*x*x + x*x - 3.0*x - 3.0;
    //double retval = 3.0*x + sin(x) - exp(x);    
    //return retval;
    //return x*x*x + 2.0*x*x - x + 5.0;
    //return (x-1.0)*(exp(x-1)-1);
    return (x+1.0)*(x+1.0)*(x+1.0);
}

double my_Newton1::df(double x)
{
    //return 3.0 + cos(x) - exp(x);
    //return 3.0*x*x + 4.0*x - 1.0;    
    return 3.0*x*x + 6.0*x + 3.0;
}

/*
 * Complex Newton single var non-linear eqn. solver
 */
class my_Complex_Newton1:public Complex_Newton1
{
    public:
    my_Complex_Newton1(double tol = 1e-08,unsigned int it = 100):Complex_Newton1(tol,it){}
protected:
    virtual complex<double> f(complex<double>);
    virtual complex<double> df(complex<double>);
};

complex<double> my_Complex_Newton1::f(complex<double> x)
{    
    
    return x*x*x + 2.0*x*x - x + 5.0;
}

complex<double> my_Complex_Newton1::df(complex<double> x)
{
    return 3.0*x*x + 4.0*x - 1.0;
}

/*
 * Muller single var non-linear eqn. solver
 */
class my_Muller1:public Muller1
{
public:
    my_Muller1(double tol = 1e-07,unsigned int it = 50):Muller1(tol,it){}
protected:
    virtual double f(double);
};

double my_Muller1::f(double x)
{
    //return 3*x+sin(x)-exp(x);
    //return (x-1.0)*(exp(x-1)-1);
    //return (x+1.0)*(x+1.0)*(x+1.0);
    double retval = 1.0+x+x*x-x*x*x;
    retval/=((1.0-x)*(1.0-x)*(1.0-x));
    retval-=0.892;
    return retval;
}

/*
 * Muller single var complex non-linear eqn. solver
 */
class my_Complex_Muller1:public Complex_Muller1
{
public:
    my_Complex_Muller1(double tol = 1e-06,unsigned int it = 100):Complex_Muller1(tol,it){}
protected:
    virtual complex<double> f(complex<double>);
};

complex<double> my_Complex_Muller1::f(complex<double> x)
{
    return x*x*x + 2.0*x*x - x + 5.0;
}

int main( int argc, char** argv )
{

    QApplication app(argc,argv);    
    vector<float> somexs(7);
    vector<float> someys(7);

    somexs[0] = 1.0;
    somexs[1] = 2.0;
    somexs[2] = 4.0;
    somexs[3] = 5.5;        
    somexs[4] = 7.5;    
    somexs[5] = 8.5;    
    somexs[6] = 9.5;    

    someys[0] = 1.5;
    someys[1] = 9.0;
    someys[2] = 9.5;
    someys[3] = 5.5;
    someys[4] = 1.0;
    someys[5] = 2.0;
    someys[6] = 7.5;
    
    Bezier1 myBezier(somexs,someys,3);
    vector<float> curvex = myBezier.getCurveXSeries();
    vector<float> curvey = myBezier.getCurveYSeries();

    ChartWidget1 *chart = new ChartWidget1();    
    chart->setNumXTicks(10);
    chart->setNumYTicks(10);
    chart->setMaxX(10.0);    
    chart->setMaxY(10.0);
    chart->addPoints(somexs,someys,"control points");    
    chart->addLine(curvex,curvey,"bezier curve");    
    chart->show();
    app.setMainWidget(chart);
    

    //Bezier1 myBezier(


    
    /*vector<float> x(9);
    vector<float> y(9);
    x[0] = 0.0;y[0] = 0.302;
    x[1] = 0.2;y[1] = 0.185;
    x[2] = 0.3;y[2] = 0.106;
    x[3] = 0.4;y[3] = 0.093;
    x[4] = 0.5;y[4] = 0.240;
    x[5] = 0.6;y[5] = 0.579;
    x[6] = 0.7;y[6] = 0.561;
    x[7] = 0.8;y[7] = 0.468;
    x[8] = 1.0;y[8] = 0.302;



    CubicSplineInterpolater csi(x,y,1);

    vector<float> interpolated;
    vector<float> xs;
    
    for(int n= 0;n<21;n++)
    {
        xs.push_back(n*0.05);
    }
    

    for(unsigned int i = 0; i<xs.size(); i++)
    {
        interpolated.push_back(csi.Interpolate(xs[i]));
        cout << "Interpolated value for " << xs[i] << " is " << csi.Interpolate(xs[i]) << endl;
    }
        
    
    ChartWidget1 *chart = new ChartWidget1();    
    chart->setNumXTicks(10);
    chart->setNumYTicks(10);
    chart->setMaxX(1.0);    
    chart->setMaxY(0.6);
    chart->addLine(x,y,"Test1");
    chart->addPoints(xs,interpolated,"TestPoints");
    app.setMainWidget(chart);
    chart->show();
    */

    return app.exec();


    /*vector<float> xs(5);
    vector<float> fxs(5);
    xs[0] = 1.0; fxs[0] = 14.2;
    xs[1] = 2.7; fxs[1] = 17.8;
    xs[2] = 3.2; fxs[2] = 22.0;
    xs[3] = 4.8; fxs[3] = 38.3;
    xs[4] = 5.6; fxs[4] = 51.7;

    float interpolated = LagrangeInterpolate(4.0,xs,fxs);

    cout << "The Lagrange interpolated value for x = 4.0 is " << interpolated << endl;
    */

    /*xs[0] = 3.2; fxs[0] = 22.0;
    xs[1] = 2.7; fxs[1] = 17.8;
    xs[2] = 1.0; fxs[2] = 14.2;
    xs[3] = 4.8; fxs[3] = 38.3;
    xs[4] = 5.6; fxs[4] = 51.7;
    
    DividedDiffInterpolater dd(xs,fxs);    
    
    
    cout << "Interpolating using divided differences for x = 4.0 gives " << dd.interpolate(4.0) << endl;

  */

    /*my_Bisection1 bisection;
    bisection.solve(0.0,1.0);
    */
    /*
    my_Secant1 secant;
    secant.solve(-1.5,0.0);
    */
    /*
    my_FalsePos1 falsepos;
    falsepos.solve(0.0,1.0);
    */
    /*
    my_Newton1 newton;
    newton.solve(-1.5);
    */
    /*
    my_Complex_Newton1 cnewton;
    complex<double> c(1.0,-1.0);
    cnewton.solve(c);
    */
    
    /*my_Muller1 muller;
    muller.solve(-0.5,0.0,-1.5);
    //muller.solve(0.5,1.5,0.0);  
    */
    /*
    my_Complex_Muller1 cmuller;
    complex<double> a(0.5,-0.5);
    complex<double> b(1.0,-1.0);
    complex<double> c(0.0,0.0);
    cmuller.solve(a,b,c);
    */
    /*Horner1 horner;
    vector<double> poly(4);
    poly[0] = 1.0;
    poly[1] = 1.0;
    poly[2] = -3.0;
    poly[3] = -3.0;    
    vector<double> Q = horner.dividePoly(poly,2.0);
    */

    /*vector<float> poly(5);
    poly[0] = 1.0;
    poly[1] = -5.7;
    poly[2] = 26.7;
    poly[3] = -42.21;
    poly[4] = 69.23;        
    */
    
    /*vector<float> poly(6);
    poly[0] = 1;
    poly[1] = -17.8;
    poly[2] = 99.41;
    poly[3] = -261.218;
    poly[4] = 352.611;        
    poly[5] = -134.105;    
    */

    /*vector<float> poly(5);
    poly[0] = 1.0;
    poly[1] = -1.0;
    poly[2] = 1.0;
    poly[3] = -1.0;
    poly[4] = 1.0;        
    */

    /*vector<float> poly(4);
    poly[0] = -1.0;
    poly[1] = -1.0;
    poly[2] = 3.0;
    poly[3] = 1.0;

    Bairstow1 bairstow;
    vector< vector<float> > factors(0);
    factors = bairstow.factorPoly(poly,0.0,0.0,1e-03,20);
    if(!factors.empty())
    {        
        cout << "Factors are:" << endl;
        for(int i =0; i < factors.size(); i++)
        {     
            cout << "(";
            for(int j = 0;j < factors[i].size(); j++)
            {
                if(!FCMP(factors[i][j],0.0))
                {
                    cout << factors[i][j] << "x^" << factors[i].size()-j-1 << " - ";
                }
            }
            cout << ")" << endl;
        }
    }*/

    //FloatMatrix mat1(4,4);
    //cout << "Mat1 after construction:" << endl << mat1 << endl;
    //mat1 = identity(mat1);
    //cout << "Mat1:" << endl << mat1 << endl;
    //FloatMatrix mat2 = mat1;    
    //cout << "mat2:" << endl << mat2 << endl;
    //FloatMatrix mat3 = mat1;
    //cout << "mat3:" << endl << mat3 << endl;    
    //FloatMatrix mat4 = mat3*2.0;
    //cout << "mat4 (=mat3*2.0)" << endl << mat4 << endl;
    //mat4*=2.0;
    //cout << "mat4 after multplied by 2" << endl;
    //cout << mat4 << endl;
    //mat4+=mat1;
    //cout << "Mat4 after adding mat1 to it:" <<endl << mat4 << endl;
    //mat4-=mat1;
    //cout <<"Mat4 after subtracting mat from it:" << endl << mat4 << endl;
    
    //FloatMatrix mat5 = mat4 + mat1;
    
    //cout << "mat5 = mat4 + mat1:" << endl;
    //cout << mat5 << endl;
    
    //cout << "mat4 after used in addition:" << endl;
    //cout << mat4 << endl;   


    //mat5+=mat5;
    //cout << mat5 << endl;
    //mat5-=mat4;
    //cout << "mat5 - mat4" << endl << mat5 << endl;
    //FloatMatrix mat6 = mat5-mat1;
    //cout << "mat5 - mat1" << endl << mat6 << endl;    
    //mat6*=mat5;
    //cout << "mat6*mat5" << endl << mat6 << endl;
    //FloatMatrix mat7 = mat6*mat5;
    //cout << "mat6*mat5 (30*6 = 180)" << endl << mat6 << endl << mat5 << endl << mat7 << endl;

    /*if(!(mat1 == mat7))
    {
        //cout << "Mat1 correctly != mat7" << endl;
    }
    else
    {
        //cout << "Mat1 incorectly thought to be eqaul to mat7!" << endl;
    }
    if(mat1 == mat2)
    {
        //cout << "Mat1 correctly == mat2" << endl;
    }
    else
    {
        //cout << "Mat1 incorrectly thought not to be eqaul to mat2!" << endl;
    }

    FloatMatrix mat8(3,3);
    for(unsigned int i =0;i<3;i++)
    {
        for(unsigned int j = 0; j<3; j++)
        {
            //mat8.contents()[mat8.cols()*i+j] = mat8.cols()*i+j;
            mat8(i,j) = mat8.cols()*i+j;
        }
    }
    //cout << "Mat8 is " << endl << mat8 << endl;
    FloatMatrix mat9 = transpose(mat8);
    //cout << "Transpose is " << endl << mat9 << endl;

    //cout << "trace of mat8 is " << trace(mat8) << " and trace of mat9 is " << trace(mat9) << endl;

    vector<float> construct(9);

    for(unsigned j = 0; j < 9; j++)
    {
        construct[j] = j*j;
    }

    FloatMatrix mat10(construct,3,3);
    FloatMatrix mat11(construct,5,2);
    
    // cout << "mat10 from vector constructor is " << endl << mat10 << endl << " and mat11 is " << endl << mat11 << endl;

    vector<float> coeffs1(16);

    coeffs1[0] = 3.0;
    coeffs1[1] = 2.0;
    coeffs1[2] = -1.0;
    coeffs1[3] = 2.0;
    coeffs1[4] = 1.0;
    coeffs1[5] = 4.0;
    coeffs1[6] = 0.0;
    coeffs1[7] = 2.0;
    coeffs1[8] = 2.0;
    coeffs1[9] = 1.0;
    coeffs1[10] = 2.0;
    coeffs1[11] = -1.0;
    coeffs1[12] = 1.0;
    coeffs1[13] = 1.0;
    coeffs1[14] = -1.0;
    coeffs1[15] = 3.0;

    vector<float> coeffs2(16);

    coeffs2[0] = 0.0;
    coeffs2[1] = 2.0;
    coeffs2[2] = 0.0;
    coeffs2[3] = 1.0;
    coeffs2[4] = 2.0;
    coeffs2[5] = 2.0;
    coeffs2[6] = 3.0;
    coeffs2[7] = 2.0;
    coeffs2[8] = 4.0;
    coeffs2[9] = -3.0;
    coeffs2[10] = 0.0;
    coeffs2[11] = 1.0;
    coeffs2[12] = 6.0;
    coeffs2[13] = 1.0;
    coeffs2[14] = -6.0;
    coeffs2[15] = -5.0;

    vector<float> coeffs3(9);
    coeffs3[0] = 4.0;
    coeffs3[1] = -2.0;
    coeffs3[2] = 1.0;
    coeffs3[3] = -3.0;
    coeffs3[4] = -1.0;
    coeffs3[5] = 4.0;
    coeffs3[6] = 1.0;
    coeffs3[7] = -1.0;
    coeffs3[8] = 3.0;
    

    FloatMatrix A(coeffs1,4,4);

    vector<float> bs1(12);    
    bs1[0] = 0.0;
    bs1[1] = -2.0;
    bs1[2] = 2.0;
    bs1[3] = 0.0;
    bs1[4] = 1.0;
    bs1[5] = 2.0;
    bs1[6] = 1.0;
    bs1[7] = 3.0;
    bs1[8] = 0.0;
    bs1[9] = 0.0;
    bs1[10] = 4.0;
    bs1[11] = 0.0;
    

    vector<float> bs2(4);
    bs2[0] = 0.0;
    bs2[1] = -2.0;
    bs2[2] = -7.0;
    bs2[3] = 6.0;
    
    vector<float> bs3(3);
    bs3[0] = 15.0;
    bs3[1] = 8.0;
    bs3[2] = 13.0;    
    

    //FloatMatrix b(bs1,3,1);
    FloatMatrix b(bs1,4,3);
    //FloatMatrix b(bs2,4,1);

    //cout << "A is " << endl << A << endl;
    //cout << "b is " << endl << b << endl;


    FloatMatrix res = GaussElimLinearSolve(A,b);

    //cout << "The results are:" << endl << res << endl;

    FloatMatrix mult = A*res;

    //cout << "A*res = " << endl << mult << endl;

    FloatMatrix x(4,3);
    //FloatMatrix x(3,1);

    GaussElimLinearSolve(A,x,b);

    //cout << "The results2 are:" << endl << x << endl;

    coeffs3[0] = 0.0;
    coeffs3[1] = 2.0;
    coeffs3[2] = 1.0;
    coeffs3[3] = 1.0;
    coeffs3[4] = 0.0;
    coeffs3[5] = 0.0;
    coeffs3[6] = 3.0;
    coeffs3[7] = 0.0;
    coeffs3[8] = 1.0;
    */
    
/*    coeffs3[0] = 3.0;
    coeffs3[1] = -1.0;
    coeffs3[2] = 2.0;
    coeffs3[3] = 1.0;
    coeffs3[4] = 2.0;
    coeffs3[5] = 3.0;
    coeffs3[6] = 2.0;
    coeffs3[7] = -2.0;
    coeffs3[8] = -1.0;    
    */

/*
    FloatMatrix LU(coeffs3,3,3);
    FloatMatrix Orig(LU);
    FloatMatrix Orig2(LU);

    cout << "Before decomposition, LU is: " << endl << LU << endl;
    
    LU.LUDecompose();
    //FloatMatrix decomposed(LU);
    cout << "LU LUDecomposed() is:" << endl << LU << endl;

    vector<float> bvec(3);
    bvec[0] = 5.0;
    bvec[1] = -1.0;
    bvec[2] = -2.0;

    FloatMatrix GausstestB(bvec,3,1);    

    //cout << " and b is " << endl << GausstestB << endl;

    LU.LUBackSubstitute(bvec);

    FloatMatrix BSolved(bvec,3,1);

    cout << "The solution vector is:" << endl << BSolved << endl;
    */

    /*cout << "LU result * LU is " << endl << Orig*BSolved << endl;
    cout << "LU decomposed * result is " << endl << decomposed*BSolved << endl;

    bvec[0] = 5.0;
    bvec[1] = -1.0;
    bvec[2] = -2.0;
    */

    //FloatMatrix GausstestA(coeffs3,3,3);    

    //FloatMatrix res2 = GaussElimLinearSolve(GausstestA,GausstestB);
    //cout << res2 << endl;   


    //Now test if ludecomposed results can recreate A
    //(it does)
    /*
    vector<float> alpha(9);
    alpha[0] = 1.0;alpha[1] = 0.0;alpha[2] = 0.0;
    alpha[3] = 0.0;alpha[4] = 1.0;alpha[5] = 0.0;
    alpha[6] = 1.0/3.0; alpha[7] = 0.0; alpha[8] = 1.0;

    vector<float> beta(9);
    beta[0] = 3.0;beta[1] = 0.0;beta[2] = 1.0;
    beta[3] = 0.0;beta[4] = 2.0;beta[5] = 1.0;
    beta[6] = 0.0; beta[7] = 0.0;beta[8] = -1.0/3.0;

    FloatMatrix Alpha(alpha,3,3);
    FloatMatrix Beta(beta,3,3);

    cout << "LUdcomposed results mutli;lied = " << endl << Alpha*Beta << endl;
    */

/*
    
    FloatMatrix inv = Orig.inverse();

    cout << "Matrix inverse is " << endl << inv << endl;

    cout << "Multplied by inverse we get idemtity:" << endl << inv*Orig2 << endl;

    //cout << "Determinant  is " << Orig2.det() << endl;

    cout << "maxcolsum norm is " << Orig2.maxcolsum_Norm() << endl;
    cout << "maxrowsum norm is " << Orig2.maxrowsum_Norm() << endl;
    cout << "Euclidean norm is " << Orig2.Euclidean_Norm() << endl;
    */

    /*Tests to check mult. works - results on p.107 of Gerald and Wheatley    
    vector<float> vA(6);
    vA[0] = 3.0;
    vA[1] = 7.0;
    vA[2] = 1.0;
    vA[3] = -2.0;
    vA[4] = 1.0;
    vA[5] = -3.0;

    vector<float> vB(6);
    vB[0] = 5.0;
    vB[1] = -2.0;
    vB[2] = 0.0;
    vB[3] = 3.0;
    vB[4] = 1.0;
    vB[5] = -1.0;

    vector<float> vx(3);
    vx[0] = -3.0;
    vx[1] = 1.0;
    vx[2] = 4.0;

    FloatMatrix mA(vA,2,3);
    FloatMatrix mB(vB,3,2);
    FloatMatrix mx(vx,3,1);

    cout << "mA:\n" << mA << endl;
    cout << "mB:\n" << mB << endl;
    cout << "mx:\n" << mx << endl << endl;

    cout << "mA*mB:\n" << mA*mB << endl;
    cout << "mA:\n" << mA << endl;
    cout << "mB:\n" << mB << endl;
    cout << "mx:\n" << mx << endl << endl;

    cout << "mB*mA:\n" << mB*mA << endl;
    cout << "mA:\n" << mA << endl;
    cout << "mB:\n" << mB << endl;
    cout << "mx:\n" << mx << endl << endl;

    cout << "mA * mx:\n" << mA*mx << endl;
    cout << "mA:\n" << mA << endl;
    cout << "mB:\n" << mB << endl;
    cout << "mx:\n" << mx << endl << endl;
    */

    
    return 0;
}




