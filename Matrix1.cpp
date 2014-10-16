#include "Matrix1.hpp"
#include "BaseUnaryNonLinearSolver.hpp"
#include <math.h>
#include <stdio.h>

Matrix::~Matrix()
{
}

FloatMatrix::~FloatMatrix()
{ 
}


FloatMatrix::FloatMatrix(unsigned long int rows,unsigned long int cols):Matrix(rows,cols),d(1.0)
{
    matrix = vector<float>(rows*cols,0.0);
    row_perm = vector<unsigned int>(cols,0);
}      

FloatMatrix::FloatMatrix(const FloatMatrix &mat):Matrix(mat.rows(),mat.cols()),row_perm(mat.row_perm),d(mat.d)
{      
    matrix=vector<float>(mat.contents());      
}

FloatMatrix::FloatMatrix(const vector<float> & valarray,unsigned long int rows,unsigned long int cols):Matrix(rows,cols),d(1.0)
{
    //initialise matrix to zeros, then assign valarray to it.
    //This ensures that if valarray is smaller than matrix,
    //additional elements of matrix are 0.
    matrix = vector<float>(rows*cols,0.0);
    matrix.assign(valarray.begin(),valarray.end());
    row_perm = vector<unsigned int>(cols,0);
}

FloatMatrix & FloatMatrix::operator=(const FloatMatrix &mat)
{       
    matrix=vector<float>(mat.contents());            
    row_perm = mat.row_perm;
    d = mat.d;
    return (*this);
}

ostream & operator<<(ostream & out, const FloatMatrix & source)
{
    for(unsigned long int i = 0; i < source.rows(); i++)
    {
        for(unsigned long int j = 0; j < source.cols(); j++)
        {
            //out << source.contents()[source.cols()*i+j] << "   ";
            out << source(i,j) << "   ";
        } 
        out << endl;
    }

    return out;
}

FloatMatrix & FloatMatrix::operator*(float num)
{
    for(unsigned long int i = 0; i < rows(); i++)
    {
        for(unsigned long int j = 0; j < cols(); j++)
        {                
            contents()[cols()*i+j]*=num;                                         
        }            
    }
    
    return *this;
}

FloatMatrix & FloatMatrix::operator*=(float num)
{    

    for(unsigned long int i = 0; i < rows(); i++)
    {
        for(unsigned long int j = 0; j < cols(); j++)
        {                
            contents()[cols()*i+j]*=num;                             
        }            
    }    
    return *this;
}

FloatMatrix operator+(const FloatMatrix &mat1,const FloatMatrix &mat2)
{       
    FloatMatrix retval = mat1;    
    retval += mat2;    

    return retval;
}

FloatMatrix operator-(const FloatMatrix &mat1,const FloatMatrix &mat2)
{
    FloatMatrix retval = mat1;
    retval -= mat2;    
    
    return retval;    
}

FloatMatrix & FloatMatrix::operator+=(const FloatMatrix &mat)
{    
    if(!mat.cols() == cols() && mat.rows() == rows())
    {
       cout << "Can't add different sized matrices!" <<  endl;           
    }
    else
    {
        for(unsigned long int i = 0; i<mat.rows(); i++)
        {
            for(unsigned long int j = 0; j<mat.cols(); j++)
            {
                contents()[cols()*i+j]+=mat(i,j);//mat.contents()[cols()*i+j];            
            }
        }
    }

    return *this;
}

FloatMatrix & FloatMatrix::operator-=(const FloatMatrix &mat)
{    
    if(!mat.cols() == cols() && mat.rows() == rows())
    {
       cout << "Can't add different sized matrices!" <<  endl;           
    }
    else
    {
        for(unsigned long int i = 0; i<mat.rows(); i++)
        {
            for(unsigned long int j = 0; j<mat.cols(); j++)
            {
                contents()[cols()*i+j]-=mat(i,j);//mat.contents()[cols()*i+j];
            }
        }
    }    
    return *this;
}


FloatMatrix & FloatMatrix::operator*=(const FloatMatrix &mult)
{
    if(cols()!=mult.rows())
    {
        cout << "Matrices nonconformable for multiplication!" <<  endl;                   
    }
    else
    {
        FloatMatrix temp(rows(),mult.cols());//to contain the result        

        for(unsigned long int i = 0; i < rows(); i++)
        {
            for(unsigned long int j = 0; j < mult.cols(); j++)
            {
                for(unsigned long int k = 0; k < cols() ; k++)
                {
                    temp(i,j) += (*this)(i,k)*mult(k,j);
                }
            }
        }    
        *this = temp;                
        this->setRows(temp.rows());
        this->setCols(temp.cols());
    }    

    return *this;
}

FloatMatrix operator*(const FloatMatrix &A,const FloatMatrix &B)
{
    FloatMatrix res = A;
    res*=B;
    return res;
}

bool operator==(const FloatMatrix &A,const FloatMatrix &B)
{
    if(A.rows()!= B.rows() || A.cols() != B.cols())
    {
        return false;
    }
    else
    {
        for(unsigned long int i=0; i<A.rows();i++)
        {
            for(unsigned long int j=0;j<A.cols();j++)
            {
                //if(A.contents()[A.cols()*i+j]!=B.contents()[A.cols()*i+j])
                if(A(i,j) != B(i,j))
                {
                    return false;
                }
            }
        }
        return true;
    }
}

float & FloatMatrix::operator()(unsigned long int row,unsigned long int col)
{
    return matrix[cols()*row + col];
}

float FloatMatrix::operator()(unsigned long int row,unsigned long int col) const
{
    return matrix[cols()*row + col];
}

void FloatMatrix::LUDecompose()
{
    //Performs LUDecomposition of this matrix, replacing it with its 
    //LU decomposition
    if(rows()!=cols())
    {
        cout << "Error - square matrix requried fro LUDecompose()." << endl;
        return;
    }
    d = 1.0;            
    vector<float> vv(cols());
    float big,temp,sum,dum;
    unsigned int imax;
    for(unsigned int i = 0; i<rows(); i++)//loop over rows to get implicit scaling information
    {
        big = 0.0;
        for(unsigned int j = 0; j<cols();j++)
        {            
            if((temp = fabs((*this)(i,j))) > big)
            {         
                big = temp;
            }
        }

        if(FCMP(big,0.0))
        {            
            return;
        }    
        
        vv[i] = 1.0/big;//save the scaling        
    }    

    for(unsigned int j = 0; j < cols();j++)//this is the loop over columsn of Crout's method
    {
        for(unsigned int i = 0; i<j; i++)//now do eqn 2.3.12 of Press et al (p44) (excpt for i = j)
        {                 
            sum = (*this)(i,j);
            for(unsigned int k =0; k<i;k++)
            {                
                sum -= (*this)(i,k) * (*this)(k,j);             
            }            
            (*this)(i,j) = sum;            
        }
        big = 0.0;//initialise the search for the largest pivot element
        
        for(unsigned int q= j; q<rows();q++)//this is i = j of eqn. 2.3.12 of Press et al (p.44)
        {                                    //and i = j+1...N of 2.3.13         
            sum = (*this)(q,j);
            for(unsigned int k = 0; k<j;k++)
            {   
                sum -= (*this)(q,k)*(*this)(k,j);             
            }
            (*this)(q,j) = sum;
            
            if( (dum = vv[q]*fabs(sum)) >= big)
            {            
                big = dum;
                imax = q;             
            }            
        }
        
        if(j != imax)//do rows need interchanging?
        {//Yes! Do it.        
            for(unsigned int k = 0; k<rows();k++)
            {                
                dum = (*this)(imax,k);             
                (*this)(imax,k) = (*this)(j,k);                
                (*this)(j,k) = dum;
            }
            
            d=-d;//row interchange - change sign of d
            
            vv[imax] = vv[j];//interchange the scale factor
            
        }        
        
        row_perm[j]= imax;

        if(  FCMP( (*this)(j,j),0.0 ) )
        {
            cout << "Error - pivot element is 0!" << endl << "Matrix is singular and only partly LUDecomposed." << endl;
            return;
        }        
        if(j != rows()-1)//finally, divide by the pivot element
        {            
            dum = 1.0/((*this)(j,j));
            for(unsigned int i = j+1; i<rows(); i++)
            {         
                (*this)(i,j)*=dum;
            }
        }             
    }    
}

void FloatMatrix::LUBackSubstitute( vector<float> &b )
{
    unsigned long int ip;
    long int ii = -1;
    float sum;

    for(unsigned long int i = 0; i<b.size(); i++)
    {
        ip = row_perm[i];        
        sum = b[ip];         
        b[ip] = b[i];  
        
        if(ii>-1)
        {
            //printf("ii not zero (=%d, with i = %d)...\n",ii,i);
            for(unsigned long int j = ii; j<=i-1;j++)
            {
                //printf("...(ii=true), subtracting a[%d][%d]*b[%d] (=%f/%f = %f) from sum (currently %f)...\n",i,j,j,(*this)(i,j),b[j],(*this)(i,j)/b[j],sum);
                sum-=(*this)(i,j)*b[j];
                //printf("..and now sum is %f\n",sum);
            }
        }
        else if(!FCMP(sum,0.0))
        {
            //printf("sum not zero (=%f), with i = %d\n",sum,i);
            ii =i;
        }
        b[i] = sum;        
    }
    for(long int k = b.size()-1;k>=0;k--)
    {        
        sum = b[k];        
        for(unsigned long int j = k+1;j<b.size();j++)
        {            
            sum-=(*this)(k,j)*b[j];
        }        
        b[k] = sum/(*this)(k,k);
    }    
}

FloatMatrix FloatMatrix::inverse()
{
    vector<float> col(cols(),0.0);
    FloatMatrix inv(rows(),cols());
    LUDecompose();
    for(unsigned long int j = 0; j < rows(); j++)
    {
        for(unsigned long int i =0; i < rows(); i++)
        {
            col[i] = 0.0;
        }
        col[j] = 1.0;     
        LUBackSubstitute(col);
        for(unsigned long int k = 0; k< rows(); k++)
        {             
            inv(k,j) = col[k];
        }
    }

    return inv;
}

float FloatMatrix::det()
{    
    LUDecompose();
    float res = d;
    for(unsigned long int i = 0; i<rows();i++)
    {
        res*=(*this)(i,i);
    }
    return res;
}

float FloatMatrix::maxcolsum_Norm()
{
    float res = 0.0;
    float max = res;
    for(unsigned long int i = 0; i< rows();i++)
    {
        res = 0.0;
        for(unsigned long int j = 0; j<cols();j++)
        {
            res+=fabs((*this)(j,i));
        }
        if(res>max)
        {
            max = res;
        }
    }
    return max;
}

float FloatMatrix::maxrowsum_Norm()
{
    float res = 0.0;
    float max = res;
    for(unsigned long int i = 0; i< rows();i++)
    {
        res = 0.0;
        for(unsigned long int j = 0; j<cols();j++)
        {
            res+=fabs((*this)(i,j));
        }
        if(res>max)
        {
            max = res;
        }
    }
    return max;
}

float FloatMatrix::Euclidean_Norm()
{
    float res = 0.0;
    for(unsigned long int i = 0; i< rows();i++)
    {     
        for(unsigned long int j = 0; j<cols();j++)
        {
            res+=((*this)(i,j)*(*this)(i,j));
        }        
    }

    return sqrt(res);
}

FloatMatrix transpose(const FloatMatrix &in)
{
    FloatMatrix temp(in.cols(),in.rows());
    for(unsigned long int i = 0; i < in.cols(); i++)
    { 
        for(unsigned long int j = 0; j < in.rows(); j++)
        {
            //temp.contents()[temp.cols()*i+j] = matrix[cols()*j+i];
            temp(i,j) = in.contents()[in.cols()*j+i];
        }
    }
    return temp;
}

FloatMatrix Identity(const FloatMatrix &in)
{
    FloatMatrix temp = in;
    for(unsigned long int i = 0; i < in.rows(); i++)
    { 
        for(unsigned long int j = 0; j < in.cols(); j++)
        {
            if(i==j)
            {                
                temp(i,j) = 1.0;
            }
            else
            {                
                temp(i,j) = 0.0;
            }
        }
    }
    return temp;
}

float trace(const FloatMatrix &in)
{
    float retval = 0.0;
    for(unsigned long int i = 0; i < in.rows(); i++)
    { 
        for(unsigned long int j = 0; j < in.cols(); j++)
        {
            if(i==j)
            {
                retval+=in(i,j);
            }
        }
    }
    return retval;
}

FloatMatrix GaussElimLinearSolve(const FloatMatrix &ain,const FloatMatrix &bin)
{     
    FloatMatrix x(bin.rows(),bin.cols());
    FloatMatrix a = ain;
    FloatMatrix b = bin;    

    if( (ain.rows() != x.rows()) || (ain.rows() != bin.rows()) || (bin.rows() != x.rows()) )
    {
        cout << "Can't solve - not all matrices in " << ain.rows() << " equations." << endl;
        return x;
    }    

    float pvt;                   
    vector<unsigned long int> pivot(ain.rows());
    unsigned long int ipvt_store;

    for(unsigned long int j = 0; j < ain.rows()-1; j++)//step over all the (columns?) from the first to the last but 1
    {
        pvt = fabs(a(j,j));//initialise pivot point to a(0,0)
        pivot[j] = j;//
        ipvt_store = j+1;

        for(unsigned long int i = j + 1; i < ain.rows(); i++)
        {//now find pivot row (if its not the first one)
         //(the loop statement above has the effect of stepping down the 
         //jth column)
            
            if(fabs(a(i,j)) > pvt)//if value at the ith row of this column is > than the current pivot value....
            {
                pvt = fabs(a(i,j));//set the pivot value to the new value
                ipvt_store = i;//?
            }
        }

        if(pivot[j] != ipvt_store)
        {
            float temp;
            pivot[j] = ipvt_store;
            pivot[ipvt_store] = j;

            for(unsigned long int k = 0; k < a.rows(); k++)
            {
                temp = a(j,k);
                a(j,k) = a(pivot[j],k);
                a(pivot[j],k) = temp;
            }

            //for multiple right hand sides... otherwise l just == 0
            for(unsigned long int l = 0; l < b.cols(); l++)
            {

                temp = b(j,l);
                b(j,l) = b(pivot[j],l);
                b(pivot[j],l) = temp;
            }            
        }

        if(FCMP(a(j,j),0.0))
        {
            cout << "Error - singular matrix encountered, no solution." << endl;
            return x;
        }

        for(unsigned long int l = j+1; l < ain.rows(); l++)
        {//store multipliers
            a(l,j)/=a(j,j);
        }

        for(unsigned long int m = j+1; m < ain.rows(); m++)
        {//create zeros below main diagonals
            for(unsigned long int n = j+1; n < ain.rows(); n++)
            {
                a(m,n)-=(a(m,j)*a(j,n));
            }

            for(unsigned long int o = 0; o < b.cols(); o++)
            {
                b(m,o)-=(a(m,j)*b(j,o));
            }            
        }        
    }
    //now for baclk substition
    for(unsigned long int n = 0; n<x.cols(); n++)
    {
        x(ain.rows()-1,n) = b(ain.rows()-1,n)/a(ain.rows()-1,ain.rows()-1);
    }    

    for(unsigned long int i = ain.rows()-2; i >= 0; i--)
    {        
        if(FCMP(a(i,i),0.0))
        {
            cout << "Error - singular matrix encountered, no solution." << endl;
            cout << "a(" << i << "," << i << ") is " << a(i,i) << ")" << endl;
            return x;
        }        

        for(unsigned  int p = 0; p < b.cols(); p++)
        {
            x(i,p) = b(i,p);

            for(unsigned int k = ain.rows()-1; k >= i+1; k--)
            {           
                x(i,p)-=(x(k,p)*a(i,k));
            }

            x(i,p)/=a(i,i);
        }        

        if(i == 0)//unsigned - dcrenmemnting after 0 just give smax _unsigned_ value, not <0!
        {
            break;
        }
    }

    /*cout << "a is:" << endl << a << endl;
    cout << "b is:" << endl << b << endl;
    cout << "x is:" << endl << x << endl;    
    */

    return x;
}

void GaussElimLinearSolve(FloatMatrix &a, FloatMatrix &x, FloatMatrix &b)
{    

    if( (a.rows() != x.rows()) || (a.rows() != b.rows()) || (b.rows() != x.rows()) )
    {
        cout << "Can't solve - not all matrices in " << a.rows() << " equations." << endl;        
    }    

    float pvt;                   
    vector<unsigned long int> pivot(a.rows());
    unsigned long int ipvt_store;

    for(unsigned long int j = 0; j < a.rows()-1; j++)//step over all the (columns?) from the first to the last but 1
    {
        pvt = fabs(a(j,j));//initialise pivot point to a(0,0)
        pivot[j] = j;
        ipvt_store = j+1;

        for(unsigned long int i = j + 1; i < a.rows(); i++)
        {//now find pivot row (if its not the first one)
         //(the loop statement above has the effect of stepping down the 
         //jth column)
            
            if(fabs(a(i,j)) > pvt)//if value at the ith row of this column is > than the current pivot value....
            {
                pvt = fabs(a(i,j));//set the pivot value to the new value
                ipvt_store = i;//?
            }
        }

        if(pivot[j] != ipvt_store)
        {
            float temp;
            pivot[j] = ipvt_store;
            pivot[ipvt_store] = j;

            for(unsigned long int k = 0; k < a.rows(); k++)
            {
                temp = a(j,k);
                a(j,k) = a(pivot[j],k);
                a(pivot[j],k) = temp;
            }

            //for multiple right hand sides... otherwise l just == 0
            for(unsigned long int l = 0; l < b.cols(); l++)
            {

                temp = b(j,l);
                b(j,l) = b(pivot[j],l);
                b(pivot[j],l) = temp;
            }      
        }

        if(FCMP(a(j,j),0.0))
        {
            cout << "Error - singular matrix encountered, no solution." << endl;            
        }

        for(unsigned long int l = j+1; l < a.rows(); l++)
        {//store multipliers
            a(l,j)/=a(j,j);
        }

        for(unsigned long int m = j+1; m < a.rows(); m++)
        {//create zeros below main diagonals
            for(unsigned long int n = j+1; n < a.rows(); n++)
            {
                a(m,n)-=(a(m,j)*a(j,n));
            }

            for(unsigned long int o = 0; o < b.cols(); o++)
            {
                b(m,o)-=(a(m,j)*b(j,o));
            }            
        }        
    }
    //now for baclk substition
    for(unsigned long int n = 0; n<x.cols(); n++)
    {
        x(a.rows()-1,n) = b(a.rows()-1,n)/a(a.rows()-1,a.rows()-1);
    }       

    for(unsigned long int i = a.rows()-2; i >= 0; i--)
    {        
        if(FCMP(a(i,i),0.0))
        {
            cout << "Error - singular matrix encountered, no solution." << endl;
            cout << "a(" << i << "," << i << ") is " << a(i,i) << ")" << endl;        
        }

        for(unsigned  int p = 0; p < b.cols(); p++)
        {
            x(i,p) = b(i,p);

            for(unsigned int k = a.rows()-1; k >= i+1; k--)
            {           
                x(i,p)-=(x(k,p)*a(i,k));
            }

            x(i,p)/=a(i,i);
        }       

        if(i == 0)//unsigned - dcrenmemnting after 0 just give smax _unsigned_ value, not <0!
        {
            break;
        }
    }

    //cout << "a is:" << endl << a << endl;
    //cout << "b is:" << endl << b << endl;
    //cout << "x is:" << endl << x << endl;        
}

vector<float> tridag(const vector<float> &a,const vector<float> &b, const vector<float> &c, const vector<float> &r)
{
    vector<float> u(b.size(),0.0);
    vector<float> gam(b.size());
    float bet;

    if(FCMP(b[0],0.0))
    {
        cout << "Error - matrix[0][0] == 0.0 - try rewriting eqns. as a set of order n-1 with u2 trivually elimianted..." << endl;
        return u;
    }

    u[0] = r[0]/(bet=b[0]);
    for(unsigned long int j = 1; j<b.size(); j++)
    {
        gam[j] = c[j-1]/bet;
        bet = b[j]-a[j]*gam[j];
        if(FCMP(bet,0.0))
        {
            cout << "Error -  cannot solve tridag system.\n" << endl;
            return u;
        }
        u[j] = (r[j]-a[j]*u[j-1])/bet;
    }
    for(unsigned long int i = b.size()-2; i>=0;i--)
    {
        u[i]-=gam[i+1]*u[i+1];
        if(i==0)
        {
            break;
        }
    }
    return u;
}
































