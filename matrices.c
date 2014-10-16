#include <math.h>
#include <stdio.h>
#include  "nr.h"
#include "matrices.h"
#include <float.h>

#ifdef FCMP
#   error "FCMP already defined"
#else
#   define FCMP(aa, bb) ((fabs((double)(aa) - (double)(bb))/\
           (FLT_EPSILON+fabs((double)(aa))+fabs((double)(bb))) < FLT_EPSILON))
#endif 

void ludcmp(float **a, int n, int *indx, float *d)
{
    int i,imax,j,k;
    float big,dum,sum,temp;
    float *vv;
    /*float ival;*/

    vv = vector(1,n);
    *d = 1.0;
    for(i = 1; i<=n;i++)
    {
        big = 0.0;
        for(j=1;j<=n;j++)
        {            
            if((temp = fabs(a[i][j])) > big)
            {                      
                big = temp;
            }
        }
        if(FCMP(big,0.0))
        {            
            return;
        }
                
        vv[i] = 1.0/big;
    }
    
    for(j=1; j<=n;j++)
    {
        for(i = 1; i<j;i++)
        {    
            sum = a[i][j];
            for(k=1;k<i;k++)
            {         
                sum -= a[i][k]*a[k][j];
            }
            a[i][j] = sum;            
        }
        big = 0.0;
        for(i=j;i<=n;i++)
        {            
            sum = a[i][j];
            for(k=1;k<j;k++)
            {
                sum -= a[i][k]*a[k][j];             
            }
            a[i][j] = sum;            
            if( (dum=vv[i]*fabs(sum)) >= big)
            {                
                big = dum;
                imax =i;             
            }
        }        
        if(j!=imax)
        {         
            for(k=1;k<=n;k++)
            {         
                dum = a[imax][k];             
                a[imax][k] = a[j][k];                 
                a[j][k] = dum;
            }            
            *d= -(*d);            
            vv[imax] = vv[j];
        }        
        indx[j] = imax;

        if(a[j][j] == 0.0)
        {
            a[j][j]=1.0e-20;
        }  

        if(j!=n)
        {        
            dum = 1.0/(a[j][j]);
            for(i=j+1;i<=n;i++)
            {         
                a[i][j]*=dum;
            }
        }        
    }
    /*free_vector(vv,1,n);*/
}

void lubksb(float **a, int n, int *indx,float b[])
{
    int i,ii=0,ip,j;
    float sum;

    for(i=1;i<=n;i++)
    {
        ip = indx[i];        
        sum = b[ip];        
        b[ip] = b[i];        
        if(ii)
        {
            /*printf("ii not zero (=%d, with i = %d)...\n",ii,i);*/
            for(j=ii;j<=i-1;j++)
            {
                /*printf("...(ii=true), subtracting a[%d][%d]*b[%d] (=%f/%f = %f) from sum (currently %f)...\n",i,j,j,a[i][j],b[j],a[i][j]/b[j],sum);*/
                sum-=a[i][j]*b[j];
                /*printf("..and now sum is %f\n",sum);*/
            }
        }
        else if (sum)
        {
            /*printf("sum not zero (=%f), with i = %d\n",sum,i);*/
            ii=i;
        }
        b[i]=sum;        
    }
    for(i=n;i>=1;i--)
    {        
        sum = b[i];
     
        for(j=i+1;j<=n;j++)
        {            
            sum-=a[i][j]*b[j];         
        }        
        b[i]=sum/a[i][i];
    }
}



            
