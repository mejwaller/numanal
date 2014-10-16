#include "matrices.h"
#include "nr.h"
#include <stdio.h>

int main( int argc, char** argv )
{
    float **A;
    int i,j;
    int *rowperms;
    float sign;
    float *b;
    float **y;

    A = matrix(1,3,1,3);

    A[1][1] = 0.0; A[1][2] = 2.0; A[1][3] = 1.0;
    A[2][1] = 1.0; A[2][2] = 0.0; A[2][3] = 0.0;
    A[3][1] = 3.0; A[3][2] = 0.0; A[3][3] = 1.0;

    rowperms = ivector(1,3);    

    printf("A is: \n");

    for(i =1; i<=3;i++)
    {
        for(j=1;j<=3;j++)
        {
            printf("%f ",A[i][j]);
        }
        printf("\n");
    }

    ludcmp(A,3,rowperms,&sign);

    printf("LUdecompsed A is:\n");

    for(i =1; i<=3;i++)
    {
        for(j=1;j<=3;j++)
        {
            printf("%f ",A[i][j]);
        }
        printf("\n");
    }

    /*Add test to see if matrix inverse created properly...*/   
    b = vector(1,3);    
    y = matrix(1,3,1,3);


    for(j=1;j<=3;j++)
    {
        for(i=1;i<=3;i++)
        {
            b[i] = 0.0;
        }
        b[j] = 1.0;
        lubksb(A,3,rowperms,b);
        for(i = 1; i<=3;i++)
        {
            y[i][j] = b[i];
        }
    }

    printf("Inverse of A is: \n");

    for(i =1; i<=3;i++)
    {
        for(j=1;j<=3;j++)
        {
            printf("%f ",y[i][j]);
        }
        printf("\n");
    }

    /*float *b;
    b = vector(1,3);
    b[1] = 5.0;
    b[2] = -1.0;
    b[3] = -2.0;

    lubksb(A,3,rowperms,b);

    printf("Solution vector is:\n");

    for(i = 1; i<=3;i++)
    {
        printf("%f\n",b[i]);
    }*/

    return (0);
}
    

