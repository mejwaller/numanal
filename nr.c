#include "nr.h"
#include <stdlib.h>
#include <stdio.h>

float *vector(long nl, long nh)
{
    float *v;
    v=(float *)malloc( (size_t) ((nh-nl+1)*sizeof(float)));
    if(!v)
    {
        printf("Allocation failure in nr vector()!\n");
        return (float *)0;
    }
    return v-nl+1;
}

int *ivector(long nl, long nh)
{
    int *v;
    v=(int *)malloc((size_t) ((nh-nl+1)*sizeof(int)));
    if(!v)
    {
        printf("Allocation failure in nr ivector()!\n");
        return (int *)0;
    }
    return v-nl+1;
}

float **matrix(long nrl, long nrh, long ncl, long nch)
{
    long i, nrow=nrh-nrl+1,ncol = nch-ncl+1;
    float **m;

    m = (float **)malloc((size_t)((nrow+1)*sizeof(float *)));

    if(!m)
    {
        printf("Allocation failure in nr matrix()!\n");
        return (float **)0;
    }
    m+=1;
    m-=nrl;

    m[nrl]=(float *)malloc((size_t)((nrow*ncol+1)*sizeof(float)));
    if(!m[nrl])
    {
        printf("Allocation failure in nr matrix()!\n");
        return (float **)0;
    }
    m[nrl]+=1;
    m[nrl]-=ncl;
    for(i=nrl+1; i<=nrh;i++)
    {
        m[i] = m[i-1]+ncol;
    }

    return m;
}

void free_vector(float *v, long nl, long nh)
{
    free((char *) (v+nl-1));
}


