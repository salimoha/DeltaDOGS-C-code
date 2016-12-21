/* Calculate the derivative of the interpolation at the given point */

#include <stdio.h>
#include <math.h>

/*
#define d 2
# define n 4
*/


void inter_grad(float* x, float* xiT, float* w, float* v, float* g)
{
    int i, j , k;
   
    for (j=0; j<d; j++)
    {
        g[j]=v[j+1];
        for (i=0; i<n; i++)
    {
       float p1=0;
        for (k=0; k<d; k++)
            p1+=pow((x[k]-xiT[i*d+k]),2);
        p1=pow(p1,0.5);
        g[j]+=3*p1*w[i]*(x[j]-xiT[i*d+j]);
    }
    }

}
