/* Calculate the gradient of search function */
/* method=0: constant K, method=1: dynamic K */

#include <stdio.h>
#include <math.h>

void search_grad(float* x, float* xc, float R2, float* xiT, float* w, float* v, int method, float* g)
{
    float p, e, ge[d];

    /* calculation of the interpoling value, and gradient */
    inter_val(x, xiT, w, v, p);
    inter_grad(x, xiT, w, v, g);
    /* calculation of the error function and its gradient*/
    e=R2;
    for (int i=0;i<d;i++)
    {
        e+=-pow((x[i]-xc[i]),2);
        ge[i]=-2*(x[i]-xc[i]);
    }
    
         if (method==0)
         for (int i=0;i<d;i++)
             g[i]+=K*ge[i];
         else
         for (int i=0;i<d;i++)
            g[i]=g[i]/e-(p-y0)*ge[i]/pow(e,2);
    
}
