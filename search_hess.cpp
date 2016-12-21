/* Calculate the Hessian of search function */
/* method=0: constant K, method=1: dynamic K */

#include <stdio.h>
#include <math.h>


void search_hess(float* x, float* xc, float R2, float* xiT, float* w, float* v, int method,float* g,float* HT)
{
    float p, g1[d], ge[d], He[d*d];
    
    /* calculation of the interpoling value, gradient and Hessian */
    inter_val(x, xiT, w, v, p);
    inter_grad(x, xiT, w, v, g1);
    inter_hess(x, xiT, w, v, HT);
    
    for (int i=0;i<d;i++) g[i]=g1[i];
    /* calculation of the error function, gradient and Hessian*/
   float  e=R2;
    for (int i=0;i<d;i++)
    {
        e+=-pow((x[i]-xc[i]),2);
        ge[i]=-2*(x[i]-xc[i]);
    }
    for (int i=0;i<d;i++)
        for (int j=0;j<d;j++)
        {
            He[i*d+j]=0;
            if (i==j) He[i*d+j]=-2;
        }
    /* calculation of the hessian f the search function*/
    
if (method==0)
         for (int i=0;i<d;i++)
             g[i]+=K*ge[i];
         else
         for (int i=0;i<d;i++)
            g[i]=g[i]/e-(p-y0)*ge[i]/pow(e,2);
      
if (method==0)
         for (int i=0;i<d*d;i++)
           HT[i]+=K*He[i];
         else
         for (int i=0;i<d;i++)
             for (int j=0; j<d;j++)
             HT[i*d+j]=HT[i*d+j]/e-(g1[i]*ge[j]+g1[j]*ge[i])/pow(e,2)+(p-y0)*(-He[i*d+j]/pow(e,2)+2*ge[i]*ge[j]/pow(e,3));
         
    
}
