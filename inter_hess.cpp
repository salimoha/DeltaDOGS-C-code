/* Calculate the derivative of the interpolation at the given point */

#include <stdio.h>
#include <math.h>

void inter_hess(float* x, float* xiT, float* w, float* v, float* HT)

{
    
    float  p1 ;/* single precision!!! */
    
    /* initialize The Hessian matrix */
    
     for (int i=0; i<d*d; i++)
         HT[i]=0;
    
    /* Calculate the Hessian matrix*/
    for (int k=0; k<n; k++)
    {
    
        /* Calculate p_1=||x-x_k|| */
        p1=0;
        for (int l=0; l<d; l++)
            p1+=pow((x[l]-xiT[k*d+l]),2);
            p1=pow(p1,0.5);
        
      if (p1>1e-4)
       {
         for (int i=0; i<d; i++)
         for (int j=0; j<d; j++)
           {
              HT[i*d+j]+=3*w[k]*(x[i]-xiT[k*d+i])*(x[j]-xiT[k*d+j])/p1;
            if (i==j) HT[i*d+j]+=3*w[k]*p1;
             
           }
        }
     }

}
