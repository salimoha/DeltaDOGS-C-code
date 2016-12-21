/* solving the matrix equation A*x=b using LAPACK */

#include <stdio.h>
#include <clapack.h>
#include <math.h>

#define d 2 /* dimension of matrix */

void  circume(float* xiT, float* xc, float &R2)
{

    float MT[d*d], xi[d][d+1], xiT[d*(d+1)];	/* single precision!!! */
    
    for (int i=0; i<d; i++)/* construct the M matrix and vector xc*/
    {
        xc[i]=0;
        for(int j=0; j<d; j++)
        {
            MT[j*d+i]=-2*(xiT[j]-xiT[(i+1)*d+j]);
            xc[i]+=pow((xiT[j]-xiT[(i+1)*d+j]),2);
        }
    
    }
    int c1=d, c2=1, pivot[d], ok;
    sgesv_(&c1, &c2, MT, &c1, pivot, xc, &c1, &ok);
    
    /* Calculate the circume radius */
    R2=0;
    for (int j=0; j<d; j++) R2+=(xc[j]-xiT[j])*(xc[j]-xiT[j]);
    
     /* Calculate the circume radius */
    
    for (int j=0; j<d; j++) printf("%e\n", xc[j]);	/* print vector x */
    
     printf("R2 is eqaul to %e\n", R2);
}

