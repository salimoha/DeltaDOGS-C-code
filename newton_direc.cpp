/* Calculate the Newton direction using modified cholesky decomposition*/
/* A:hessian, b:gradient */

#include <stdio.h>
#include <math.h>

void newton_direc(float* HT, float* g)
{
    int i, j, k;
    
    float alpha, beta, theta, A[d][d], L[d][d], c[d][d], D[d];/* single precision!!! */
    
    /*  The parameters for the hessian modification */
    alpha=0.5;beta=10;
    
    for (i=0;i<d;i++)
        for (j=0;j<d;j++)
            A[i][j]=HT[d*i+j];
    /*
    for (i=0;i<d*d;i++)
        printf("%e/n",HT[i]);
    */
    /* Define L=I */
     for (i=0; i<d; i++)
     for (j=0; j<d; j++)
     {
         L[i][j]=0;
         if (i==j) L[i][j]=1;
     }
    
    /* D(0)=max(alpha,A(0,0))*/
/* Calculate the modified LDL decomposition */
    D[0]=A[0][0];
    if (D[0]<alpha) D[0]=alpha;
    /* c(:,1)=A(:,1) */
    for (i=0; i<d; i++)
        c[i][0]=A[i][0];
    /* L(2:n,1)=c(2:n,1)/D(1) */
    for (i=1; i<d; i++)
        L[i][0]=c[i][0]/D[0];
    
    for (j=1; j<d; j++)
    {
        for (i=j;i<d;i++)
    {
            c[i][j]=A[i][j];
            for (k=0; k<j; k++)
               c[i][j]+=-1*L[j][k]*L[i][k]*D[k];
    }
        theta=c[j+1][j];
        for (k=j+2; k<d; k++)
            if (c[k][j]>theta) theta=c[k][j];
        
        D[j]=alpha;
         if (pow((theta/beta),2)>D[j]) D[j]=pow((theta/beta),2);
         if (c[j][j]>D[j]) D[j]=c[j][j];
         if (-c[j][j]>D[j]) D[j]=-c[j][j];
        
        for (i=j+1;i<d;i++)
           L[i][j]=c[i][j]/D[j];
    }
/* Solve the linear system */

    for (i=0;i<d;i++)
    {
        for (j=0;j<i;j++)
            g[i]+=-1*(g[j]*L[i][j]);
        g[i]=g[i]/L[i][i];
    }
    for (i=0;i<d;i++)
        g[i]=g[i]/D[i];
    
    for (i=d-1;i>-1;i--)
    {
        for (j=d-1;j>i;j--)
            g[i]+=-1*(g[j]*L[j][i]);
        g[i]=g[i]/L[i][i];
    }
    
    
}
