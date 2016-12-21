/* solving the matrix equation A*x=b using LAPACK */

#include <stdio.h>
#include <math.h>


#define d 2				/* dimension of matrix */
void newton_direc(float* HT, float* g);
int main()
{
    float A[d][d], b[d],AT[d*d];	/* single precision!!! */
    
    
    A[0][0]=0.5466;  A[0][1]=-0.6019;  //A[0][2]=1; A[0][3]=1;
    
    A[1][0]=A[0][1];  A[1][1]=0.8867; // A[1][2]=2.8284; A[1][3]=2;
    
    /*
    A[2][0]=A[0][2];  A[2][1]=A[1][2];  A[2][2]=1; A[2][3]=3;
    
    A[3][0]=A[0][3];  A[3][1]=A[1][3];  A[3][2]=A[2][3]; A[3][3]=0;
    */
    b[0]=0.2147;			/* if you define b as a matrix then you */
    b[1]=0.5516;			/* can solve multiple equations with */
   // b[2]=1;
   // b[3]=1;
    /* the same A but different b */
    
    for (int i=0;i<d;i++)
        for (int j=0;j<d;j++)
            AT[d*i+j]=A[i][j];
    
    
    
    
    newton_direc(AT, b);
   for (int i=0;i<d;i++)
    printf("%e\n",b[i]);
    
   }
void newton_direc(float* HT, float* g)
{
    
    float alpha, beta, theta, A[d][d], L[d][d], c[d][d], D[d];/* single precision!!! */
    
    /*  The parameters for the hessian modification */
    alpha=0.5;beta=2;
    
    for (int i=0;i<d;i++)
        for (int j=0;j<d;j++)
            A[i][j]=HT[d*i+j];

            /*
     for (i=0;i<d*d;i++)
     printf("%e/n",HT[i]);
     */
    /* Define L=I */
    for (int i=0; i<d; i++)
    {
        for (int j=0; j<d; j++)
            L[i][j]=0;
        L[i][i]=1;
    }
    
    /* D(0)=max(alpha,A(0,0))*/
    /* Calculate the modified LDL decomposition */
    D[0]=A[0][0];
    if (D[0]<-A[0][0]) D[0]=-A[0][0];
    if (D[0]<alpha) D[0]=alpha;
    
    for (int j=0; j<d; j++)
    {
        for (int i=j;i<d;i++)
        {
            c[i][j]=A[i][j];
            for (int k=0; k<j; k++)
                c[i][j]-=L[j][k]*L[i][k]*D[k];
        }
        if (j<d-1)
        {
        theta=c[j+1][j];
        for (int k=j+2; k<d; k++)
            if (c[k][j]>theta) theta=c[k][j];
        D[j]=pow((theta/beta),2);
        }
        else
            D[j]=0;
        if (alpha>D[j]) D[j]=alpha;
        if (c[j][j]>D[j]) D[j]=c[j][j];
        if (-c[j][j]>D[j]) D[j]=-c[j][j];
        
        for (int i=j+1;i<d;i++)
            L[i][j]=c[i][j]/D[j];
    }
  //  for (int i=0;i<d;i++)
    //    printf("%e\n",L[i][1]);
    /* Solve the linear system */
    
    for (int i=0;i<d;i++)
    {
        for (int j=0;j<i;j++)
            g[i]-=g[j]*L[i][j];
        g[i]=g[i]/L[i][i];
    }
    for (int i=0;i<d;i++)
        g[i]/=D[i];

    for (int i=d-1;i>-1;i--)
        for (int j=d-1;j>i;j--)
            g[i]-=g[j]*L[j][i];
    for (int i=0;i<d;i++)
        printf("%e\n",g[i]);
    
}