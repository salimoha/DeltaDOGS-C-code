/* solving the matrix equation A*x=b using LAPACK */

#include <stdio.h>
#include <clapack.h>


#define size 7				/* dimension of matrix */

int main()
{
    int i, j , c1, c2, pivot[size], ok;
    float A[size][size], b[size], AT[size*size];	/* single precision!!! */
    
    
    A[0][0]=0;  A[0][1]=1;  A[0][2]=1;
    A[0][3]=0.3536; A[0][4]=1; A[0][5]=0; A[0][6]=0;
    
    A[1][0]=1;  A[1][1]=0;  A[1][2]=2.8284;
    A[1][3]=0.3536; A[1][4]=1; A[1][5]=1.0; A[1][6]=0;
    
    A[2][0]=1;  A[0][1]=1;  A[0][2]=1;
    A[2][3]=0.3536; A[0][4]=1; A[0][5]=0; A[0][6]=1;
    
    A[3][0]=0.3536;  A[3][1]=1;  A[3][2]=1;
    A[3][3]=0.3536; A[3][4]=1; A[3][5]=0; A[3][6]=0.5;
    
    A[4][0]=1;  A[4][1]=1;  A[4][2]=1;
    A[4][3]=0.3536; A[4][4]=1; A[4][5]=0; A[4][6]=0;
    
    A[5][0]=0;  A[5][1]=1;  A[5][2]=1;
    A[5][3]=0.3536; A[5][4]=1; A[5][5]=0; A[5][6]=0;
    
    A[6][0]=0;  A[6][1]=1;  A[6][2]=1;
    A[6][3]=0.3536; A[6][4]=1; A[6][5]=0; A[6][6]=0;
    
    b[0]=1;			/* if you define b as a matrix then you */
    b[1]=1;			/* can solve multiple equations with */
    b[2]=1;			/* the same A but different b */
    
    for (i=0; i<size; i++)		/* to call a Fortran routine from C we */
    {				/* have to transform the matrix */
        for(j=0; j<size; j++) AT[j+size*i]=A[j][i];
    }
    
    c1=size;			/* and put all numbers we want to pass */
    c2=1;    			/* to the routine in variables */
    
    /* find solution using LAPACK routine SGESV, all the arguments have to */
    /* be pointers and you have to add an underscore to the routine name */
    sgesv_(&c1, &c2, AT, &c1, pivot, b, &c1, &ok);
    
    /*
     parameters in the order as they appear in the function call
     order of matrix A, number of right hand sides (b), matrix A,
     leading dimension of A, array that records pivoting,
     result vector b on entry, x on exit, leading dimension of b
     return value */ 
    
    for (j=0; j<size; j++) printf("%e\n", b[j]);	/* print vector x */
    return 1;
}
