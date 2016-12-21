/* solving the matrix equation A*x=b using LAPACK */

#include <stdio.h>
#include <clapack.h>
#include <vector>
#include <map.h>


#define size 7				/* dimension of matrix */

int main()
{
    int i, j , c1, c2, pivot[size], ok;
    float A[size][size] ;	/* single precision!!! */

    vector<double> aa;
    
    float AT[49] = {  0  ,  1.0000  ,  1.0000  ,  0.3536  ,  1.0000  ,       0  ,       0  ,   1.0000  ,       0  ,  2.8284  ,  0.3536  ,  1.0000  ,  1.0000  ,       0, 1.0000  ,  2.8284  ,       0  ,  0.3536  ,  1.0000  ,       0  ,  1.0000, 0.3536  ,  0.3536  ,  0.3536  ,       0  ,  1.0000  ,  0.5000  ,  0.5000, 1.0000  ,  1.0000  ,  1.0000  ,  1.0000  ,       0  ,       0  ,       0, 0  ,  1.0000  ,       0  ,  0.5000  ,       0  ,       0  ,       0, 0  ,       0  ,  1.0000  ,  0.5000  ,       0  ,       0  ,       0 };
    float b[7] = { 0, 1, 1, 0.5, 0, 0, 0};
   
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


