/* Calculate the coefficient of the polyharmonic spline interpolation*/
// The first d+1 points are not used in the interpolation process
void inter_coeff(point* xi, float* v, int n)
{
    float* MT;
    float* b;
    int* pivot;
     MT= new float[(n+d+1)*(n+d+1)];
     b= new float[n+d+1];
    pivot= new int[n+d+1];
    
        /* construct the F part in M*/
    for (int i=0; i<n; i++)/* construct the M matrix and vector b*/
    {
        for(int j=0; j<n; j++)
        {
		MT[i+(n+d+1)*j]=0;
            for (int k=0; k<d; k++)
                MT[i+(n+d+1)*j]+=pow((xi[i+d+1].x[k]-xi[j+d+1].x[k]),2);
            MT[i+(n+d+1)*j]=pow(MT[i+(n+d+1)*j],1.5);
        }
    }
    /* construct the V part in M*/
    for (int i=0; i<n; i++)/* construct the M matrix and vector b*/
    {
        MT[n+(n+d+1)*i]=1;  MT[i+(n+d+1)*n]=1;
        for(int j=0; j<d; j++)
        {
            MT[i+(n+d+1)*(n+j+1)]=xi[i+d+1].x[j];
            MT[n+j+1+(n+d+1)*i]=xi[i+d+1].x[j];
        }
    }
   

    /* construct the zero part in M*/
    for (int i=n; i<n+d+1; i++)/* construct the M matrix and vector b*/
     for(int j=n; j<n+d+1; j++)
            MT[i+(n+d+1)*j]=0;
        /* constrcut the b */
    
    for (int i=0; i<n; i++)
        b[i]=xi[d+1+i].y;
    for (int i=n; i<n+d+1; i++)
        b[i]=0;
   /* calculate the linear system */
   int  c1=n+d+1, c2=1, ok;    			/* to the routine in variables */

    sgesv_(& c1,& c2, MT, &c1, pivot, b, &c1, &ok);
    
        for (int i=0; i<n; i++)
        xi[i+d+1].w=b[i];
    for (int i=0; i<d+1; i++)
        v[i]=b[i+n];
    
}
