/* Calculate the Hessian, gradient and value of the polyharmonic spline interpolation */
void Interpolate(float* x,inter_info &In)
{
    /* calculation of the interpoling value, gradient and Hessian */
    In.p=v[0];
    for (int i=0; i<d; i++)
    {
        In.p+=x[i]*v[i+1];
        In.g[i]=v[i+1];
        for (int j=0; j<d; j++)
            In.H[i][j]=0;
    }

    for (int k=0; k<Ni; k++)
    {
        
        /* Calculate p_1=||x-x_k|| */
        float  p1=0;
        for (int l=0; l<d; l++)
            p1+=pow((x[l]-xi[k+d+1].x[l]),2);
        p1=sqrt(p1);
        
        if (p1>1e-4)
        {
            for (int i=0; i<d; i++)
            {
                for (int j=0; j<d; j++)
                {
                    In.H[i][j]+=3*xi[k+d+1].w*(x[i]-xi[k+d+1].x[i])*(x[j]-xi[k+d+1].x[j])/p1;
                }
                In.g[i]+=3*p1*xi[k+d+1].w*(x[i]-xi[k+d+1].x[i]);
                In.H[i][i]+=3*xi[k+d+1].w*p1;

            }
        In.p+=pow(p1,3)*xi[k+d+1].w;
        }
    }
    
}
