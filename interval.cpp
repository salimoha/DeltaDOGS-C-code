/* Calculate the value of the polyharmonic spline interpolation */
void interval(float* x,float &p)
{
    /* calculation of the interpoling value, gradient and Hessian */
    p=v[0];
    for (int i=0; i<d; i++)
        p+=x[i]*v[i+1];
    for (int k=0; k<Ni; k++)
        
    {     /* Calculate p_1=||x-x_k|| */
        float  p1=0;
        for (int l=0; l<d; l++)
            p1+=pow((x[l]-xi[k+d+1].x[l]),2);
            p1=sqrt(p1);
        if (p1>1e-4) p+=pow(p1,3)*xi[k+d+1].w;
    }
    
}
