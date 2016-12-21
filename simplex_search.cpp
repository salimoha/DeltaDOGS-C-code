/* Calculate the Newton direction for the search function */
/* method=0: constant K, method=1: dynamic K */
#include "interpolate.cpp"
#include "interval.cpp"
void newton_direc(inter_info In, float* g);
void Newton_search_direction(float* x, float &val, float* P, float &r, int method, float* xc, float R2);
void feasible_const(float* x, float* xc,float R2, float* p, float &alpha);
void search_val(float* x, float R2, float* xc, int method, float& p);
void simplex_search(float* x, float &p, float R2, float* xc, int &method)
{
    float alpha, etas, rho, xn[d], pn, r, P[d];
    etas=0.01; rho=0.9;
    method=0;
    while (1)
    {
        Newton_search_direction(x,p,P,r,method,xc,R2);
        feasible_const(x,xc,R2,P,alpha);
        /* Breaking conditions on the boundary or derivative zero*/
        float err=0; for (int i=0;i<d;i++) err+=P[i]*P[i];
        if (err<1e-4) break;
        float c=0, dd=0;
        for (int i=0;i<d;i++) {c+=pow(x[i]-xc[i],2); dd+=(x[i]-xc[i])*P[i];}
        if (c>0.9*R2 && dd>0) break;
        c=0, dd=0;
        for (int i=0;i<d;i++) {c+=pow(x[i],2); dd+=(x[i])*P[i];}
        if (c>0.96 && dd>0) break;
        
        
        while (1)
        {
            for (int i=0;i<d;i++) xn[i]=x[i]+alpha*P[i];
            search_val(xn,R2,xc,method,pn);
            if ((pn-p)/(r*alpha)>etas)
            {
                for(int i=0;i<d;i++) x[i]=xn[i]; p=pn;
                break;
            }
            else alpha*=rho;
        }
        if (method<2 & p<0) method++;
    }

    
}
void Newton_search_direction(float* x, float &val, float* P, float &r, int method, float* xc, float R2)
{
    inter_info In;
    Interpolate(x,In);
    for (int i=0;i<d;i++) P[i]=-In.g[i];
    /* calculation of the error function, its gradient, and its  */
if (method<2)
{
   /* calculation of the error function, gradient and Hessian*/
    float  e, ge[d], He[d][d];
    e=R2;
    for (int i=0;i<d;i++)
        for (int j=0;j<d;j++)
            He[i][j]=0;
    for (int i=0;i<d;i++)
    {
        e+=-pow((x[i]-xc[i]),2);
        ge[i]=-2*(x[i]-xc[i]);
        He[i][i]=-2;
    }
    
    /* Calculate the hessian of the search function*/
    
if (method==0)
{
         for (int i=0;i<d;i++)
         {
             In.g[i]-=K_max*ge[i];
             In.H[i][i]-=K_max*He[i][i];
         }
    In.p-=K_max*e+y_0;
 for (int i=0;i<d;i++) P[i]=-In.g[i];
}
else
{
         for (int i=0;i<d;i++)
         {
            P[i]=-(In.g[i]/e-(In.p-y_0)*ge[i]/(e*e));
            for (int j=0; j<d;j++)
            {In.H[i][j]/=e; In.H[i][j]+=-(In.g[i]*ge[j]+In.g[j]*ge[i])/(e*e)+(In.p-y_0)*(-He[i][j]/(e*e)+2*ge[i]*ge[j]/(e*e*e));}
         }
    In.p=(In.p-y_0)/e;
    for (int i=0;i<d;i++) In.g[i]=-P[i];
}
}
    newton_direc(In,P);
    /* The directional derivative over the P direction */
    r=0; for (int i=0;i<d;i++) r+=P[i]*In.g[i];
    val=In.p;
}
/* Restrict the line search to the feasible domain */
 void feasible_const(float* x, float* xc,float R2, float* p, float &alpha)
 {
 float a=0, b=0, c=-0.95*R2, b1=0, alpha1;
 float c1=-0.95;
 for (int i=0; i<d; i++)
 {
 a+=pow(p[i],2);
 b+=2*(x[i]-xc[i])*p[i];
 c+=pow((x[i]-xc[i]),2);
 b1+=2*(x[i])*p[i];
 c1+=pow((x[i]),2);
 }
 alpha=(sqrt(b*b-4*a*c)-b)/(2*a);
 alpha1=(sqrt(b1*b1-4*a*c1)-b1)/(2*a);
 if (alpha>1)
 alpha=1;
 if (alpha>alpha1)
 alpha=alpha1;
 
 }
/* Calcualate the value of the search function */
void search_val(float* x, float R2, float* xc, int method, float& p)

{
    /* Calculate the interpolating value */
    interval(x,p);
    /* Calculate the error function*/
    if (method<2)
    {
        float e=R2;
        for (int i=0; i<d; i++) e-=pow((x[i]-xc[i]),2);
        if (method==0)
            p-=K_max*e+y_0;
        else
            p=(p-y_0)/e;
    }
    
    
}

/* Calculate the newton direction using the hessian and the derivative */
void newton_direc(inter_info In, float* g)
{
 
    float alpha, beta, theta, A[d][d], L[d][d], c[d][d], D[d];/* single precision!!! */
    
    /*  The parameters for the hessian modification */
    alpha=0.01;beta=100;
    
    for (int i=0;i<d;i++)
        for (int j=0;j<d;j++)
            A[i][j]=In.H[i][j];
    
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
    /* Solve the linear system */
    
    for (int i=0;i<d;i++)
    {
        for (int j=0;j<i;j++)
            g[i]-=g[j]*L[i][j];
        g[i]=g[i]/L[i][i];
    }
    for (int i=0;i<d;i++)
        g[i]=g[i]/D[i];
    
    for (int i=d-1;i>-1;i--)
    {
        for (int j=d-1;j>i;j--)
            g[i]-=g[j]*L[j][i];
        g[i]=g[i]/L[i][i];
    }
    
    
}