/* Find a local minimum for the search function in a simplex */
/* We have considered the simple sphere constraints ||x-x0||< R2*/
/* method=0: constant K, method=1: dynamic K */

#include <stdio.h>
#include <math.h>
#include <clapack.h>

# define d 2	/* dimension of matrix */
# define n 4
# define K 1
# define y0 -2

#include "search_hess1.cpp"


void  inter_val(float* x, float* xiT, float* w, float* v, float R2, float* xc, int method, float& p);
void  search(float* x, float* xiT, float* w, float* v, float R2, float* xc, float &p, int method);
void  circume(float* xiT, float* xc, float &R2);
void feasible_const(float* x, float* xc,float R2, float* p, float &alpha);

int main()
{
    float xi[d][n], xiT[d*n], p, xir[d*(d+1)], xc[d], R2;
    /* Initialize the points */
    
    xi[0][0]=0;  xi[0][1]=1;  xi[0][2]=0;	/* matrix xi=[x_0,x_1,x_2] */
    xi[1][0]=0;  xi[1][1]=0.5; xi[1][2]=10;
    
    for (int i=0; i<d; i++)		/* put the points in a vector */
        for(int j=0; j<d+1; j++) xir[d*j+i]=xi[i][j];
    
    
    circume(xir,xc,R2);
    
    for (int j=0; j<d; j++) printf("%e\n", xc[j]);	/* print vector x */
    
    printf("R2 is equal to %e\n", R2);
    
    xi[0][0]=0;  xi[0][1]=1;  xi[0][2]=0; xi[0][3]=0.5;
    xi[1][0]=0;  xi[1][1]=0; xi[1][2]=1; xi[1][3]=0.5;
    
    for (int i=0; i<d; i++)
       for (int j=0; j<n; j++)
        xiT[j*d+i]=xi[i][j];
    

    
    /* Initialize the coefficients */
    
    float w[n]={0,0.3536,0.3536,-0.7071};
    
    float v[d+1]={-0.4571,0.7071,0.7071};
    
    /* Calculation point */
    float  x[d];
    for (int i=0;i<d; i++) x[i]=xc[i];
    
    
   /* calculation of the Newton direction */
    int  method=1;
    
   /* calculate until we found a negative point */
    if (method==1)
    {
        search(x,xiT,w,v,R2,xc,p,method);
    if (p<0)
    {
        method=3;
        search(x,xiT,w,v,R2,xc,p,method);
    }
    }
    else
        search(x,xiT,w,v,R2,xc,p,method);
     for (int i=0;i<d;i++) printf("%e\n",x[i]);
     printf("%e\n",p);
     printf("%e\n",1.0*method);
}

void feasible_const(float* x, float* xc,float R2, float* p, float &alpha)
{
    float a=0, b=0, c=-0.95*R2, b1=0, alpha1;
    float x0[d];
    for (int i=0;i<d;i++) x0[i]=xc[i];
    float R20=6, c1=-0.95*R20;
    for (int i=0; i<d; i++)
    {
        a+=pow(p[i],2);           
        b+=2*(x[i]-xc[i])*p[i];
        c+=pow((x[i]-xc[i]),2);
        b1+=2*(x[i]-x0[i])*p[i];
        c1+=pow((x[i]-x0[i]),2);
    }
    alpha=(sqrt(b*b-4*a*c)-b)/(2*a);
    alpha1=(sqrt(b1*b1-4*a*c1)-b1)/(2*a);
       if (alpha>1)
        alpha=1;
    if (alpha>alpha1)
        alpha=alpha1;
    
}
void inter_val(float* x, float* xiT, float* w, float* v, float R2, float* xc, int method, float& p)

{
/* Calculate the interpolating value */
    p=v[0];
    for (int i=0; i<d; i++)
        p+=x[i]*v[i+1];
    
    for (int k=0; k<n; k++)
    {
     float   p1=0;
        for (int l=0; l<d; l++)
            p1+=pow((x[l]-xiT[k*d+l]),2);
        p+=w[k]*pow(p1,1.5);
    }
/* Calculate the error function*/
if (method<2)
{
    float e=R2;
    for (int i=0; i<d; i++) e-=pow((x[i]-xc[i]),2);
 if (method==0)
     p-=K*e;
 else
    p=(p-y0)/e;
}


}

void search(float* x, float* xiT, float* w, float* v, float R2, float* xc, float &p, int method)
{
float g[d],xn[d],alpha, pn, P[d], etas=0.01, rho=0.9, tol=1e-4;

    while (1)
 {
    search_hess1(x,xc,R2,xiT,w,v,method,p,g,P);
    feasible_const(x,xc,R2,P,alpha);
     float c=0, dd=0;
     for (int i=0;i<d;i++) {c+=pow(x[i]-xc[i],2); dd+=(x[i]-xc[i])*P[i];}
     //printf("dd value\n");
     //printf("%e\n",dd);
     if (c>0.9*R2 && dd>0) break;
    while (1)
    {
        for (int i=0;i<d;i++) xn[i]=x[i]+alpha*P[i];
            inter_val(xn,xiT,w,v,R2,xc,method,pn);
            float ff=0;
            for (int i=0; i<d; i++) ff+=g[i]*P[i]*alpha;
                if ((pn-p)/ff>etas)
                {
                    for(int i=0;i<d;i++) x[i]=xn[i];
                    p=pn;
                        break;
                   // for(int i=0;i<d;i++) printf("%e\n",x[i]);
                }
                else alpha*=rho;
    }
     float err=0;
    for(int i=0;i<d;i++) err+=pow(alpha*P[i],2);
           if (err<tol)
              break;
           if (method==1 && p<0)
              break;
 }
}
void  circume(float* xiT, float* xc, float &R2)
{
    
    float MT[d*d];	/* single precision!!! */
    
    for (int i=0; i<d; i++)/* construct the M matrix and vector b*/
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
}




