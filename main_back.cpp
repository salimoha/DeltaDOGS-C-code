/* Implementing the delaunay triangulation */
#include <stdio.h>
#include <math.h>
#include <clapack.h>
#include <set>


# define d 2	/* dimension of matrix */

struct simplex {
    int p[d+1];    /* The pointers of vertices of the simplex*/
    int ad[d+1];  /* The pointers of neighbers of the simplex*/
    int stat; /* stat=0: killed simplex */
    int ill; /* the illposed simplex */
    float xc[d]; /* circume center of the simplex */
    float xm[d]; /* the best point of the simplex */
    float R2; /* square of circumeradius of the simplex */
    float pm; /* Value of the search function */
};
// defining the structure for simplex
struct point {
    float x[d];/* The pointers of vertices of the simplex*/
    float y;
    float w;
};

struct inter_info {
    float p;
    float g[d];/* The pointers of vertices of the simplex*/
    float H[d][d];
};
float v[d+1];
inter_info In;
int method;
float y_0, x0[d], R20, K_max;
/*  The boundary circle */
int Ni, N_tri, N_new;
simplex tri[20];
point xi[20];
void  circume(point* xir, float* xc, float &R2);
void  test_func(point x, float& y);
void min_weight(float* x, float &w, int &ind, simplex tri);
void initial_simplex( simplex &tri);
void  simplex_find(int &ind, float* xn);
#include "inter_coeff.cpp"
#include "simplex_search.cpp"


int main()
{
    /* Initialization of the points and simplics */
// Calculate an uniform simplex around origin
    for (int i=2;i<d+1; i++) x0[i]=0;
    R20=1; y_0=0; K_max=10;
    method=1;
    xi[0].x[0]=-1; xi[1].x[0]=1;
    for (int i=2;i<d+1; i++)
    {
        for (int j=0;j<i; j++){for(int k=0;k<(i-1); k++) xi[j].x[k]*=sqrt(1-1/pow(i,2)); xi[j].x[i-1]=-1.0/i; }
            for(int k=0;k<(i-2); k++) xi[i].x[k]=0; xi[i].x[i-1]=1;
    }
// Calculate the points
   for (int i=0;i<d+1; i++) for (int j=0; j<d; j++) { xi[d+1+i].x[j]=xi[i].x[j];  xi[i].x[j]*=-(d+2.0)/d;}
// Calculate the simplices
   for (int i=0;i<d+1; i++)
   {
       point xir[d+1];
       for (int j=0; j<d; j++)  {for (int k=0; k<d+1;k++) xir[k].x[j]=xi[d+1+k].x[j];xir[i].x[j]=xi[i].x[j];}
       circume(xir, tri[i].xc, tri[i].R2);
       for (int k=0; k<d+1;k++) {tri[i].ad[k]=-1; tri[i].p[k]=k+d+1;}
       tri[i].p[i]=i; tri[i].ad[i]=d+1;
       tri[i].stat=1; tri[i].ill=1;
   }
    for (int k=0; k<d+1;k++) {tri[d+1].ad[k]=k; tri[d+1].p[k]=k+d+1;}
    for (int j=0; j<d; j++) {tri[d+1].xc[j]=0;tri[d+1].xm[j]=0;}
    tri[d+1].stat=1; tri[d+1].ill=0; tri[d+1].R2=1;
    
    Ni=d+1; N_tri=d+2;
    /* calculater the cost function at initial points */
    for (int i=0; i<d+1; i++) {float y; test_func(xi[i+d+1],y); xi[i+d+1].y=y;}
    /* Calculate the coefficient of the interpolation */
    inter_coeff(xi,v, Ni);
    /* Calculate the minimum of the search function */
   
    
    float p_min=K_max;
    float p1, xn[d];
    int ind;
    for (int ii=0;ii<N_tri;ii++)
        if (tri[ii].ill==0 && tri[ii].stat==1)
        {
            initial_simplex(tri[ii]);
            simplex_search(tri[ii].xm,p1,tri[ii].R2, tri[ii].xc, method);
            if (method>1)
            {
                for (int j; j<d; j++) xn[j]=tri[ii].xm[j]; ind=ii;
            break;
            }
            else
                if (p1<p_min)
                {for (int j=0; j<d; j++) xn[j]=tri[ii].xm[j]; ind=ii; p_min=p1;}
          
        }
     if (method<2)
     {float w=0; int index; min_weight(xn,w,index,tri[ind]);
         if (w<-1e-2)
         {
             float dis=0;
             for(int i=0; i<d;i++) dis+=pow((xn[i]-x0[i]),2);
             for(int i=0; i<d;i++) xn[i]=x0[i]+(xn[i]-x0[i])*sqrt(R20/dis);
       /* Check if the new point is in circume sphere of the original simplex */
             dis=0 ;for(int i=0; i<d;i++) dis+=pow((xn[i]-tri[ind].xc[i]),2);
             if (dis>tri[ind].R2)            simplex_find(ind,xn);
         }
     }
}
/* Find the simplex which includes the new point */
void simplex_find(int &ind, float* xn)
{
float w, dis; 
int index;
while (1)
{
min_weight(xn, w, index, tri[ind]);
ind=tri[ind].ad[index]; 
dis=0;for(int i=0; i<d;i++) dis+=pow((xn[i]-tri[ind].xc[i]),2);
if (dis<tri[ind].R2)
break;
}
}
void  circume(point* xir, float* xc, float &R2)
/* Calculation of the circume center and circume radius */
{
    float MT[d*d];	/* single precision!!! */
    for (int i=0; i<d; i++)/* construct the M matrix and vector b*/
    {
        xc[i]=0;
        for(int j=0; j<d; j++)
        {
            MT[j*d+i]=2*(xir[0].x[j]-xir[i+1].x[j]);
            xc[i]+=pow(xir[0].x[j],2)-pow(xir[i+1].x[j],2);
        }
        
    }
    int c1=d, c2=1, pivot[d], ok;
    sgesv_(&c1, &c2, MT, &c1, pivot, xc, &c1, &ok);
    
    /* Calculate the circume radius */
    R2=0;
    for (int j=0; j<d; j++) R2+=pow((xc[j]-xir[0].x[j]),2);
    
}
/* Calculae the weight of a point in a simplex */
void min_weight(float* x, float &w,int &ind, simplex tri)
{
    float w1[d+1], MT[(d+1)*(d+1)];
    for (int i=0; i<d; i++)
    {w1[i]=x[i]; for (int j=0; j<d+1;j++) MT[j*(d+1)+i]=xi[tri.p[j]].x[i];}
    for (int j=0; j<d+1;j++) MT[j*(d+1)+d]=1; w1[d]=1;
    int c1=d+1, c2=1, pivot[d+1], ok;
        sgesv_(&c1, &c2, MT, &c1, pivot, w1, &c1, &ok);
    for(int j=0; j<d+1; j++) if(w1[j]<w) {w=w1[j]; ind=j;}
}
/* Calculate the initial point of the simplex */
void initial_simplex( simplex &tri)
{
    float w=0;
    int  index;
    min_weight(tri.xc,w,index, tri);
    float alpha=1.0/(1.0-w);
    for(int i=0; i<d; i++) {tri.xm[i]=alpha*tri.xc[i]; for (int k=0; k<d+1;k++) tri.xm[i]+=(1-alpha)*xi[tri.p[k]].x[i]/(d+1); }
    
}

/* Calculation of the test function */
void  test_func(point x, float &y)
{
    /* the parabola test function */
    float x0[d], x00=0.3;
    for (int i=0; i<d; i++) x0[i]=x00;
    y=0;
    for (int i=0; i<d; i++) y+=pow((x.x[i]-x0[i]),2);
}


