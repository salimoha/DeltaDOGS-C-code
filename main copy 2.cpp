/* Implementing the delaunay triangulation */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <clapack.h>
#include <set>


# define d 5	/* dimension of matrix */

struct simplex {
    int p[d+1];    /* The pointers of vertices of the simplex*/
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
    int simplex[100];
    int simplex_length;
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
simplex tri[1000];
int killed[1000];
int k_all;
point xi[20];
void  circume(point* xir, float* xc, float &R2);
void  test_func(point x, float& y);
void min_weight(float* x, float &w, int &ind, simplex tri);
void initial_simplex( simplex &tri);
void  simplex_find(int &ind, float* xn);
void  face_find(simplex* tri, point* xi,float* xn, int ind, int* killed, int &k_all);
void  new_simplex(simplex* tri, point* xi, int Ni, int N_new, int k_all);
#include "inter_coeff.cpp"
#include "simplex_search.cpp"


int main()
{
/* Put the initial points to the data set */
    FILE *file;
    /* Read the xi elements */
    file = fopen("xi.txt", "r");
    Ni=22;
    for (int i=0;i<Ni;i++) for (int j=0; j<d;j++) fscanf(file, "%f", &xi[i].x[j]);
    fclose(file);
    /* Read the tri elements */
    file = fopen("tri.txt", "r");
    N_tri=363;
    for (int i=0;i<N_tri;i++) for (int j=0; j<d+1;j++) {fscanf(file, "%d", &tri[i].p[j]); tri[i].p[j]--; }
    fclose(file);
    /* Read the tri elements */
    file = fopen("adj.txt", "r");
    for (int i=0;i<N_tri;i++) for (int j=0; j<d+1;j++) {fscanf(file, "%d", &tri[i].ad[j]); tri[i].ad[j]--;}
    fclose(file);
    // Calculate the circume center of all simplices
    point xir[d+1];
     for (int i=0; i<N_tri;i++)
        {
     for (int j=0; j<d+1; j++) for (int k=0; k<d;k++) xir[j].x[k]=xi[tri[i].p[j]].x[k];
         circume(xir, tri[i].xc, tri[i].R2);
            tri[i].stat=1;
        }
    // Calculate the new point
    int ind=0;
    float xn[d];
    for (int k=0; k<d; k++)
    {xn[k]=0;
        for (int j=0; j<d+1; j++)
            xn[k]+=xi[tri[ind].p[j]].x[k];
    xn[k]/=(d+1);}
    // Update the delaunay triangulation
    face_find(tri,xi,xn,ind,killed,k_all);
    new_simplex(tri,xi,Ni,N_new,k_all);
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



void  face_find(simplex* tri, point* xi, float* xn, int ind, int* killed, int &k_all)
{
    int k_read=0;
    k_all=1;
    killed[k_read]=ind;
    N_new=0;
    while (k_read<k_all)
    {
        ind=killed[k_read];
        k_read++;
        tri[ind].stat=0;
        for (int i=0;i<d+1;i++)
        {
            // find the i-the neighber of this simplex //
            int k=tri[ind].ad[i];
            if (k==-1)
            {       // constrcut the simplex
                    // calculate its point
                    for (int l=0;l<(d+1);l++) tri[N_tri+N_new].p[l]=tri[ind].p[l]; tri[N_tri+N_new].p[i]=Ni;
                    // calculate known neighber
                    tri[N_tri+N_new].ad[i]=k;
                    // calculate the circume center
                    point xir[d+1];
                    for (int l=0;l<d+1;l++) for (int t=0;t<d;t++) xir[l].x[t]=xi[tri[N_tri+N_new].p[l]].x[t];
                    circume(xir,tri[N_tri+N_new].xc,tri[N_tri+N_new].R2);
                    // calculate stat and ill
                    tri[N_tri+N_new].stat=1;
                    // Check the illness of the new simplex
                    tri[N_tri+N_new].ill=0;
                //for (int ss=0; ss<d+1; ss++) if (tri[N_tri+N_new])
                    
                    // change the index
                    N_new++;
            }
            else
                if (tri[k].stat==1)
                {
                    float p=0;
                    for (int l=0;l<d;l++) p+=pow((tri[k].xc[l]-xn[l]),2);
                    if (p<tri[k].R2) {tri[k].stat=2; killed[k_all]=k; k_all++;}
                    else
                    {
                        // constrcut the simplex
                        // calculate its point
                        for (int l=0;l<d+1;l++) tri[N_tri+N_new].p[l]=tri[ind].p[l]; tri[N_tri+N_new].p[i]=Ni;
                        // calculate known neighber
                        tri[N_tri+N_new].ad[i]=k;
                        // calculate the circume center
                        point xir[d+1];
                        for (int l=0;l<d+1;l++) for (int t=0;t<d;t++) xir[l].x[t]=xi[tri[N_tri+N_new].p[l]].x[t]; circume(xir,tri[N_tri+N_new].xc,tri[N_tri+N_new].R2);
                        for (int l=0;l<d+1;l++) if (tri[k].ad[l]==ind) {tri[k].ad[l]=N_tri+N_new; break; }
                        // calculate stat and ill
                        tri[N_tri+N_new].stat=1; tri[N_tri+N_new].ill=0;
                        // change the index
                        N_new++;
                        
                    }
                    
                }
        }
    }
    printf("%i\n",k_all);
}
/* Find the  adjancey relationship for new simplex*/
void  new_simplex(simplex* tri, point* xi, int Ni, int N_new, int k_all)
{
    int N1=N_tri+N_new-1;
    for (int k=0; k<N_new-1; k++)
        for (int i=0; i<d+1; i++)
            if(tri[N1-k].p[i]!=Ni)
            {
                int c1[d+1]; for(int j=0; j<d+1; j++) c1[j]=tri[N1-k].p[j]; c1[i]=-1;
                for (int k1=k+1; k1<N_new; k1++)
                    for (int i1=0;i1<d+1; i1++)
                    { int c2[d+1]; for(int j=0; j<d+1; j++) c2[j]=tri[N1-k1].p[j]; c2[i1]=-1;
                        std::set<int> s1 (c1,c1+(d+1));
                        std::set<int> s2(c2,c2+(d+1));
                        if (s1==s2)
                        {tri[N1-k].ad[i]=N1-k1; tri[N1-k1].ad[i1]=N1-k;
                            if (k<k_all) tri[N1-k1].ad[i1]=killed[k];
                            if (k1<k_all) tri[N1-k].ad[i]=killed[k1];
                        }
                    }
 
                }
             else
             {
                 int ss=tri[N_tri+k].ad[i];
                 if (ss!=-1)
            for (int l=0;l<d+1;l++) if (tri[ss].ad[l]==N1-k) {tri[ss].ad[l]=killed[k]; break; }
                     
             }
    for (int k=0; k<k_all; k++) tri[killed[k]]=tri[N1-k];
}




