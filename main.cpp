/* Implementing the Delta_Dogs Algorithm  */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <clapack.h>
#include <set>
#include <vector>
#include <algorithm>


# define d 4	/* dimension of matrix */
# define cash 100
# define max_tri 1000
# define max_num 100

struct simplex {
    int p[d+1];    /* The pointers of vertices of the simplex*/
    int stat; /* stat=0: killed simplex */
    int ill; /* the illposed simplex */
    float xc[d]; /* circume center of the simplex */
    float xm[d]; /* the best point of the simplex */
    float R2; /* square of circumeradius of the simplex */
    float pm; /* Value of the search function */
    int method; /* the method that used for this simplex */
};
// defining the structure for simplex
struct point {
    float x[d];/* The pointers of vertices of the simplex*/
    float y;
    float w;
    int simplex[cash];
    int simplex_length;
};

struct inter_info {
    float p;
    float g[d];/* The pointers of vertices of the simplex*/
    float H[d][d];
};
/* Define the parameters for search */
float K_max=0.5, y_0=0;
int method=0;
/* Define variables for interpolation */
float v[d+1];
inter_info In;
/*  The boundary circle */
int Ni, N_tri, N_new;
/* Define values for the delaunay trianglation */
simplex tri[max_tri];
point xi[max_num];
int killed[cash], k_all;
/* Definition of functions */
void  circume(point* xir, float* xc, float &R2);
void initial_simplex( simplex &tri);
void neighber(int ind, int i, int &j);
void delaunay_update(simplex* tri, point* xi, float* xn, int ind);
void simplex_find(int &ind, float* xn);
void  insimplex_check(int &ind, float* xn, int &check);
void initial_simplex( simplex &tri);

/* include cost function evaluation, simplex.search and interpolation functions */
#include "test_func.cpp"
#include "inter_coeff.cpp"
#include "simplex_search.cpp"

/* The main code */

int main()
{
    /* Initialization of the points and simplics */
    // Calculate an uniform simplex around origin
    // Calculate the uniform simplex
    xi[0].x[0]=-1; xi[1].x[0]=1;
    for (int i=2;i<d+1; i++)
    {
        for (int j=0;j<i; j++){for(int k=0;k<(i-1); k++) xi[j].x[k]*=sqrt(1-1/pow(i,2)); xi[j].x[i-1]=-1.0/i; }
        for(int k=0;k<(i-2); k++) xi[i].x[k]=0; xi[i].x[i-1]=1;
    }
    // Calculate the scale factor for the neighber points
    float p=0; for (int i=0; i<d;i++) {p+=xi[0].x[i]*xi[1].x[i]; xi[i].simplex[0]=0; xi[i].simplex_length=1;}
    point xic[d+1];
    // scaling the points
    for (int i=0;i<d+1; i++) for (int j=0; j<d; j++) {xic[i].x[j]=xi[i].x[j]; xi[i].x[j]/=(1-i*0.01)*p;}
    // Calculate the initial simplex
    for (int i=0;i<d+1; i++) tri[0].p[i]=i;
    tri[0].stat=1; tri[0].ill=1;
    int ind=0;
    float xn[d];
    Ni=d+1; N_tri=1;
    for (int i=0;i<d+1; i++)
    {
        ind=N_tri-1;
    for (int j=0;j<d; j++) xn[j]=xic[i].x[j];
        simplex_find(ind,xn);
    delaunay_update(tri,xi,xn,ind);
    }
    /* Calculate the cost function at initial points */
    for (int i=0; i<d+1;i++) test_func(xi[i+d+1]);
    /* Calculate the coefficients of the interpolation */
    for (int kk=0;kk<8;kk++)
    {
    inter_coeff(xi,v,Ni-d-1);
    /* perform the search in all simplices */
    float p1;
    float p_min=xi[d+2].y-y_0;
    method=0;
    for (int ii=0;ii<N_tri;ii++)
        if (tri[ii].ill==0 && tri[ii].stat==1)
        {
            
            initial_simplex(tri[ii]);
            simplex_search(tri[ii].xm,tri[ii].pm,tri[ii].R2, tri[ii].xc, tri[ii].method);
            if (tri[ii].method>1) { for (int j; j<d; j++) xn[j]=tri[ii].xm[j]; ind=ii;break;}
            else
                if (tri[ii].pm<p_min && tri[ii].method==method)
                {for (int j=0; j<d; j++) xn[j]=tri[ii].xm[j]; ind=ii; p_min=tri[ii].pm;}
                if (tri[ii].method>method)
                {for (int j=0; j<d; j++) xn[j]=tri[ii].xm[j]; ind=ii; p_min=tri[ii].pm; method++;}
            
        }
        /*
        printf("iteration %i\n",kk);
        for (int i=0;i<d;i++) printf("%e\n",xn[i]);
        printf("index %i\n",ind);
        printf("method %i\n",method);
        */
        
    /* Reflect the point on the boundary if it is not in righ simplices */
    if (method<2)
    { int check=1;
      insimplex_check(ind,xn,check);
        if (check==0)
        {
            /*reflect the point on the boundary */
            float dis=0;
            simplex_find(ind,xn);
            if (tri[ind].ill==1)
            {
            for(int i=0; i<d;i++) dis+=pow(xn[i],2); dis=sqrt(dis);
            for(int i=0; i<d;i++) xn[i]/=dis;
            }
            /* Check if the new point is in circume sphere of the original simplex */
            dis=0 ;for(int i=0; i<d;i++) dis+=pow((xn[i]-tri[ind].xc[i]),2);
            if (dis>tri[ind].R2)  simplex_find(ind,xn);
        }
    }
        
         printf("iteration %i\n",kk);
         for (int i=0;i<d;i++) printf("%e\n",xn[i]);
         printf("index %i\n",ind);
         printf("method %i\n",method);
        
    /* Update the Delaunay trinagulation */
    delaunay_update(tri,xi,xn,ind);
        test_func(xi[Ni-1]);
    }
    
    
}

/* Functions definitions */
/* Calculation of the circume center and circume radius */
void  circume(point* xir, float* xc, float &R2)
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
/* Calculate of the i-the neighber of the simplex ind. j is the negihber */
void  neighber(int ind, int i, int &j)
{
if (i!=0)
{
    std::set<int> s1 (xi[tri[ind].p[0]].simplex,xi[tri[ind].p[0]].simplex+xi[tri[ind].p[0]].simplex_length);
    for (int ii=1;ii<d+1;ii++)
     if (ii!=i)
     {
         std::set<int> s3;
         std::set<int> s2 (xi[tri[ind].p[ii]].simplex,xi[tri[ind].p[ii]].simplex+xi[tri[ind].p[ii]].simplex_length);
         std::set_intersection(s1.begin(),s1.end(),s2.begin(),s2.end(),std::inserter(s3, s3.end()));
         s1=s3;
     }
    for (std::set<int>::iterator ii=s1.begin(); ii!=s1.end(); ii++)
    { int p=*ii; if (p!=ind && p<N_tri) if (tri[p].stat>0) {j=p; break; } }

}
 else
{
    std::set<int> s1 (xi[tri[ind].p[d]].simplex,xi[tri[ind].p[d]].simplex+xi[tri[ind].p[d]].simplex_length);
    for (int ii=1;ii<d;ii++)
 {
     std::set<int> s3;
     std::set<int> s2 (xi[tri[ind].p[ii]].simplex,xi[tri[ind].p[ii]].simplex+xi[tri[ind].p[ii]].simplex_length);
     std:: set_intersection(s1.begin(),s1.end(),s2.begin(),s2.end(),std::inserter(s3, s3.end()));
     s1=s3;
 }
    for (std::set<int>::iterator ii=s1.begin(); ii!=s1.end(); ii++)
    { int p=*ii; if (p!=ind && p<N_tri) if (tri[p].stat>0) {j=p; break; } }

}
}


/* Update the delaunay triangulation */
void  delaunay_update(simplex* tri, point* xi, float* xn, int ind)
{
    for (int i=0;i<d;i++)  xi[Ni].x[i]=xn[i];
    int k_read=0, k_all=1;
    killed[k_read]=ind;
    N_new=0;
    while (k_read<k_all)
    {
        ind=killed[k_read];
        k_read++;
        tri[ind].stat=2;
        for (int i=0;i<d+1;i++)
        {
            // find the i-the neighber of this simplex //
            int k=-1; neighber(ind,i,k);
            if (k==-1)
            {
                    // constrcut the simplex
                    // calculate its point
                    for (int l=0;l<d+1;l++) tri[N_tri+N_new].p[l]=tri[ind].p[l]; tri[N_tri+N_new].p[i]=Ni;
                    // calculate the circume center
                    point xir[d+1];
                    for (int l=0;l<d+1;l++) for (int t=0;t<d;t++) xir[l].x[t]=xi[tri[N_tri+N_new].p[l]].x[t];
                    circume(xir,tri[N_tri+N_new].xc,tri[N_tri+N_new].R2);
                    // calculate stat and ill
                    tri[N_tri+N_new].stat=1; tri[N_tri+N_new].ill=0;
                   for (int l=0;l<d+1;l++) if (tri[N_tri+N_new].p[l]<d+1) {tri[N_tri+N_new].ill=1; break;}

                    //update the data points
                    for (int ii=0;ii<d+1;ii++) if (ii!=i)
                    {xi[tri[ind].p[ii]].simplex[xi[tri[ind].p[ii]].simplex_length]=N_tri+N_new;
                    xi[tri[ind].p[ii]].simplex_length++;}
                    // change the index
                    N_new++;
            }
            else
                if (tri[k].stat==1)
                {
                    float p=0;
                    for (int l=0;l<d;l++) p+=pow((tri[k].xc[l]-xn[l]),2);
                    if (p<tri[k].R2)
                    {tri[k].stat=2; killed[k_all]=k; k_all++;}
                    else
                    {
                        // constrcut the simplex
                        // calculate its point
                        for (int l=0;l<d+1;l++) tri[N_tri+N_new].p[l]=tri[ind].p[l]; tri[N_tri+N_new].p[i]=Ni;
                        // calculate the circume center
                        point xir[d+1];
                        for (int l=0;l<d+1;l++) for (int t=0;t<d;t++) xir[l].x[t]=xi[tri[N_tri+N_new].p[l]].x[t]; circume(xir,tri[N_tri+N_new].xc,tri[N_tri+N_new].R2);
                        // calculate stat and ill
                        tri[N_tri+N_new].stat=1; tri[N_tri+N_new].ill=0;
                        for (int l=0;l<d+1;l++) if (tri[N_tri+N_new].p[l]<d+1) {tri[N_tri+N_new].ill=1; break;}
                        //update the data points
                        for (int ii=0;ii<d+1;ii++) if (ii!=i)
                        {xi[tri[ind].p[ii]].simplex[xi[tri[ind].p[ii]].simplex_length]=N_tri+N_new;
                            xi[tri[ind].p[ii]].simplex_length++;}
                        // change the index
                        N_new++;
                        
                    }
                    
                }
        }
    }
    for (int i=0;i<k_all;i++) tri[killed[i]].stat=0;
    xi[Ni].simplex_length=N_new;
    for (int i=0;i<N_new;i++)  {tri[N_tri+i].stat=1; xi[Ni].simplex[i]=N_tri+i;}
    N_tri+=N_new; Ni++;
}

/* Find the simplex which includes the point */
void  simplex_find(int &ind, float* xn)
{
float w[d+1], MT[(d+1)*(d+1)];
while (1)
{
for( int t=0; t<d; t++) w[t]=xn[t]; w[d]=1;
for( int t=0; t<d+1; t++) {for (int l=0; l<d; l++) MT[t*(d+1)+l]=xi[tri[ind].p[t]].x[l]; MT[t*(d+1)+d]=1;}
int c1=d+1, c2=1, pivot[d+1], ok, change;
sgesv_(&c1, &c2, MT, &c1, pivot, w, &c1, &ok);
change=1;
    for( int t=0; t<d+1; t++) if (w[t]<-1e-3) {int k=-1; neighber(ind,t,k); ind=k; change=0; break;}
if (change) break;
}

}
/* Check if the points is in a simplex*/
void  insimplex_check(int &ind, float* xn, int &check)
{
    float w[d+1], MT[(d+1)*(d+1)];
        for( int t=0; t<d; t++) w[t]=xn[t]; w[d]=1;
        for( int t=0; t<d+1; t++) {for (int l=0; l<d; l++) MT[t*(d+1)+l]=xi[tri[ind].p[t]].x[l]; MT[t*(d+1)+d]=1;}
        int c1=d+1, c2=1, pivot[d+1], ok, change;
        sgesv_(&c1, &c2, MT, &c1, pivot, w, &c1, &ok);
    for( int t=0; t<d+1; t++) if (w[t]<-1e-3) {check=0; break;}
}
/* Calculate the initial point of the simplex */
void initial_simplex( simplex &tri)
{
    float w[d+1], MT[(d+1)*(d+1)],w1;
    for( int t=0; t<d; t++) w[t]=tri.xc[t]; w[d]=1;
    for( int t=0; t<d+1; t++) {for (int l=0; l<d; l++) MT[t*(d+1)+l]=xi[tri.p[t]].x[l]; MT[t*(d+1)+d]=1;}
    int c1=d+1, c2=1, pivot[d+1], ok, change;
    sgesv_(&c1, &c2, MT, &c1, pivot, w, &c1, &ok);
    w1=0;  for( int t=0; t<d+1; t++) if (w[t]<w1) w1=w[t];
    float alpha=1.0/(1-(d+1)*w1);
    for(int i=0; i<d; i++)
    {
    tri.xm[i]=tri.xc[i]*alpha;
    for (int k=0; k<d+1;k++) tri.xm[i]+=(1-alpha)*xi[tri.p[k]].x[i]/(d+1);
    }
    
}

