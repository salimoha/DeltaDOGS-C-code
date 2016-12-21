/* Implementing the delaunay triangulation */
#include <stdio.h>
#include <math.h>
#include <clapack.h>
#include <set>

# define d 2	/* dimension of matrix */
# define n 3
# define K 1
# define y0 -279
struct simplex {
    int p[d+1];    /* The pointers of vertices of the simplex*/
    int ad[d+1];  /* The pointers of neighbers of the simplex*/
    int stat; /* stat=0: killed simplex */
    int ill; /* the illposed simplex */
    float xc[d]; /* circume center of the simplex */
    float xm[d]; /* the best point of the simplex */
    float R2; /* square of circumeradius of the simplex */
};
// defining the structure for simplex
struct point {
    float x[d];/* The pointers of vertices of the simplex*/
    float y;
};

int Ni, N_tri, N_new;
simplex tri[10];
point xi[10];

void circume(point* xir, float* xc, float &R2);
void  face_find(simplex* tri, point* xi,float* xn, int ind);
void  new_simplex(simplex* tri, point* xi, int Ni, int N_new);

int main()
{
    // defining the structure for simplex
        // Initialization of the points
    Ni=5; N_tri=4;
    // initilaizatin of the points
    xi[0].x[0]=0; xi[0].x[1]=0;
    xi[1].x[0]=1; xi[1].x[1]=0;
    xi[2].x[0]=0; xi[2].x[1]=1;
    xi[3].x[0]=1; xi[3].x[1]=1;
    xi[4].x[0]=0.5; xi[4].x[1]=0.5;
    // initilizatiation of the  triangulation
    tri[0].p[0]=0;    tri[0].p[1]=1;   tri[0].p[2]=4;
    tri[1].p[0]=0;    tri[1].p[1]=2;   tri[1].p[2]=4;
    tri[2].p[0]=3;    tri[2].p[1]=1;   tri[2].p[2]=4;
    tri[3].p[0]=3;    tri[3].p[1]=2;   tri[3].p[2]=4;
    // initilizatiation of the adjancy matrix
    tri[0].ad[0]=2;    tri[0].ad[1]=1;   tri[0].ad[2]=-1;
    tri[1].ad[0]=3;    tri[1].ad[1]=0;   tri[1].ad[2]=-1;
    tri[2].ad[0]=0;    tri[2].ad[1]=3;   tri[2].ad[2]=-1;
    tri[3].ad[0]=1;    tri[3].ad[1]=2;   tri[3].ad[2]=-1;
    for (int k=0; k<4; k++)
    {
        
        point xir[d+1];
        for (int j=0;j<d+1;j++)
            for (int i=0; i<d; i++)
            xir[j].x[i]=xi[tri[k].p[j]].x[i];
    circume(xir,tri[k].xc,tri[k].R2);
    for(int i=0;i<d;i++) tri[k].xm[i]=tri[k].xc[i];
        tri[k].stat=1; tri[k].ill=0;
    }
    // the new point
    float xn[d]={0.25, 0.75}; int ind=1;
    for (int t=0;t<d;t++)xi[Ni].x[t]=xn[t];
    face_find(tri,xi,xn,ind);
    new_simplex(tri,xi,Ni,N_new);
}
/* Calculation of the circume center and circume sphere */
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

/* find all faces of the new simplices */
void  face_find(simplex* tri, point* xi, float* xn, int ind)
{
int killed[10], k_read=0, k_all=1;
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
        {
            float w[d+1], MT[(d+1)*(d+1)];
            for( int t=0; t<d; t++) w[t]=xn[t]; w[d]=1;
            for( int t=0; t<d+1; t++) {for (int l=0; l<d; l++) MT[t*(d+1)+l]=xi[tri[ind].p[t]].x[l]; MT[t*(d+1)+d]=1;}
            int c1=d+1, c2=1, pivot[d+1], ok;
            sgesv_(&c1, &c2, MT, &c1, pivot, w, &c1, &ok);
            for( int t=0; t<d+1; t++) printf("%e\n",w[t]);
            if (w[i]>0)
            {
             // constrcut the simplex
                // calculate its point
                for (int l=0;l<d+1;l++) tri[N_tri+N_new].p[l]=tri[ind].p[l]; tri[N_tri+N_new].p[i]=Ni;
                // calculate known neighber
                tri[N_tri+N_new].ad[i]=k;
                // calculate the circume center
                    point xir[d+1];
                    for (int l=0;l<d+1;l++) for (int t=0;t<d;t++) xir[l].x[t]=xi[tri[N_tri+N_new].p[l]].x[t];
                    circume(xir,tri[N_tri+N_new].xc,tri[N_tri+N_new].R2);
                // calculate stat and ill
                tri[N_tri+N_new].stat=1; tri[N_tri+N_new].ill=0;

                // change the index
                    N_new++;
            }
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
}
/* Find the  adjancey relationship for new simplex*/
void  new_simplex(simplex* tri, point* xi, int Ni, int N_new)
{
    for (int k=0; k<N_new-1; k++)
        for (int i=0; i<d+1; i++)
            if(tri[N_tri+k].p[i]!=Ni)
            {
                int c1[d+1]; for(int j=0; j<d+1; j++) c1[j]=tri[N_tri+k].p[j]; c1[i]=-1;
                for (int k1=k+1; k1<N_new; k1++)
                    for (int i1=0;i1<d+1; i1++)
                    { int c2[d+1]; for(int j=0; j<d+1; j++) c2[j]=tri[N_tri+k1].p[j]; c2[i1]=-1;
                        std::set<int> s1 (c1,c1+(d+1));
                        std::set<int> s2(c2,c2+(d+1));
                        if (s1==s2)
                        {tri[N_tri+k].ad[i]=N_tri+k1; tri[N_tri+k1].ad[i1]=N_tri+k;}
                    }
            }
}




