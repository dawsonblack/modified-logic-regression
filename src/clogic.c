#include <math.h>
#include <stdio.h>
#include "R.h"
#include "clogic.h"
#define Salloc(n, t)  (t *)R_alloc((long)(n), (int)sizeof(t))

void F77_NAME(slogreg)(int *, int *, int *, int *, 
   float *, float *, int *, int *, float *,
   float *, int *, int *, float *, float *,
   int *, int *, int *, int *,
   float *, int *, int*, 
   double *, double *, float *, int *, int *, int *, int *, int *, int *);

static int *ilvector(int l);
static float *flvector(int l);
static double *dlvector(int l);

void clogreg(int *n1,int *n2,int *nsep,int *intpars,float *rpars,float *seps,
	int *dcph,int* orders,float *resp,float *weight,
  int *datri,int *iotrees,float *iocoef,float *ioscores,int *rd4)
{
   int *ntrx,*nknx,*storage,*storage2,*storage3,bmax,jmax,tmax;
   float *storage4,*wur1;
   double *wud1,*wud2;
   int length,ip4,*wui1,*wui2,*wui3,i;

   bmax=55;
   ip4 = 2*intpars[3]+1;
   ntrx = ilvector(intpars[5]);
   nknx = ilvector(ip4);
   length = 2*intpars[5]*ip4*n1[0];
   storage = ilvector(length);
   length = 7*intpars[5]*(ip4+1)*n2[0]*4;
   storage3 = ilvector(length);
   storage4 = flvector(length);
   storage2 = ilvector(n1[0]*n2[0]); 
   wud1 = dlvector(n1[0]*(bmax*bmax+bmax+6));
   wud2 = dlvector(16384*2);
   wur1 = flvector(n1[0]*(2*bmax+8));
   wui1 = ilvector(n1[0]*(3*bmax+8));
   jmax = 2*intpars[3];
   if(jmax<2)jmax=2;
   wui2 = ilvector((n1[0]+2)*jmax);
   tmax = 1;
   for(i=0;i<intpars[5];i++)tmax=tmax*2;
   tmax = tmax+1;
   wui3 = ilvector(n1[0]*tmax);

   F77_CALL(slogreg)(n1,n2,nsep,intpars,
      rpars,seps,dcph,orders,resp,
      weight,datri,iotrees,iocoef,ioscores,
      ntrx,nknx,storage,storage3,
      storage4,storage2,rd4,wud1,wud2,
      wur1,wui1,&bmax,wui2,&jmax,wui3,&tmax);
   return;
}
/******************************************************************************/

static int *ilvector(int l)
/* allocate an int vector with subscript range v[0...(l-1)] */
{
   int *v,i;
   v=(int *)Salloc(l,int);
   for(i=0;i<l;i++)v[i]=0;
   return v;
}
/******************************************************************************/
static float *flvector(int l)
/* allocate a double vector with subscript range v[0...(l-1)] */
{
   float *v;
   int i;
   v=(float *)Salloc(l,float);
   for(i=0;i<l;i++) v[i]=0.;
   return v;
}
/******************************************************************************/
static double *dlvector(int l)
/* allocate a double vector with subscript range v[0...(l-1)] */
{
   double *v;
   int i;
   v=(double *)Salloc(l,double);
   for(i=0;i<l;i++) v[i]=0.;
   return v;
}
/******************************************************************************/
