#include <stdio.h>
void cwrite_(float *lvisit,int *visit, int *k)
{
   FILE *ff;
   int i;
   ff=fopen("slogiclisting.tmp","a");
   fprintf(ff,"%f %f ",lvisit[0],lvisit[1]);
   for(i=1;i<(*k);i++)fprintf(ff,"%d ",visit[i]);
   fprintf(ff,"\n");
   fclose(ff);
}
void ciwrite_(void)
{
   FILE *ff;
   ff=fopen("slogiclisting.tmp","w");
   fclose(ff);
}


