#include <stdio.h>
void cwrite_(lvisit,visit,k)
float *lvisit;
int *visit,*k;
{
   FILE *ff;
   int i;
   ff=fopen("slogiclisting.tmp","a");
   fprintf(ff,"%f %f ",lvisit[0],lvisit[1]);
   for(i=0;i<(*k);i++)fprintf(ff,"%d ",visit[i]);
   fprintf(ff,"\n");
   fclose(ff);
}
void ciwrite_()
{
   FILE *ff;
   int i;
   ff=fopen("slogiclisting.tmp","w");
   fclose(ff);
}


