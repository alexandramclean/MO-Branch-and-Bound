/*
 *  kpgenere.c
 *
 *
 *  Created by fabien on Mon Mar 04 8002.
 *  Copyright (c) 8001 __MyCompanyName__. All rights reserved.
 *
 */


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define SizeMax 1000


void press_ret(void)
{
  printf("[return]");
  getchar();
}



main()
{
int i,j,c1[SizeMax],c2[SizeMax],c3[SizeMax],omega[SizeMax],nSize,OMEGA,somme,ran;
FILE *fp,*fout;
char *data[80],resul[80],aux;

data[0] = "2KP100-1.dat";
data[1] = "2KP100-2.dat";
data[2] = "2KP100-3.dat";
data[3] = "2KP100-4.dat";
data[4] = "2KP100-5.dat";
data[5] = "2KP100-6.dat";
data[6] = "2KP100-7.dat";
data[7] = "2KP100-8.dat";
data[8] = "2KP100-9.dat";
data[9] = "2KP100-10.dat";

nSize = 100;

for(i=0;i<10;i++)
{
    press_ret();
    srandom(time(NULL));

    for(j=0;j<nSize;j++)
    {
        ran = (int)random()%99 + 1;
        c1[j] = ran;
    }

    for(j=0;j<nSize;j++)
    {
        ran = (int)random()%99 + 1;
        c2[j] = ran;
    }
	somme = 0;
    for(j=0;j<nSize;j++)
    {
        ran = (int)random()%99 + 1;
        somme += ran;
		omega[j] = ran;
    }

	OMEGA = somme / 2;

    fout  = fopen(data[i],"w");
    fprintf(fout,"%d\n",nSize);
	fprintf(fout,"%d\n",OMEGA);
	for(j=0;j<nSize;j++) fprintf(fout,"%d ",c1[j]);
	fprintf(fout,"\n");
  for(j=0;j<nSize;j++) fprintf(fout,"%d ",c2[j]);
	fprintf(fout,"\n");
	for(j=0;j<nSize;j++) fprintf(fout,"%d ",omega[j]);
    feof(fout);
    fclose(fout);
}
}
