#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include "tp.h"

void TP(double *T, char *name, double *cp, double *h, double *s, double *g)
{
	FILE *fp = NULL;
	char tmp[200] = "";
	strcat(tmp,"./thermo_glenn/");
	strcat(tmp,name);
	strcat(tmp,".csv");
	fp = fopen(tmp,"r");
	double t1,t2,a1,a2,a3,a4,a5,a6,a7,b1,b2;
	//&t1,&t2,&a1,&a2,&a3,&a4,&a5,&a6,&a7,&b1,&b2

	while(fscanf(fp,"%lf,%lf,%le,%le,%le,%le,%le,%le,%le,%le,%le",
				&t1,&t2,&a1,&a2,&a3,&a4,&a5,&a6,&a7,&b1,&b2) != EOF)
	{
		//printf("%lf,%lf,%le,%le,%le,%le,%le,%le,%le,%le,%le\n",
		//		t1,t2,a1,a2,a3,a4,a5,a6,a7,b1,b2);
		if((*T)>t1 && (*T)<=t2)
		{
			*cp = a1*pow(*T,-2) +
				a2*pow(*T,-1) + 
				a3*1 +
				a4*(*T) +
				a5*pow(*T,2) + 
				a6*pow(*T,3) +
				a7*pow(*T,4);

			*h = - a1*pow(*T,-2) +
				a2*log(*T)/(*T) +
				a3 +
				a4*(*T)/2 + 
				a5*pow(*T,2)/3 +
				a6*pow(*T,3)/4 +
				a7*pow(*T,4)/5 +
				b1/(*T);

			*s = -a1*pow(*T,-2)/2 -
				a2*pow(*T,-1) +
				a3*log(*T) +
				a4*(*T) +
				a5*pow(*T,2)/2 +
				a6*pow(*T,3)/3 +
				a7*pow(*T,4)/4 +
				b2;

			*g = *h - *s;
		}
	}
	fclose(fp);
	return;
}

double KP(double *T,int nprod, int *sprod, char **prod,
		int nreag, int *sreag, char **reag)
{
	//printf("%s, %s\n", prod[0], reag[0]);
	double CP, H, S, G, GP=0, GR=0, GT=0;
	int sp = 0, sr = 0;
	for(int i = 0; i<nprod; ++i)
	{
		TP(T,prod[i],&CP,&H,&S,&G);
		GP += (G*sprod[i]);
		sp += sprod[i];
	}

	for(int i = 0; i<nreag; ++i)
	{
		TP(T,reag[i],&CP,&H,&S,&G);
		GR += (G*sreag[i]);
		sr += sreag[i];
	}
	
	GT = GP - GR;
	double pref = 1e5;
	double b = (double)(sp - sr); 
	double A = pow(pref,b);
	return A*exp(-GT);
}
