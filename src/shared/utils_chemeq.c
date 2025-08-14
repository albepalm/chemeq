#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stdint.h>
#include <math.h>
#include "utils_chemeq.h"

void CHEMEQ_SPLASH()
{
	CHEMEQ_CLS();
	printf("=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=\n");
	fflush(stdout);
	printf("-=-=-=- CHEMICAL  EQUILIBRIUM -=-=-=-\n");
	fflush(stdout);
	printf("= for homogeneous gasseous mixtures =\n");
	fflush(stdout);
	printf("-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-\n");
	fflush(stdout);
	return;
}

void CHEMEQ_CLS()
{
	//printf("\e[1;1H\e[2J");
	if(system("clear"))
	{
		fprintf(stderr,"shell cmd error\nabort\n");
		exit(1);
	}
	return;
}

int CHEMEQ_USR_INT()
{
	char line[256];
	int i;
	if (fgets(line, sizeof(line), stdin)) {
	    if (1 != sscanf(line, "%d", &i)) {
			fprintf(stderr,"INPUT ERROR! Retry\n");
			i = CHEMEQ_USR_INT();
	    }
	}
	return i;
}

double CHEMEQ_USR_DBL()
{
	char line[256];
	double i;
	if (fgets(line, sizeof(line), stdin)) {
	    if (1 != sscanf(line, "%lf", &i)) {
			fprintf(stderr,"INPUT ERROR! Retry\n");
			i = CHEMEQ_USR_DBL();
	    }
	}
	return i;
}

int CHEMEQ_MAIN()
{
	CHEMEQ_SPLASH();
	printf("\n\n 1\tPure Nitrogen\n 2\tPure Oxygen\n 3\tOxygen - Nitrogen Mixture\n");
	printf("Choose Elements: ");
	
	int i = 0;
	i = CHEMEQ_USR_INT();
	while(i != 1 && i != 2 && i != 3)
	{
		printf("No valid option chosen!\nChoose Elements: ");
		i =  CHEMEQ_USR_INT();
	}
	return i;
}

int CHEMEQ_SPEC()
{
	CHEMEQ_SPLASH();
	printf("\n\n 1\tNeutral Species Only\n 2\tNeutral Species and Major Iones\n 3\tNeutral Species and All Iones\n");
	printf("Choose Species: ");

	int i = 0;
	i =  CHEMEQ_USR_INT();
	while(i != 1 && i != 2 && i != 3)
	{
		printf("No valid option chosen!\nChoose Species: ");
		i =  CHEMEQ_USR_INT();
	}
	return i;
}

double CHEMEQ_O2N_RAT()
{
	double ret = -2;
	printf("Oxygen to Nitrogen ratio (>0): ");
	ret = CHEMEQ_USR_DBL();
	while(ret <= 0)
	{
		printf("Invalid Ratio!\nOxygen to Nitrogen ratio (>0): ");
		ret =  CHEMEQ_USR_INT();
	}
	return ret;
}

void CHEMEQ_PRES(double *min, double *max, int *n)
{
	CHEMEQ_SPLASH();
	printf("Minimum pressure (Pascal): ");
	(*min) = CHEMEQ_USR_DBL();
	while((*min)<=0)
	{
		printf("Invalid value!\nMinimum pressure (Pascal): ");
		(*min) = CHEMEQ_USR_DBL();
	}
	printf("Maximum pressure (Pascal): ");
	(*max) = CHEMEQ_USR_DBL();
	while((*max)<(*min))
	{
		printf("Invalid value!\nMaximum pressure (Pascal): ");
		(*max) = CHEMEQ_USR_DBL();
	}
	if(!((*min)==(*max)))
	{
		printf("Number of values: ");
		(*n) = CHEMEQ_USR_INT();
		while((*n)<=1)
		{
			printf("Invalid value!\nNumber of values: ");
			(*n) = CHEMEQ_USR_INT();
		}
	}
	else
	{
		(*n) = 1;
	}
	return;
}

void CHEMEQ_TEMP(double *min, double *max, int *n)
{
	CHEMEQ_SPLASH();
	printf("Minimum temperature (Kelvin): ");
	(*min) = CHEMEQ_USR_DBL();
	while((*min)<=0)
	{
		printf("Invalid value!\nMinimum temperature (Kelvin): ");
		(*min) = CHEMEQ_USR_DBL();
	}
	printf("Maximum temperature (Kelvin): ");
	(*max) = CHEMEQ_USR_DBL();
	while((*max)<(*min))
	{
		printf("Invalid value!\nMaximum temperature (Kelvin): ");
		(*max) = CHEMEQ_USR_DBL();
	}
	if(!((*min)==(*max)))
	{
		printf("Number of values: ");
		(*n) = CHEMEQ_USR_INT();
		while((*n)<=1)
		{
			printf("Invalid value!\nNumber of values: ");
			(*n) = CHEMEQ_USR_INT();
		}
	}
	else
	{
		(*n) = 1;
	}
	return;
}

void CHEMEQ_LINSPACE(double **lin, double min, double max, uint64_t n)
{
	(*lin) = (double *) malloc(n*sizeof(double));
	if (n == 1)
	{
		(*lin)[0] = min;
	}
	else
	{
		for(uint64_t i = 0; i < n; ++i)
		{
			(*lin)[i] = min + i*(max-min)/(n-1);
		}
	}
	return;
}

void CHEMEQ_LOG10SPACE(double **lin, double min, double max, uint64_t n)
{
	(*lin) = (double *) malloc(n*sizeof(double));
	if (n == 1)
	{
		(*lin)[0] = min;
	}
	else
	{
		for(uint64_t i = 0; i < n; ++i)
		{
			(*lin)[i] = min*pow((max/min), (double)i/(n-1));
		}
	}
	return;
}
