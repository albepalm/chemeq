#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "linalg.h"

void LINALG_BACKSUBS(LINALG_VEC *x, LINALG_MAT *U, LINALG_VEC *b)
{
	uint64_t n = (*U).rows;
	if(n == (*U).cols && n==(*b).len)
	{
		double sum;
		uint64_t ii;
		for(uint64_t i = n; i > 0; --i)
		{
			sum = 0;
			ii = i-1;
			for(uint64_t j = ii+1; j < n; ++j)
			{
				sum += ((*U).pt[ii*n+j])*((*x).pt[j]);
			}
			if((*U).pt[ii*n+ii]!=0)
			{
				(*x).pt[ii] = ((*b).pt[ii]-sum)/(*U).pt[ii*n+ii];
			}
			else
			{
				printf("Division by 0: singular matrix or not triangular!\n");
			}
		}
	}
	else
	{
		printf("The backsubsitution can not be performed!\n");
	}
	return;
}
