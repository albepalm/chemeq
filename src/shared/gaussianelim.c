#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "linalg.h"

bool LINALG_GAUSSELIM(LINALG_MAT *A, LINALG_VEC *b,
		LINALG_MAT *eA, LINALG_VEC *eb)
{
	bool SOLVER_FAIL = false;
	uint64_t n = (*A).rows;
	if(n != (*A).cols)
	{
		fprintf(stderr,"NON SQUARE MATRIX!\n");
		exit(1);
	}
	//if(n != LINALG_RANK(A))
	//{
	//	fprintf(stderr,"SINGULAR MATRIX!\n");
	//	exit(1);
	//}
	LINALG_COPY_MAT(eA,A);
	LINALG_COPY_VEC(eb,b);
	for(uint64_t i = 0; i < n; ++i)
	{
		//if((*A).pt[i*(n+1)] != 0)
		//{
			double pivmax = 0;
			uint64_t pospivmax = 0;
			for(uint64_t j = i; j < n; ++j)
			{
				if(fabs((*eA).pt[j*n+i]) > fabs(pivmax))
				{
					pivmax = (*eA).pt[j*n+i];
					pospivmax = j;
				}
			}
			double tmp = 0;
			for(uint64_t j = i; j < n; ++j)
			{
				tmp = (*eA).pt[i*n+j];
				(*eA).pt[i*n+j] = (*eA).pt[pospivmax*n+j];
				(*eA).pt[pospivmax*n+j] = tmp;
			}
			tmp = (*eb).pt[i];
			(*eb).pt[i] = (*eb).pt[pospivmax];
			(*eb).pt[pospivmax] = tmp;
			if(pivmax == 0.0)
			{
				printf("WARNING! Matrix is singular.\n");
				fflush(stdout);
				SOLVER_FAIL = true;
				//getchar();
				return SOLVER_FAIL;
			}
			else
			{
				for(uint64_t j = i+1; j < n; ++j)
				{
					tmp = (*eA).pt[j*n+i]/pivmax;
					for(uint64_t k = i; k < n; ++k)
					{
						(*eA).pt[j*n+k] -= tmp*(*eA).pt[i*n+k];
					}
					(*eb).pt[j] -= tmp*(*eb).pt[i];
				}
			}
			//PRINTM(*eA);
			//PRINTV(*eb);
			//getchar();
		//}
	}
	return SOLVER_FAIL;
}
