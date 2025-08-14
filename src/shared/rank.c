#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "linalg.h"

uint64_t LINALG_RANK(LINALG_MAT *A)
{
	uint64_t R = (*A).rows;
	uint64_t C = (*A).cols;
	uint64_t rank = C;
	LINALG_MAT B;
	LINALG_INIT_MAT(&B);
	LINALG_COPY_MAT(&B,A);
	//if(c>r) {LINALG_TRS_MAT(&B,&B);}

	for(uint64_t row = 0; row < rank; ++row)
	{
		if(B.pt[row*C+row] != 0)
		{
			for(uint64_t col = 0; col < R; ++col)
			{
				if(col != row)
				{
					double mult = (double) B.pt[col*C+row]/B.pt[row*C+row];
					for(uint64_t i = 0; i < rank; ++i)
					{
						B.pt[col*C+i] -= mult*B.pt[row*C+i];
					}
				}
			}
		}
		else
		{
			bool reduce = true;
			for(uint64_t i = row+1; i < R; ++i)
			{
				if(B.pt[i*C+row] != 0)
				{
					for(uint64_t ii = 0; ii < rank; ++ii)
					{
						double temp = B.pt[row*C+ii];
						B.pt[row*C+ii] = B.pt[i*C+ii];
						B.pt[i*C+ii] = temp;
					}
					reduce = false;
					break;
				}
			}
			if(reduce)
			{
				rank--;
				for(uint64_t i = 0; i < R; ++i)
				{
					B.pt[i*C+row] = B.pt[i*C+rank];
				}
			}
			row--;
		}
	}
	LINALG_DEALLOC_MAT(&B);
	return rank;
}
