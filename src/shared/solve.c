#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "linalg.h"

bool LINALG_SOLVE(LINALG_VEC *x, LINALG_MAT *A, LINALG_VEC *b, uint64_t d)
{
	uint64_t n = (*A).rows;
	bool SOLVER_FAIL = false;
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
	LINALG_MAT eA;
	LINALG_VEC eb;
	LINALG_INIT_MAT(&eA);
	LINALG_ZERO_MAT(&eA,d,d);
	LINALG_INIT_VEC(&eb);
	LINALG_ZERO_VEC(&eb,d);
	SOLVER_FAIL = LINALG_GAUSSELIM(A,b,&eA,&eb);
	if(!SOLVER_FAIL)
	{
		LINALG_BACKSUBS(x,&eA,&eb);
	}
	LINALG_DEALLOC_MAT(&eA);
	LINALG_DEALLOC_VEC(&eb);
	return SOLVER_FAIL;
}

void LINALG_SOLVEWITHNORM(LINALG_VEC *x, LINALG_MAT *A, LINALG_VEC *b,
		uint64_t d)
{
	LINALG_MAT tmpm;
	LINALG_VEC tmpv;
	LINALG_INIT_MAT(&tmpm);
	LINALG_ZERO_MAT(&tmpm,d,1);
	LINALG_INIT_VEC(&tmpv);
	LINALG_ZERO_VEC(&tmpv,d);

	LINALG_SOLVE(x,A,b,d);

	LINALG_V2MC(&tmpm,x);
	LINALG_MUL_MAT(&tmpm,A,&tmpm);
	LINALG_MC2V(&tmpv,&tmpm);
	LINALG_NEG_VEC(&tmpv,&tmpv);
	LINALG_SUM_VEC(&tmpv,b,&tmpv);
	double res = -1;
	res = LINALG_NORMMAX_VEC(&tmpv);
	printf("res: %7.5e\n",res);
	LINALG_DEALLOC_MAT(&tmpm);
	LINALG_DEALLOC_VEC(&tmpv);
	return;
}
