#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "linalg.h"

void LINALG_INIT_MAT(LINALG_MAT *m)
{
	(*m).rows = 0;
	(*m).cols = 0;
	(*m).pt = NULL;
	return;
}

void LINALG_ASSIGN_MAT(LINALG_MAT *m, double *usrm, uint64_t r, uint64_t c)
{
	if((*m).rows == r && (*m).cols == c && (*m).pt != NULL)
	{
		for(uint64_t i = 0; i < r; ++i)
		{
			for(uint64_t j = 0; j < c; ++j)
			{
				(*m).pt[i*c+j] = usrm[i*c+j];
			}
		}
	} else
		ASSIGN_ERROR;
	return;
}

bool LINALG_MALLOC_MAT(LINALG_MAT *m, uint64_t r, uint64_t c)
{
	bool ret;
	(*m).rows = r;
	(*m).cols = c;
	(*m).pt = (double *) malloc(r*c*sizeof(double));
	if ((*m).pt == NULL) {
		MALLOC_ERROR;
		ret = true;
	}
	else {
		ret = false;
	}
	return ret;
}

bool LINALG_REALLOC_MAT(LINALG_MAT *m, uint64_t r, uint64_t c)
{
	bool ret;
	(*m).rows = r;
	(*m).cols = c;
	(*m).pt = (double *) realloc((*m).pt,r*c*sizeof(double));
	if ((*m).pt == NULL) {
		MALLOC_ERROR;
		ret = true;
	}
	else {
		ret = false;
	}
	return ret;
}

void LINALG_DEF_MAT(LINALG_MAT *m, double *usrm, uint64_t r, uint64_t c)
{
	if((*m).rows == 0 && (*m).cols == 0 && (*m).pt == NULL) {
		if(!LINALG_MALLOC_MAT(m,r,c))
		{
			LINALG_ASSIGN_MAT(m, usrm, r, c);
			//printf("New matrix defined!\n");
		}
		else
			MALLOC_ERROR;
	} else
		INIT_ERROR;
	return;
}

void LINALG_REDEF_MAT(LINALG_MAT *m, double *usrm, uint64_t r, uint64_t c)
{
	if (!LINALG_REALLOC_MAT(m,r,c)) {
		LINALG_ASSIGN_MAT(m,usrm,r,c);
		//printf("Matrix redefined!\n");
	}
	else
		MALLOC_ERROR;
	return;
}

void LINALG_ZERO_MAT(LINALG_MAT *m, uint64_t r, uint64_t c)
{
	if((*m).pt == NULL && (*m).rows == 0 && (*m).cols == 0)
	{
		LINALG_MALLOC_MAT(m,r,c);
	} else
	{
		LINALG_REALLOC_MAT(m,r,c);
	}
	for(uint64_t i = 0; i < r*c; ++i)
		(*m).pt[i] = 0;
	return;
}

void LINALG_ID_MAT(LINALG_MAT *m, uint64_t r)
{
	if((*m).pt == NULL && (*m).rows == 0 && (*m).cols == 0)
	{
		LINALG_MALLOC_MAT(m,r,r);
	} else
	{
		LINALG_REALLOC_MAT(m,r,r);
	}
	for(uint64_t i = 0; i < r; ++i)
		(*m).pt[i*(r+1)] = 1;
	return;
}

void LINALG_ONES_MAT(LINALG_MAT *m, uint64_t r, uint64_t c)
{
	if((*m).pt == NULL && (*m).rows == 0 && (*m).cols == 0)
	{
		LINALG_MALLOC_MAT(m,r,c);
	} else
	{
		LINALG_REALLOC_MAT(m,r,c);
	}
	for(uint64_t i = 0; i < r*c; ++i)
		(*m).pt[i] = 1;
	return;
}

void LINALG_RAND_MAT(LINALG_MAT *m, uint64_t r, uint64_t c)
{
	srand(time(NULL));
	if((*m).pt == NULL && (*m).rows == 0 && (*m).cols == 0)
	{
		LINALG_MALLOC_MAT(m,r,c);
	} else
	{
		LINALG_REALLOC_MAT(m,r,c);
	}
	for(uint64_t i = 0; i < r*c; ++i)
		(*m).pt[i] = ((double) rand())/((double) RAND_MAX);
	return;
}

void LINALG_DEALLOC_MAT(LINALG_MAT *m)
{
	free((*m).pt);
	LINALG_INIT_MAT(m);
}

void LINALG_PRINT_MAT(LINALG_MAT *m, char *name)
{
	printf("%s = [\n", name);
	for(uint64_t i = 0; i < (*m).rows; ++i)
	{
		for(uint64_t j = 0; j < (*m).cols; ++j)
			printf("\t%7.5e", (*m).pt[i*(*m).cols+j]);
			//printf("\t%7.5f", (*m).pt[i*(*m).cols+j]);
		printf("\n");
	}
	printf("]\tmatrix(%lu x %lu)\n\n",(*m).rows,(*m).cols);
	return;
}

void LINALG_COPY_MAT(LINALG_MAT *d, LINALG_MAT *s)
{
	if((*s).pt != NULL && (*s).rows != 0 && (*s).cols != 0)
	{
		uint64_t r = (*s).rows;
		uint64_t c = (*s).cols;
		uint64_t l = r*c;
		double *dn = NULL;
		if(!LINALG_MALLOC_DOUBLE(&dn, l))
		{
			for(uint64_t i = 0; i < l; ++i)
				dn[i] = (*s).pt[i];
			if((*d).pt == NULL)
				LINALG_DEF_MAT(d, dn, r, c);
			else
				LINALG_REDEF_MAT(d, dn, r, c);
			free(dn);
		} else
		{
			MALLOC_ERROR;
		}
	} else
	{
		COPY_ERROR;
	}
	return;
}

void LINALG_SUM_MAT(LINALG_MAT *s, LINALG_MAT *m1, LINALG_MAT *m2)
{
	if((*m1).rows ==(*m2).rows
			&& (*m1).rows != 0
			&& (*m1).cols ==(*m2).cols
			&& (*m1).cols != 0
			&& (*m1).pt != NULL
			&& (*m2).pt != NULL)
	{
		uint64_t r = (*m1).rows;
		uint64_t c = (*m1).cols;
		uint64_t l = r*c;
		double *d = NULL;
		if(!LINALG_MALLOC_DOUBLE(&d, l)) {
			for(uint64_t i = 0; i < l; ++i )
				d[i] = (*m1).pt[i] + (*m2).pt[i];
			if((*s).pt == NULL)
			{
				LINALG_DEF_MAT(s, d, r, c);
			}
			else
				LINALG_REDEF_MAT(s, d, r, c);
			free(d);
		} else {
			MALLOC_ERROR;
		}
	} else {
		SUM_ERROR;
	}
	return;
}

void LINALG_NEG_MAT(LINALG_MAT *d, LINALG_MAT *s) {
	if((*s).pt != NULL && (*s).rows != 0 && (*s).cols != 0) {
		uint64_t l = (*s).rows * (*s).cols;
		LINALG_COPY_MAT(d,s);
		for(uint64_t i = 0; i < l; ++i ) {
			(*d).pt[i] *= -1;
		} 
	} else
		OP_ERROR;
	return;
}

void LINALG_TIMES_MAT(LINALG_MAT *d, LINALG_MAT *s, double *a)
{
	if((*s).pt != NULL && (*s).rows != 0 && (*s).cols != 0) {
		uint64_t l = (*s).rows * (*s).cols;
		LINALG_COPY_MAT(d,s);
		for(uint64_t i = 0; i < l; ++i ) {
			(*d).pt[i] *= (*a);
		} 
	} else
		OP_ERROR;
	return;
}

void LINALG_MUL_MAT(LINALG_MAT *p, LINALG_MAT *m1, LINALG_MAT *m2)
{
	if((*m1).cols == (*m2).rows
			&& (*m1).rows != 0
			&& (*m1).cols != 0
			&& (*m2).rows != 0
			&& (*m2).cols != 0
			&& (*m1).pt != NULL
			&& (*m2).pt != NULL)
	{
		uint64_t r = (*m1).rows;
		uint64_t _p = (*m1).cols;
		uint64_t c = (*m2).cols;
		uint64_t l = r*c;
		double *d = NULL;
		if(!LINALG_MALLOC_DOUBLE(&d, l)) {
			for(uint64_t i = 0; i < r; ++i )
			{
				//d[i] = (*m1).pt[i] + (*m2).pt[i];
				for(uint64_t j = 0; j < c; ++j)
				{
					double sum = 0;
					for(uint64_t k = 0; k < _p; ++k)
						sum += (*m1).pt[i*_p+k] * (*m2).pt[k*c+j];
					d[i*c+j] = sum;
				}
			}
			if((*p).pt == NULL)
			{
				LINALG_DEF_MAT(p, d, r, c);
			}
			else
				LINALG_REDEF_MAT(p, d, r, c);
			free(d);
		} else {
			MALLOC_ERROR;
		}
	} else {
		OP_ERROR;
	}
	return;
}

void LINALG_TRS_MAT(LINALG_MAT *t, LINALG_MAT *m)
{
	if((*m).rows != 0
			&& (*m).cols != 0
			&& (*m).pt != NULL)
	{
		uint64_t r = (*m).rows;
		uint64_t c = (*m).cols;
		uint64_t tr = c;
		uint64_t tc = r;
		uint64_t l = r*c;
		double *d = NULL;
		if(!LINALG_MALLOC_DOUBLE(&d, l)) {
			for(uint64_t i = 0; i < tr; ++i)
			{
				for(uint64_t j = 0; j < tc; ++j)
				{
					d[i*tc+j] = (*m).pt[j*c+i];
				}
			}
			if((*t).pt == NULL)
			{
				LINALG_DEF_MAT(t, d, tr, tc);
			}
			else
				LINALG_REDEF_MAT(t, d, tr, tc);
			free(d);
		} else {
			MALLOC_ERROR;
		}
	} else {
		OP_ERROR;
	}
	return;
}


void LINALG_MC2V(LINALG_VEC *v, LINALG_MAT *m)
{
	uint64_t l = (*m).rows;
	for(uint64_t i = 0; i < l; ++i)
	{
		(*v).pt[i] = (*m).pt[i];
	}
	return;
}

void LINALG_MR2V(LINALG_VEC *v, LINALG_MAT *m)
{
	uint64_t l = (*m).cols;
	for(uint64_t i = 0; i < l; ++i)
	{
		(*v).pt[i] = (*m).pt[i];
	}
	return;
}
