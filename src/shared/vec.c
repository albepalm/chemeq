#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "linalg.h"

bool LINALG_MALLOC_DOUBLE(double **d, uint64_t l) {
	bool res;
	*d = (double *) malloc(l*sizeof(double));
	//d = NULL;
	if(d == NULL) {
		res = true;
		//printf("double");
		MALLOC_ERROR;
	}
	else
		res = false;
	return res;
}
void LINALG_INIT_VEC(LINALG_VEC *v) {
	(*v).len = 0;
	(*v).pt = NULL;
	return;
}

void LINALG_ASSIGN_VEC(LINALG_VEC *v, double *usrv, uint64_t l) {
	for(uint64_t i = 0; i < l; ++i)
		(*v).pt[i] = usrv[i];
	return;
}

bool LINALG_MALLOC_VEC(LINALG_VEC *v, uint64_t l) {
	bool ret;
	(*v).len = l;
	(*v).pt = (double *) malloc(l*sizeof(double));
	if ((*v).pt == NULL) {
		MALLOC_ERROR;
		ret = true;
	}
	else {
		ret = false;
	}
	return ret;
}

bool LINALG_REALLOC_VEC(LINALG_VEC *v, uint64_t l) {
	bool ret;
	(*v).len = l;
	(*v).pt = (double *) realloc((*v).pt, l*sizeof(double));
	if ((*v).pt == NULL) {
		MALLOC_ERROR;
		ret = true;
	}
	else {
		ret = false;
	}
	return ret;
}

void LINALG_DEF_VEC(LINALG_VEC *v, double *usrv, uint64_t l) {
	//LINALG_INIT_VEC(v);
	if((*v).len == 0 && (*v).pt == NULL) {
		if (!LINALG_MALLOC_VEC(v,l))
			LINALG_ASSIGN_VEC(v, usrv, l);
		else
			MALLOC_ERROR;
	} else
		INIT_ERROR;
	//printf("New vector defined!\n");
	return;
}

void LINALG_REDEF_VEC(LINALG_VEC *v, double *usrv, uint64_t l) {
	//if((*v).len != 0 && (*v).pt != NULL) {
		if (!LINALG_REALLOC_VEC(v,l)) {
			LINALG_ASSIGN_VEC(v, usrv, l);
			//printf("Vector redefined!\n");
		}
		else
			MALLOC_ERROR;
	//} else
	//		INIT_ERROR;
	return;
}

void LINALG_ZERO_VEC(LINALG_VEC *v, uint64_t l)
{
	if((*v).pt == 0 && (*v).len == 0)
	{
		LINALG_MALLOC_VEC(v, l);
	} else
	{
		LINALG_REALLOC_VEC(v, l);
	}
	for(uint64_t i = 0; i < l; ++i)
	{
		(*v).pt[i] = 0;
	}
	return;
}

void LINALG_ONES_VEC(LINALG_VEC *v, uint64_t l)
{
	if((*v).pt == 0 && (*v).len == 0)
	{
		LINALG_MALLOC_VEC(v, l);
	} else
	{
		LINALG_REALLOC_VEC(v, l);
	}
	for(uint64_t i = 0; i < l; ++i)
	{
		(*v).pt[i] = 1;
	}
	return;
}

void LINALG_RAND_VEC(LINALG_VEC *v, uint64_t l)
{
	if((*v).pt == 0 && (*v).len == 0)
	{
		LINALG_MALLOC_VEC(v, l);
	} else
	{
		LINALG_REALLOC_VEC(v, l);
	}
	srand(time(NULL));
	for(uint64_t i = 0; i < l; ++i)
	{
		(*v).pt[i] = ((double) rand())/((double) RAND_MAX);
	}
	return;
}

void LINALG_DEALLOC_VEC(LINALG_VEC *v) {
	free((*v).pt);
	LINALG_INIT_VEC(v);
}

void LINALG_PRINT_VEC(LINALG_VEC *v, char *name) {
	printf("%s = [\n", name);
	for(uint64_t i = 0; i < (*v).len; ++i)
		printf("\t%7.5e\n", (*v).pt[i]);
		//printf("\t%7.5f\n", (*v).pt[i]);
	printf("]\tvector(%lu)\n\n",(*v).len);
	return;
}

void LINALG_COPY_VEC(LINALG_VEC *d, LINALG_VEC *s) {
	if((*s).pt != NULL && (*s).len != 0) {
		uint64_t l = (*s).len;
		double *dn = NULL;
		if(!LINALG_MALLOC_DOUBLE(&dn, l)) {
			for(uint64_t i = 0; i < l; ++i)
				dn[i] = (*s).pt[i];
			if((*d).pt == NULL)
				LINALG_DEF_VEC(d, dn, l);
			else
				LINALG_REDEF_VEC(d, dn, l);
			free(dn);
		} else {
			MALLOC_ERROR;
		}
	} else {
		COPY_ERROR;
	}
	return;
}

void LINALG_SUM_VEC(LINALG_VEC *s, LINALG_VEC *v1, LINALG_VEC *v2) {
	if((*v1).len==(*v2).len 
			&& (*v1).len != 0
			&& (*v1).pt != NULL
			&& (*v2).pt != NULL) {
		uint64_t l = (*v1).len;
		double *d = NULL;
		if(!LINALG_MALLOC_DOUBLE(&d, l)) {
			for(uint64_t i = 0; i < l; ++i )
				d[i] = (*v1).pt[i] + (*v2).pt[i];
			if((*s).pt == NULL)
				LINALG_DEF_VEC(s, d, l);
			else
				LINALG_REDEF_VEC(s, d, l);
			free(d);
		} else {
			MALLOC_ERROR;
		}
	} else {
		SUM_ERROR;
	}
	return;
}

void LINALG_NEG_VEC(LINALG_VEC *d, LINALG_VEC *s) {
	if((*s).pt != NULL && (*s).len != 0) {
		uint64_t l = (*s).len;
		LINALG_COPY_VEC(d,s);
		for(uint64_t i = 0; i < l; ++i ) {
			(*d).pt[i] *= -1;
		} 
	} else
		OP_ERROR;
	return;
}

void LINALG_TIMES_VEC(LINALG_VEC *d, LINALG_VEC *s, double *a) {
	if((*s).pt != NULL && (*s).len != 0) {
		uint64_t l = (*s).len;
		LINALG_COPY_VEC(d,s);
		for(uint64_t i = 0; i < l; ++i ) {
			(*d).pt[i] *= (*a);
		} 
	} else
		OP_ERROR;
	return;
}

double LINALG_DOT_VEC(LINALG_VEC *v1, LINALG_VEC *v2) {
	double res = 0;
	if((*v1).len==(*v2).len 
			&& (*v1).len != 0
			&& (*v1).pt != NULL
			&& (*v2).pt != NULL) {
		uint64_t l = (*v1).len;
		for(uint64_t i = 0; i < l; ++i )
			res += ((*v1).pt[i]) * ((*v2).pt[i]);
	} else {
		DOT_ERROR;
	}
	return res;
}

double LINALG_NORMMAX_VEC(LINALG_VEC *v)
{
	double res = 0;
	if((*v).len != 0
			&& (*v).pt != NULL)
	{
		uint64_t l = (*v).len;
		for(uint64_t i = 0; i < l; ++i )
		{
			if(fabs((*v).pt[i]) > res)
			{
				res = fabs((*v).pt[i]);
			}
		}
	} else {
		OP_ERROR;
	}
	return res;
}

void LINALG_V2MC(LINALG_MAT *m, LINALG_VEC *v)
{
	uint64_t r = (*v).len;
	for(uint64_t i = 0; i < r; ++i)
	{
		(*m).pt[i] = (*v).pt[i];
	}
	return;
}

void LINALG_V2MR(LINALG_MAT *m, LINALG_VEC *v)
{
	uint64_t c = (*v).len;
	for(uint64_t i = 0; i < c; ++i)
	{
		(*m).pt[i] = (*v).pt[i];
	}
	return;
}
