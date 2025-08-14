#ifndef LINALG_H
#define LINALG_H

#include <stdint.h>
#include <stdbool.h>

#define GET_NAME(VAR) #VAR
#define LEN_VEC(X) \
	( sizeof( (X) ) / sizeof( double ) )
//////////////////////////////////////////////////
//VECTOR HANDLING MACRS
#define INITV(V) LINALG_INIT_VEC(&(V))
#define ZEROV(V,L) \
	do{ LINALG_INIT_VEC(&(V)); \
		LINALG_ZERO_VEC(&(V),L); \
	} while(0)
#define ONESV(V,L) \
	do{ LINALG_INIT_VEC(&(V)); \
		LINALG_ONES_VEC(&(V),L); \
	} while(0)
#define RANDV(V,L) \
	do{ LINALG_INIT_VEC(&(V)); \
		LINALG_RAND_VEC(&(V),L); \
	} while(0)
#define DEFV(A,AUSR) \
	LINALG_DEF_VEC(&(A), (AUSR), LEN_VEC(AUSR))
#define REDEFV(A,AUSR) \
	LINALG_REDEF_VEC(&(A), (AUSR), LEN_VEC(AUSR))
#define PRINTV(A) \
	LINALG_PRINT_VEC(&(A),GET_NAME(A))
#define DEALLOCV(V) LINALG_DEALLOC_VEC(&(V))
#define SUMV(S,A,B) LINALG_SUM_VEC(&(S),&(A),&(B))
#define COPYV(D,S) LINALG_COPY_VEC(&(D),&(S))
#define NEGV(D,S) LINALG_NEG_VEC(&(D),&(S))
#define TIMESV(D,S,A) LINALG_TIMES_VEC(&(D),&(S),&(A))
#define DOTV(A,B) LINALG_DOT_VEC(&(A),&(B))
#define SOLVE(X,A,B,D) \
	LINALG_SOLVE(&(X),&(A),&(B),D)
#define SOLVEWITHNORM(X,A,B,D) \
	LINALG_SOLVEWITHNORM(&(X),&(A),&(B),D)
//////////////////////////////////////////////////
//MATRIX HANDLING MACRS
#define INITM(M) LINALG_INIT_MAT(&(M))
#define ZEROM(V,R,C) \
	do{ LINALG_INIT_MAT(&(V)); \
		LINALG_ZERO_MAT(&(V),R,C); \
	} while(0)
#define ONESM(V,R,C) \
	do{ LINALG_INIT_MAT(&(V)); \
		LINALG_ONES_MAT(&(V),R,C); \
	} while(0)
#define RANDM(V,R,C) \
	do{ LINALG_INIT_MAT(&(V)); \
		LINALG_RAND_MAT(&(V),R,C); \
	} while(0)
#define IDM(V,R) \
	do{ LINALG_INIT_MAT(&(V)); \
		LINALG_ID_MAT(&(V),R); \
	} while(0)
#define DEFM(A,AUSR,R,C) \
	LINALG_DEF_MAT(&(A), *(AUSR), (R), (C))
#define REDEFM(A,AUSR,R,C) \
	LINALG_REDEF_MAT(&(A), *(AUSR), (R), (C))
#define PRINTM(A) \
	LINALG_PRINT_MAT(&(A),GET_NAME(A))
#define DEALLOCM(M) LINALG_DEALLOC_MAT(&(M))
#define SUMM(S,A,B) LINALG_SUM_MAT(&(S),&(A),&(B))
#define COPYM(D,S) LINALG_COPY_MAT(&(D),&(S))
#define NEGM(D,S) LINALG_NEG_MAT(&(D),&(S))
#define TIMESM(D,S,A) LINALG_TIMES_MAT(&(D),&(S),(A))
#define MULM(P,A,B) LINALG_MUL_MAT(&(P),&(A),&(B))
//////////////////////////////////////////////////
//ERROR MACRO
#define INIT_ERROR \
	do { printf("Init failure\n"); } while(0) 
#define ASSIGN_ERROR \
	do { printf("Assign failure\n"); } while(0) 
#define MALLOC_ERROR \
	do { printf("Malloc failure\n"); } while(0) 
#define REALLOC_ERROR \
	do { printf("Realloc failure\n"); } while(0) 
#define COPY_ERROR \
	do { printf("Copy failure\n"); } while(0)
#define DOT_ERROR \
	do { printf("Dot product failure\n"); } while(0) 
#define SUM_ERROR \
	do { printf("Sum failure\n"); } while(0) 
#define OP_ERROR \
	do { printf("Operation failure\n"); } while(0) 

// SOLUTION STRUCT - simple structure of pointers
//typedef struct {}
// FUNCTION HANDLER:
typedef double		(*LINALG_FUNCTION_HANDLE)(double *, double *);

// DOUBLE MALLOC
bool LINALG_MALLOC_DOUBLE(double **d, uint64_t l);

// VECTOR HANDLER
typedef struct {
	uint64_t	len;
	double		*pt;
} LINALG_VEC;

// VECTOR INITIALIZATION (0 & NULL)
void LINALG_INIT_VEC(LINALG_VEC *v);

// VECTORT ASSIGNMENT
void LINALG_ASSIGN_VEC(LINALG_VEC *d, double *s, uint64_t l);

// VECTOR MEM ALLOCATION
bool LINALG_MALLOC_VEC(LINALG_VEC *v, uint64_t l);

// VECTOR MEM REALLOCATION
bool LINALG_REALLOC_VEC(LINALG_VEC *v, uint64_t l);

// VECTOR VALUE DEFINITION
void LINALG_DEF_VEC(LINALG_VEC *v, double *usrv, uint64_t l);

// BASIC VECTORS
void LINALG_ZERO_VEC(LINALG_VEC *v, uint64_t l);
void LINALG_ONES_VEC(LINALG_VEC *v, uint64_t l);
void LINALG_RAND_VEC(LINALG_VEC *v, uint64_t l);

// VECTOR VALUE REDEFINITION
void LINALG_REDEF_VEC(LINALG_VEC *v, double *usrv, uint64_t l);

// VECTOR DEALLOCATION
void LINALG_DEALLOC_VEC(LINALG_VEC *v);

// VECTOR PRINT
void LINALG_PRINT_VEC(LINALG_VEC *v, char *name);

// VECTOR OPERATIONS
void LINALG_COPY_VEC(LINALG_VEC *d, LINALG_VEC *s);
void LINALG_SUM_VEC(LINALG_VEC *s, LINALG_VEC *v1, LINALG_VEC *v2);
void LINALG_NEG_VEC(LINALG_VEC *d, LINALG_VEC *s);
void LINALG_TIMES_VEC(LINALG_VEC *d, LINALG_VEC *s, double *a);
double LINALG_DOT_VEC(LINALG_VEC *v1, LINALG_VEC *v2);
double LINALG_NORMMAX_VEC(LINALG_VEC *v);

// MATRIX HANDLER
typedef struct {
	uint64_t	rows;
	uint64_t	cols;
	double		*pt;
} LINALG_MAT;

// MATRIX INITIALIZATION
void LINALG_INIT_MAT(LINALG_MAT *m);

// MATRIX ASSIGNMENT
void LINALG_ASSIGN_MAT(LINALG_MAT *m, double *usrm, uint64_t r, uint64_t c);

// MATRIX ALLOCATION
bool LINALG_MALLOC_MAT(LINALG_MAT *m, uint64_t r, uint64_t c); 

// MATRIX ALLOCATION
bool LINALG_REALLOC_MAT(LINALG_MAT *m, uint64_t r, uint64_t c); 

// DEFINE MATRIX
void LINALG_DEF_MAT(LINALG_MAT *m, double *usrm, uint64_t r, uint64_t c);

// REDEFINE MATRIX
void LINALG_REDEF_MAT(LINALG_MAT *m, double *usrm, uint64_t r, uint64_t c); 

// MATRIX DEALLOCATION	
void LINALG_DEALLOC_MAT(LINALG_MAT *m); 

// MATRIX PRINT
void LINALG_PRINT_MAT(LINALG_MAT *m, char *name); 

// BASIC MATRIX
void LINALG_ZERO_MAT(LINALG_MAT *m, uint64_t r, uint64_t c);
void LINALG_ID_MAT(LINALG_MAT *m, uint64_t r);
void LINALG_ONES_MAT(LINALG_MAT *m, uint64_t r, uint64_t c);
void LINALG_RAND_MAT(LINALG_MAT *m, uint64_t r, uint64_t c);

// MATRIX OPERATIONS
void LINALG_COPY_MAT(LINALG_MAT *d, LINALG_MAT *s);
void LINALG_SUM_MAT(LINALG_MAT *s, LINALG_MAT *m1, LINALG_MAT *m2);
void LINALG_NEG_MAT(LINALG_MAT *d, LINALG_MAT *s);
void LINALG_TIMES_MAT(LINALG_MAT *d, LINALG_MAT *s, double *a);
void LINALG_MUL_MAT(LINALG_MAT *p, LINALG_MAT *m1, LINALG_MAT *m2);
void LINALG_TRS_MAT(LINALG_MAT *t, LINALG_MAT *m);
void LINALG_MC2V(LINALG_VEC *v, LINALG_MAT *m);
void LINALG_MR2V(LINALG_VEC *v, LINALG_MAT *m);
void LINALG_V2MR(LINALG_MAT *m, LINALG_VEC *v);
void LINALG_V2MC(LINALG_MAT *m, LINALG_VEC *v);

// LINEAR SYSTEM
uint64_t LINALG_RANK(LINALG_MAT *A);
bool LINALG_GAUSSELIM(LINALG_MAT *A, LINALG_VEC *b,
		LINALG_MAT *eA, LINALG_VEC *eb);
void LINALG_BACKSUBS(LINALG_VEC *x, LINALG_MAT *U, LINALG_VEC *b);
bool LINALG_SOLVE(LINALG_VEC *x, LINALG_MAT *A, LINALG_VEC *b, uint64_t d);
void LINALG_SOLVEWITHNORM(LINALG_VEC *x, LINALG_MAT *A, LINALG_VEC *b, uint64_t d);

// NEWTON SOLVER
bool LINALG_NEWTON(LINALG_FUNCTION_HANDLE *F,
		LINALG_FUNCTION_HANDLE *DF, uint64_t d,
		LINALG_VEC *x, LINALG_VEC *x0, LINALG_VEC *fpar,
		uint64_t max_iter, double abs_tol, 
		double lim_inf, double lim_sup,
		int output);
void LINALG_QNEWTON(LINALG_FUNCTION_HANDLE *F,
		double epsJ, uint64_t d,
		LINALG_VEC *x, LINALG_VEC *x0, LINALG_VEC *fpar,
		uint64_t max_iter, double abs_tol, 
		double lim_inf, double lim_sup,
		int output);
#endif //LINALG_H
