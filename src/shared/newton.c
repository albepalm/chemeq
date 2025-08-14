#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <float.h>
#include "linalg.h"

bool LINALG_NEWTON(LINALG_FUNCTION_HANDLE *F,
		LINALG_FUNCTION_HANDLE *DF, uint64_t d,
		LINALG_VEC *x, LINALG_VEC *x0, LINALG_VEC *fpar,
		uint64_t max_iter, double abs_tol, 
		double lim_inf, double lim_sup,
		int output)
{
	if(lim_inf>lim_sup)
	{
		fprintf(stderr,"INVALID SOLUTION LIMITS!\n");
		exit(1);
	}
	bool SOLVER_FAIL = false;
	LINALG_VEC neg_system;
	LINALG_INIT_VEC(&neg_system);
	LINALG_VEC y;
	LINALG_INIT_VEC(&y);
	LINALG_MAT jacobian;
	LINALG_INIT_MAT(&jacobian);

	uint64_t i = 0;
	bool is_maxit_reached = false;
	bool is_tol_reached = false;
	double res = 0;

	do
	{
		LINALG_ZERO_VEC(&neg_system,d);
		LINALG_ZERO_MAT(&jacobian,d,d);
		LINALG_ZERO_VEC(&y,d);

		if(i == 0) {LINALG_COPY_VEC(x,x0);}

		// DOMAIN CHECK
		for(uint64_t ii = 0; ii < d; ++ii)
		{
			// to be completed...
			if(isnan((*x).pt[ii]))
			{
				(*x).pt[ii] = lim_inf;
			}
			(*x).pt[ii] = fabs((*x).pt[ii]);
			if((*x).pt[ii] < lim_inf)
			{
				(*x).pt[ii] = lim_inf;
			}
			if((*x).pt[ii] > lim_sup)
			{
				(*x).pt[ii] = lim_sup;
			}
		}
		for(uint64_t j = 0; j < d; ++j)
		{
			neg_system.pt[j] = (F[j]((*x).pt, (*fpar).pt))*(-1.0);
			for(uint64_t k = 0; k < d; ++k)
			{
				jacobian.pt[j*d+k] = (DF[j*d+k]((*x).pt, (*fpar).pt));
			}
		}
		//LINALG_PRINT_VEC(x,"Solution");
		//LINALG_PRINT_VEC(&neg_system,"Function");
		//LINALG_PRINT_MAT(&jacobian,"Jacobian");
		//getchar();
		res = LINALG_NORMMAX_VEC(&neg_system);
		if(i != 0 && (output >= 2))
		{
			printf("Iteration: %6lu\tResidual: %7.5e\n",i,res);
			if(output >= 3) {LINALG_PRINT_VEC(x,"x_i");}
			if(output >= 4) {getchar();}
		}
		if(i>=max_iter)
		{
			is_maxit_reached = true;
			if(fabs(res)>1e-5)
			{
				SOLVER_FAIL = true;
				printf("Convergenge not reached\n");
				fflush(stdout);
				return SOLVER_FAIL;
			}
		}
		if(fabs(res) < abs_tol)
		{
			is_tol_reached = true;
		}
		if(is_maxit_reached || is_tol_reached)
		{
			LINALG_NEG_VEC(&neg_system,&neg_system);
			if(output >= 1 )
			{
				printf("-=-=-=-=-=-=-=-=- NEWTON RAPHSON SOLVER -=-=-=-=-=-=-=-=- \n");
				printf("Iterations:\t%lu\n", i);
				LINALG_PRINT_VEC(x,"Solution");
				LINALG_PRINT_VEC(&neg_system,"Function");
				LINALG_PRINT_MAT(&jacobian,"Jacobian");
				printf("%c[%c] Max Iterations Reached (max_iter = %lu)\n",
						0x20,is_maxit_reached?0x78:0x20,max_iter);
				printf("%c[%c] Tol Abs (tol_abs = %7.5e)\n",
						0x20,is_tol_reached?0x78:0x20,abs_tol);
			}
		}
		if(!is_maxit_reached && !is_tol_reached)
		{
			++i;
			SOLVER_FAIL = LINALG_SOLVE(&y,&jacobian,&neg_system,d);
			if(SOLVER_FAIL)
			{
				return SOLVER_FAIL;
			}
			LINALG_SUM_VEC(x,&y,x);
		}
		//LINALG_PRINT_VEC(x,"Solution");
		//LINALG_PRINT_VEC(&neg_system,"Function");
		//LINALG_PRINT_MAT(&jacobian,"Jacobian");
		//getchar();
	} while(!is_maxit_reached && !is_tol_reached);

	LINALG_DEALLOC_VEC(&y);
	LINALG_DEALLOC_VEC(&neg_system);
	LINALG_DEALLOC_MAT(&jacobian);
	return SOLVER_FAIL;
}
