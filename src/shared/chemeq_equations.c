#include <math.h>
#include "chemeq_definitions.h"
// PRESSURE
double FP(double *x, double *par)
{
	double res = 0;
	for(int i = 0; i < SD; ++i)
	{
		res += x[i];
	}
	return res-par[0];
}

double FNEUTR(double *x, double *par)
{
	double res = 0*par[0];
	if (e == 1)
		res -= x[SD-1];
	if (N == 3)
	{
		if (sN[2] == '1' || sN[3] == '1')
			res += x[2];
	}
	if (N == 4)
	{
		res += x[2] + x[3];
	}
	if (O == 3)
	{
		if (sO[2] == '1' || sO[3] == '1')
			res += x[N+2];
	}
	if (O == 4)
	{
		res += x[N+2] + x[N+3];
	}
	if (NO == 2)
		res += x[N+O+1];
	if (Ar == 2)
		res += x[SD-2];
	return res;
}

double FN2O(double *x, double *par)
{
	double res = 0;
	if (N >= 1)
		res += 2*x[0];
	if (N >= 2)
		res += x[1];
	if (N == 3)
	{
		if (sN[2] == '1')
			res += 2*x[2];
		if (sN[3] == '1')
			res += x[2];
	}
	if (N == 4)
		res += 2*x[2] + x[3];

	if (O >= 1)
		res -= 2*x[N]*par[eq_r1];
	if (O >= 2)
		res -= x[N+1]*par[eq_r1];
	if (O == 3)
	{
		if (sO[2] == '1')
			res -= 2*x[N+2]*par[eq_r1];
		if (sN[3] == '1')
			res -= x[N+2]*par[eq_r1];
	}
	if (O == 4)
		res -= (2*x[N+2] + x[N+3])*par[eq_r1];

	if (NO >= 1)
	{
		res += x[N+O]*(1-par[eq_r1]);
		if (sNO[1] == '1')
		{
			res += x[N+O+1]*(1-par[eq_r1]);
		}
	}
	return res;
}

double JRN1(double *x, double *par)
{
	return 1+0*x[0]*par[0];
}

double JRN2(double *x, double *par)
{
	return 2+0*x[0]*par[0];
}

double JRO1(double *x, double *par)
{
	return -1*par[eq_r1]+0*x[0]*par[0];
}

double JRO2(double *x, double *par)
{
	return -2*par[eq_r1]+0*x[0]*par[0];
}

double JRNO(double *x, double *par)
{
	return (1-par[eq_r1])+0*x[0]*par[0];
}

double FN2Ar(double *x, double *par)
{
	double res = 0;
	if (N >= 1)
		res += 2*x[0];
	if (N >= 2)
		res += x[1];
	if (N == 3)
	{
		if (sN[2] == '1')
			res += 2*x[2];
		if (sN[3] == '1')
			res += x[2];
	}
	if (N == 4)
		res += 2*x[2] + x[3];
	
	if (sAr[0] == '1')
		res -= x[N+O+NO]*par[eq_r2];
	if (sAr[1] == '1')
		res -= x[N+O+NO+1]*par[eq_r2];
	return res;
}

double JRAr(double *x, double *par)
{
	return -1*par[eq_r2]+0*x[0]*par[0];
}

double FDN(double *x, double *par)
{
	double res = pow(x[1],2) - par[offset+1]*x[0];
	if (N == 0)
		res = 0;
	return res;
}

double JDNp(double *x, double *par)
{
	double res = 2*x[1] + 0*par[0];
	if (N == 0)
		res = 0;
	return res;
}

double JDNr(double *x, double *par)
{
	double res = -par[offset+1] + 0*x[0];
	if (N == 0)
		res = 0;
	return res;
}

double FIN2(double *x, double *par)
{
	double res = x[2]*x[SD-1] - par[offset+2]*x[0];
	if (N == 0)
		res = 0;
	return res;
}

double JIp(double *x, double *par)
{
	double res = x[SD-1] + par[0]*0;
	return res;
}

double JIN2e(double *x, double *par)
{
	double res = x[2] + par[0]*0;
	if (N == 0)
		res = 0;
	return res;
}

double JIN2r(double *x, double *par)
{
	double res = - par[offset+2] + 0*x[0];
	if (N == 0)
		res = 0;
	return res;
}

double FIN1_1(double *x, double *par)
{
	double res = x[2]*x[SD-1] - par[offset+2]*x[1];
	if (N == 0)
		res = 0;
	return res;
}

double JIN1_1e(double *x, double *par)
{
	double res = x[2] - par[0]*0;
	if (N == 0)
		res = 0;
	return res;
}

double JIN1_1r(double *x, double *par)
{
	double res = - par[offset+2]+0*x[0];
	if (N == 0)
		res = 0;
	return res;
}

double FIN1_2(double *x, double *par)
{
	double res = x[3]*x[SD-1] - par[offset+3]*x[1];
	if (N == 0)
		res = 0;
	return res;
}

double JIN1_2e(double *x, double *par)
{
	double res = x[3] + par[0]*0;
	if (N == 0)
		res = 0;
	return res;
}

double JIN1_2r(double *x, double *par)
{
	double res = - par[offset+3] + 0*x[0];
	if (N == 0)
		res = 0;
	return res;
}

double FDO(double *x, double *par)
{
	if (N != 0) idx = offset+N;
	else idx = offset+N+1;
	double res = pow(x[N+1],2) - par[idx]*x[N];
	if (O == 0)
		res = 0;
	return res;
}

double JDOp(double *x, double *par)
{
	if (N != 0) idx = offset+N;
	else idx = offset+N+1;
	double res = 2*x[N+1] - par[0]*0;
	if (O == 0)
		res = 0;
	return res;
}

double JDOr(double *x, double *par)
{
	if (N != 0) idx = offset+N;
	else idx = offset+N+1;
	double res = - par[idx] + 0*x[0];
	if (O == 0)
		res = 0;
	return res;
}

double FIO2(double *x, double *par)
{
	if (N != 0) idx = offset+N+1;
	else idx = offset+N+2;
	double res = x[N+2]*x[SD-1] - par[idx]*x[N];
	if (O == 0)
		res = 0;
	return res;
}

double JIO2e(double *x, double *par)
{
	if (N != 0) idx = offset+N+1;
	else idx = offset+N+2;
	double res = x[N+2] - par[0]*0;
	if (O == 0)
		res = 0;
	return res;
}

double JIO2r(double *x, double *par)
{
	if (N != 0) idx = offset+N+1;
	else idx = offset+N+2;
	double res = - par[idx]+0*x[0];
	if (O == 0)
		res = 0;
	return res;
}

double FIO1_1(double *x, double *par)
{
	if (N != 0) idx = offset+N+1;
	else idx = offset+N+2;
	double res = x[N+2]*x[SD-1] - par[idx]*x[N+1];
	if (O == 0)
		res = 0;
	return res;
}

double JIO1_1e(double *x, double *par)
{
	if (N != 0) idx = offset+N+1;
	else idx = offset+N+2;
	double res = x[N+2] - par[0]*0;
	if (O == 0)
		res = 0;
	return res;
}

double JIO1_1r(double *x, double *par)
{
	if (N != 0) idx = offset+N+1;
	else idx = offset+N+2;
	double res = - par[idx]+0*x[0];
	if (O == 0)
		res = 0;
	return res;
}

double FIO1_2(double *x, double *par)
{
	if (N != 0) idx = offset+N+1;
	else idx = offset+N+2;
	double res = x[N+3]*x[SD-1] - par[idx+1]*x[N+1];
	if (O == 0)
		res = 0;
	return res;
}

double JIO1_2e(double *x, double *par)
{
	if (N != 0) idx = offset+N+1;
	else idx = offset+N+2;
	double res = x[N+3] - par[0]*0;
	if (O == 0)
		res = 0;
	return res;
}

double JIO1_2r(double *x, double *par)
{
	if (N != 0) idx = offset+N+1;
	else idx = offset+N+2;
	double res = - par[idx+1]+0*x[0];
	if (O == 0)
		res = 0;
	return res;
}

double FNO(double *x, double *par)
{
	double res = pow(x[N+O],2) - par[offset+N+O-1]*x[0]*x[N];
	if (NO == 0)
		res = 0;
	return res;
}

double JNO_N2(double *x, double *par)
{
	double res = - par[offset+N+O-1]*x[N];
	if (NO == 0)
		res = 0;
	return res;
}

double JNO_O2(double *x, double *par)
{
	double res = - par[offset+N+O-1]*x[0];
	if (NO == 0)
		res = 0;
	return res;
}

double JNO_2NO(double *x, double *par)
{
	double res = 2*x[N+O] + par[offset+N+O-1]*0;
	if (NO == 0)
		res = 0;
	return res;
}

double FINO(double *x, double *par)
{
	double res = x[N+O+1]*x[SD-1] - par[offset+N+O]*x[N+O];
	if (NO != 2)
		res = 0;
	return res;
}

double JINOr(double *x, double *par)
{
	double res = - par[offset+N+O]+0*x[N+O];
	if (NO != 2)
		res = 0;
	return res;
}

double JINOe(double *x, double *par)
{
	double res = x[N+O+1] - par[0]*0;
	if (NO != 2)
		res = 0;
	return res;
}

double FIAr(double *x, double *par)
{
	double res = x[SD-2]*x[SD-1] - par[SD-1]*x[SD-3];
	if (Ar != 2)
	{
		res = 0;
	}
	return res;
}

double JIArr(double *x, double *par)
{
	double res = 0*x[SD-2]*x[SD-1] - par[SD-1];
	if (Ar != 2)
	{
		res = 0;
	}
	return res;
}

double JIAre(double *x, double *par)
{
	double res = x[SD-2] - par[SD-1]*0;
	if (Ar != 2)
	{
		res = 0;
	}
	return res;
}

double J0(double *x, double *par)
{
	return 0*x[0]*par[0];
}

double Jp1(double *x, double *par)
{
	return 1+0*x[0]*par[0];
}

double Jm1(double *x, double *par)
{
	return -1+0*x[0]*par[0];
}
/*//
double 1F3(double *x, double *par)
{
	return x[2]*x[3] - par[2]*x[0];
}

double 1F4(double *x, double *par)
{
	return x[2]*x[3] - par[2]*x[1];
}

double 1F5(double *x, double *par)
{
	return x[2]*x[3] - par[2]*x[1];
}
*/
