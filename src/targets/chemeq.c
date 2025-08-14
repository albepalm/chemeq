#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <unistd.h>

#include "linalg.h"
#include "tp.h"
#include "utils_chemeq.h"
#include "chemeq_definitions.h"
#include "chemeq_equations.h"
#define STRLEN 256

//GLOBAL VARIABLES
int SD;
int offset;
int eq_r1;
int eq_r2;
int idx;
int N;
int O;
int NO;
int Ar;
int e;
char sN[5];
char sO[5];
char sNO[3];
char sAr[3];
double rN;
double rO;
double rAr;

int main(int argc, char *argv[])
{
	CHEMEQ_SPLASH();
	char output_path[STRLEN] = "output/";
	char output_pp_path[STRLEN] = "output/";
	char config_file[STRLEN] = "chemeq.config";
	char input_path[STRLEN];
	if (argc != 2)
		exit(1);
	strcpy(input_path,argv[1]);
	char mod;
	uint64_t max_iter = 1*1000;
	double abs_tol = 1e-15;
	int OUTFORM = 1;
	bool SOLVER_FAIL = true;
	int SOLVER_FAIL_COUNTER = 0;
	double p_min;
	double p_max;
	double t_min;
	double t_max;
	uint64_t p_step;
	uint64_t t_step;
	char specs[STRLEN];
	FILE *fp_input = fopen(input_path,"r");
	if (fscanf(fp_input,
				"%c %lf %lf %lu %lf %lf %lu %s %d %s %d %s %d %s %d %s %d %lf %lf %lf\n",
			&mod,
			&p_min, &p_max, &p_step,
			&t_min, &t_max, &t_step,
			specs,
			&N, sN,
			&O, sO,
			&NO, sNO,
			&Ar, sAr,
			&e,
			&rN, &rO, &rAr) != 20)
		exit(1);
	fclose(fp_input);
	FILE *fp_config = fopen(config_file, "r");
	if (fp_config == NULL){exit(1);}
	if (fscanf(fp_config,"%lu %lf %d\n", &max_iter, &abs_tol, &OUTFORM) != 3)
	{
		//printf("%lu %e %d\n", max_iter, abs_tol, OUTFORM);
		exit(1);
	}
	fclose(fp_config);
	strcat(output_path,specs);
	strcat(output_path,".txt");
	strcat(output_pp_path,specs);
	strcat(output_pp_path,"_pp.txt");
	double *p_vec = NULL;
	double *t_vec = NULL;
	CHEMEQ_LINSPACE(&t_vec, t_min, t_max, t_step);
	CHEMEQ_LOG10SPACE(&p_vec, p_min, p_max, p_step);
	
	SD = N+O+NO+Ar+e;
//	printf("%d\n",SD);
//	printf("%d %s\n",N,sN);
//	printf("%d %s\n",O,sO);
//	printf("%d %s\n",NO,sNO);
//	printf("%d %s\n",Ar,sAr);
//	printf("%d\n",e);
//	printf("%lf %lf %lf\n",rN,rO,rAr);
//	printf("%lf %lf %lu\n",p_min,p_max,p_step);
//	printf("%lf %lf %lu\n",t_min,t_max,t_step);
//	getchar();
	
	// SYSTEM VECTORS
	LINALG_VEC x, x0, fpar;
	LINALG_INIT_VEC(&x);
	LINALG_INIT_VEC(&x0);
	LINALG_INIT_VEC(&fpar);
	LINALG_ZERO_VEC(&x,SD);
	LINALG_RAND_VEC(&x0,SD);
	LINALG_ZERO_VEC(&fpar,SD);
	LINALG_FUNCTION_HANDLE *F = 
		(LINALG_FUNCTION_HANDLE *) malloc(SD*sizeof(LINALG_FUNCTION_HANDLE));
	/*
	LINALG_FUNCTION_HANDLE **J =
		(LINALG_FUNCTION_HANDLE **) malloc(SD*sizeof(LINALG_FUNCTION_HANDLE *));
	for(int i = 0; i < SD; ++i)
	{
		J[i] =
		(LINALG_FUNCTION_HANDLE *) malloc(SD*sizeof(LINALG_FUNCTION_HANDLE));
	}
	*/
	LINALG_FUNCTION_HANDLE *J =
		(LINALG_FUNCTION_HANDLE *) malloc(SD*SD*sizeof(LINALG_FUNCTION_HANDLE));
	for(int i = 0; i < SD; ++i)
	{
		for(int j = 0; j < SD; ++j)
		{
			J[i*SD+j] = J0;
			//printf("%lf\n",J[i][j](x.pt,fpar.pt));
		}
	}

    //SYSTEM MATRIX GEN AND SOLVE
	FILE *fp_output = fopen(output_path,"w");
	FILE *fp_output_pp = fopen(output_pp_path,"w");
	int first_pressure_iteration = 1;
	for(uint64_t i = 0; i < p_step; ++i)
	{
		if (first_pressure_iteration == 1)
		{
			first_pressure_iteration = 1;
			offset = 0;
			fpar.pt[0] = p_vec[i];
			F[0] = FP;
			for(int ij = 0; ij < SD; ++ij)
			{
				J[ij] = Jp1;
			}
			if (e != 0)
			{
				offset += 1;
				F[offset] = FNEUTR;
				J[offset*SD+SD-1] = Jm1;
				if (N == 3) {J[offset*SD+2]=Jp1;}
				if (N == 4) {J[offset*SD+2]=Jp1;J[offset*SD+3]=Jp1;}
				if (O == 3) {J[offset*SD+N+2]=Jp1;}
				if (O == 4) {J[offset*SD+N+2]=Jp1;J[offset*SD+N+3]=Jp1;}
				if (NO == 2) {J[offset*SD+N+O+1]=Jp1;}
				if (Ar == 2) {J[offset*SD+N+O+NO+1]=Jp1;}

			}
			if (N != 0 && O != 0)
			{
				offset += 1;
				eq_r1 = offset;
				F[offset] = FN2O;
				fpar.pt[offset] = rN/rO;

				if (N >= 1) {J[offset*SD+0]=JRN2;}
				if (N >= 2) {J[offset*SD+1]=JRN1;}
				if (N == 3)
				{
					if (sN[2] == '1') {J[offset*SD+2]=JRN2;}
					if (sN[3] == '1') {J[offset*SD+2]=JRN1;}
				}
				if (N == 4) {J[offset*SD+2]=JRN2; J[offset*SD+3]=JRN1;}

				if (O >= 1) {J[offset*SD+N]=JRO2;}
				//if (O >= 1) {J[offset*SD+N]=Jm1;}
				if (O >= 2) {J[offset*SD+N+1]=JRO1;}
				//if (O >= 2) {J[offset*SD+N+1]=Jm1;}
				if (O == 3)
				{
					if (sO[2] == '1') {J[offset*SD+N+2]=JRO2;}
					if (sO[3] == '1') {J[offset*SD+N+2]=JRO1;}
				}
				if (O == 4) {J[offset*SD+N+2]=JRO2; J[offset*SD+N+3]=JRO1;}

				if (NO >= 1) {J[offset*SD+N+O]=JRNO;}
				if (NO == 2) {J[offset*SD+N+O+1]=JRNO;}

				if (Ar != 0)
				{
					offset += 1;
					eq_r2 = offset;
					F[offset] = FN2Ar;
					fpar.pt[offset] = rN/rAr;

					if (N >= 1) {J[offset*SD+0]=JRN2;}
					if (N >= 2) {J[offset*SD+1]=JRN1;}
					if (N == 3)
					{
						if (sN[2] == '1') {J[offset*SD+2]=JRN2;}
						if (sN[3] == '1') {J[offset*SD+2]=JRN1;}
					}
					if (N == 4) {J[offset*SD+2]=JRN2; J[offset*SD+3]=JRN1;}

					if (Ar >= 1) {J[offset*SD+N+O+NO]=JRAr;}
					if (Ar == 2) {J[offset*SD+N+O+NO+1]=JRAr;}
				}
			}
		}
		for(uint64_t j = 0; j < t_step; ++j)
		{
			if (N>=2)
			{
				fpar.pt[1+offset] = KP((t_vec+j),npdn2,spdn2,pdn2,
						nrdn2,srdn2,rdn2);
				F[1+offset] = FDN;
				J[(1+offset)*SD+0] = JDNr;
				J[(1+offset)*SD+1] = JDNp;
			}
			if (N==3)
			{
				if (sN[2] == '1')
				{
					fpar.pt[2+offset] = KP((t_vec+j),npin2,spin2,pin2,
						nrin2,srin2,rin2);
					F[2+offset] = FIN2;
					J[(2+offset)*SD+0] = JIN2r;
					J[(2+offset)*SD+2] = JIp;
					J[(2+offset)*SD+SD-1] = JIN2e;
				}
				if (sN[3] == '1')
				{
					fpar.pt[2+offset] = KP((t_vec+j),npin,spin,pin,
						nrin,srin,rin);
					F[2+offset] = FIN1_1;
					J[(2+offset)*SD+1] = JIN1_1r;
					J[(2+offset)*SD+2] = JIp;
					J[(2+offset)*SD+SD-1] = JIN1_1e;
				}
			}
			if (N==4)
			{
				fpar.pt[2+offset] = KP((t_vec+j),npin2,spin2,pin2,
					nrin2,srin2,rin2);
				fpar.pt[3+offset] = KP((t_vec+j),npin,spin,pin,
					nrin,srin,rin);
				F[2+offset] = FIN2;
				J[(2+offset)*SD+0] = JIN2r;
				J[(2+offset)*SD+2] = JIp;
				J[(2+offset)*SD+SD-1] = JIN2e;
				F[3+offset] = FIN1_2;
				J[(3+offset)*SD+1] = JIN1_2r;
				J[(3+offset)*SD+3] = JIp;
				J[(3+offset)*SD+SD-1] = JIN1_2e;
			}
			if (O>=2)
			{
				if (N != 0) idx = offset+N;
				else idx = offset+1;
				fpar.pt[idx] = KP((t_vec+j),npdo2,spdo2,pdo2,
						nrdo2,srdo2,rdo2);
				F[idx] = FDO;
				J[(idx)*SD+N] = JDOr;
				J[(idx)*SD+N+1] = JDOp;
			}
			if (O==3)
			{
				if (N != 0) idx = offset+N+1;
				else idx = offset+2;
				if (sO[2] == '1')
				{
					fpar.pt[idx] = KP((t_vec+j),npio2,spio2,pio2,
						nrio2,srio2,rio2);
					F[idx] = FIO2;
					J[(idx)*SD+N] = JIO2r;
					J[(idx)*SD+N+2] = JIp;
					J[(idx)*SD+SD-1] = JIO2e;
				}
				if (sO[3] == '1')
				{
					fpar.pt[idx] = KP((t_vec+j),npio,spio,pio,
						nrio,srio,rio);
					F[idx] = FIO1_1;
					J[(idx)*SD+N+1] = JIO1_1r;
					J[(idx)*SD+N+2] = JIp;
					J[(idx)*SD+SD-1] = JIO1_1e;
				}
			}
			if (O==4)
			{
				if (N != 0) idx = offset+N+1;
				else idx = offset+2;
				fpar.pt[idx] = KP((t_vec+j),npio2,spio2,pio2,
					nrio2,srio2,rio2);
				fpar.pt[idx+1] = KP((t_vec+j),npio,spio,pio,
					nrio,srio,rio);
				F[idx] = FIO2;
				J[(idx)*SD+N] = JIO2r;
				J[(idx)*SD+N+2] = JIp;
				J[(idx)*SD+SD-1] = JIO2e;
				F[idx+1] = FIO1_2;
				J[(idx+1)*SD+N+1] = JIO1_2r;
				J[(idx+1)*SD+N+3] = JIp;
				J[(idx+1)*SD+SD-1] = JIO1_2e;
			}
			if (NO>=1)
			{
				if (sNO[0] == '1')
				{
					fpar.pt[offset+N+O-1] = KP((t_vec+j),npno,spno,pno,
						nrno,srno,rno);
					F[offset+N+O-1] = FNO;
					J[(offset+N+O-1)*SD+0] = JNO_N2;
					J[(offset+N+O-1)*SD+N] = JNO_O2;
					J[(offset+N+O-1)*SD+N+O] = JNO_2NO;
				}
				if (sNO[1] == '1')
				{
					fpar.pt[offset+N+O] = KP((t_vec+j),npino,spino,pino,
						nrino,srino,rino);
					F[offset+N+O] = FINO;
					J[(offset+N+O)*SD+N+O] = JINOr;
					J[(offset+N+O)*SD+N+O+1] = JIp;
					J[(offset+N+O)*SD+SD-1] = JINOe;
				}	
			}
			if (Ar == 2)
			{
				fpar.pt[SD-1] = KP((t_vec+j),npiar,spiar,piar,
					nriar,sriar,riar);
				F[SD-1] = FIAr;
				J[(SD-1)*SD+SD-3] = JIArr;
				J[(SD-1)*SD+SD-2] = JIp;
				J[(SD-1)*SD+SD-1] = JIAre;
			}
			
			while(SOLVER_FAIL)
			{
				SOLVER_FAIL = LINALG_NEWTON(F,J,SD,&x,&x0,&fpar,
					max_iter,abs_tol,DBL_MIN,DBL_MAX,OUTFORM);
				if(SOLVER_FAIL)
				{
					SOLVER_FAIL_COUNTER += 1;
					if(SOLVER_FAIL_COUNTER > 10)
					{
						SOLVER_FAIL = false;
						//printf("Something is wrong...\n\n");
						//exit(1);
					}
					printf("Solver fails! New x0 vector (fail count: %d)\n",
							SOLVER_FAIL_COUNTER);
					fflush(stdout);
					LINALG_RAND_VEC(&x0,SD);
				}
			}
			SOLVER_FAIL_COUNTER = 0;
			SOLVER_FAIL = true;
			LINALG_COPY_VEC(&x0,&x);

			if (OUTFORM == 0)
			{
				printf("STATUS: Current Pressure:%12.5e (%lu of %lu)\t"
						"Current Temperature:%12.5e (%lu of %lu)\r",
						p_vec[i],i+1,p_step,t_vec[j],j+1,t_step);
				fflush(stdout);
			}
			if (mod == 't')
			{
				fprintf(fp_output,"%7.5e",t_vec[j]);
				fprintf(fp_output_pp,"%7.5e",t_vec[j]);
				for(int k = 0; k < SD; ++k)
				{
					fprintf(fp_output,"\t%7.5e",x.pt[k]/p_vec[i]);
					fprintf(fp_output_pp,"\t%7.5e",x.pt[k]);
				}
				fprintf(fp_output,"\n");
				fprintf(fp_output_pp,"\n");
			}
			if (mod == 'm')
			{
				fprintf(fp_output,"%7.5e",t_vec[j]);
				fprintf(fp_output,"\t%7.5e",p_vec[i]);
				fprintf(fp_output_pp,"%7.5e",t_vec[j]);
				fprintf(fp_output_pp,"\t%7.5e",p_vec[i]);
				for(int k = 0; k < SD; ++k)
				{
					fprintf(fp_output,"\t%7.5e",x.pt[k]/p_vec[i]);
					fprintf(fp_output_pp,"\t%7.5e",x.pt[k]);
				}
				fprintf(fp_output,"\n");
				fprintf(fp_output_pp,"\n");
			}
		}
	}
	printf("\n\n\n");
	fclose(fp_output);
	fclose(fp_output_pp);
	LINALG_DEALLOC_VEC(&x);
	LINALG_DEALLOC_VEC(&x0);
	LINALG_DEALLOC_VEC(&fpar);
	free(F);
	free(J);
	return 0;
}
