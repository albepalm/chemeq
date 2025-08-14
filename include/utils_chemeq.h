#ifndef UTILS_CHEMEQ_H
#define UTILS_CHEMEQ_H

#include <stdint.h>

void CHEMEQ_SPLASH();
void CHEMEQ_CLS();
int CHEMEQ_MAIN();
int CHEMEQ_USR_INT();
double CHEMEQ_USR_DBL();
int CHEMEQ_SPEC();
double CHEMEQ_O2N_RAT();
void CHEMEQ_PRES(double *min, double *max, int *n);
void CHEMEQ_TEMP(double *min, double *max, int *n);
void CHEMEQ_ATOMS_CHECK(char **input, int *input_len, int *atoms_array,
		int *atoms_number_list);
void CHEMEQ_SPECIES_CHECK(char **input, int *input_len, int *species_arrayi,
		int *species_number_list, int* true_species_number_list);

void CHEMEQ_LINSPACE(double **lin, double min, double max, uint64_t n);
void CHEMEQ_LOG10SPACE(double **lin, double min, double max, uint64_t n);
#endif	// UTILS_CHEMEQ_H
