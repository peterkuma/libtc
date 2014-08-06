/*
 * misc.h
 *
 * Miscellaneous functions.
 *
 */

#include <stddef.h>
#include <stdio.h>

#define debug(...) fprintf(stderr, __VA_ARGS__)

#define MAX(x, y) ((x)>=(y)?(x):(y))
#define MIN(x, y) ((x)<=(y)?(x):(y))

double log_beta(double a, double b);

double frand(void);

double frand1(void);

size_t sample(size_t n, const double p[]);

void init_gsl(void);

double rtnorm(double mean, double sd, double a, double b);
