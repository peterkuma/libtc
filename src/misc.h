/*
 * misc.h
 *
 * Miscellaneous functions.
 *
 */

#include <stddef.h>
#include <stdio.h>

#ifdef DEBUG
#define debug(...) fprintf(stderr, __VA_ARGS__)
#else
#define debug(...)
#endif /* DEBUG */

#define MAX(x, y) ((x)>=(y)?(x):(y))
#define MIN(x, y) ((x)<=(y)?(x):(y))

double log_beta(double a, double b);

double frand(void);

double frand1(void);

size_t sample(size_t n, const double p[]);

void init_gsl(void);

void deinit_gsl(void);

double rtnorm(double mean, double sd, double a, double b);

void *array_insert(
    const void *array,
    size_t n,
    const void *element,
    size_t i,
    size_t size
);

void *array_remove(const void *array, size_t n, size_t i, size_t size);
