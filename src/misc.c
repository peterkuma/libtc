/*
 * misc.c
 *
 * Miscellaneous functions.
 *
 */

#define _BSD_SOURCE

#include <math.h>
#include <assert.h>
#include <strings.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include "misc.h"

/* GNU Scientific Library RNG. */
gsl_rng *rng = NULL;

/*
 * Log-Beta function (Stirling approx.).
 */
double
log_beta(double a, double b)
{
    return 0.5*log(2*M_PI) +
        (a - 0.5)*log(a) +
        (b - 0.5)*log(b) -
        (a + b - 0.5)*log(a + b);
}

/*
 * Generate a floating-point pseudorandom number from the interval <0, 1).
 */
double
frand(void)
{
    return rand()/(RAND_MAX + 1.0);
}

/*
 * Generate a floating-point pseudorandom number from the interval (0, 1).
 */
double
frand1(void)
{
    return (rand() + 1.)/(RAND_MAX + 2.);
}

/*
 * Sample from `n` possible outcomes with probabilities `p`. If `p` is NULL,
 * all outcomes are assumed to have equal probability. `n` has to be greater
 * than 0.
 */
size_t
sample(size_t n, const double p[])
{
    double sum = 0;
    size_t i = 0;

    assert(n > 0);

    if (p != NULL) {
        for (i = 0; i < n; i++)
            sum += p[i];
    }

    double f = 0, r = frand();
    if (p == NULL) return r*n;
    for (i = 0; i < n; i++) {
        f += p[i]/sum;
        if (r < f) return i;
    }
    assert(0);
}

/*
 * Initialize the GNU Scientific Library.
 */
void
init_gsl(void)
{
    const gsl_rng_type *T;
    gsl_rng_env_setup();
    T = gsl_rng_default;
    rng = gsl_rng_alloc(T);
}

/*
 * Deinitialize the GNU Scientific Library.
 */
void
deinit_gsl(void)
{
   gsl_rng_free(rng);
}


/*
 * Generate a random number from the truncated normal distribution defined
 * by mean `mean`, standard deviation `sd`, and bounds `a`, `b`.
 * The generated number is in the open interval (a, b).
 *
 * Note: This is a trivial implementation of rejection sampling,
 * which could be performing very poorly depending on the choice of bounds.
 *
 * TODO: Should be replaced by a better implementation.
 */
double
rtnorm(double mean, double sd, double a, double b)
{
    double x = 0;
    x = mean + gsl_ran_gaussian(rng, sd);
    while (!(x > a && x < b))
        x = mean + gsl_ran_gaussian(rng, sd);
    return x;
}

/*
 * Insert `element` into `array`. `n` is the length of array, `i`
 * is the insert position, and `size` is the size of array elements in bytes.
 * Returns a newly-allocated array with the element inserted,
 * or NULL on failure. It is the responsibility of the callee to free the array.
 */
void *
array_insert(
    const void *array,
    size_t n,
    const void *element,
    size_t i,
    size_t size
) {
    void *new_array = NULL;
    void *p = NULL;
    new_array = calloc(n + 1, size);
    if (new_array == NULL) return NULL;
    p = new_array;
    bcopy(array, p, size*i);
    p += size*i;
    bcopy(element, p, size);
    p += size;
    bcopy(array + size*i, p, size*(n - i));
    return new_array;
}

void *
array_remove(const void *array, size_t n, size_t i, size_t size)
{
    void *new_array = NULL;
    void *p = NULL;
    assert(n > 0);
    assert(i >= 0);
    assert(i < n);
    new_array = calloc(n - 1, size);
    if (new_array == NULL) return NULL;
    p = new_array;
    bcopy(array, p, size*i);
    p += size*i;
    bcopy(array + size*(i+1), p, size*(n - i - 1));
    return new_array;
}
