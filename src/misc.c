/*
 * misc.c
 *
 * Miscellaneous functions.
 *
 */

#define _BSD_SOURCE

#include <math.h>
#include <assert.h>
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
 * all outcomes are assumed to have equal probability.
 */
size_t
sample(size_t n, const double p[])
{
    double sum = 0;
    size_t i = 0;

    if (p != NULL) {
        for (i = 0; i < n; i++)
            sum += p[i];
        assert(sum == 1);
    }

    double f = 0, r = frand();
    if (p == NULL) return r*n;
    for (i = 0; i < n; i++) {
        f += p[i];
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
