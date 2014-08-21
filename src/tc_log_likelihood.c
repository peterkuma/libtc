/*
 * tc_log_likelihood.c
 *
 * tc_log_likelihood implementation.
 *
 */

#include <math.h>
#include <stdlib.h>

#include "misc.h"
#include "tree.h"
#include "tc.h"

double
tc_log_likelihood(
    const struct tc_tree *tree,
    const void *ds[],
    size_t N
) {
    double l1 = 0, l2 = 0, l = 0; /* Likelihood. */
    size_t s = 0;
    size_t S = 0; /* Number of segments. */
    struct tc_segment *segment;
    struct tc_segment *segments = NULL;

    segments = tc_segments(tree, ds, N, &S);
    if (segments == NULL)
        return NAN;

    /*
     * Calculate the log-likelihood.
     * Part I: contribution of volumes.
     */
    l1 = 0;
    for (s = 0; s < S; s++) {
        segment = &segments[s];
        if (segment->NX == 0)
            continue;
        if (segment->V == 0)
            continue;
        l1 -= segment->NX*log(segment->V);
    }
    // debug("l1 = %lf\n", l1);

    /*
     * Calculate the log-likelihood
     * Part II: contribution of numbers of elements.
     */
    double a = 0, b = 0;
    l2 = 0;
    for (s = S-1; s != -1; s--) {
        a = segments[s].NX + 1;
        b += (s + 1 < S ? segments[s+1].NX : 0) + 1;
        // debug("a = %lf, b = %lf\n", a, b);
        l2 += log_beta(a, b);
    }
    // debug("l2 = %lf\n", l2);

    l = l1 + l2;

cleanup:
    if (segments != NULL) {
        tc_free_segments(segments, S);
        free(segments);
        segments = NULL;
    }
    return l;
}
