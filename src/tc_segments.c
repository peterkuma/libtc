/*
 * tc_segments.c
 *
 * tc_segments implementation.
 *
 */

#include <stdlib.h>
#include <assert.h>
#include <errno.h>
#include <math.h>

#include "misc.h"
#include "tc.h"
#include "tree.h"

struct tc_segment *
tc_segments(
    const struct tc_tree *tree,
    const void *ds[],
    size_t N,
    size_t *S
) {
    size_t n = 0, s = 0, k = 0, i = 0;
    struct tc_node *node = NULL;
    const struct tc_param_def *pd = NULL;
    union tc_valuep data;
    struct tc_range *range = NULL;
    struct tc_segment *segments = NULL;
    struct tc_segment *segment = NULL;
    double *p = NULL;
    double V = 0;

    *S = count_segments(tree);
    segments = calloc(*S, sizeof(struct tc_segment));
    if (segments == NULL) {
        errno = ENOMEM;
        return NULL;
    }

    for (s = 0; s < *S; s++)
        init_segment(&segments[s], tree->K);

    s = 0;
    for (node = tree->first; node != NULL; node = node->next) {
        if (!is_segment(node))
            continue;
        segment = &segments[s];
        node->_aux = segment;
        /* Determine volume of the segment. */
        segment->V = 1.;
        for (k = 0; k < tree->K; k++) {
            range = &segment->ranges[k];
            pd = &tree->param_def[k];
            node_range(node, k, range);
            V = range->max - range->min;
            segment->V *= V > 0 ? V : 1;
        }
        s++;
    }

    for (n = 0; n < N; n++) {
        /*
         * Find element in tree.
         */
        node = tree->root;
        while (node != NULL) {
            if (is_segment(node)) {
                segment = (struct tc_segment *) node->_aux;
                segment->NX++;
                break;
            }

            pd = &tree->param_def[node->param];
            data.buf = ds[node->param];

            if (pd->type == TC_METRIC) {
                if (isnan(data.float64[n])) {
                    /* Assign element randomly according to size of children. */
                    p = calloc(node->nchildren, sizeof(double));
                    if (p == NULL) {
                        errno = ENOMEM;
                        return NULL;
                    }
                    for (i = 0; i < node->nchildren; i++) {
                        struct tc_range range;
                        node_range(node->children[i], node->param, &range);
                        p[i] = range.max - range.min;
                        free_range(&range);
                    }
                    i = sample(node->nchildren, p);
                    free(p);
                    p = NULL;
                } else {
                    for (i = 0; i  < node->ncuts; i++) {
                        if (data.float64[n] <= node->cuts[i])
                            break;
                    }
                }
            } else if (pd->type == TC_NOMINAL) {
                i = node->categories[data.int64[n]];
            } else {
                assert(0);
            }
            node = node->children[i];
        }
    }

    return segments;
}
