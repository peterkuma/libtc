/*
 * tc_segments.c
 *
 * tc_segments implementation.
 *
 */

#include <stdlib.h>
#include <assert.h>

#include "misc.h"
#include "tc.h"
#include "tree.h"

#define IS_INT64(pd) ((pd)->size == TC_INT64)
#define IS_FLOAT64(pd) ((pd)->size == TC_FLOAT64)

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

    *S = count_segments(tree);
    segments = calloc(*S, sizeof(struct tc_segment));
    for (s = 0; s < *S; s++)
        init_segment(&segments[s], tree->K);

    s = 0;
    for (node = tree->first; node != NULL; node = node->next) {
        if (node->nchildren != 0) continue; /* Not a segment. */
        segment = &segments[s];
        node->_aux = segment;
        /* Determine volume of the segment. */
        segment->V = 1.;
        for (k = 0; k < tree->K; k++) {
            range = &segment->ranges[k];
            pd = &tree->param_def[k];
            node_range(tree, node, k, range);
            if (IS_INT64(pd))
                segment->V *= range->max.int64 - range->min.int64;
            else if (IS_FLOAT64(pd))
                segment->V *= range->max.float64 - range->min.float64;
            else assert(0);
        }
        s++;
    }

    for (n = 0; n < N; n++) {
        /*
         * Find element in tree.
         */
        node = tree->root;
        while (node != NULL) {
            pd = &tree->param_def[node->param];
            data.buf = ds[node->param];

            if (node->nchildren == 0) {
                segment = (struct tc_segment *) node->_aux;
                segment->NX++;
                break;
            }

            if (pd->type == TC_METRIC) {
                size_t ncuts = node->nchildren - 1;
                for (i = 0; i < ncuts; i++) {
                    if (IS_INT64(pd)) {
                        if (data.int64[n] <= node->part.int64[i]) break;
                    } else if (IS_FLOAT64(pd)) {
                        if (data.float64[n] <= node->part.float64[i]) break;
                    } else assert(0);
                }
            } else if (pd->type == TC_NOMINAL) {
                i = node->part.int64[data.int64[n]];
            } else {
                assert(0);
            }
            node = node->children[i];
        }
    }

    return segments;
}
