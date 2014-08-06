/*
 * tc_clustering.c
 *
 * tc_clustering implementation.
 *
 */

#include <stdlib.h>
#include <stdbool.h>
#include <assert.h>
#include <err.h>
#include <math.h>

#include "misc.h"
#include "tree.h"
#include "tc.h"

struct tc_opts tc_default_opts = {
    .burnin = 1000,
    .nsamples = 1000,
    .split1_p = 0.1,
    .split2_p = 0.1,
    .merge_p = 0.1,
    .move_p = 0.7,
    .move_sd_frac = 0.1
};

enum action {
    SPLIT1,
    SPLIT2,
    MERGE,
    MOVE
};

/*
 * Check validity of options. Returns true if correct, false if incorrect.
 */
static bool
check_opts(const struct tc_opts *opts)
{
    bool cond;

    cond = opts->merge_p + opts->split1_p + opts->split2_p + opts->move_p == 1;
    if (!cond) return false;

    return true;
}

int
tc_clustering(
    const void *ds[],
    size_t N,
    const struct tc_param_def param_def[],
    size_t K,
    tc_clustering_cb cb,
    void *cb_data,
    const struct tc_opts *opts
) {
    struct tc_tree *tree = NULL;
    enum action action = 0; /* Action to take. */
    double l = 0; /* Log-likelihood. */
    double lx = 0; /* Proposal log-likelihood. */
    double p = 0; /* Acceptance probability. */
    bool accept = false; /* Accept proposal? */
    size_t nsamples = 0;

    if (!check_opts(opts))
        return 1;

    init_gsl();

    tree = tc_new_tree(10024, param_def, K);
    if (tree == NULL) return 1;

    /* DEBUG: Create a dummy node. */
    // node = tc_new_node(tree, 0, 2, (double[]) { 2 });
    // tc_replace_node(tree, tree->root, node);
    // node = NULL;

    // tc_dump_tree_simple(tree, NULL);

    l = tc_log_likelihood(tree, ds, N);
    // tc_dump_segments_json(tree);

    cb(tree, ds, N, cb_data);

    nsamples = 0;
    while (nsamples < opts->burnin + opts->nsamples) {
        action = sample(4, (double[]){
            opts->merge_p,
            opts->split1_p,
            opts->split2_p,
            opts->move_p
        });
        switch (action) {
        case SPLIT1:
            // debug("SPLIT1\n");
            // SS = count_supersegments(tree);
            // ss = sample(SS);
            // node = select_supersegment(tree, ss);
            // i0 = sample(node->nchildren);
            // n0 = node->children[i0];
            // node_range(tree, n0, n0->param, &range);
            // part = rand_part(&range);
            // free_range(&range);
            // if (part == NULL) continue;
            // part =
            // new_node = tc_new_node(tree, node->param, node->nchildren+1, part);
            // if (new_node == NULL) continue;
            // tc_replace_node(node, new_node);
            // lx = tc_log_likelihood(tree, ds, N);
            // p = fmin(1, exp(lx - l));
            // accept = sample(2, (double[]){1-p, p});
            // if (!accept) {
            //     tc_replace_node(new_node, node);
            // }

//          pd = param_def[node->param];
//
//          bcopy(node->part, part, sz*k0);
//          if (IS_INT64(pd)) part.int64[k] = cut.int64;
//          else if (IS_FLOAT64(pd)) part.float64[k] = cut.float64;
//          else assert(0);
//          bcopy(node->part, part.buf + sz*(k0+1), sz*(node->nchildren - k0 - 1));
//
//          /* Perform the split. */
//          size_t sz = TC_SIZE[pd->size];
//          part.buf = calloc(node->nchildren, sz);
//          if (part.buf == NULL) return 1;
//          bcopy(node->part, part, sz*k0);
//          if (IS_INT64(pd)) part.int64[k] = cut.int64;
//          else if (IS_FLOAT64(pd)) part.float64[k] = cut.float64;
//          else assert(0);
//          bcopy(node->part, part.buf + sz*(k0+1), sz*(node->nchildren - k0 - 1));
//
//          new_node = tc_new_node(tree, node->param, node->nchildren+1, part);
//          bcopy(
//              node->children,
//              new_node->children,
//              sizeof(struct node *)*k0
//          );
//          node->children[k0] = tc_new_node(tree, 0, 0, NULL);
//          node->children[k0].parent = node;
//          bcopy(
//              node->children + k0,
//              new_node->children + k0 + 1,
//              new_node->nchildren - k0 - 1
//          );
            break;
        case SPLIT2:
            ;
            size_t S = 0;
            size_t s = 0;
            size_t k = 0;
            const struct tc_param_def *pd = NULL;
            struct tc_range range;
            struct tc_node *node = NULL, *new_node = NULL;
            union tc_valuep part;

            if (K == 1) continue;
            S = count_segments(tree);
            s = sample(S, NULL);
            node = select_segment(tree, s);
            if (node->parent) {
                k = sample(K-1, NULL);
                if (k >= node->parent->param) k++;
            } else{
                k = sample(K, NULL);
            }
            node_range(tree, node, k, &range);
            pd = &tree->param_def[node->param];
            // debug("%lf, %lf\n", range.min.float64, range.max.float64);
            part = rand_part(&range, pd);
            free_range(&range);
            if (part.buf == NULL) continue;
            new_node = tc_new_node(tree, k, 2, part.buf);
            free(part.buf);
            if (new_node == NULL) continue;
            part.buf = NULL;
            tc_replace_node(tree, node, new_node);
            lx = tc_log_likelihood(tree, ds, N);
            p = fmin(1, exp(lx - l));
            // debug("l = %lf, lx = %lf, p = %lf\n", l, lx, p);
            accept = sample(2, (double[]){1-p, p});
            if (accept) {
                // debug("SPLIT2\n");
                if (++nsamples > opts->burnin)
                    cb(tree, ds, N, cb_data);
                // size_t S = 0;
                // struct tc_segment *segments = tc_segments(tree, ds, N, &S);
                // for (size_t s = 0; s < S; s++) {
                //     debug("s = %zu, NX[s] = %zu, V[s] = %lf, ((%lf, %lf),(%lf, %lf))\n",
                //         s,
                //         segments[s].NX,
                //         segments[s].V,
                //         segments[s].ranges[0].min.float64,
                //         segments[s].ranges[0].max.float64,
                //         segments[s].ranges[1].min.float64,
                //         segments[s].ranges[1].max.float64
                //     );
                // }
                // tc_free_segments(segments, S);
                // free(segments);
                // segments = NULL;

                // // tc_dump_tree_simple(tree, NULL);
                // tc_dump_segments_json(tree);
                l = lx;
            } else {
                tc_replace_node(tree, new_node, node);
                // tree_free_node(new_node);
            }
            break;
        case MERGE:
            // debug("MERGE\n");
            // SS = count_supersegments(tree);
            // ss = sample(SS);
            // node = select_supersegment(tree, ss);
            // if (pd->type == TC_METRIC)
            //     assert(node->nchildren > 0);
            //     i = sample(node->nchildren - 1);
            //     part =
            //     new_node = tc_new_node(tree, node->param, node->nchildren - 1, part);
            //     if (new_node == NULL) continue;
            //     tc_replace_node(node, new_node);
            //     lx = tc_log_likelihood(tree, ds, N);
            //     p = fmin(1, exp(lx - l));
            //     accept = sample(2, (double[]){1-p, p});
            //     if (!accept) {
            //         tc_replace_node(new_node, node);
            //         tree_free_node(new_node);
            //     }
            // } else if (pd->type == TC_NOMINAL) {
            //     /* Not implemented. */
            // } else assert(0);
            break;
        case MOVE:
            ;
            double w1 = 0, w2 = 0;
            union tc_value cut;
            union tc_value new_cut;
            double lx = 0;
            double p = 0;
            size_t SS = 0, ss = 0;
            size_t i = 0;

            SS = count_supersegments(tree);
            if (SS == 0) continue;
            // debug("SS = %zu\n", SS);
            ss = sample(SS, NULL);
            node = select_supersegment(tree, ss);
            pd = &param_def[node->param];
            if (pd->type == TC_METRIC) {
                assert(node->nchildren > 0);
                i = sample(node->nchildren - 1, NULL);
                if (IS_INT64(pd))
                    cut.int64 = node->part.int64[i];
                else if(IS_FLOAT64(pd))
                    cut.float64 = node->part.float64[i];
                else assert(0);
                w1 = node_width(tree, node->children[i]);
                w2 = node_width(tree, node->children[i+1]);
                new_cut.float64 = cut.float64 + rtnorm(0, (w1 + w2)*opts->move_sd_frac, -w1, w2);
                // debug("cut = %lf, new_cut = %lf\n", cut.float64, new_cut.float64);
                if (IS_INT64(pd)) /* Not implemented. */
                    node->part.int64[i] = new_cut.float64;
                else if (IS_FLOAT64(pd))
                    node->part.float64[i] = new_cut.float64;
                else assert(0);
                lx = tc_log_likelihood(tree, ds, N);

                p = fmin(1, exp(lx - l));
                // debug("l = %lf, lx = %lf, p = %lf\n", l, lx, p);
                accept = sample(2, (double[]){1-p, p});
                if (accept) {
                    // debug("MOVE\n");
                    // debug("%lf\n", lx);
                    if (++nsamples > opts->burnin)
                        cb(tree, ds, N, cb_data);

                    // size_t S = 0;
                    // struct tc_segment *segments = tc_segments(tree, ds, N, &S);
                    // for (size_t s = 0; s < S; s++) {
                    //     debug("s = %zu, NX[s] = %zu, V[s] = %lf, ((%lf, %lf),(%lf, %lf))\n",
                    //         s,
                    //         segments[s].NX,
                    //         segments[s].V,
                    //         segments[s].ranges[0].min.float64,
                    //         segments[s].ranges[0].max.float64,
                    //         segments[s].ranges[1].min.float64,
                    //         segments[s].ranges[1].max.float64
                    //     );
                    // }
                    // tc_free_segments(segments, S);
                    // free(segments);
                    // segments = NULL;

                    // tc_dump_segments_json(tree);

                    // tc_dump_tree_simple(tree, NULL);
                    // debug("accept\n");
                    l = lx;
                } else {
                    if (IS_INT64(pd)) /* Not implemented. */
                        node->part.int64[i] = cut.float64;
                    else if (IS_FLOAT64(pd))
                        node->part.float64[i] = cut.float64;
                    else assert(0);
                }
            } else if (pd->type == TC_NOMINAL) {
                /* Not implemented. */
            } else assert(0);
            break;
        default: assert(0);
        }
        // tc_dump_tree_simple(tree, NULL);
    }

    return 0;

error:
    if (tree != NULL) free(tree);
    return 0;
}
