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
#include <errno.h>

#include "misc.h"
#include "tree.h"
#include "tc.h"

struct tc_opts tc_default_opts = {
    .burnin = 0,
    .nsamples = 10,
    .split1_p = 0.1,
    .split2_p = 0,
    .merge_p = 0.1,
    .move_p = 0.8,
    .move_sd_frac = 0.1
};

enum action {
    MOVE,
    SPLIT1,
    SPLIT2,
    MERGE
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

    if (!check_opts(opts)) {
        errno = EINVAL;
        return -1;
    }

    init_gsl();

    tree = tc_new_tree(1000024, param_def, K);
    if (tree == NULL) {
        errno = ENOMEM;
        return -1;
    }

    /* DEBUG: Create a dummy node. */
    // node = tc_new_node(tree, 0, 2, (double[]) { 2 });
    // tc_replace_node(tree, tree->root, node);
    // node = NULL;

    // tc_dump_tree_simple(tree, NULL);

    l = tc_log_likelihood(tree, ds, N);
    // tc_dump_segments_json(tree);

    cb(tree, l, ds, N, cb_data);
    tc_dump_tree_simple(tree, NULL);

    nsamples = 0;
    while (nsamples < opts->burnin + opts->nsamples) {
        size_t i = 0, j = 0, k = 0;
        size_t S = 0, s = 0;
        size_t SS = 0, ss = 0;
        size_t C = 0, c = 0;
        double w1 = 0, w2 = 0;
        struct tc_node *node = NULL, *parent = NULL;
        struct tc_node *old_node = NULL, *new_node = NULL;
        const struct tc_param_def *pd = NULL;
        struct tc_range range;
        union tc_valuep part;
        union tc_value cut;
        union tc_value new_cut;

        assert(check_tree(tree));

        tc_dump_tree_simple(tree, NULL);

        action = sample(4, (double[]){
            opts->move_p,
            opts->split1_p,
            opts->split2_p,
            opts->merge_p
        });

        switch (action) {
        case SPLIT1:
            S = count_segments(tree);
            s = sample(S, NULL);
            node = select_segment(tree, s);
            k = sample(K, NULL);
            parent = node->parent;

            debug("k = %zu\n", k);
            debug("s = %zu\n", s);

            if (parent != NULL && k == parent->param) {
                debug("split1-1\n");
                pd = &param_def[k];
                i = find_child(parent, node);
                node_range(tree, node, k, &range);
                cut.float64 = range.min.float64 + frand1()*(range.max.float64 - range.min.float64);
                free_range(&range);
                part.buf = array_insert(
                    parent->part.buf,
                    parent->nchildren - 1,
                    &cut.float64,
                    i,
                    TC_SIZE[pd->size]
                );
                if (part.buf == NULL) {
                    errno = ENOMEM;
                    return -1;
                }
                new_node = tc_new_node(
                    tree,
                    parent->param,
                    parent->nchildren + 1,
                    part.buf
                );
                if (new_node == NULL) {
                    errno = ENOMEM;
                    return -1;
                }
                tc_replace_node(tree, parent, new_node);
                old_node = parent;
                for (j = 0; j < i; j++)
                    tc_replace_node(tree, new_node->children[j], parent->children[j]);
                for (j = i; j < node->nchildren; j++)
                    tc_replace_node(tree, new_node->children[j+1], parent->children[j]);
                assert(check_tree(tree));
            } else {
                debug("split1-2\n");
                node_range(tree, node, k, &range);
                pd = &tree->param_def[k];
                // debug("%lf, %lf\n", range.min.float64, range.max.float64);
                part = rand_part(&range, pd);
                free_range(&range);
                if (part.buf == NULL) {
                    errno = ENOMEM;
                    return -1;
                }
                new_node = tc_new_node(tree, k, 2, part.buf);
                free(part.buf);
                part.buf = NULL;
                if (new_node == NULL) {
                    errno = ENOMEM;
                    return -1;
                }
                tc_replace_node(tree, node, new_node);
                old_node = node;
                assert(check_tree(tree));
            }

            lx = tc_log_likelihood(tree, ds, N);
            p = fmin(1, exp(lx - l));
            accept = sample(2, (double[]){1-p, p});
            if (accept) {
                debug("SPLIT1\n");
                l = lx;
                if (++nsamples > opts->burnin) {
                    cb(tree, l, ds, N, cb_data);
                    tc_dump_tree_simple(tree, NULL);

                    size_t S = 0;
                    struct tc_segment *segments = tc_segments(tree, ds, N, &S);
                    for (size_t s = 0; s < S; s++) {
                        printf("s = %zu, NX[s] = %zu, V[s] = %lf, ((%lf, %lf),(%lf, %lf))\n",
                            s,
                            segments[s].NX,
                            segments[s].V,
                            segments[s].ranges[0].min.float64,
                            segments[s].ranges[0].max.float64,
                            segments[s].ranges[1].min.float64,
                            segments[s].ranges[1].max.float64
                        );
                    }
                    tc_free_segments(segments, S);
                    free(segments);
                    segments = NULL;

                }
            } else {
                tc_replace_node(tree, new_node, old_node);
                for (i = 0; i < old_node->nchildren; i++)
                    old_node->children[i]->parent = old_node;
                assert(check_tree(tree));
                // tree_free_node(new_node);
            }

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
            debug("SPLIT2\n");

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
                l = lx;
                if (++nsamples > opts->burnin) {
                    cb(tree, l, ds, N, cb_data);
                    tc_dump_tree_simple(tree, NULL);
                }
            } else {
                tc_replace_node(tree, new_node, node);
                // tree_free_node(new_node);
            }
            break;
        case MERGE:
            debug("MERGE\n");
            SS = count_supersegments(tree);
            if (SS == 0) continue;
            ss = sample(SS, NULL);
            node = select_supersegment(tree, ss);
            pd = &tree->param_def[node->param];
            if (pd->type == TC_METRIC) {
                C = count_movable_cuts(tree, node);
                c = sample(C, NULL);
                i = select_movable_cut(tree, node, c);
                part.float64 = array_remove(
                    node->part.float64,
                    node->nchildren,
                    i,
                    TC_SIZE[pd->size]
                );
                if (part.float64 == NULL) {
                    errno = ENOMEM;
                    return -1;
                }
                new_node = tc_new_node(tree, node->param, node->nchildren - 1, part.buf);
                tc_replace_node(tree, node, new_node);
                old_node = node;
                for (j = 0; j < new_node->nchildren; j++) {
                    tc_replace_node(
                        tree,
                        new_node->children[i],
                        node->children[j < i ? j : j + 1]
                    );
                }
                lx = tc_log_likelihood(tree, ds, N);
                p = fmin(1, exp(lx - l));
                accept = sample(2, (double[]){1-p, p});
                tc_dump_tree_simple(tree, NULL);
                debug("l = %lf, lx = %lf, p = %lf\n", l, lx, p);
                if (accept) {
                    debug("accept\n");
                    l = lx;
                    if (++nsamples > opts->burnin) {
                        cb(tree, l, ds, N, cb_data);
                        tc_dump_tree_simple(tree, NULL);
                    }
                } else {
                    tc_replace_node(tree, new_node, old_node);
                    for (j = 0; j < node->nchildren; j++)
                        old_node->children[j]->parent = old_node;
                    // tree_free_node(new_node);
                }
            } else if (pd->type == TC_NOMINAL) {
                /* Not implemented. */
            } else assert(0);
            break;
        case MOVE:
            SS = count_supersegments(tree);
            if (SS == 0) continue;
            // debug("SS = %zu\n", SS);
            ss = sample(SS, NULL);
            node = select_supersegment(tree, ss);
            pd = &param_def[node->param];
            if (pd->type == TC_METRIC) {
                C = count_movable_cuts(tree, node);
                c = sample(C, NULL);
                i = select_movable_cut(tree, node, c);
                cut.float64 = node->part.float64[i];
                node_range(tree, node->children[i], node->param, &range);
                w1 = range.max.float64 - range.min.float64;
                free_range(&range);
                node_range(tree, node->children[i+1], node->param, &range);
                w2 = range.max.float64 - range.min.float64;
                free_range(&range);

                debug("k = %zu\n", node->param);
                debug("w1 = %lf, w2 = %lf\n", w1, w2);
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
                    debug("MOVE\n");
                    l = lx;
                    if (++nsamples > opts->burnin) {
                        cb(tree, l, ds, N, cb_data);
                        tc_dump_tree_simple(tree, NULL);
                    }
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
    }

    return 0;

error:
    if (tree != NULL) free(tree);
    return 0;
}
