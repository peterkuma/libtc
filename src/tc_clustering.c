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
#include <mcheck.h>

#include "misc.h"
#include "tree.h"
#include "tc.h"

struct tc_opts tc_default_opts = {
    .burnin = 0,
    .nsamples = 10,
    .split_p = 0.1,
    .merge_p = 0.1,
    .move_p = 0.8,
    .move_sd_frac = 0.1
};

enum action {
    MOVE,
    SPLIT,
    MERGE
};

/*
 * Check validity of options. Returns true if correct, false if incorrect.
 */
static bool
check_opts(const struct tc_opts *opts)
{
    bool cond;

    cond = opts->merge_p + opts->split_p + opts->move_p == 1;
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
    size_t i = 0, j = 0, k = 0;
    size_t S = 0, s = 0;
    size_t SS = 0, ss = 0;
    size_t C = 0, c = 0;
    double w1 = 0, w2 = 0;
    struct tc_node *node = NULL, *parent = NULL;
    struct tc_node *old_node = NULL, *new_node = NULL;
    const struct tc_param_def *pd = NULL;
    struct tc_range range;
    double *cuts = NULL;
    double cut;
    double new_cut;

    errno = 0;
    mtrace();

    if (!check_opts(opts)) {
        errno = EINVAL;
        goto cleanup;
    }

    init_gsl();

    tree = tc_new_tree(1000024, param_def, K);
    if (tree == NULL) {
        errno = ENOMEM;
        goto cleanup;
    }

    /* DEBUG: Create a dummy node. */
    // node = tc_new_node(tree, 0, 2, (double[]) { 2 });
    // tc_replace_node(tree->root, node);
    // node = NULL;

    // tc_dump_tree_simple(tree, NULL);

    l = tc_log_likelihood(tree, ds, N);
    // tc_dump_segments_json(tree);

    cb(tree, l, ds, N, cb_data);
    tc_dump_tree_simple(tree, NULL);

    nsamples = 0;
    while (nsamples < opts->burnin + opts->nsamples) {
        assert(check_tree(tree));

        tc_dump_tree_simple(tree, NULL);

        action = sample(3, (double[]){
            opts->move_p,
            opts->split_p,
            opts->merge_p
        });

        switch (action) {
        case SPLIT:
            S = count_segments(tree);
            s = sample(S, NULL);
            node = select_segment(tree, s);
            k = sample(K, NULL);
            parent = node->parent;

            debug("k = %zu\n", k);
            debug("s = %zu\n", s);

            if (parent != NULL && k == parent->param) {
                debug("split-1\n");
                pd = &param_def[k];
                i = find_child(parent, node);
                node_range(node, k, &range);
                cut = range.min + frand1()*(range.max - range.min);
                free_range(&range);
                cuts = array_insert(
                    parent->cuts,
                    parent->ncuts,
                    &cut,
                    i,
                    sizeof(cut)
                );
                if (cuts == NULL) {
                    errno = ENOMEM;
                    goto cleanup;
                }
                new_node = tc_new_node(
                    tree,
                    parent->param,
                    parent->nchildren + 1,
                    cuts
                );
                free(cuts);
                cuts = NULL;
                if (new_node == NULL) {
                    errno = ENOMEM;
                    goto cleanup;
                }
                tc_replace_node(parent, new_node);
                old_node = parent;
                for (j = 0; j < i; j++)
                    tc_replace_node(new_node->children[j], parent->children[j]);
                for (j = i; j < node->nchildren; j++)
                    tc_replace_node(new_node->children[j+1], parent->children[j]);
                assert(check_tree(tree));
            } else {
                debug("split-2\n");
                node_range(node, k, &range);
                cut = range.min + frand1()*(range.max - range.min);
                new_node = tc_new_node(tree, k, 2, (double[]){ cut });
                if (new_node == NULL) {
                    errno = ENOMEM;
                    goto cleanup;
                }
                tc_replace_node(node, new_node);
                old_node = node;
                assert(check_tree(tree));
            }

            lx = tc_log_likelihood(tree, ds, N);
            p = fmin(1, exp(lx - l));
            accept = sample(2, (double[]){1-p, p});
            if (accept) {
                debug("SPLIT\n");
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
                            segments[s].ranges[0].min,
                            segments[s].ranges[0].max,
                            segments[s].ranges[1].min,
                            segments[s].ranges[1].max
                        );
                    }
                    tc_free_segments(segments, S);
                    free(segments);
                    segments = NULL;

                }
            } else {
                tc_replace_node(new_node, old_node);
                for (i = 0; i < old_node->nchildren; i++)
                    old_node->children[i]->parent = old_node;
                assert(check_tree(tree));
                // tree_free_node(new_node);
            }
            break;
        case MERGE:
            SS = count_supersegments(tree);
            if (SS == 0) continue;
            ss = sample(SS, NULL);
            node = select_supersegment(tree, ss);
            pd = &tree->param_def[node->param];
            if (pd->type == TC_METRIC) {
                C = count_movable_cuts(node);
                c = sample(C, NULL);
                i = select_movable_cut(node, c);
                cuts = array_remove(node->cuts, node->ncuts, i, sizeof(cut));
                if (cuts == NULL) {
                    errno = ENOMEM;
                    goto cleanup;
                }
                new_node = tc_new_node(tree, node->param, node->nchildren - 1, cuts);
                if (new_node == NULL) {
                    errno = ENOMEM;
                    goto cleanup;
                }
                free(cuts);
                cuts = NULL;
                tc_replace_node(node, new_node);
                old_node = node;
                for (j = 0; j < new_node->nchildren; j++) {
                    tc_replace_node(
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
                    debug("MERGE\n");
                    l = lx;
                    if (++nsamples > opts->burnin) {
                        cb(tree, l, ds, N, cb_data);
                        tc_dump_tree_simple(tree, NULL);
                    }
                } else {
                    tc_replace_node(new_node, old_node);
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
                C = count_movable_cuts(node);
                c = sample(C, NULL);
                i = select_movable_cut(node, c);
                cut = node->cuts[i];
                node_range(node->children[i], node->param, &range);
                w1 = range.max - range.min;
                free_range(&range);
                node_range(node->children[i+1], node->param, &range);
                w2 = range.max - range.min;
                free_range(&range);

                debug("k = %zu\n", node->param);
                debug("w1 = %lf, w2 = %lf\n", w1, w2);
                new_cut = cut + rtnorm(0, (w1 + w2)*opts->move_sd_frac, -w1, w2);
                node->cuts[i] = new_cut;
                lx = tc_log_likelihood(tree, ds, N);
                p = fmin(1, exp(lx - l));
                accept = sample(2, (double[]){1-p, p});
                if (accept) {
                    debug("MOVE\n");
                    l = lx;
                    if (++nsamples > opts->burnin) {
                        cb(tree, l, ds, N, cb_data);
                        tc_dump_tree_simple(tree, NULL);
                    }
                } else {
                    node->cuts[i] = cut;
                }
            } else if (pd->type == TC_NOMINAL) {
                /* Not implemented. */
            } else assert(0);
            break;
        default: assert(0);
        }
    }

cleanup:
    if (tree != NULL) free(tree);
    if (cuts != NULL) free(cuts);
    deinit_gsl();
    muntrace();
    return errno != 0 ? -1 : 0;
}
