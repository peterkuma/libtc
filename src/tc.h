/*
 * tc.h
 *
 * Tree Clustering.
 *
 * Implementation of Cluster Analysis using Decision Trees
 * by the Metropolis-Hastings algorithm.
 *
 * See README.md for API reference.
 *
 */

#ifndef TC_H
#define TC_H

#include <stddef.h>
#include <inttypes.h>

extern size_t TC_SIZE[];

enum tc_param_type {
    TC_METRIC,
    TC_NOMINAL
};

enum tc_param_size {
    TC_FLOAT64,
    TC_INT64
};

union tc_value {
    int64_t int64;
    double float64;
};

union tc_valuep {
    int64_t *int64;
    double *float64;
    uint8_t *buf;
};

struct tc_param_def {
    enum tc_param_type type; /* Parameter types. */
    enum tc_param_size size; /* Parameter sizes. */
    union tc_value min; /* Minimum parameter value. */
    union tc_value max; /* Maximum parameter value. */
};

struct tc_tree {
    const struct tc_param_def *param_def; /* Parameter definitions. */
    size_t K; /* Number of parameters. */
    size_t size; /* Size of buffer in bytes. */
    uint8_t *p; /* Pointer to the end of the buffer. */
    struct tc_node *first; /* First node (for sequential traversing). */
    struct tc_node *last; /* Last node (for sequential traversing). */
    struct tc_node *root; /* Root node. */
    uint8_t buf[]; /* Buffer in which nodes are stored. */
};

struct tc_node {
    struct tc_node *parent; /* Parent node. */
    size_t nchildren; /* Number of child nodes. */
    struct tc_tree *tree; /* Reference to node's tree. */
    size_t param; /* Parameter number. */
    double *cuts; /* Cuts. */
    size_t ncuts; /* Number of cuts. */
    int64_t *categories; /* Categories partitioning. */
    size_t ncategories; /* Number of categories. */
    struct tc_node **children; /* Child nodes. */
    struct tc_node *next; /* Next node (for sequential traversing). */
    struct tc_node *prev; /* Previous node (for sequential traversing). */
    void *_aux; /* Auxillary data for any purpose. */
};

struct tc_opts {
    size_t burnin; /* Number of burn-in samples. */
    size_t nsamples; /* Number of samples to generate (excl. burn-in). */
    double split_p; /* Probability of split. */
    double merge_p; /* Probability of merge. */
    double move_p; /* Probability of move. */
    double move_sd_frac; /* Move standard deviation as a fraction. */
};

extern struct tc_opts tc_default_opts;

struct tc_range {
    double min;
    double max;
    int64_t *categories;
    size_t ncategories;
};

struct tc_segment {
    size_t NX;
    size_t K;
    double V;
    struct tc_range *ranges;
};

typedef void tc_clustering_cb(
    const struct tc_tree *tree,
    double l,
    const void **ds,
    size_t N,
    void *data
);

void
tc_param_def_init(
    struct tc_param_def *pd,
    const void *data,
    size_t N
);

struct tc_tree *
tc_new_tree(size_t size, const struct tc_param_def *param_def, size_t K);

struct tc_node *tc_new_leaf(struct tc_tree *tree);

struct tc_node *
tc_new_node(
    struct tc_tree *tree,
    size_t param,
    size_t nchildren,
    void *partitioning
);

int
tc_replace_node(
    struct tc_node *orig,
    struct tc_node *node
);

void tc_free_segments(struct tc_segment *segments, size_t S);

int
tc_clustering(
    const void *ds[],
    size_t N,
    const struct tc_param_def param_def[],
    size_t K,
    tc_clustering_cb cb,
    void *data,
    const struct tc_opts *opts
);

double
tc_log_likelihood(
    const struct tc_tree *tree,
    const void *ds[],
    size_t N
);

struct tc_segment *
tc_segments(
    const struct tc_tree *tree,
    const void *ds[],
    size_t N,
    size_t *S
);

void
tc_dump_tree_simple(const struct tc_tree *tree, const struct tc_node *node);

void tc_dump_segments_json(const struct tc_tree *tree);

#endif /* TC_H */
