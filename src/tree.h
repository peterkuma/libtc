/*
 * tree.c
 *
 * Tree structures and functions.
 *
 */

#include <stdbool.h>

#include "tc.h"

#define IS_INT64(pd) ((pd)->size == TC_INT64)
#define IS_FLOAT64(pd) ((pd)->size == TC_FLOAT64)

union tc_value min(union tc_valuep data, size_t N, enum tc_param_size size);

union tc_value max(union tc_valuep data, size_t N, enum tc_param_size size);

void free_range(struct tc_range *range);

void free_segment(struct tc_segment *segment);

void *tree_alloc(struct tc_tree *tree, size_t size);

struct tc_node *
tree_alloc_node(
    struct tc_tree *tree,
    size_t part_size,
    size_t nchildren
);

void tree_attach_node(struct tc_node *node);

void tree_detach_node(struct tc_node *node);

struct tc_node *copy_node(const struct tc_node *node, struct tc_tree *tree);

int compact_tree(struct tc_tree *new, const struct tc_tree *old);

void init_segment(struct tc_segment *segment, size_t K);

void
node_range(
    const struct tc_node *node,
    size_t param,
    struct tc_range *range
);

double node_width(const struct tc_node *node);

size_t find_child(const struct tc_node *node, const struct tc_node *child);

union tc_valuep rand_part(
    const struct tc_range *range,
    const struct tc_param_def *pd
);

bool check_tree(const struct tc_tree *tree);

bool check_subtree(const struct tc_tree *tree, const struct tc_node *node);

bool is_segment(const struct tc_node *node);

size_t count_segments(const struct tc_tree *tree);

struct tc_node *select_segment(const struct tc_tree *tree, size_t s);

bool is_supersegment(const struct tc_node *node);

size_t count_supersegments(const struct tc_tree *tree);

struct tc_node *select_supersegment(const struct tc_tree *tree, size_t ss);

bool is_movable_cut(const struct tc_node *node, size_t i);

size_t count_movable_cuts(const struct tc_node *node);

size_t select_movable_cut(const struct tc_node *node, size_t c);
