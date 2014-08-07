/*
 * tc.c
 *
 * Tree Clustering.
 *
 * Implementation of Cluster Analysis using Decision Trees
 * by the Metropolis-Hastings algorithm.
 *
 * See README.md for API reference.
 *
 */

#define _BSD_SOURCE

#include <math.h>
#include <stdlib.h>
#include <assert.h>
#include <stddef.h>
#include <stdint.h>
#include <strings.h>
#include <inttypes.h>
#include <stdbool.h>

#include "misc.h"
#include "tree.h"
#include "tc.h"

/* tc_param_size to number of bytes mapping. */
size_t TC_SIZE[] = {
    8,
    8
};

void
tc_param_def_init(struct tc_param_def *pd, const void *data, size_t N)
{
	union tc_valuep data_;
	data_.buf = data;
	pd->min = min(data_, N, pd->size);
	pd->max = max(data_, N, pd->size);
}

struct tc_tree *
tc_new_tree(size_t size, const struct tc_param_def *param_def, size_t K)
{
	struct tc_tree *tree = NULL;
	tree = calloc(1, sizeof(tree) + size);
	if (tree == NULL) goto error;
	tree->size = size;
	tree->last = NULL;
	tree->p = tree->buf;
	tree->param_def = param_def;
	tree->K = K;
	tree->root = tc_new_leaf(tree);
	if (tree->root == NULL) goto error;
	tree_attach_node(tree, tree->root);
	return tree;
error:
	if (tree != NULL) free(tree);
	return NULL;
}

struct tc_node *
tc_new_leaf(struct tc_tree *tree)
{
	return tc_new_node(tree, 0, 0, NULL);
}

struct tc_node *
tc_new_node(
	struct tc_tree *tree,
	size_t param,
	size_t nchildren,
	void *part
) {
	size_t i = 0;
	struct tc_node *node = NULL;
	const struct tc_param_def *pd;
	size_t part_size;

	pd = &tree->param_def[param];

	/* Determine part_size. */
	if (nchildren == 0) {
		part_size = 0;
	} else {
		if (pd->type == TC_METRIC) {
			part_size = (nchildren - 1)*TC_SIZE[pd->size];
		} else if (pd->type == TC_NOMINAL) {
			part_size = (pd->max.int64 - pd->min.int64 + 1)*sizeof(size_t);
		} else assert(0);
	}

	node = tree_alloc_node(tree, part_size, nchildren);
	if (node == NULL) return NULL;
	node->param = param;
	bcopy(part, node->part.buf, part_size);
	for (i = 0; i < node->nchildren; i++) {
		node->children[i] = tc_new_leaf(tree);
		if (node->children[i] == NULL) return NULL;
		node->children[i]->parent = node;
	}
	return node;
}

void
tc_replace_node(
	struct tc_tree *tree,
	struct tc_node *orig,
	struct tc_node *node
) {
	size_t i = 0;
	if (orig->parent != NULL) {
		i = find_child(orig->parent, orig);
		orig->parent->children[i] = node;
		node->parent = orig->parent;
		orig->parent = NULL;
	} else {
		tree->root = node;
	}
	tree_detach_node(tree, orig);
	tree_attach_node(tree, node);
}

void
tc_free_segments(struct tc_segment *segments, size_t S)
{
	size_t s = 0;
	for (s = 0; s < S; s++)
		free_segment(&segments[s]);
}

void
tc_dump_tree_simple(const struct tc_tree *tree, const struct tc_node *node)
{
	size_t i = 0;
	const struct tc_param_def *pd = NULL;
	if (node == NULL) node = tree->root;
	if (node == NULL) {
		printf("()\n");
		return;
	}
	pd = &tree->param_def[node->param];
	printf("(%zu, [", node->param);
	if (pd->type == TC_METRIC) {
		for (i = 0; i + 1 < node->nchildren; i++) {
			if (IS_INT64(pd)) printf("%zu", node->part.int64[i]);
			else if (IS_FLOAT64(pd)) printf("%lf", node->part.float64[i]);
			else assert(0);
			if (i < node->nchildren - 2) printf(", ");
		}
	} else if (pd->type == TC_NOMINAL) {
		/* Not implemented. */
	} else assert(0);
	printf("], [");
	for (i = 0; i < node->nchildren; i++) {
		tc_dump_tree_simple(tree, node->children[i]);
		if (i + 1 != node->nchildren)
			printf(", ");
	}
	printf("])");
	if (node == tree->root)
		printf("\n");
}

void
tc_dump_segments_json(const struct tc_tree *tree)
{
	size_t k = 0;
	struct tc_range range;
	const struct tc_node *node = NULL;
	const struct tc_param_def *pd = NULL;
	bool is_first = true;

	printf("{\"segments\":[");
	for (node = tree->first; node != NULL; node = node->next) {
		if (node->nchildren != 0)
			continue;
		if (!is_first) printf(",");
		is_first = false;
		pd = &tree->param_def[node->param];
		printf("[");
		for (k = 0; k < tree->K; k++) {
			if (pd->type == TC_METRIC) {
				node_range(tree, node, k, &range);
				printf("[%lf,%lf]", range.min.float64, range.max.float64);
				free_range(&range);
				if (k != tree->K - 1) printf(",");
			} else if (pd->type == TC_NOMINAL) {
				/* Not implemented. */
			} else assert(0);
		}
		printf("]");
	}
	printf("]}\n");
}

// void
// tc_dump_segments_json(const struct tc_tree *tree)
// {
//     size_t S = 0;
//     struct tc_segment *segments = NULL;
//     const struct tc_param_def *pd;
//     segments = tc_segments(tree, ds, N, &S);
//     struct tc_range range;
//     printf("{\"segments\":[");
//     for (size_t s = 0; s < S; s++) {
// 		if (s != 0) printf(",");
// 		printf("[");
// 		for (k = 0; k < tree->K; k++) {
// 			if (k != 0) printf(", ");
// 			pd = tree->param_def[k];
// 			if (pd->type == TC_METRIC) {
// 				range = segments[s].ranges[k];
// 				if (IS_FLOAT64(pd))
// 					printf("[%lf,%lf]", range.min.float64, range.max.float64);
// 				else if (IS_INT64(pd))
// 					/* Not implemented. */
// 				else assert(0);
// 			} else if (pd->type == TC_NOMINAL) {
// 				/* Not implemented. */
// 			} else assert(0);
// 		}
// 		printf("]");
//     }
//     printf("]}\n");
//     tc_free_segments(segments, S);
//     free(segments);
//     segments = NULL;
// }
