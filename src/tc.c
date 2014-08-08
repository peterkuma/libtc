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
#include <errno.h>

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
	if (tree == NULL) {
		errno = ENOMEM;
		goto error;
	}
	tree->size = size;
	tree->last = NULL;
	tree->p = tree->buf;
	tree->param_def = param_def;
	tree->K = K;
	tree->root = tc_new_leaf(tree);
	if (tree->root == NULL) {
		errno = ENOMEM;
		goto error;
	}
	tree_attach_node(tree->root);
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
	void *partitioning
) {
	size_t i = 0;
	struct tc_node *node = NULL;

    node = tree_alloc(tree, sizeof(struct tc_node));
    if (node == NULL) {
    	errno = ENOMEM;
    	return NULL;
    }

    node->nchildren = nchildren;
    node->children = tree_alloc(tree, nchildren*sizeof(node->children));
    if (node->children == NULL) {
    	errno = ENOMEM;
    	return NULL;
    }

	node->tree = tree;
	node->param = param;
	node->ncuts = 0;
	node->ncategories = 0;

	if (IS_METRIC(node)) {
		node->ncuts = nchildren > 0 ? nchildren - 1 : 0;
		node->cuts = tree_alloc(tree, node->ncuts*sizeof(node->cuts));
		if (node->cuts == NULL) {
			errno = ENOMEM;
			return NULL;
		}
		bcopy(
			partitioning,
			node->cuts,
			node->ncuts*sizeof(node->cuts)
		);
	} else if (IS_NOMINAL(node)) {
		node->ncategories = PD(node)->max.int64 - PD(node)->min.int64 + 1;
		node->categories = tree_alloc(tree, node->ncategories*sizeof(node->categories));
		if (node->categories == NULL) {
			errno = ENOMEM;
			return NULL;
		}
		bcopy(
			partitioning,
			node->categories,
			node->ncategories*sizeof(node->categories)
		);
	} else assert(0);

	for (i = 0; i < node->nchildren; i++) {
		node->children[i] = tc_new_leaf(tree);
		if (node->children[i] == NULL) {
			errno = ENOMEM;
			return NULL;
		}
		node->children[i]->parent = node;
	}
	return node;
}

int
tc_replace_node(
	struct tc_node *orig,
	struct tc_node *node
) {
	size_t i = 0;
	struct tc_tree *tree = NULL;
	if (orig->tree != node->tree) {
		errno = EINVAL;
		return -1;
	}
	tree = orig->tree;
	if (orig->parent != NULL) {
		i = find_child(orig->parent, orig);
		orig->parent->children[i] = node;
		node->parent = orig->parent;
		orig->parent = NULL;
	} else {
		tree->root = node;
	}
	tree_detach_node(orig);
	tree_attach_node(node);
	return 0;
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
		for (i = 0; i < node->ncuts; i++) {
			printf("%lf", node->cuts[i]);
			if (i < node->ncuts - 1) printf(", ");
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
				node_range(node, k, &range);
				printf("[%lf,%lf]", range.min, range.max);
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
