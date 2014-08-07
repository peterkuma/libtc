/*
 * tree.c
 *
 * Tree structures and functions.
 *
 */

#include <stdlib.h>
#include <strings.h>
#include <assert.h>
#include <math.h>

#include "misc.h"
#include "tree.h"
#include "tc.h"

union tc_value
min(union tc_valuep data, size_t N, enum tc_param_size size)
{
    size_t n = 0;
    union tc_value value;
    if (size == TC_FLOAT64) {
        value.float64 = INFINITY;
        for (n = 0; n < N; n++)
            if (data.float64[n] < value.float64)
                value.float64 = data.float64[n];
    } else if (size == TC_INT64) {
        value.int64 = INT64_MAX;
        for (n = 0; n < N; n++)
            if (data.int64[n] < value.int64)
                value.int64 = data.int64[n];
    } else assert(0);
    return value;
}

union tc_value
max(union tc_valuep data, size_t N, enum tc_param_size size)
{
    size_t n = 0;
    union tc_value value;
    if (size == TC_FLOAT64) {
        value.float64 = -INFINITY;
        for (n = 0; n < N; n++)
            if (data.float64[n] > value.float64)
                value.float64 = data.float64[n];
    } else if (size == TC_INT64) {
        value.int64 = INT64_MAX;
        for (n = 0; n < N; n++)
            if (data.int64[n] > value.int64)
                value.int64 = data.int64[n];
    } else assert(0);
    return value;
}

/*
 * Free range `range`. The range needs to be deallocated in addition
 * if allocated dynamically.
 */
void
free_range(struct tc_range *range)
{
    if (range->categories != NULL) free(range->categories);
    range->categories = NULL;
    range->ncategories = 0;
}

/*
 * Free segment `segment`. The segment needs to be deallocated in addition
 * if allocated dynamically.
 */
void
free_segment(struct tc_segment *segment)
{
    size_t k = 0;
    for (k = 0; k < segment->K; k++) {
        free_range(&segment->ranges[k]);
    }
    segment->K = 0;
    segment->NX = 0;
}

/*
 * Allocate buffer of size `size` on the tree `tree`. Returns NULL on failure.
 */
void *
tree_alloc(struct tc_tree *tree, size_t size)
{
    void *obj = tree->p;
    if ((tree->size - (tree->p - tree->buf)) < size) {
        return NULL; /* No space left in buf. */
    }
    tree->p += size;
    return obj;
}

struct tc_node *
tree_alloc_node(struct tc_tree *tree, size_t part_size, size_t nchildren)
{
    struct tc_node *node = NULL;
    node = tree_alloc(tree, sizeof(struct tc_node));
    if (node == NULL) return NULL;
    node->part_size = part_size;
    node->nchildren = nchildren;
    node->part.buf = tree_alloc(tree, part_size);
    if (node->part.buf == NULL) return NULL;
    node->children = tree_alloc(tree, nchildren*sizeof(struct tc_node *));
    if (node->children == NULL) return NULL;
    return node;
}

/*
 * Attach `node` to `tree`. The node is added to the sequential traversing.
 * The node needs to be added to the tree structure in addition.
 */
void
tree_attach_node(struct tc_node *node)
{
    size_t i = 0;
    struct tc_tree *tree = NULL;
    tree = node->tree;
    node->next = NULL;
    node->prev = tree->last;
    if (tree->last != NULL)
        tree->last->next = node;
    tree->last = node;
    if (tree->first == NULL)
        tree->first = node;
    for (i = 0; i < node->nchildren; i++)
        tree_attach_node(node->children[i]);
}

/*
 * Detach `node` from `tree`. The node is removed from sequential traversing.
 * The node needs to be removed from the tree structure in addition.
 */
void
tree_detach_node(struct tc_node *node)
{
    size_t i = 0;
    struct tc_tree *tree = NULL;
    tree = node->tree;
    for (i = 0; i < node->nchildren; i++)
        tree_detach_node(node->children[i]);
    if (node->prev) node->prev->next = node->next;
    if (node->next) node->next->prev = node->prev;
    if (tree->last == node) tree->last = node->prev;
    if (tree->first == node) tree->first = node->next;
    node->prev = NULL;
    node->next = NULL;
}

/*
 * Copy node `node` to the tree `tree`. Child nodes are preserved
 * (point to the old tree). Returns pointer to the new node or NULL on failure.
 */
struct tc_node *
copy_node(const struct tc_node *node, struct tc_tree *tree)
{
    struct tc_node *copy;
    copy = tree_alloc_node(tree, node->part_size, node->nchildren);
    if (copy == NULL) return NULL;
    *copy = *node;
    bcopy(
        node->part.buf,
        copy->part.buf,
        node->part_size
    );
    bcopy(
        node->children,
        copy->children,
        node->nchildren
    );
    return copy;
}


/*
 * Compact tree `old` into tree `new`. Returns 0 on success, -1 on failure.
 */
int
compact_tree(struct tc_tree *new, const struct tc_tree *old)
{
    size_t i = 0;
    const struct tc_node *node = NULL;
    copy_node(old->root, new);
    if (node == NULL) return -1;
    while (node != NULL) {
        /* Copy children to the new tree. */
        for (i = 0; i < node->nchildren; i++) {
            node->children[i] = copy_node(node->children[i], new);
            if (node->children[i] == NULL) return -1;
        }
        node = node->next;
    }
    return 0;
}

/*
 * Initialize segment `segment`. `K` is the number of child nodes.
 * Initialized segment needs to be freed with free_segment.
 */
void
init_segment(struct tc_segment *segment, size_t K)
{
    segment->K = K;
    segment->NX = 0;
    segment->V = 0;
    segment->ranges = calloc(K, sizeof(struct tc_range));
}

/*
 * Determine range of node `node` in parameter `param`. The result is saved
 * to `range. The callee should call free_range on `range` after use.
 */
void
node_range(
    const struct tc_node *node,
    size_t param,
    struct tc_range *range
) {
    size_t i = 0;
    const struct tc_node *n = NULL, *child = NULL;
    const struct tc_param_def *pd = NULL;

    n = node;
    pd = &node->tree->param_def[param];
    range->min = pd->min;
    range->max = pd->max;
    range->categories = NULL;
    range->ncategories = 0;

    //debug("pd->min,max = (%" PRId64 ", %" PRId64 ")\n", range->min.int64, range->max.int64);

    for (n = node->parent, child = node; n != NULL; child = n, n = n->parent) {
        if (n->param != param) continue;

        if (pd->type == TC_METRIC) {
            for (i = 0; i < n->nchildren; i++) {
                if (n->children[i] == child) break; /* Found. */
            }
            assert(i != n->nchildren);
            //debug("%zu", i);
            if (i != 0) {
                if (IS_FLOAT64(pd))
                    range->min.float64 = MAX(range->min.float64, n->part.float64[i-1]);
                else if (IS_INT64(pd))
                    range->min.int64 = MAX(range->min.int64, n->part.int64[i-1]);
            }
            if (i + 1 != n->nchildren) {
                if (IS_FLOAT64(pd))
                    range->max.float64 = MIN(range->max.float64, n->part.float64[i]);
                else if (IS_INT64(pd))
                    range->max.int64 = MIN(range->max.int64, n->part.int64[i]);
            }
        } else if(pd->type == TC_NOMINAL) {
            /* Not implemented. */
        } else assert(0);
    }
}

/*
 * Return width of `node` (range size, resp. the number of categories).
 */
double
node_width(const struct tc_node *node)
{
    double width = 0;
    struct tc_range range;
    const struct tc_param_def *pd = NULL;
    pd = &node->tree->param_def[node->param];
    node_range(node, node->param, &range);
    if (pd->type == TC_NOMINAL) {
        width = range.ncategories;
    } else if (pd->type == TC_METRIC) {
        if (IS_INT64(pd))
            width = range.max.int64 - range.min.int64;
        else if (IS_FLOAT64(pd))
            width = range.max.float64 - range.min.float64;
        else assert(0);
    } else {
        assert(0);
    }
    free_range(&range);
    return width;
}

/*
 * Find child node `child` of node `node`. Returns the index of the child node,
 * or -1 if not found.
 */
size_t
find_child(const struct tc_node *node, const struct tc_node *child)
{
    size_t i = 0;
    for (i = 0; i < node->nchildren; i++) {
        if (node->children[i] == child) return i;
    }
    return -1;
}

/*
 * Generate a random partitioning of `range` into two partitions.
 * Returns the partitioning (allocated dynamically), or NULL on failure.
 * The partitioning needs to be freed with `free(part.buf)`.
 */
union tc_valuep
rand_part(const struct tc_range *range, const struct tc_param_def *pd)
{
    union tc_valuep part;
    union tc_value cut;
    part.buf = NULL;
    if (pd->type == TC_METRIC) {
        part.buf = calloc(1, TC_SIZE[pd->size]);
        if (part.buf == NULL) return part;
        if (IS_INT64(pd)) {
            if (range->max.int64 - range->min.int64 <= 1)
                return part; /* We can't split. */
            cut.int64 = rand() % (range->max.int64 - range->min.int64 - 1) + 1;
            part.int64[0] = cut.int64;
        } else if(IS_FLOAT64(pd)) {
            cut.float64 = range->min.float64 +
                frand1()*(range->max.float64 - range->min.float64);
            part.float64[0] = cut.float64;
        } else assert(0);
    } else if (pd->type == TC_NOMINAL) {
        /* Not implemented. */
    } else assert(0);
    return part;
}

/*
 * Check tree `tree` for correctness.
 * Returns true when correct, false otherwise.
 */
bool
check_tree(const struct tc_tree *tree)
{
    if (tree->param_def == NULL) return false;
    if (tree->K < 0) return false;
    if (tree->root == NULL) return false;
    if (tree->first == NULL) return false;
    if (tree->last == NULL) return false;
    return check_subtree(tree, tree->root);
}

/*
 * Check subtree of tree `tree` at node `node`.
 * Returns true when correct, false otherwise.
 */
bool
check_subtree(const struct tc_tree *tree, const struct tc_node *node)
{
    size_t i = 0;
    const struct tc_param_def *pd = NULL;
    pd = &tree->param_def[node->param];
    if (node->nchildren < 0)
        return false;
    for (i = 0; i < node->nchildren; i++) {
        if (node->children[i]->parent != node)
            return false;
        if (!check_subtree(tree, node->children[i]))
            return false;
    }
    if (pd->type == TC_METRIC) {
        for (i = 1; i + 1 < node->nchildren; i++)
            if (node->part.float64[i] < node->part.float64[i-1])
                return false;
    }
    return true;
}

/*
 * Returns true if `node` is a segment, or false otherwise.
 * Segment is a node with no child nodes.
 */
bool
is_segment(const struct tc_node *node)
{
    return node->nchildren == 0;
}

/*
 * Count the number of segments (left nodes) in tree `tree`.
 */
size_t
count_segments(const struct tc_tree *tree)
{
    size_t n = 0;
    struct tc_node *node = NULL;
    node = tree->first;
    while (node != NULL) {
        if (is_segment(node)) n++;
        node = node->next;
    }
    return n;
}

/*
 * Return the `s`-th segment of tree `tree` or NULL if s is greater than
 * the number of segments.
 */
struct tc_node *
select_segment(const struct tc_tree *tree, size_t s)
{
    size_t n = 0;
    struct tc_node *node = NULL;
    node = tree->first;
    while (node != NULL) {
        if (is_segment(node)) {
            if (n == s) return node;
            n++;
        }
        node = node->next;
    }
    return NULL;
}

/* Return true if a `node` is a supersegment, or false otherwise.
 * Supersegment is a node with two adjacent segments (metric parameter),
 * resp. two or more segments (nominal parameter).
 */
bool
is_supersegment(const struct tc_node *node)
{
    size_t i = 0;
    size_t S = 0;
    const struct tc_param_def *pd = NULL;
    pd = &node->tree->param_def[node->param];
    if (pd->type == TC_METRIC) {
        /* Need two adjacent segments. */
        for (i = 1; i < node->nchildren; i++) {
            if (is_segment(node->children[i-1]) &&
                is_segment(node->children[i]))
            {
                return true;
            }
        }
    } else if(pd->type == TC_NOMINAL) {
        /* Need two segments. */
        S = 0;
        for (i = 0; i < node->nchildren; i++)
            if (is_segment(node->children[i])) S++;
        return S >=2;
    } else assert(0);
    return false;
}

/*
 * Count the number of supersegments in `tree`.
 */
size_t
count_supersegments(const struct tc_tree *tree)
{
    size_t n = 0;
    struct tc_node *node = NULL;
    node = tree->first;
    while (node != NULL) {
        if (is_supersegment(node))
            n++;
        node = node->next;
    }
    return n;
}

/*
 * Return the `ss`-th supersegment of `tree`, or NULL if `ss` is greater
 * than the number of supersegments.
 */
struct tc_node *
select_supersegment(const struct tc_tree *tree, size_t ss)
{
    size_t i = 0;
    struct tc_node *node = NULL;
    node = tree->first;
    while (node != NULL) {
        if (is_supersegment(node))
            if (i++ == ss) return node;
        node = node->next;
    }
    return NULL;
}

/*
 * Return true if cut `i` of a metric node `node` is movable, false otherwise.
 * Movable cut is a cut between two segments.
 */
bool
is_movable_cut(const struct tc_node *node, size_t i)
{
    const struct tc_param_def *pd = NULL;
    pd = &node->tree->param_def[node->param];
    assert(pd->type == TC_METRIC);
    return is_segment(node->children[i]) &&
        is_segment(node->children[i+1]);
}

/*
 * Return count of movable cuts is `node`.
 */
size_t
count_movable_cuts(const struct tc_node *node)
{
    size_t i = 0;
    size_t n = 0;
    for (i = 0; i + 1 < node->nchildren; i++) {
        if (is_movable_cut(node, i))
            n++;
    }
    return n;
}

/*
 * Return the index of `c`-th movable cut of `node`.
 */
size_t
select_movable_cut(const struct tc_node *node, size_t c)
{
    size_t i = 0;
    size_t n = 0;
    for (i = 0; i + 1 < node->nchildren; i++)
        if (is_movable_cut(node, i))
            if (n++ == c) return i;
    return -1;
}
