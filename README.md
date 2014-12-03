Tree Clustering Library
=======================

The tree clustering library performs clustering of elements using trees.
Elements are defined by a set of numerical (metric) parameters.
A tree approximates the probability density function from which elements
have been drawn by partitioning the parameter space hierarchically
into a number of rectangular segments, on which the probability
density function is assumed to have a uniform distribution.
Clustering trees are generated using a Metropolis-Hastings sampler.

**Note:** The package [rtc](https://github.com/peterkuma/rtc) provides an
easy-to-use binding for the R language.

Motivation
----------

libtc is a library which implements an MCMC sampler to produce a clustering
(here in the sense of partitioning a parameter space into rectangular areas)
of a number of elements. In a single dimension, this is equivalent to producing
a histogram with adaptable bounds. The unique aspect about this technique
is that it works with arbitrary number of dimensions, metric and nominal
parameters and is statistically objective.

Elements are defined by their parameters, i.e. they are points in the
parameter space. The partitioning of the parameter space is rectangular,
which means it is divided into a finite set of non-overlapping rectangles
spanning the whole parameter space. This allows future elements to be
attributed into a partition.

libtc was developed as an algorithm for automatically producing segmentation
of players of a free-to-play computer game into a relatively small number of
segments according to their differing parameters, e.g. age, country, level,
etc.

The partitioning is formulated as a tree, where each node defines partitioning
of the parent rectangle into a number of child rectangles by means of straight
lines (surfaces). New partitioning can be added, moved, or removed
in the MCMC process, eventually leading to a relatively optimal partitioning.

The goodness of partitioning, i.e. its likelihood (posterior probability
in Bayesian terms) is determined by the application of the Bayes theorem.
Elements are thought of as having been drawn from a probability density function
of such a form where the density is constant over each rectangle of the
partitioning, but otherwise unknown. This is a valid situation in the Bayesian
statistics, and we can infer what the posterior probability of our concrete
partitioning is, given this "observation" of elements. It turns out that
this probability can be derived analytically, and is similar to the formula
for entropy, but in a more extended form. Crucially, the likelihood
is the result of balance between volume of rectangles and how populated by
elements they are. Favored are large partitions with few elements,
small partitions with many elements, and highly uneven populations
among partitions. It turns out that these aspects
are naturally in opposition, and give rise to strong gradients,
allowing us to discern between better or worse partitionings,
and yield practically meaningful results.

The formula for likelihood/posterior probability of a partitioning takes
the form:

	1/V1^N1 1/V2^N2 ... 1/Vn^Nn N1!N2!...Nn! / (N + n)!

where V1, ..., Vn are volumes of partitions and N1, ..., Nn are populations
of elements in partitions. It is the result of integrating over all possible
constants defining the probability density function on the individual
partitions (can be performed analytically using the Beta function).
The prior probability for the constants is assumed to be uniform on the
in interval [0, 1].

Installation
------------

Requirements:

* POSIX-compatible environment such as Linux
* gcc compiler supporting C99
* [GNU Scientific Library](https://www.gnu.org/software/gsl/)
* [SCons](http://www.scons.org/)

You can build the library with the following command:

	scons

Install with:

	scons install [prefix=<path>]

The installation path can be chosen with the optional argument `prefix`.
If omitted, the library is installed under `/usr/local/`.

Example
-------

```C
#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <strings.h>
#include <err.h>

#include <tc.h>

tc_clustering_cb cb;

int
main(int argc, char *argv[])
{
    size_t N = 8; /* Number of elements. */
    size_t K = 2; /* Number of parameters. */

    const void *ds[2]; /* Dataset. */
    ds[0] = (double[]) { 1, 2, 1, 2, 4, 5, 4, 5 };
    ds[1] = (double[]) { 1, 1, 2, 2, 4, 4, 5, 5 };

    struct tc_param_def param_def[2];
    param_def[0] = (struct tc_param_def) {
        .type = TC_METRIC,
        .size = TC_FLOAT64,
        .fragment_size = 1
    };
    tc_param_def_init(&param_def[0], ds[0], N);
    param_def[1] = (struct tc_param_def) {
        .type = TC_METRIC,
        .size = TC_FLOAT64,
        .fragment_size = 1
    };
    tc_param_def_init(&param_def[1], ds[1], N);

    struct tc_opts opts = tc_default_opts;
    int res = tc_clustering(ds, N, param_def, K, cb, NULL, &opts);
    if (res != 0) {
        err(1, "Clustering failed");
    }
    return 0;
}

bool
cb(const struct tc_tree *tree, double l, const void **ds, size_t N, void *data)
{
    size_t S = 0;
    struct tc_segment *segments = tc_segments(tree, ds, N, &S);
    printf("l = %lf\n", l);
    for (size_t s = 0; s < S; s++) {
        printf("%zu: NX = %zu, V = %lf, ((%lf, %lf),(%lf, %lf))\n",
            s,
            segments[s].NX,
            segments[s].V,
            segments[s].ranges[0].min,
            segments[s].ranges[0].max,
            segments[s].ranges[1].min,
            segments[s].ranges[1].max
        );
    }
    printf("\n");
    tc_free_segments(segments, S);
    free(segments);
    segments = NULL;
    return true;
}
```

Compile with:

	gcc -std=c99 -Wall -o tc-example tc-example.c -ltc

Run:

	./tc-example

API Reference
-------------

### Dataset definition

Dataset is a number of elements defined by a set of parameters. Elements
have the same number parameters, all of which are required.

Dataset is defined by arrays of values, one array for each parameter.
Only metric parameters of type `double` are supported.
E.g. a dataset with two metric parameters and 8 elements can be defined as:

```C
void **ds[2];
ds[0] = (double[]) { 1, 2, 1, 2, 4, 5, 4, 5 };
ds[1] = (double[]) { 1, 1, 2, 2, 4, 4, 5, 5 };
```

### Parameters definition

Parameters are defined by an array of `tc_param_def` instances:

```C
struct tc_param_def {
	enum tc_param_type type;
	enum tc_param_size size;
    union tc_value min;
    union tc_value max;
	double fragment_size;
};
```

`type` is the type of parameter:

* **TC_METRIC** – Metric parameter (real valued).
* **TC_NOMINAL** – Nominal parameter (categorical).

`size` is the size of parameter data:

* **TC_FLOAT64** – values of type `double`.
* **TC_INT64** – values of type `int64_t`.

Nominal parameters are only compatible with size TC_INT64.

`min`, `max` are the minimum and maximum parameter values, which
define the range of the parameter space.

`fragment_size` is the size of fragment, i.e. the smallest unit by which
partitioning is performed.

### Functions

#### Main functions

##### tc_param_def_init

```C
void tc_param_def_init(struct tc_param_def *pd,	const void *data, size_t N)
```

Initialize parameter definition `pd`.
`data` is an array of data values in a given parameter
(i.e. a subset of a dataset). `N` is the number of elements in data.
Determines ranges of `data` necessary for subsequent computations.

##### tc_new_tree

```C
struct tc_tree *tc_new_tree(size_t size, const struct tc_param_def *param_def, size_t K)
```

Create a new tree. `size` is the size of tree memory buffer in bytes,
`param_def` are the parameter definitions and `K` is the number of
parameters.

Returns a pointer to the new tree or NULL on failure.
The tree should be deallocated with `free`.

##### tc_clustering

```C
int tc_clustering(
    const void *ds[],
    size_t N,
    const struct tc_param_def param_def[],
    size_t K,
    tc_clustering_cb cb,
    void *cb_data,
    const struct tc_opts *opts
)
```

Perform clustering of dataset `ds`. This function implements a
Metropolis-Hastings sampler to generate clustering trees.

`ds` is the dataset, `N` is the number of elements in dataset,
`param_def` are parameter definitions, `K` is the number of parameters.

`opts` is a stucture containing additional options:

```C
struct tc_opts {
    size_t nsamples; /* Number of samples to generate (excl. burn-in). */
    size_t maxiter; /* Maximum number of iterations. */
    double split_p; /* Probability of split. */
    double merge_p; /* Probability of merge. */
    double move_p; /* Probability of move. */
    double move_sd_frac; /* Move standard deviation as a fraction. */
    size_t max_segments; /* Maximum number of segments. */
};
```

The callback function `cb` is called for every sample accepted by the
sampler. It has the following form:

```C
bool cb(
    const struct tc_tree *tree,
    double l,
    const void **ds,
    size_t N,
    void *data
);
```

where `tree` is the tree, `l` is the log-likelihood of dataset being
drawn from the segmentation, `ds` and `N` are as in `tc_clustering`,
and `data` is an arbitrary user-supplied pointer passed to `tc_clustering`
as `cb_data`.

##### tc_segments

```C
struct tc_segment *tc_segments(
	const struct tc_tree *tree
	const union value *ds[],
	size_t N,
	size_t *S
)
```

Return segments defined by tree `tree`.
Segments correspend to leaf nodes of the tree.
`ds` is the data set, and `N` is the number of elements in data set.
The total number of segments is stored in `S`.

Returns an array of segments, which the callee should
free with `tc_free_segments` and `free`.

Segment is an instance of `tc_segment`:

```C
struct tc_segment {
	size_t NX;
	double V;
	struct tc_range ranges[];
};
```

where `NX` is the number of elements in segment, `V` is the volume of
segment (product of ranges), and `ranges` is an array
of ranges in each parameter.

`tc_range` is a structure defining a metric or nominal range:

```C
struct tc_range {
	double min;
	double max;
	int64_t *categories;
	size_t ncategories;
}
```

In the case of TC_METRIC parameter, `min` and `max` is the range
of segment in the respective parameter.

In the case of TC_NOMINAL parameter, categories is an array of `categories`
belonging to the segment, and `ncategories` is the number of categories.

##### tc_free_segments

```C
void tc_free_segments(struct tc_segment *segments, size_t S)
```

Free array of segments `segments`. `S` is the number of segments.
This only frees the internal structures. If allocated
dynamically, the array itself needs to be freed with `free`.

#### Miscellaneous functions

##### tc_new_node

```C
struct tc_node *tc_new_node(
	struct tc_tree *tree,
	size_t param,
	size_t nchildren,
	void *partitioning
)
```

Create a new node in `tree`. `param` is the parameter number over which
node splits the parameter space, `nchildren` is the number of child nodes,
and `partitioning` is the definition of partitioning. For metric parameters,
partitioning is an array of splits of type `double`.

Returns a pointer to the new node or NULL on failure. The node does
not need to be freed (it is allocated in the tree buffer).

##### tc_new_leaf

```C
struct tc_node *c_new_leaf(struct tc_tree *tree)
```

Create a new leaf node in `tree`.

Returns a pointer to the new node or NULL on failure. The node does
not need to be freed (it is allocated in the tree buffer).

##### tc_reaplce_node

```C
int tc_replace_node(struct tc_node *orig, struct tc_node *node)
```

Replace node `orig` with `node`. The original node is removed from the
tree structure, but remains allocated in the tree buffer.

Returns 0 on success, -1 on failure.

##### tc_log_likelihood

```C
double tc_log_likelihood(
	const struct tc_tree *tree,
	const union tc_value *ds[],
	size_t N
)
```

Calculate the log-likelihood of drawing data `ds` from tree `tree`.
`N` is the number of elements in `ds`.

Thanks
------

This library was developed thanks to the support of PIXEL FEDERATION, s.r.o.

License
-------

This software is released under the terms of the MIT License (see `LICENSE`).
