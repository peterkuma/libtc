Tree Clustering
===============

**NOT FINISHED**

The tree clustering library performs clustering of data using decision trees.
A tree approximates the probability density function of data.
Optimal decision tree is determined by the maximum likelihood method.

API Reference
-------------

### Parameters definition

Parametric space is defined by an array of `tc_param_def` instances:

	struct tc_param_def {
		enum tc_param_type type; /* Parameter types. */
		enum tc_param_size size; /* Parameter sizes. */
		/* ... */
	};

`type` is the type of parameter:

* **TC_METRIC** – Metric parameter (real of integer valued).
* **TC_NOMINAL** – Nominal parameter (categorical).

`size` is the size of parameter data:

* **TC_FLOAT64** – values of type `double`.
* **TC_INT64** – values of type `int64_t`.

Nominal parameters are only compatible with size TC_INT64.

### Dataset definition

Dataset is defined by an array of pointers to data values of type
`union value`, e.g. a dataset with two metric parameters of type TC_INT64
and 8 elements can be defined as:

	union tc_value *ds[2];
	ds[0] = (union tc_value *)(int64_t[]) { 1, 2, 1, 2, 4, 5, 4, 5 };
	ds[1] = (union tc_value *)(int64_t[]) { 1, 1, 2, 2, 4, 4, 5, 5 };

### Functions

* void **tc_param_def_init**(struct tc_param_def **pd*, const union tc_value *data*[], size_t *N*)

	Initialize parameter definition `pd`.
	`data` is an array of data values in a given parameter
	(i.e. a subset of a dataset).
	Determines ranges of `data` necessary
	for subsequent computations. `N` is the number of elements in data.

* struct tc_tree ***tc_new_tree**(size_t *size*, const struct tc_param_def **param_def*, size_t *K*)

	Create a new tree. `size` is the size of tree memory buffer in bytes,
	`param_def` are the parameter definitions and `K` is the number of
	parameters.

	Returns a pointer to the new tree or NULL on failure.
	The tree should be deallocated with `free`.

* struct tc_node ***tc_new_node**(
	struct tc_tree **tree*,
	size_t *param*,
	size_t *nchildren*,
	void **part*
)

	Create a new node in `tree`. `param` is the parameter number over which
	node splits the parameter space, `nchildren` is the number of child nodes,
	and `part` is the definition of partitioning (split).

	Returns a pointer to the new node or NULL on failure. The node does
	not need to be freed (it is allocated in the tree buffer).

* void **tc_normalize_tree**(struct tc_tree **tree*)

	Normalize tree `tree`. Creates empty leaf nodes where missing and sets
	parent node references. This function should be called after manually
	constructing a tree.

* struct tc_segment ***tc_segments**(
	const struct tc_tree **tree*
	const union value **ds*[],
	size_t *N*,
	size_t **S*
)

	Determine segments of `tree`, and the number of elements in each segment.
	Segments correspend to leaf nodes of the tree. `ds` is the data set,
	and `N` is the number of elements in data set. The number of segments is
	saved in `S`. Returns an array of segments, which the callee should
	free with `tc_free_segments` and `free`.

	Segment is an instance of `tc_segment`:

		struct tc_segment {
			size_t NX;
			double V;
			struct tc_range ranges[];
		};

	where `NX` is the number of elements in segment, `V` is the volume of
	segment (product of ranges), and `ranges` is an array
	of ranges in each parameter. `tc_range` is a structure defining a metric
	or nominal range:

		struct tc_range {
			union tc_value min;
			union tc_value max;
			int64_t *categories;
			size_t ncategories;
		}

	In case of TC_METRIC parameter, `min` and `max` is the range
	of the segment in a parameter.
	In case of TC_NOMINAL parameter, categories is an array of `categories`
	belonging to the segment, and `ncategories` is the lengths of the array.

* void **tc_free_segments**(struct tc_segment **segments*, size_t *S*)

	Free array of segments `segments`. `S` is the number of segments.
	This only frees the internal structures, the array needs to
	be freed with `free` in addition (if allocated dynamically).

* double **tc_log_likelihood**(
	const struct tc_tree **tree*,
	const union tc_value **ds*[],
	size_t *N*
)

	Calculate the log-likelihood of drawing data `ds` from tree `tree`.
	`N` is the number of elements in `ds`.
