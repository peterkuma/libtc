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
