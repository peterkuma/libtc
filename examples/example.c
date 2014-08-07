#include <stdlib.h>
#include <stdio.h>
#include <err.h>

#include <tc.h>

#define debug(...) fprintf(stderr, __VA_ARGS__)

/*
 * Usage example.
 */

// #define N 5;
// #define K 2;

// enum country {
//  UK,
//  IT,
//  BR
// };

// size_t nsegments = 10;

// int64_t age[N] = { 22, 30, 55, 43, 15 };
// enum country country[N] = { UK, UK, IT, BR, IT };
// void *data[N];
// data[0] = age;
// data[1] = country;

// struct tc_param_def param_def[K];
// param_def[0].type = TC_METRIC;
// param_def[0].size = TC_INT64;
// param_def[1].type = TC_NOMINAL;
// param_def[1].size = TC_INT64;

// struct tc_opts opts;
// tc_opts_init(opts);

// struct tc_tree *tree = tc_new_tree(1024);
// int res;
// res = tree_clustering(tree, data, N, K, &param_def, nsegments, opts);

// /* Do something with tree... */

// free(tree);

tc_clustering_cb cb;

void
cb(const struct tc_tree *tree, double l, const void **ds, size_t N, void *data)
{
    size_t S = 0;

    struct tc_segment *segments = tc_segments(tree, ds, N, &S);
    printf("l = %lf\n", l);
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
    printf("\n");
    tc_free_segments(segments, S);
    free(segments);
    segments = NULL;
}

int
main(int argc, char *argv[])
{
    size_t N = 0;
    size_t K = 2;
    struct tc_tree *tree;
    struct tc_node *root, *node;
    int res;

    struct tc_param_def param_def[K];
    param_def[0].type = TC_METRIC;
    param_def[0].size = TC_FLOAT64;
    param_def[0].min.float64 = 1;
    param_def[0].max.float64 = 25;
    param_def[0].max.float64 = 5;
    param_def[1].type = TC_METRIC;
    param_def[1].size = TC_FLOAT64;
    param_def[1].min.float64 = 1;
    param_def[1].max.float64 = 25;
    // param_def[1].max.float64 = 5;

    const void *ds[2]; /* Data set. */
    // N = 8;
    // ds[0] = (double[]) { 1, 2, 1, 2, 4, 5, 4, 5 };
    // ds[1] = (double[]) { 1, 1, 2, 2, 4, 4, 5, 5 };
    // ds[0] = (union tc_value *)(int64_t[N]) { 1, 2, 1, 2, 1, 2, 1, 2 };
    //ds[1] = (union tc_value *)(int64_t[N]) { 1, 1, 2, 2, 1, 2, 1, 2 };

    N = 84;
    ds[0] = (double[]) {1,2,3,4,5,6,7,8,9,10,11,12,1,2,3,4,5,6,7,8,9,10,11,12,1,2,3,4,5,6,7,8,9,10,11,12,1,2,3,4,5,6,7,8,9,10,11,12,1,2,3,4,5,6,7,8,9,10,11,12,1,2,3,4,5,6,7,8,9,10,11,12,1,2,3,4,5,6,7,8,9,10,11,1};
    ds[1] = (double[]) {1,1,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,2,2,2,3,3,3,3,3,3,3,3,3,3,3,3,4,4,4,4,4,4,4,4,4,4,4,4,5,5,5,5,5,5,5,5,5,5,5,5,6,6,6,6,6,6,6,6,6,6,6,6,7,7,7,7,7,7,7,7,7,7,7,7};

    // ds[0] = (double[]) {8 ,9 ,10 ,11 ,12 ,13 ,14 ,15 ,16 ,17 ,18 ,19 ,8 ,9 ,10 ,11 ,12 ,13 ,14 ,15 ,16 ,17 ,18 ,19 ,8 ,9 ,10 ,11 ,12 ,13 ,14 ,15 ,16 ,17 ,18 ,19 ,8 ,9 ,10 ,11 ,12 ,13 ,14 ,15 ,16 ,17 ,18 ,19 ,8 ,9 ,10 ,11 ,12 ,13 ,14 ,15 ,16 ,17 ,18 ,19 ,8 ,9 ,10 ,11 ,12 ,13 ,14 ,15 ,16 ,17 ,18 ,19 ,8 ,9 ,10 ,11 ,12 ,13 ,14 ,15 ,16 ,17 ,18 ,19};
    // ds[1] = (double[]) {11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17};

    // ds[0] = (double[]) {2 ,3 ,4 ,5 ,6 ,7 ,8 ,9 ,10 ,11 ,12 ,13 ,2 ,3 ,4 ,5 ,6 ,7 ,8 ,9 ,10 ,11 ,12 ,13 ,2 ,3 ,4 ,5 ,6 ,7 ,8 ,9 ,10 ,11 ,12 ,13 ,2 ,3 ,4 ,5 ,6 ,7 ,8 ,9 ,10 ,11 ,12 ,13 ,2 ,3,4 ,5 ,6 ,7 ,8 ,9 ,10 ,11 ,12 ,13 ,2 ,3 ,4 ,5 ,6 ,7 ,8 ,9 ,10 ,11 ,12 ,13 ,2 ,3 ,4 ,5 ,6 ,7 ,8 ,9 ,10 ,11 ,12 ,13 ,14 ,15 ,16 ,17 ,18 ,19 ,20 ,21 ,22 ,23 ,24 ,25 ,14,15 ,16 ,17 ,18 ,19 ,20 ,21 ,22 ,23 ,24 ,25 ,14 ,15 ,16 ,17 ,18 ,19 ,20 ,21 ,22 ,23 ,24 ,25 ,14 ,15 ,16 ,17 ,18 ,19 ,20 ,21 ,22 ,23 ,24 ,25 ,14 ,15 ,16 ,17 ,18 ,19,20 ,21 ,22 ,23 ,24 ,25 ,14 ,15 ,16 ,17 ,18 ,19 ,20 ,21 ,22 ,23 ,24 ,25 ,14 ,15 ,16 ,17 ,18 ,19 ,20 ,21 ,22 ,23 ,24 ,25};
    // ds[1] = (double[]) {3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 18, 18, 18, 18, 18, 18, 18, 18, 18, 18, 18, 18, 19, 19, 19, 19, 19, 19, 19, 19, 19, 19, 19, 19, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 22, 22, 22, 22, 22, 22, 22, 22, 22, 22, 22, 22, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23};

    // tc_param_def_init(&param_def[0], ds[0], N);
    // tc_param_def_init(&param_def[1], ds[1], N);

    // tree = tc_new_tree(1024, param_def, K);
    // root = tc_new_node(tree, 0, 2, (double[]) { 3 });
    // tc_replace_node(tree, tree->root, root);
    // node = tc_new_node(tree, 1, 2, (double[]) { 3 });
    // tc_replace_node(tree, root->children[0], node);
    // node = tc_new_node(tree, 1, 2, (double[]) { 3 });
    // tc_replace_node(tree, root->children[1], node);

    //root->children[0] = tc_new_node(tree, 1, 2, (int64_t[]) { 3 });
    //root->children[1] = tc_new_node(tree, 1, 2, (int64_t[]) { 3 });
    //root->children[1]->children[1] = tc_new_node(tree, 1, 2, (int64_t[]) { 5 });
    //tc_normalize_tree(tree);

    // tc_dump_tree_simple(tree, NULL);

    // size_t S = 0;
    // struct tc_segment *segments = tc_segments(tree, ds, N, &S);
    // for (size_t s = 0; s < S; s++) {
    //     debug("s = %zu, NX[s] = %zu, V[s] = %lf, ((%lf, %lf),(%lf, %lf))\n",
    //         s,
    //         segments[s].NX,
    //         segments[s].V,
    //         segments[s].ranges[0].min.float64,
    //         segments[s].ranges[0].max.float64,
    //         segments[s].ranges[1].min.float64,
    //         segments[s].ranges[1].max.float64
    //     );
    // }
    // tc_free_segments(segments, S);
    // free(segments);
    // segments = NULL;

    // double l = tc_log_likelihood(tree, ds, N);
    // printf("l = %lf\n", l);

    // free(tree);

    struct tc_opts opts = tc_default_opts;
    opts.burnin = 0;
    opts.nsamples = 30;
    res = tc_clustering(ds, N, param_def, K, cb, NULL, &opts);
    if (res != 0) {
        err(1, "Clustering failed");
    }
    return 0;
}
