#include <stdlib.h>
#include <stdio.h>
#include <err.h>
#include <getopt.h>
#include <string.h>
#include <errno.h>

#include <tc.h>

#define debug(...) fprintf(stderr, __VA_ARGS__)

const char *program_name = NULL;

void
usage(void)
{
    fprintf(stderr, "Usage: %s FILE\n\n", program_name);
    fprintf(stderr, "Try `%s --help` for help.\n", program_name);
}

void
help(void)
{
    fprintf(stderr, "Usage: %s FILE\n\n", program_name);
    fprintf(stderr, "Perform segmentation of data using tree clustering.\n\n");
    fprintf(stderr, "Optional arguments:\n");
    fprintf(stderr, "  -h,--help        print this help information and exit\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "FILE is a file containing tab-separated data.\n");
}

void **
read_dataset(const char *filename, size_t *N, size_t *K)
{
    size_t n = 0;
    char *line = NULL;
    char *s = NULL;
    const char *p = NULL;
    int off = 0;
    FILE *fp = NULL;
    void **ds = NULL;
    bool first = true;
    size_t N_ = 1024;

    fp = fopen(filename, "r");
    if (fp == NULL)
        err(1, "%s", filename);

    while (1) {
        line = NULL;
        errno = 0;
        n = fscanf(fp, "%m[^\n]\n", &line);
        if (n != 1) {
            if (errno != 0)
                err(1, "%s", filename);
            break;
        }
        if (line[0] == '#') goto cleanup;

        if (first) {
            /* Count the number of fields. */
            p = line;
            while ((n = sscanf(p, "%ms%n", &s, &off)) == 1) {
                *K++;
                free(s);
                s = NULL;
                p += off;
            }
            ds = calloc(*K, sizeof(void *));
            if (ds == NULL)
                err(1, "");
            for (k = 0; k < *K; k++) {
                ds[k] = calloc(N_, sizeof)

            }
            first = false;
        }
        p = line;
        while ((n = sscanf(p, "%ms%n", &s, &off)) == 1) {
            printf("FIELD=%s\n", s);
            free(s);
            s = NULL;
            p += off;
        }


        printf("\n");

    cleanup:
        if (line != NULL) {
            free(line);
            line = NULL;
        }
    }

    return data;
}

int
main(int argc, char *argv[])
{
    const char *filename = NULL;
    int c = 0;
    int option_index = 0;
    static struct option long_options[] = {
        {"help", no_argument, 0, 0},
        {NULL, 0, NULL, 0}
    };

    program_name = argv[0];

    while ((c = getopt_long(
        argc,
        argv,
        "h",
        long_options,
        &option_index)
    ) != -1) {
        switch (c) {
        case 0:
            if (strcmp(long_options[option_index].name, "help") == 0) {
                help();
                exit(0);
            }
        case 'h':
            help();
            exit(0);
            break;
        }
    }

    if (argc - optind != 1) {
        usage();
        exit(1);
    }

    filename = argv[optind];

    void **ds = NULL;
    size_t N = 0;
    size_t K = 0;
    ds = read_dataset(filename, &N, &K);
    return 0;
}
