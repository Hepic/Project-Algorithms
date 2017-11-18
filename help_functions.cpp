#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <algorithm>
#include "help_functions.h"

char *get_arguments(const char *argv[], int len, const char *flag) {
    for (int i = 0; i < len; ++i) {
        if (!strcmp(argv[i], flag)) {
            return (char*)argv[i+1]; 
        }
    }

    return (char*)"";
}

void perror_exit(const char *msg) {
    perror(msg);
    exit(EXIT_FAILURE);
}

double uniform_distribution(double A, double B) {
    double zero_to_one = (double)rand() / ((double)RAND_MAX + 1.0);
    double ret = A + zero_to_one * (B - A);

    return ret;
}
