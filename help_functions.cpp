#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <algorithm>
#include <math.h>
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

double normal_distribution(double mean, double stddev) {
    double U, V, S;
    static double X, Y;
    static bool Y_time = false;
    
    if(Y_time)
    {
        Y_time = false;
        return (Y * stddev + mean);
    }
    
    do {
        U = uniform_distribution(-1, 1);
        V = uniform_distribution(-1, 1);
        S = U*U + V*V;
    } while (S >= 1.0);
    
    X = U * sqrt((-2.0 * log(S)) / S);
    Y = V * sqrt((-2.0 * log(S)) / S);
    Y_time = true;
    
    return (X*stddev + mean);
}

double dot_product(const vector<double> &vec_1, const vector<double> &vec_2) {
    double dot = 0;
    
    for (int i = 0; i < (int)vec_1.size(); ++i) {
        dot += vec_1[i] * vec_2[i];
    }
    
    return dot;
}

void insert_curves_into_hashtables(vector<HashTable> &hashtables, int L, double delta, int k, const char *hash_function) {
    vector<Curve> concat_curves(input_curves.size(), Curve());
    
    for (int i = 0; i < L; ++i) { 
        multiple_grids(concat_curves, delta, k);

        for (int j = 0; j < (int)concat_curves.size(); ++j) {
            hashtables[i].insert(input_curves[j], concat_curves[j], hash_function);
            concat_curves[j].clear_curve();
        }
    }
}

void search_curves_from_hashtables(const vector<HashTable> &hashtables, int L, double delta, int k, double R, const char *hash_function, const char *dist_function, vector<vector<Curve> > &R_closest_curves, const vector<bool> &grid_curves_found, vector<Curve> &search_curves, bool check) {
    vector<double> min_curve_dist((int)search_curves.size(), -1);
    
    for (int i = 0; i < (int)search_curves.size(); ++i) {
        if (grid_curves_found[i] && R == 0) {
            continue;
        }
        
        for (int j = 0; j < L; ++j) {
            Curve concat_curve = multiple_grids_curve(delta, k, hashtables[j].get_id(), search_curves[i]); // retrieve concatenated grid_curve
            
            vector<Curve> closer_curves;
            closer_curves = hashtables[j].search(search_curves[i], concat_curve, hash_function, dist_function, R, check);
            
            if (!closer_curves.empty()) {
                    for (int k = 0; k < (int)closer_curves.size(); ++k) {
                        R_closest_curves[i].push_back(closer_curves[k]);
                    }
            }
        }
    }
}

void general_search(const vector<HashTable> &hashtables, int L, double delta, int k, double R, const char *hash_function, const char *dist_function, vector<vector<Curve> > &R_closest_curves, vector<Curve> &search_curves, vector<bool> &grid_curves_found) {
    // check first if grid_curve is same
    search_curves_from_hashtables(hashtables, L, delta, k, R, hash_function, dist_function, R_closest_curves, grid_curves_found, search_curves); 
    
    // keep curves that grid_curve was found
    for (int i = 0; i < (int)search_curves.size(); ++i) { 
        if (!R_closest_curves[i].empty()) {
            grid_curves_found[i] = true;
        }
    }

    // search without checking grid_curve for curves that we did not find anything before
    search_curves_from_hashtables(hashtables, L, delta, k, R, hash_function, dist_function, R_closest_curves, grid_curves_found, search_curves, false); 
}
