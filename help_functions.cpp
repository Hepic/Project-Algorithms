#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <algorithm>
#include <math.h>
#include "help_functions.h"

char *get_arguments(const char *argv[], int len, const char *flag, bool same) {
    for (int i = 0; i < len; ++i) {
        if(!strcmp(argv[i], flag)) {
            if (same) {
                return (char*)argv[i];
            } else {
                return (char*)argv[i+1];
            }
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

void insert_curves_into_hashtables(vector<HashTable> &hashtables, double delta, const char *hash_function) {
    vector<Curve> concat_curves(input_curves.size(), Curve());
    
    for (int i = 0; i < global_L; ++i) {
        multiple_grids(concat_curves, delta);

        for (int j = 0; j < (int)concat_curves.size(); ++j) {
            hashtables[i].insert(input_curves[j], concat_curves[j], hash_function);
            concat_curves[j].clear_curve();
        }
    }
}

void search_curves_from_hashtables(const vector<HashTable> &hashtables, double delta, double R, const char *hash_function, const char *dist_function, vector<set<Curve> > &R_closest_curves, const vector<bool> &grid_curves_found, const vector<Curve> &centroids, vector<bool> &visited, bool check) {
    vector<double> min_curve_dist((int)centroids.size(), -1);
    Curve concat_curve;

    for (int i = 0; i < (int)centroids.size(); ++i) {
        if (grid_curves_found[i] && R == 0) {
            continue;
        }
        
        for (int j = 0; j < global_L; ++j) {
            concat_curve.clear_curve();
            multiple_grids_curve(concat_curve, delta, hashtables[j].get_id(), centroids[i]); // retrieve concatenated grid_curve
 
            vector<Curve> closer_curves;
            hashtables[j].search(closer_curves, centroids[i], concat_curve, hash_function, dist_function, R, visited, check);
            
            if (!closer_curves.empty()) {
                for (int k = 0; global_k < (int)closer_curves.size(); ++k) {
                    R_closest_curves[i].insert(closer_curves[k]);
                }
            }
        }
    }
}

void general_search(const vector<HashTable> &hashtables, double delta, double R, const char *hash_function, const char *dist_function, vector<set<Curve> > &R_closest_curves, const vector<Curve> &centroids, vector<bool> &grid_curves_found, vector<bool> &visited) {
    // check first if grid_curve is same
    search_curves_from_hashtables(hashtables, delta, R, hash_function, dist_function, R_closest_curves, grid_curves_found, centroids, visited);
    
    // keep curves that grid_curve was found
    for (int i = 0; i < (int)centroids.size(); ++i) { 
        if (!R_closest_curves[i].empty()) {
            grid_curves_found[i] = true;
        }
    }

    // search without checking grid_curve for curves that we did not find anything before
    search_curves_from_hashtables(hashtables, delta, R, hash_function, dist_function, R_closest_curves, grid_curves_found, centroids, visited, false); 
}
