#include <iostream>
#include <cstdlib>
#include "cluster.h"
#include "curve.h"
#include "distances.h"
#include "help_functions.h"

vector<int> k_random_selection(int len, int k) {
    vector<int> all_random_points(len);
    vector<int> pick_random_points(k);
    
    for (int i = 0; i < len; ++i) {
        all_random_points[i] = i;
    }

    for (int i = 0; i < k; ++i) {
        int pos = rand() % (len - i) + i;
        swap(all_random_points[pos], all_random_points[i]);
        
        pick_random_points[i] = all_random_points[i];
    }

    return pick_random_points;
}

vector<int> k_means_pp(int len, int k, const char *dist_func) {
    vector<int> centroids;
    vector<double> min_distance(len, -1);
    vector<bool> is_centroid(len, false);
    int pos;
    double max_sum = 0;

    centroids.reserve(k);

    for (int t = 0; ; ++t) {
        if (!t) {
            pos = rand() % len;

            is_centroid[pos] = true; 
            centroids.push_back(pos);
        }
        
        if (t == k - 1) {
            break;
        }

        for (int i = 0; i < len; ++i) {
            if (is_centroid[i]) {
                continue;
            }

            double dist = compute_distance(input_curves[i], input_curves[pos], dist_func);
 
            if (min_distance[i] == -1 || dist < min_distance[i]) {
                if (min_distance[i] != -1) {
                    max_sum -= min_distance[i] * min_distance[i];
                }

                max_sum += dist * dist;
                min_distance[i] = dist;
            }
        }
        
        double value = uniform_distribution(0, max_sum);
        double curr = 0;

        for (int i = 0; i < len; ++i) {
            if (is_centroid[i]) {
                continue;
            }

            curr += min_distance[i] * min_distance[i];

            if (curr >= value) {
                pos = i;
                break;
            }
        }
        
        is_centroid[pos] = true;
        centroids.push_back(pos);
        
        max_sum -= min_distance[pos] * min_distance[pos];
    }

    return centroids;
}

vector<int> loyd_assignment(vector<int> centroids) {
    vector<int> assignment(input_curves.size());
    double minim, dist;
    int num_curve;
    
    for (int i = 0; i < (int)input_curves.size(); ++i) {
        minim = -1;
        for (int j=0; j < (int)centroids.size(); ++j) {
            if (centroids[j] == i) {
                continue;
            }
            dist = compute_distance(input_curves[i], input_curves[centroids[j]],"DFT");
            
            if (minim == -1 || dist < minim) {
                minim = dist;
                num_curve = centroids[j];
            }
        }
        
        assignment[i] = num_curve;
    }
    
    return assignment;
}
