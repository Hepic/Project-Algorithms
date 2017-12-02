#include <iostream>
#include <cstdlib>
#include "cluster.h"
#include "curve.h"
#include "distances.h"
#include "help_functions.h"
#include "binary_mean_tree.h"
#include "update.h"

double swap_update_centroid(const vector<const Curve*> &centroids, int old_centr, int new_centr, const vector<int> &assign, const vector<double> &close_dist, const vector<double> &close_dist_sec, const char *metric) {
    double value = 0;
    
    for (int i = 0; i < (int)assign.size(); ++i) {
        if ((i != old_centr && (centroids[assign[i]]->get_int_id() == i)) || i == new_centr) {
            continue;
        }
        
        double dist = compute_distance(input_curves[i], input_curves[new_centr], metric);
        
        if (centroids[assign[i]]->get_int_id() != old_centr) {
            value += min(dist, close_dist[i]);
        } else {
            value += min(dist, close_dist_sec[i]);
        }
    }
    
    return value;
}

bool PAM_update(vector<const Curve*> &centroids, double value, const vector<vector<int> > &clusters, const char *metric) {
    vector<double> close_dist((int)input_curves.size(), -1), close_dist_sec((int)input_curves.size(), -1);
    vector<int> assign((int)input_curves.size());
    
    double min_value = value;
    int p_new_cent, p_clust;
    
    for (int i = 0; i < (int)input_curves.size(); ++i) {
        for (int j = 0; j < (int)centroids.size(); ++j) {
            double dist = compute_distance(input_curves[i], *centroids[j], metric);
            
            if (close_dist[i] == -1 || dist < close_dist[i]) {
                close_dist_sec[i] = close_dist[i];
                close_dist[i] = dist;
                assign[i] = j;
            } else if (close_dist_sec[i] == -1 || close_dist_sec[i] < dist) {
                close_dist_sec[i] = dist;
            }
        }
    }
    
    for (int i = 0; i < (int)clusters.size(); ++i) {
        for (int j = 0; j < (int)clusters[i].size(); ++j) {
            if (clusters[i][j] == centroids[i]->get_int_id()) {
                continue;
            }
            
            double new_value = swap_update_centroid(centroids, centroids[i]->get_int_id(), clusters[i][j], assign, close_dist, close_dist_sec, metric);
            
            if (new_value < min_value) {
                min_value = new_value;
                p_new_cent = clusters[i][j];
                p_clust = i;
            }
        }
    }
    
    if (value > min_value) {
        centroids[p_clust] = &input_curves[p_new_cent];
        return true;
    }
    
    return false;
}

bool mean_frechet_update(vector<const Curve*> &centroids, const vector<vector<int> > &clusters) {
    bool check = false;
    
    for (int i = 0; i < (int)clusters.size(); ++i) {
        BinaryMeanTree tree(clusters[i]);
        const Curve *mean_curve = tree.get_mean();
        
        mean_curve->print_curve();

        if (!centroids[i]->equal_curves(*mean_curve)) {
            centroids[i] = mean_curve;
            check = true;
        }
    }
    
    return check;
}

