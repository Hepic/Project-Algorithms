#include <iostream>
#include <cstdlib>
#include "cluster.h"
#include "curve.h"
#include "distances.h"
#include "help_functions.h"
#include "binary_mean_tree.h"
#include "initialization.h"
#include "assignment.h"
#include "update.h"

void clustering(const vector<HashTable> &hashtables, double delta, vector<double> &silhouette_cluster, vector<const Curve*> &centroids, vector<vector<int> > &clusters) {
    vector<int> assignment(input_curves.size());
    bool check;
    double value;

    if (method_init == 1) {
        k_means_pp(centroids, input_curves.size(), "DFT");
    } else if (method_init == 2) {
        k_random_selection(centroids, input_curves.size());
    }
    cout << "initialization ended" << endl;
    
    do {
        value = loyd_assignment(centroids, clusters);
        if (method_update == 1) {
            check = mean_frechet_update(centroids, clusters);
        } else if (method_update == 2) {
            check = PAM_update(centroids, value, clusters);
        }
    } while(check);
    
     silhouette(centroids, clusters, silhouette_cluster);
    
    cout << "value = " << value << endl;

    for (int i = 0; i < (int)input_curves.size(); ++i) {
        cout << assignment[i] << " ";
    }

    cout << endl;
}

void silhouette(const vector<const Curve*> &centroids, vector<vector<int> > &clusters, vector<double> &silhouette_cluster) {
    double min_s = -1, max_s = -1;
    vector<double> close_dist((int)input_curves.size(), -1), close_dist_sec((int)input_curves.size(), -1);
    
    for (int i = 0; i < (int)input_curves.size(); ++i) {
        for (int j = 0; j < (int)centroids.size(); ++j) {
            double dist = compute_distance(input_curves[i], *centroids[j], "DFT");
            
            if (close_dist[i] == -1 || dist < close_dist[i]) {
                close_dist_sec[i] = close_dist[i];
                close_dist[i] = dist;
            } else if (close_dist_sec[i] == -1 || close_dist_sec[i] < dist) {
                close_dist_sec[i] = dist;
            }
        }
    }
    
    for (int i = 0; i < num_of_clusters; ++i) {
        double diss_a = 0, diss_b = 0, res = 0;
        
        for (int j = 0; j < (int)clusters[i].size(); ++j) {
            diss_a += close_dist[clusters[i][j]];
            diss_b += close_dist_sec[clusters[i][j]];
            
            double s_elem = (diss_b - diss_a) / max(diss_a, diss_b);
            res += s_elem;
        }
        
        res = (double)(res/clusters[i].size());
        cout << "res:" << res << endl;
        silhouette_cluster[i] = res;
        min_s = (min_s == -1 ? res : min(min_s, res));
        max_s = (max_s == -1 ? res : max(max_s, res));
    }
}
