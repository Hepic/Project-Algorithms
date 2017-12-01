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

void clustering(const vector<HashTable> &hashtables, int L, int k, double delta) {
    vector<const Curve*> centroids;
    vector<int> assignment(input_curves.size());
    vector<double> close_dist(input_curves.size()), close_dist_sec(input_curves.size());
    bool check;
    double value;
    int num_of_clusters = 2;

    k_means_pp(centroids, input_curves.size(), num_of_clusters, "DFT");
    cout << "initialization ended" << endl;
    
    vector<vector<int> > clusters(num_of_clusters);
    
    do {
        value = loyd_assignment(centroids, clusters);
        //check = PAM_update(centroids, value, clusters); 
        check = mean_frechet_update(centroids, clusters);
    } while(check);

    double min_s = -1, max_s = -1;

    for (int i = 0; i < num_of_clusters; ++i) {
        double diss_a = 0, diss_b = 0, res = 0;

        for (int j = 0; j < (int)clusters[i].size(); ++j) {
            diss_a += close_dist[clusters[i][j]];
            diss_b += close_dist_sec[clusters[i][j]];
            
            double s_elem = (diss_b - diss_a) / max(diss_a, diss_b);
            res += s_elem;
        }
        
        res /= clusters[i].size();
        min_s = (min_s == -1 ? res : min(min_s, res));
        max_s = (max_s == -1 ? res : max(max_s, res));
    }

    cout << "value = " << value << endl;

    for (int i = 0; i < (int)input_curves.size(); ++i) {
        cout << assignment[i] << " ";
    }

    cout << endl;
}
