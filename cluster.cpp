#include <iostream>
#include <cstdlib>
#include "cluster.h"
#include "curve.h"
#include "distances.h"
#include "help_functions.h"
#include "binary_mean_tree.h"

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

void k_means_pp(vector<int> &centroids_ind, vector<const Curve*> &centroids, int len, int k, const char *dist_func) {
    vector<double> min_distance(len, -1);
    vector<bool> is_centroid(len, false);
    int pos;
    double max_sum = 0;

    centroids_ind.reserve(k);
    centroids.reserve(k);

    for (int t = 0; ; ++t) {
        if (!t) {
            pos = rand() % len;

            is_centroid[pos] = true; 
            centroids_ind.push_back(pos);
            centroids.push_back(&input_curves[pos]);
        }
        
        if (t == k - 1) {
            break;
        }

        for (int i = 0; i < len; ++i) {
            if (is_centroid[i]) {
                continue;
            }

            double dist = compute_distance(i, pos, dist_func);
 
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
        centroids_ind.push_back(pos);
        centroids.push_back(&input_curves[pos]);
        
        max_sum -= min_distance[pos] * min_distance[pos];
    }
}

double loyd_assignment(const vector<int> &centroids, vector<int> &assign, vector<double> &close_dist, vector<double> &close_dist_sec, vector<vector<int> > &clusters, char *metric) {
    double min_dist, min_dist_sec, dist, value = 0;
    int p_centr;
    
    for (int i = 0; i < (int)clusters.size(); ++i) {
        clusters[i].clear();
    }

    for (int i = 0; i < (int)input_curves.size(); ++i) {
        min_dist = -1;
        
        for (int j = 0; j < (int)centroids.size(); ++j) {
            dist = compute_distance(i, centroids[j], metric);

            if (min_dist == -1 || dist < min_dist) {
                min_dist_sec = min_dist;
                min_dist = dist;
                p_centr = j;
            } else if (min_dist_sec == -1 || dist < min_dist_sec) {
                min_dist_sec = dist;
            }
        }

        clusters[p_centr].push_back(i);
        assign[i] = centroids[p_centr];

        close_dist[i] = min_dist;
        close_dist_sec[i] = min_dist_sec;

        value += min_dist;
    }

    return value;
}

double loyd_assignment(const vector<const Curve*> &centroids, vector<vector<int> > &clusters) {
    double min_dist, dist, value = 0;
    int p_centr;
    
    for (int i = 0; i < (int)clusters.size(); ++i) {
        clusters[i].clear();
    }

    for (int i = 0; i < (int)input_curves.size(); ++i) {
        min_dist = -1;
        
        for (int j = 0; j < (int)centroids.size(); ++j) {
            dist = discrete_frechet_distance(input_curves[i], *centroids[j]);

            if (min_dist == -1 || dist < min_dist) {
                min_dist = dist;
                p_centr = j;
            }
        }

        clusters[p_centr].push_back(i);
        value += min_dist;
    }

    return value;
}

double swap_update_centroid(int old_centr, int new_centr, const vector<int> &assign, const vector<double> &close_dist, const vector<double> &close_dist_sec, char *metric) {
    double value = 0;

    for (int i = 0; i < (int)assign.size(); ++i) {
        if ((i != old_centr && (assign[i] == i)) || i == new_centr) {
            continue;
        }

        double dist = compute_distance(i, new_centr, metric);

        if (assign[i] != old_centr) {
            value += min(dist, close_dist[i]);
        } else {
            value += min(dist, close_dist_sec[i]);
        }
    }

    return value;
}

bool PAM_update(vector<int> &centroids, const vector<int> &assign, const vector<double> &close_dist, const vector<double> &close_dist_sec, double value, const vector<int> &cluster, int p_clust, char *metric) {
    double min_value = value; 
    int new_cent = -1;
    
    for (int i = 0; i < (int)cluster.size(); ++i) {
        if (cluster[i] == centroids[p_clust]) {
            continue;
        }
        
        double new_value = swap_update_centroid(centroids[p_clust], cluster[i], assign, close_dist, close_dist_sec, metric);

        if (new_value < min_value) {
            min_value = new_value;
            new_cent = cluster[i];
        }
    }

    if (value > min_value) {
        centroids[p_clust] = new_cent;
        return true;
    }

    return false;
}

bool mean_frechet_update(vector<const Curve*> &centroids, const vector<vector<int> > &clusters) {
    bool check = false;

    for (int i = 0; i < (int)clusters.size(); ++i) {
        BinaryMeanTree tree(clusters[i]); 
        const Curve *mean_curve = tree.get_mean();

        if (!centroids[i]->equal_curves(*mean_curve)) {
            centroids[i] = mean_curve; 
            check = true;    
        }
    }

    return check;
}

void clustering(char *metric) {
    vector<const Curve*> centroids;
    vector<int> assignment(input_curves.size()), centroids_ind;
    vector<double> close_dist(input_curves.size()), close_dist_sec(input_curves.size());
    bool check;
    double value;
    int num_of_clusters = 2;
    
    k_means_pp(centroids_ind, centroids, input_curves.size(), num_of_clusters, metric);
    cout << "initialization ended" << endl;
    
    vector<vector<int> > clusters(num_of_clusters);
    
    do {
        for (int i = 0; i < num_of_clusters; ++i) {
            //value = loyd_assignment(centroids_ind, assignment, close_dist, close_dist_sec, clusters);
            value = loyd_assignment(centroids, clusters);
            //check = PAM_update(centroids_ind, assignment, close_dist, close_dist_sec, value, clusters[i], i, metric); 
            check = mean_frechet_update(centroids, clusters);
        }
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
