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
        centroids.push_back(pos);
        
        max_sum -= min_distance[pos] * min_distance[pos];
    }

    return centroids;
}

double loyd_assignment(const vector<int>& centroids, vector<int> &assign, vector<double> &close_dist, vector<double> &close_dist_sec, vector<vector<int> > &clusters) {
    double min_dist, min_dist_sec, dist, value = 0;
    int p_centr;
    
    for (int i = 0; i < (int)clusters.size(); ++i) {
        clusters[i].clear();
    }

    for (int i = 0; i < (int)input_curves.size(); ++i) {
        min_dist = -1;
        
        for (int j = 0; j < (int)centroids.size(); ++j) {
            dist = compute_distance(i, centroids[j], "DFT");

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

double swap_update_centroid(int old_centr, int new_centr, const vector<int> &assign, const vector<double> close_dist, const vector<double> close_dist_sec) {
    double value = 0;

    for (int i = 0; i < (int)assign.size(); ++i) {
        if ((i != old_centr && (assign[i] == i)) || i == new_centr) {
            continue;
        }
        
        double dist = compute_distance(i, new_centr, "DFT");
        
        if (assign[i] != old_centr) {
            value += min(dist, close_dist[i]);
        } else {
            value += min(dist, close_dist_sec[i]);
        }
    }

    return value;
}

bool PAM_update(vector<int> &centroids, const vector<int> &assign, const vector<double> &close_dist, const vector<double> &close_dist_sec, double value, const vector<int> &cluster, int p_clust) {
    double min_value = value; 
    int new_cent = -1;
    
    for (int i = 0; i < (int)cluster.size(); ++i) {
        if (cluster[i] == centroids[p_clust]) {
            continue;
        }
        
        double new_value = swap_update_centroid(centroids[p_clust], cluster[i], assign, close_dist, close_dist_sec); 

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

Curve mean_curve_cluster(const vector<int> &cluster) {
    int len = cluster.size();
    vector<Curve> curves(len);

    for (int i = 0; i < len; ++i) {
        curves[i] = input_curves[cluster[i]];
    }

    while(len) {
        Curve mean_curve;
        int pos = 0;

        for (int i = 1; i < len; i += 2) {
            discrete_frechet_distance(curves[i - 1], curves[i], mean_curve, true);
            curves[pos++] = mean_curve;
        }

        if (len & 1) {
            curves[pos++] = curves[len - 1];
        }

        len >>= 1;
    }

    return curves[0];
}

double mean_frechet(vector<Curve> &centroids, const vector<int> &assign) {
    vector<vector<int> > clusters(centroids.size());

    for (int i = 0; i < (int)assign.size(); ++i) {
        clusters[assign[i]].push_back(i);
    }
    
    for (int i = 0; i < (int)centroids.size(); ++i) {
        centroids[i] = mean_curve_cluster(clusters[i]);
    }
    
    return 0;
}

void clustering() {
    vector<int> assignment(input_curves.size());
    vector<double> close_dist(input_curves.size()), close_dist_sec(input_curves.size());
    bool check;
    double value;
    
    for (int num_of_clusters = 2; num_of_clusters <= (int)input_curves.size(); ++num_of_clusters) {
        vector<int> centroids = k_means_pp(input_curves.size(), num_of_clusters, "DFT");
        cout << "initialization ended" << endl;
        
        vector<double> dissimilarity;
        vector<vector<int> > clusters;

        clusters.resize(num_of_clusters + 5);
        dissimilarity.resize(num_of_clusters + 5);

        do {
            for (int i = 0; i < num_of_clusters; ++i) {
                value = loyd_assignment(centroids, assignment, close_dist, close_dist_sec, clusters);
                check = PAM_update(centroids, assignment, close_dist, close_dist_sec, value, clusters[i], i); 
            }
        } while(check);
        
        double min_s = -1, max_s = -1;
        
        for (int i = 0; i < num_of_clusters; ++i) {
            double diss_a = 0, diss_b = 0;
            
            for (int j = 0; j < (int)clusters[i].size(); ++j) {
                diss_a += close_dist[clusters[i][j]];
                diss_b += close_dist_sec[clusters[i][j]];
            }
            
            double res = (diss_b - diss_a) / max(diss_a, diss_b);
            min_s = (min_s == -1 ? res : min(min_s, res));
            max_s = (max_s == -1 ? res : max(max_s, res));
        }
        
        if (max_s - min_s < 0.1) {
            cout << "k = " << num_of_clusters << endl;
            break;
        }
    } 
    
    cout << value << endl;

    for (int i = 0; i < (int)input_curves.size(); ++i) {
        cout << assignment[i] << " ";
    }

    cout << endl;
}
