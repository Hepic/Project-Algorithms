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

void k_means_pp(vector<const Curve*> &centroids, int len, int k, const char *dist_func) {
    vector<double> min_distance(len, -1);
    vector<bool> is_centroid(len, false);
    int pos;
    double max_sum = 0;

    centroids.reserve(k);

    for (int t = 0; ; ++t) {
        if (!t) {
            pos = rand() % len;

            is_centroid[pos] = true; 
            centroids.push_back(&input_curves[pos]);
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
        centroids.push_back(&input_curves[pos]);
        
        max_sum -= min_distance[pos] * min_distance[pos];
    }
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
            dist = compute_distance(input_curves[i], *centroids[j], "DFT");

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

vector<int> range_search(const vector<HashTable> &hashtables, const vector<Curve> &centroids, int dim, int L, int k, double delta) {
    vector<int> assignment((int)input_curves.size(), -1);
    vector<set<Curve> > R_closest_curves((int)centroids.size());
    vector<bool> grid_curves_found((int)centroids.size(), false), visited((int)input_curves.size(), false);
    double minim = -1;
    
    for (int i = 0; i < (int)centroids.size(); ++i) {
        centroids[i].print_curve();

        for (int j = i + 1; j < (int)centroids.size(); ++j) {
            double dist = compute_distance(centroids[i], centroids[j], "DFT");
            
            if (minim == -1 || dist < minim) {
                minim = dist;
            }
        }
    }

    double R = minim / 2.0;
    cout << "R = " << R << endl;
    
    while (1) {
        general_search(hashtables, L, delta, k, R, "classic", "DFT", R_closest_curves, centroids, grid_curves_found, visited);
        bool found = false;

        for (int i = 0; i < (int)R_closest_curves.size(); ++i) {
            for (set<Curve>::const_iterator itr = R_closest_curves[i].begin(); itr != R_closest_curves[i].end(); ++itr) {
                int id = itr->get_int_id();
                
                itr->print_curve();

                if (assignment[id] == -1) {
                    assignment[id] = i;
                    found = true;
                } else {
                    int dist_1 = compute_distance(*itr, centroids[i], "DFT");
                    int dist_2 = compute_distance(*itr, centroids[assignment[id]], "DFT"); 
                    
                    if (dist_1 < dist_2) {
                        assignment[id] = i;
                        found = true;
                    }
                }
            }
        }
        
        if (!found) {
            break;
        }

        R *= 2;
    }
    
    // assignment for long points
    for (int i = 0; i < (int)assignment.size(); ++i) {
        if (assignment[i] == -1) {
            double minim = -1;
            
            for (int j = 0; j < (int)centroids.size(); ++j) {
                double dist = compute_distance(input_curves[i], centroids[j], "DFT");
                
                if (minim == -1 || dist < minim) {
                    minim = dist;
                    assignment[i] = j;
                }
            }
        }
    }
    
    return assignment;
}

double swap_update_centroid(const vector<const Curve*> &centroids, int old_centr, int new_centr, const vector<int> &assign, const vector<double> &close_dist, const vector<double> &close_dist_sec) {
    double value = 0;

    for (int i = 0; i < (int)assign.size(); ++i) {
        if ((i != old_centr && (centroids[assign[i]]->get_int_id() == i)) || i == new_centr) {
            continue;
        }

        double dist = compute_distance(input_curves[i], input_curves[new_centr], "DFT");

        if (centroids[assign[i]]->get_int_id() != old_centr) {
            value += min(dist, close_dist[i]);
        } else {
            value += min(dist, close_dist_sec[i]);
        }
    }

    return value;
}

bool PAM_update(vector<const Curve*> &centroids, double value, const vector<vector<int> > &clusters) { 
    vector<double> close_dist((int)input_curves.size(), -1), close_dist_sec((int)input_curves.size(), -1);
    vector<int> assign((int)input_curves.size());

    double min_value = value; 
    int p_new_cent, p_clust;
    
    for (int i = 0; i < (int)input_curves.size(); ++i) {
        for (int j = 0; j < (int)centroids.size(); ++j) {
            double dist = compute_distance(input_curves[i], *centroids[j], "DFT");

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
            
            double new_value = swap_update_centroid(centroids, centroids[i]->get_int_id(), clusters[i][j], assign, close_dist, close_dist_sec); 

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

        if (!centroids[i]->equal_curves(*mean_curve)) {
            centroids[i] = mean_curve;
            check = true;    
        }
    }

    return check;
}

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

