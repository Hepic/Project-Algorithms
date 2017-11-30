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

double loyd_assignment(const vector<int> &centroids, vector<int> &assign, vector<double> &close_dist, vector<double> &close_dist_sec, vector<vector<int> > &clusters) {
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

vector<int> range_search(const vector<int>& centroids, int dim) {
    vector<int> assignment(input_curves.size(),0);
    vector<int> prev_assignment(input_curves.size());
    vector<Curve> search_curves(centroids.size());
    vector<vector<Curve> > R_closest_curves(centroids.size());
    vector<bool> grid_curves_found(search_curves.size(), false);
    int table_size, max_points = 0, min_points = -1;
    string choice;
    double delta = 2.0, minim, dist;
    int k = 2, L = 3;
    int value;
    
    for (int i = 0; i < centroids.size(); ++i) {
        search_curves.push_back(input_curves[centroids[i]]);
    }
    
    for (int i = 0; i < input_curves.size(); ++i) {
        value = input_curves[i].get_length();
        max_points = max(max_points, value);
        min_points = (min_points == -1 ? value : min(min_points, value));
    }
    
    for (int i = 0; i < max_points * k * dim; ++i) {
        vec_r.push_back(rand() % MAX_R);
    }
    
    // Create L hashTables
    table_size = (int)input_curves.size() / 4 + 1;
    vector<HashTable> hashtables(L + 1, HashTable(table_size));
    
    for (int i = 0; i < L; ++i) {
        hashtables[i].set_id(i);
    }
    
    insert_curves_into_hashtables(hashtables, L, delta, k, "DFT"); // Insert input_curves into hashtables
    
    for (int i = 0; i < search_curves.size()-1; ++i) {
        for (int j = i+1; j < search_curves.size(); ++j) {
            dist = compute_distance(search_curves[i], search_curves[j], "DFT");
            
            if (minim == -1 || dist < minim) {
                minim = dist;
            }
        }
    }
    
    double R = (double)(minim/2);
    int first_time = 1;
    int id;
    
    while (1) {
        if (first_time == 0) {
            R = (double)(R*2);
        }
        
        general_search(hashtables, L, delta, k, R, "classic", "DFT", R_closest_curves, search_curves, grid_curves_found);
        
        for (int i = 0; i < R_closest_curves.size(); ++i) {
            for (int j = 0; j < R_closest_curves[i].size(); ++j) {
                id = stoi(R_closest_curves[i][j].get_id());
                
                if (assignment[id] == 0) {
                    assignment[id] = centroids[i];
                } else {
                    int dist_1 = compute_distance(R_closest_curves[i][j], search_curves[i], "DFT");
                    int dist_2 = compute_distance(R_closest_curves[i][j], input_curves[assignment[id]], "DFT"); //prev 
                    
                    if (dist_1 < dist_2) {
                        assignment[id] = centroids[id];
                    }
                }
            }
        }
        
        if (first_time == 0 && assignment == prev_assignment) {
            break;
        } else if (first_time == 1) {
            first_time = 0;
        }
    }
    
    // assignment for long points
    for (int i = 0; i < assignment.size(); ++i) {
        if (assignment[i] == 0) {
            double minim = -1;
            int centr_id; sσ
            
            for (int j = 0; j < search_curves.size(); ++j) {
                dist = compute_distance(input_curves[i], search_curves[j], "DFT");
                
                if (minim == -1 || dist < minim) {
                    minim = dist;
                    centr_id = centroids[j];
                }
            }
            
            assignment[i] = centr_id;
        }
    }
    
    return assignment;
}

double swap_update_centroid(int old_centr, int new_centr, const vector<int> &assign, const vector<double> &close_dist, const vector<double> &close_dist_sec) {
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

void clustering() {
    vector<const Curve*> centroids;
    vector<int> assignment(input_curves.size()), centroids_ind;
    vector<double> close_dist(input_curves.size()), close_dist_sec(input_curves.size());
    bool check;
    double value;
    int num_of_clusters = 2;
    
    k_means_pp(centroids_ind, centroids, input_curves.size(), num_of_clusters, "DFT");
    cout << "initialization ended" << endl;
    
    vector<vector<int> > clusters(num_of_clusters);
    
    do {
        for (int i = 0; i < num_of_clusters; ++i) {
            //value = loyd_assignment(centroids_ind, assignment, close_dist, close_dist_sec, clusters);
            value = loyd_assignment(centroids, clusters);
            //check = PAM_update(centroids_ind, assignment, close_dist, close_dist_sec, value, clusters[i], i); 
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