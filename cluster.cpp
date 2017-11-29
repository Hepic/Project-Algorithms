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

double loyd_assignment(const vector<int>& centroids, vector<int> &assign, vector<double> &close_dist, vector<double> &close_dist_sec) {
    double min_dist, min_dist_sec, dist, value = 0;
    int centr;
    
    for (int i = 0; i < (int)input_curves.size(); ++i) {
        min_dist = -1;
        
        for (int j = 0; j < (int)centroids.size(); ++j) {
            dist = compute_distance(input_curves[i], input_curves[centroids[j]], "DFT");

            if (min_dist == -1 || dist < min_dist) {
                min_dist_sec = min_dist;
                min_dist = dist;
                centr = centroids[j];
            } else if (min_dist_sec == -1 || dist < min_dist_sec) {
                min_dist_sec = dist;
            }
        }

        assign[i] = centr;
        close_dist[i] = min_dist;
        close_dist_sec[i] = min_dist_sec;

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
            int centr_id;
            
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

double swap_update_centroid(int old_centr, int new_centr, double value, const vector<int> &assign, const vector<double> close_dist, const vector<double> close_dist_sec) {
    for (int i = 0; i < (int)assign.size(); ++i) {
        if (i != old_centr && (assign[i] == i)) {
            continue;
        }
        
        double dist = compute_distance(input_curves[i], input_curves[new_centr], "DFT");
        
        if (assign[i] != old_centr) {
            value += min(dist - close_dist[i], 0.0);
        } else {
            value += min(dist, close_dist_sec[i]) - close_dist[i];
        }
    }

    return value;
}

bool PAM_update(vector<int> &centroids, const vector<int> &assign, const vector<double> &close_dist, const vector<double> &close_dist_sec, double value) {
    double min_value = value; 
    int pos_old_cent = -1, new_cent = -1;

    for (int i = 0; i < (int)centroids.size(); ++i) {
        for (int j = 0; j < (int)assign.size(); ++j) {
            if (centroids[i] == j || centroids[i] != assign[j]) {
                continue;
            }
            
            double new_value = swap_update_centroid(centroids[i], j, value, assign, close_dist, close_dist_sec); 

            if (new_value < min_value) {
                min_value = new_value;
                pos_old_cent = i;
                new_cent = j;
            }
        }
    }

    if (value > min_value) {
        centroids[pos_old_cent] = new_cent;
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
    
    vector<int> centroids = k_means_pp(input_curves.size(), 2, "DFT");
    cout << "initialization ended" << endl;
    
    do {
        value = loyd_assignment(centroids, assignment, close_dist, close_dist_sec);
        check = PAM_update(centroids, assignment, close_dist, close_dist_sec, value);  
    } while(check);
   
    cout << value << endl;

    for (int i = 0; i < (int)input_curves.size(); ++i) {
        cout << assignment[i] << " ";
    }
    
    cout << endl;
}
