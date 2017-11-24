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

double loyd_assignment(const vector<int>& centroids, vector<int> &assignment, vector<int> &assignment_second) {
    double min_dist, min_dist_sec, dist, value = 0;
    int centr, centr_sec;
    
    for (int i = 0; i < (int)input_curves.size(); ++i) {
        min_dist = -1;
        
        for (int j = 0; j < (int)centroids.size(); ++j) {
            dist = compute_distance(input_curves[i], input_curves[centroids[j]], "DFT");

            if (min_dist == -1 || dist < min_dist) {
                min_dist_sec = min_dist;
                centr_sec = centr;
            
                min_dist = dist;
                centr = centroids[j];
            } else if (min_dist_sec == -1 || dist < min_dist_sec) {
                min_dist_sec = dist;
                centr_sec = centroids[j];
            }
        }

        assignment[i] = centr;
        assignment_second[i] = centr_sec;

        value += min_dist;
    }

    return value;
}

double swap_update_centroid(int old_centr, int new_centr, double value, const vector<int> &assignment, const vector<int> &assignment_second) {
    for (int i = 0; i < (int)assignment.size(); ++i) {
        if (i != old_centr && (assignment[i] == i)) {
            continue;
        }
        
        if (assignment[i] != old_centr) {
            double dist_1 = compute_distance(input_curves[i], input_curves[new_centr], "DFT");
            double dist_2 = compute_distance(input_curves[i], input_curves[assignment[i]], "DFT");

            if (dist_1 < dist_2) {
                value += (dist_1 - dist_2);
            }
        } else {
            double dist_1 = compute_distance(input_curves[i], input_curves[new_centr], "DFT");
            double dist_2 = compute_distance(input_curves[i], input_curves[assignment_second[i]], "DFT"); 
            double dist_3 = compute_distance(input_curves[i], input_curves[assignment[i]], "DFT");                    
            
            if (dist_1 < dist_2) {
                value += (dist_1 - dist_3);
            } else {
                value += (dist_2 - dist_3);
            }
        }
    }

    return value;
}

double PAM_update(const vector<int> &centroids, const vector<int> &assign, const vector<int> &assign_sec, double value, int &p_old_cent, int &new_cent) {
    double min_value = value; 
    p_old_cent = new_cent = -1;

    for (int i = 0; i < (int)centroids.size(); ++i) {
        for (int j = 0; j < (int)assign.size(); ++j) {
            if (centroids[i] == j || centroids[i] != assign[j]) {
                continue;
            }
            
            double new_value = swap_update_centroid(centroids[i], j, value, assign, assign_sec); 

            if (new_value < min_value) {
                min_value = new_value;
                p_old_cent = i;
                new_cent = j;
            }
        }
    }

    return min_value;
}

void clustering(vector<int> &centroids) {
    vector<int> assignment(input_curves.size()), assignment_second(input_curves.size());
    int pos_old_centr, new_centr;

    while (1) {
        double value = loyd_assignment(centroids, assignment, assignment_second);
        double new_value = PAM_update(centroids, assignment, assignment_second, value, pos_old_centr, new_centr);    

        if (value > new_value) {
            centroids[pos_old_centr] = new_centr; 
        } else {
            cout << value << endl;
            break;
        }
    }
   
    for (int i = 0; i < (int)input_curves.size(); ++i) {
        cout << assignment[i] << " ";
    }
    
    cout << endl;
}
