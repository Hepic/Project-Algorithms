#include <iostream>
#include <cstdlib>
#include "cluster.h"
#include "curve.h"
#include "distances.h"
#include "help_functions.h"
#include "assignment.h"

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

vector<int> range_search(const vector<HashTable> &hashtables, const vector<Curve> &centroids, int dim, double delta) {
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
        general_search(hashtables, delta, R, "classic", "DFT", R_closest_curves, centroids, grid_curves_found, visited);
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
