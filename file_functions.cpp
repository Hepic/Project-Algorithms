#include <iostream>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <string>
#include <vector>
#include <string>
#include "file_functions.h"

using namespace std;

int num_of_clusters;
int global_k;
int global_L;

void read_file(const char *file_name, int &dim) {
    vector<double> point;
    ifstream file(file_name);
    string str, id;
    int num_points, int_id = 0;
    char chr;
    bool read_id = false;
    
    file >> str;

    if (str == "dimension:") {
        file >> dim;
    } else {
        read_id = true;
        id = str;
    }

    point.resize(dim);

    while (1) {
        if (!read_id) {
            file >> id;
        } else {
            read_id = false;
        }
         
        if (file.eof()) {
            break;
        }
        
        file >> num_points;
        Curve curve(id, int_id++, dim);

	    for (int i = 0; i < num_points; i++) {
            for(int j = 0; j < dim; ++j) {
                file >> chr >> point[j];
            }
            
            file >> chr;
            
            if (i < num_points - 1) {
                file >> chr;
            }

            if (!curve.is_empty() && curve.get_last_point() == point) { // remove duplicates
                continue;
            }

            curve.insert_point(point);
        }
        
        input_curves.push_back(curve);
    }
    file.close();
}

void read_configuration_file(const char *file_name) {
    ifstream file(file_name);
    string str;
    
    while(1) {
        if (file.eof()) {
            break;
        }
        
        file >> str;
        if (str == "number_of_clusters:") {
            file >> num_of_clusters;
        } else if (str == "number_of_grid_curves:") {
            file >> global_k;
        } else if (str == "number_of_hash_tables:") {
            file >> global_L;
        }
    }
    file.close();
}

void print_file(const char *file_name,const char *metric, vector<double> &silhouette_cluster, double time, vector<const Curve*> &centroids, vector<vector<int> > &clusters, int dim, const char *complete) {
    ofstream out_file(file_name);
    
    out_file << "I" << method_init << "A" << method_assign << "U" << method_update << "\n";
    out_file << "Metric: " << metric << "\n";
    
    for (int i = 0; i < clusters.size(); ++i) {
        out_file << "CLUSTER-" << i << " {size: " << clusters[i].size() << ", centroid: ";
        if (method_update == 2) {
            out_file << centroids[i]->get_id() << "}" << "\n";
        } else {
            out_file << "[";
            for (int j = 0; j < centroids[i]->get_length(); ++j) {
                out_file << "(";
                for (int k = 0; k < dim - 1; ++k) {
                    out_file << centroids[i]->get_coord_point(k,j) << ", ";
                }
                out_file << centroids[i]->get_coord_point(dim-1,j);
                if (j == centroids[i]->get_length() - 1) {
                    out_file << ")";
                } else {
                    out_file << "), ";
                }
            }
            out_file << "]}\n";
        }
    }
    
    out_file << "clustering_time: " << time << "\n";
    out_file << "Silhouette: [";
    double avg_s = 0;
    for (int i = 0; i < silhouette_cluster.size(); ++i) {
        avg_s += silhouette_cluster[i];
        out_file << silhouette_cluster[i] << ",";
    }
    avg_s /= silhouette_cluster.size();
    out_file << avg_s << "]" << "\n";
    
    if (strlen(complete)) {
        for (int i = 0; i < clusters.size(); ++i) {
            out_file << "CLUSTER-" << i << " {";
            for (int j = 0; j < clusters[i].size() - 1; ++j) {
                out_file << clusters[i][j] << ",";
            }
            out_file << clusters[i][clusters[i].size()-1] << "}\n";
        }
    }
    out_file.close();
}
