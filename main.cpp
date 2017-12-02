#include <iostream>
#include <cstdlib>
#include "file_functions.h"
#include "help_functions.h"
#include "cluster.h"
#include "distances.h"
#include "hashtable.h"
#include "hash_functions.h"
#include "help_functions.h"

using namespace std;

int method_init;
int method_assign;
int method_update;


int main(int argc, const char *argv[]) {
    ios::sync_with_stdio(false);
    srand(time(NULL));
    
    char *input_file = get_arguments(argv, argc, "-i");
    char *conf_file = get_arguments(argv, argc, "-c");
    char *output_file = get_arguments(argv, argc, "-o");
    char *metric = get_arguments(argv, argc, "-d");
    char *complete = get_arguments(argv, argc, "-complete", true);
    method_init = atoi(get_arguments(argv, argc, "-init"));
    method_assign = atoi(get_arguments(argv, argc, "-assign"));
    method_update = atoi(get_arguments(argv, argc, "-update"));
    
    num_of_clusters = 3;
    global_k = 2;
    global_L = 3;
    
    int dim = 2; 
    double delta = 0.2;
    read_file(input_file, dim);
    read_configuration_file(conf_file);

    mem_distance = new double*[(int)input_curves.size()];

    for (int i = 0; i < (int)input_curves.size(); ++i) {
        mem_distance[i] = new double[(int)input_curves.size()];

        for (int j = 0; j < (int)input_curves.size(); ++j) {
            mem_distance[i][j] = -1;
        }
    }
    
    int max_points = 0;

    for (int i = 0; i < (int)input_curves.size(); ++i) {
        int value = input_curves[i].get_length();
        max_points = max(max_points, value);
    }
    
    for (int i = 0; i < max_points * global_k * dim; ++i) {
        vec_r.push_back(rand() % MAX_R);
    }
    
    int table_size = (int)input_curves.size() / 4 + 1;
    vector<HashTable> hashtables(global_L + 1, HashTable(table_size));
    
    for (int i = 0; i < global_L ; ++i) {
        hashtables[i].set_id(i);
    }
    
    vector<const Curve*> centroids;
    vector<vector<int> > clusters(num_of_clusters);
    vector<double> silhouette_cluster(num_of_clusters);
    
    if (method_assign == 2) {
        insert_curves_into_hashtables(hashtables, delta, "classic");
    }

    clock_t begin = clock();
    clustering(hashtables, delta, silhouette_cluster, centroids, clusters, metric);
    clock_t end = clock();

    double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
    print_file(output_file, metric, silhouette_cluster, elapsed_secs, centroids, clusters, dim, complete);

    return 0;
}
