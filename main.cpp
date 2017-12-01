#include <iostream>
#include <cstdlib>
#include "file_functions.h"
#include "help_functions.h"
#include "cluster.h"
#include "distances.h"

using namespace std;

int main(int argc, const char *argv[]) {
    ios::sync_with_stdio(false);
    srand(time(NULL));
    
    char *input_file = get_arguments(argv, argc, "-i");
    char *conf_file = get_arguments(argv, argc, "-c");
    char *output_file = get_arguments(argv, argc, "-o");
    char *metric = get_arguments(argv, argc, "-d");
    
    int dim = 2, num_of_clusters, k = 2, L = 3;
    read_file(input_file, dim);
    read_configuration_file(conf_file,num_of_clusters,k,L);

    clock_t begin = clock();

    mem_distance = new double*[(int)input_curves.size()];

    for (int i = 0; i < (int)input_curves.size(); ++i) {
        mem_distance[i] = new double[(int)input_curves.size()];

        for (int j = 0; j < (int)input_curves.size(); ++j) {
            mem_distance[i][j] = -1;
        }
    }

    clustering(metric);

    clock_t end = clock();
    double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
    
    cout << elapsed_secs << endl;
    
    return 0;
}
