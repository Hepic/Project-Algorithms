#include <iostream>
#include <cstdlib>
#include "file_functions.h"
#include "cluster.h"
#include "distances.h"
#include "hashtable.h"
#include "hash_functions.h"
#include "help_functions.h"

using namespace std;

int main(int argc, const char *argv[]) {
    ios::sync_with_stdio(false);
    srand(time(NULL));

    clock_t begin = clock();

    int dim = 2, k = 2, L = 2;
    double delta = 0.2;
    read_file("test.txt", dim); 

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
    
    for (int i = 0; i < max_points * k * dim; ++i) {
        vec_r.push_back(rand() % MAX_R);
    }
    
    int table_size = (int)input_curves.size() / 4 + 1;
    vector<HashTable> hashtables(L + 1, HashTable(table_size));
    
    for (int i = 0; i < L; ++i) {
        hashtables[i].set_id(i);
    }
    
    insert_curves_into_hashtables(hashtables, L, delta, k, "classic"); 
    clustering(hashtables, L, k, delta); 

    clock_t end = clock();
    double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
    
    cout << elapsed_secs << endl;

    return 0;
}
