#include <iostream>
#include <cstdlib>
#include "file_functions.h"
#include "cluster.h"
#include "distances.h"

using namespace std;

int main(int argc, const char *argv[]) {
    ios::sync_with_stdio(false);
    srand(time(NULL));
    
    clock_t begin = clock();
    
    int dim = 2;
    read_file("test_1.txt", dim); 

    mem_distance = new double*[(int)input_curves.size()];

    for (int i = 0; i < (int)input_curves.size(); ++i) {
        mem_distance[i] = new double[(int)input_curves.size()];

        for (int j = 0; j < (int)input_curves.size(); ++j) {
            mem_distance[i][j] = -1;
        }
    }
    
    clustering();
    
    clock_t end = clock();
    double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
    
    cout << elapsed_secs << endl;

    return 0;
}
