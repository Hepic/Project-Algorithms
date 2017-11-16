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
    read_file("trajectories_dataset", dim);
    
    vector<int> ret = k_means_pp(input_curves.size(), 3, "DFT");
    
    for (int i = 0; i < ret.size(); ++i) {
        cout << ret[i] << endl;
    }    
    
    clock_t end = clock();
    double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
    
    cout << elapsed_secs << endl;

    return 0;
}
