#include <iostream>
#include <cstdlib>
#include "file_functions.h"

using namespace std;

int main(int argc, const char *argv[]) {
    srand(time(NULL));

    int dim = 2;
    read_file("trajectories_dataset", dim);
    
    //input_curves[0].print_curve();

    return 0;
}
