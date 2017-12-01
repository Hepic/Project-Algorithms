#include <iostream>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <string>
#include <vector>
#include <string>
#include "file_functions.h"

using namespace std;

void read_file(const char *file_name, int &dim) {
    vector<double> point;
    ifstream file(file_name);
    string str, id;
    int num_points;
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
        Curve curve(id, dim);

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
