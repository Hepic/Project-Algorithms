#ifndef FILE_FUNCTIONS_H
#define FILE_FUNCTIONS_H

#include <vector>
#include "curve.h"

void read_file(const char*, int&);
void read_configuration_file(const char*);
void print_file(const char*, const char*, vector<double>&, double, vector<const Curve*> &centroids, vector<vector<int> >&, int, const char*);

#endif
