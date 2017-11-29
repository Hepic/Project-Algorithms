#ifndef HASH_FUNCTIONS_H
#define HASH_FUNCTIONS_H

#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <vector>
#include "curve.h"

using namespace std;

vector<Curve> grid_hashing(double, const vector<double>&);
Curve grid_hashing_curve(double, const vector<double>&, const Curve&);
vector<Curve> multiple_grids(double, int, int = -1);
Curve multiple_grids_curve(double, int, int, const Curve&);
vector<double> prob_lsh_euclid_hashing(const vector<double>&, int, int, int);
long long vector_hashing(const vector<double>&);

extern vector<int> vec_r;
extern vector<vector<vector<double> > > t_shift, vec_line;

#endif
