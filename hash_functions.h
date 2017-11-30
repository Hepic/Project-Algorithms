#ifndef HASH_FUNCTIONS_H
#define HASH_FUNCTIONS_H

#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <vector>
#include "curve.h"

using namespace std;

void grid_hashing(vector<Curve>&, double, const vector<double>&);
Curve grid_hashing_curve(double, const vector<double>&, const Curve&);
void multiple_grids(vector<Curve>&, double, int);
Curve multiple_grids_curve(double, int, int, const Curve&);
vector<double> prob_lsh_euclid_hashing(const vector<double>&, int, int, int);
long long vector_hashing(const vector<double>&);

extern vector<int> vec_r;
extern vector<vector<vector<double> > > t_shift, vec_line;

#endif
