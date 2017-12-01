#ifndef HASH_FUNCTIONS_H
#define HASH_FUNCTIONS_H

#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <vector>
#include "curve.h"

using namespace std;

void grid_hashing(vector<Curve>&, double, const vector<double>&);
void grid_hashing_curve(Curve&, double, const vector<double>&, const Curve&);
void multiple_grids(vector<Curve>&, double, int);
void multiple_grids_curve(Curve&, double, int, int, const Curve&);
vector<double> prob_lsh_euclid_hashing(const vector<double>&, int, int, int);
long long vector_hashing(const vector<double>&);

extern vector<int> vec_r;
extern vector<vector<vector<double> > > t_shift, vec_line;

#endif
