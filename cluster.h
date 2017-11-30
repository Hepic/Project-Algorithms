#ifndef CLUSTER_H
#define CLUSTER_H

#include "curve.h"
#include <vector>

using namespace std;

vector<int> k_random_selection(int, int);
void k_means_pp(vector<int>&, vector<const Curve*>&, int, int, const char*);
double loyd_assignment(const vector<int>&, vector<int>&, vector<double>&, vector<double>&, vector<vector<int> >&);
double loyd_assignment(const vector<const Curve*>&, vector<vector<int> >&);
vector<int> range_search(const vector<int>&, int);
bool PAM_update(vector<int>&, const vector<int>&, const vector<double>&, const vector<double>&, double, const vector<int>&, int);
bool mean_frechet_update(vector<const Curve*>&, const vector<vector<int> >&);
void clustering();

#endif
