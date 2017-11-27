#ifndef CLUSTER_H
#define CLUSTER_H

#include "curve.h"
#include <vector>

using namespace std;

vector<int> k_random_selection(int, int);
vector<int> k_means_pp(int, int, const char*);
double loyd_assignment(const vector<int>&, vector<int>&, vector<double>&, vector<double>&, vector<vector<int> >&);
bool PAM_update(vector<int>&, const vector<int>&, const vector<double>&, const vector<double>&, double, const vector<int>&, int);
double mean_frechet(vector<Curve>&, const vector<int>&);
void clustering();

#endif
