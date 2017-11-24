#ifndef CLUSTER_H
#define CLUSTER_H

#include "curve.h"
#include "hashtable.h"
#include <vector>

using namespace std;

vector<int> k_random_selection(int, int);
void k_means_pp(vector<const Curve*>&, int, int, const char*);
double loyd_assignment(const vector<const Curve*>&, vector<vector<int> >&);
vector<int> range_search(const vector<HashTable> &, const vector<Curve>&, int, int, int, double);
bool PAM_update(vector<int>&, double, const vector<int>&);
bool mean_frechet_update(vector<const Curve*>&, const vector<vector<int> >&);
void clustering(const vector<HashTable> &hashtables, int, int, double);

#endif
