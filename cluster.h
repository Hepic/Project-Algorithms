#ifndef CLUSTER_H
#define CLUSTER_H

#include <vector>

using namespace std;

vector<int> k_random_selection(int, int);
vector<int> k_means_pp(int, int, const char*);
vector<int> loyd_assignment(const vector<int>&);

#endif
