#ifndef CLUSTER_H
#define CLUSTER_H

#include <vector>

using namespace std;

vector<int> k_random_selection(int, int);
vector<int> k_means_pp(int, int, const char*);
double loyd_assignment(const vector<int>&, vector<int>&, vector<int>&);
double PAM_update(const vector<int>&, const vector<int>&, const vector<int>&, double, int&, int&);
void clustering(vector<int>&);

#endif
