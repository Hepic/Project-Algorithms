#ifndef CLUSTER_H
#define CLUSTER_H

#include "curve.h"
#include "hashtable.h"
#include <vector>

using namespace std;

void clustering(const vector<HashTable> &hashtables, double, vector<double>&, vector<const Curve*> &centroids, vector<vector<int> >&, const char*);
void silhouette(const vector<const Curve*> &centroids, vector<vector<int> >&, vector<double>&, const char*);

#endif
