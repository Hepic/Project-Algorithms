#ifndef CLUSTER_H
#define CLUSTER_H

#include "curve.h"
#include "hashtable.h"
#include <vector>

using namespace std;

void clustering(const vector<HashTable> &hashtables, int, int, double);

#endif
