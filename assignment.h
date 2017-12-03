#ifndef ASSIGNMENT_H
#define ASSIGNMENT_H

#include "curve.h"
#include "hashtable.h"
#include <vector>

using namespace std;

double loyd_assignment(const vector<const Curve*>&, vector<vector<int> >&, const char*);
double range_search(const vector<HashTable>&, const vector<const Curve*>&, vector<vector<int> >&, double, const char*);

#endif
