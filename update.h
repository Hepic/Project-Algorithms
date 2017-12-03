#ifndef UPDATE_H
#define UPDATE_H

#include "curve.h"
#include "hashtable.h"
#include <vector>

using namespace std;

bool PAM_update(vector<const Curve*>&, double, const vector<vector<int> >&, char*);
bool mean_frechet_update(vector<const Curve*>&, const vector<vector<int> >&);

#endif
