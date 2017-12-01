#ifndef INITIALIZATION_H
#define INITIALIZATION_H

#include "curve.h"
#include "hashtable.h"
#include <vector>

using namespace std;

vector<int> k_random_selection(int, int);
void k_means_pp(vector<const Curve*>&, int, int, const char*);

#endif
