#ifndef INITIALIZATION_H
#define INITIALIZATION_H

#include "curve.h"
#include "hashtable.h"
#include <vector>

using namespace std;

void k_random_selection(vector<const Curve*>&, int);
void k_means_pp(vector<const Curve*>&, int, const char*);

#endif
