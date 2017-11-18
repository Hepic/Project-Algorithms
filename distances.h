#ifndef DISTANCES_H
#define DISTANCES_H

#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <vector>
#include "curve.h"

using namespace std;

double euclidean_distance(const vector<double>&, const vector<double>&);
double discrete_frechet_distance(const Curve&, const Curve&);
double dynamic_time_wrapping(const Curve&, const Curve&);
double compute_distance(const Curve&, const Curve&, const char*);

#endif
