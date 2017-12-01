#ifndef DISTANCES_H
#define DISTANCES_H

#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <vector>
#include "curve.h"

using namespace std;

static Curve default_curve;
extern double **mem_distance;

vector<double> find_closest_point(const vector<double>&, const vector<double>&, double);
double euclidean_distance(const vector<double>&, const vector<double>&);
vector<double> mean_point(const vector<double>&, const vector<double>&);
double discrete_frechet_distance(const Curve&, const Curve&, Curve& = default_curve, bool = false);
double dynamic_time_wrapping(const Curve&, const Curve&);
double compute_distance(int, int, const char*);

#endif
