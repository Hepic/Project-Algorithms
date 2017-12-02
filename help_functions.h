#ifndef HELP_FUNCTIONS_H
#define HELP_FUNCTIONS_H

#include <iostream>
#include <string>
#include <vector>
#include <set>
#include "curve.h"
#include "hashtable.h"
#include "hash_functions.h"
#include "distances.h"

using namespace std;

#define MAX_R 100

char* get_arguments(const char*[], int, const char*, bool = false);
void perror_exit(const char*);
double uniform_distribution(double, double);
double normal_distribution(double, double);
double dot_product(const vector<double>&, const vector<double>&);
void insert_curves_into_hashtables(vector<HashTable>&, double, const char*);
void search_curves_from_hashtables(const vector<HashTable>&, double, double, const char*, const char* , vector<set<Curve> >&, const vector<bool>&, const vector<Curve>&, vector<bool>&, bool = true);
void general_search(const vector<HashTable>&, double, double, const char*, const char*, vector<set<Curve> >&, const vector<Curve>&, vector<bool>&, vector<bool>&);

#endif
