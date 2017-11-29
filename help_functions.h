#ifndef HELP_FUNCTIONS_H
#define HELP_FUNCTIONS_H

#include <iostream>
#include <string>
#include "curve.h"
#include "hashtable.h"
#include "hash_functions.h"
#include "distances.h"
#include <vector>

using namespace std;

#define MAX_R 100

char *get_arguments(const char *[], int, const char*);
void perror_exit(const char*);
double uniform_distribution(double, double);
double normal_distribution(double, double);
double dot_product(const vector<double>&, const vector<double>&);
void insert_curves_into_hashtables(vector<HashTable>&, int, double, int, const char*);
void search_curves_from_hashtables(const vector<HashTable>&, int, double, int, double, const char*, const char* , vector<vector<Curve> >&, const vector<bool>&, vector<Curve>&, bool = true);
void general_search(const vector<HashTable>&, int, double, int, double, const char*, const char*, vector<vector<Curve> >&, vector<Curve>&, vector<bool>&);

#endif
