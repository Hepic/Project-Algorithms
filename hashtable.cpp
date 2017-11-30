#include <cstdio>
#include <cstdlib>
#include <cstring>
#include "hashtable.h"
#include "hash_functions.h"
#include "help_functions.h"
#include "curve.h"

HashTable::HashTable(int len) {
    this->id = -1;
    this->len = len;
    this->table = new List[len+5];
}

HashTable::HashTable(const HashTable &obj) {
    this->id = obj.id;
    this->len = obj.len;
    this->table = new List[len+5];
}

HashTable::~HashTable() {
    delete[] table;
}

void HashTable::set_id(int id) {
    this->id = id;
}

int HashTable::get_id() const {
    return id;
}

long long HashTable::get_index(const Curve &curve, const char *hash_function) const {
    vector<double> vec = curve.get_convert_vector();
    vector<double> vec_lsh, line;
    vector<vector<double> > table_lines;
    int k_vec = 3;
    long long index;

    if (!strcmp(hash_function, "classic")) {
        index = vector_hashing(vec) % len;
    } else if(!strcmp(hash_function, "probabilistic")) {
        if ((int)vec_line.size() <= id) { // build k_vec lines for that hashtable once
            for (int i = 0; i < k_vec; ++i) {
                line.clear();

                for (int j = 0; j < (int)vec.size(); ++j) {
                    double val = normal_distribution(0, 1);
                    line.push_back(val);
                }

                table_lines.push_back(line);
            }
            
            vec_line.push_back(table_lines);
        }
        
        vec_lsh = prob_lsh_euclid_hashing(vec, len, k_vec, id);
        index = (long long)vector_hashing(vec_lsh) % (long long)len;
    }
    
    return index;
}

void HashTable::insert(const Curve &curve, const Curve &grid_curve, const char *hash_function) {
    int index = get_index(grid_curve, hash_function);
    table[index].insert(curve, grid_curve);
}

vector<Curve> HashTable::search(const Curve &curve, const Curve &grid_curve, const char *hash_function, const char *dist_function, double R, bool check) const {
    int index = get_index(grid_curve, hash_function);
    vector<Curve> closer_curve = table[index].search(curve, grid_curve, hash_function, dist_function, R, check);

    return closer_curve;
}

void HashTable::print_bucket(int index) const {
    table[index].print_list();
}
