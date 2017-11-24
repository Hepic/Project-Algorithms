#ifndef HASHTABLE_H
#define HASHTABLE_H

#include "list.h"

class HashTable {
    int id, len;
    List *table;

    public:
        HashTable(int);
        HashTable(const HashTable&);
        ~HashTable();
        void set_id(int);
        int get_id() const;
        long long get_index(const Curve &curve, const char*) const;
        void insert(const Curve&, const Curve&, const char*);
        void search(vector<Curve>&, const Curve&, const Curve&, const char*, const char*, double, vector<bool>&, bool = true) const;
        void print_bucket(int) const;
};

#endif
