#ifndef LIST_H
#define LIST_H

#include "curve.h"

class List {
    class Node {
        public:
            Curve curve, grid_curve;
            Node *next;

            Node(const Curve&, const Curve&);
    };

    Node *head, *tail;

    public:
        List();
        ~List();
        void insert(const Curve&, const Curve&);
        void print_list() const;
        void search(vector<Curve>&, const Curve&, const Curve&, const char*, const char*, double, vector<bool>&, bool = true) const;
};

#endif
