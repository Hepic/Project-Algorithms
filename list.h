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
        vector<Curve> search(const Curve&, const Curve&, const char*, const char*, double, bool = true) const;
};

#endif
