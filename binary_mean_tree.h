#ifndef BINARY_MEAN_TREE_H
#define BINARY_MEAN_TREE_H

#include <math.h>
#include "curve.h"

class BinaryMeanTree {
    class Node {
        public:
            Curve curve;
            Node *left, *right;
    };

    Node *head;
    int pos;

    public:
        BinaryMeanTree(const vector<int>&);
        void build_tree(Node*, int, int, const vector<int>&);
        const Curve* get_mean() const;
};

#endif
