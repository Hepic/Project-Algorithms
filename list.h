#ifndef LIST_H
#define LIST_H

#include "curve.h"

class List {
    class Node {
        public:
            int val;
            Node *next;

            Node(int val);
    };

    Node *head, *tail;

    public:
        List();
        ~List();
        void insert(int val);
        void print_list() const;
};

#endif
