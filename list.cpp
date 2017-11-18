#include <cstdio>
#include <cstdlib>
#include <cstring>
#include "list.h"

List::Node::Node(int val) {
    this->val = val;
    this->next = NULL;
}

List::List() {
    head = tail = NULL;
}

List::~List() {
    Node *node = head, *tmp;

    while (node != NULL) {
        tmp = node;
        node = node->next;

        delete tmp;
    }
}

void List::insert(int val) {
    Node *new_node = new Node(val);

    if (head == NULL) {
        head = new_node;
        tail = head;
    } else {
        tail->next = new_node;
        tail = new_node;
    }
}

void List::print_list() const {
    Node *node = head;

    while (node != NULL) {
        cout << "Value: " << val << "\n";
        node = node->next;
    }

    cout << "\n";
}
