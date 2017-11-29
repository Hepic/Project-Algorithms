#include <cstdio>
#include <cstdlib>
#include <cstring>
#include "list.h"
#include "distances.h"

List::Node::Node(const Curve &curve, const Curve &grid_curve) {
    this->curve = curve;
    this->grid_curve = grid_curve;
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

void List::insert(const Curve &curve, const Curve &grid_curve) {
    Node *new_node = new Node(curve, grid_curve);

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
        cout << "Curve: " << endl;
        node->curve.print_curve();
        
        cout << "Grid Curve: " << endl;
        node->grid_curve.print_curve();

        node = node->next;
    }

    cout << endl;
}

vector<Curve> List::search(const Curve &curve, const Curve &grid_curve, const char *hash_function, const char *dist_function, double R, bool check) const {
    Node *node = head;
    vector<Curve> min_curve;
    double min_dist = -1, dist;
    
    while (node != NULL) {
        bool cmp = false;
        
        if (check) {
            cmp = grid_curve.equal_curves(node->grid_curve);
        }
        
        if (cmp || !check) {
            dist = compute_distance(node->curve, curve, dist_function);
            
            if ((R && dist < R) || (!R && (min_dist == -1 || dist < min_dist))) {
                min_dist = dist;

                if (!R) {
                    if (!min_curve.empty()) {
                        min_curve.front() = node->curve;
                    } else {
                        min_curve.push_back(node->curve);
                    }
                    
                    min_dist = dist;
                } else {
                    min_curve.push_back(node->curve);
                }
            }
        }
        
        node = node->next;
    }

    return min_curve;
}
