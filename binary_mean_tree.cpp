#include "distances.h"
#include "binary_mean_tree.h"

BinaryMeanTree::BinaryMeanTree(const vector<int> &cluster) {
    this->pos = 0;
    this->head = new Node;

    int height = ceil(log2(cluster.size()));
    build_tree(head, 0, height, cluster);
}

void BinaryMeanTree::build_tree(Node *node, int curr, int height, const vector<int> &cluster) { 
    if (curr < height) {
        node->left = new Node;
        node->right = new Node;

        build_tree(node->left, curr + 1, height, cluster);
        
        if (this->pos < (int)cluster.size()) {
            build_tree(node->right, curr + 1, height, cluster);
            discrete_frechet_distance(node->left->curve, node->right->curve, node->curve, true);   
        } else {
            node->curve = node->left->curve;
        }

    } else if (curr == height) {
        node->curve = input_curves[cluster[this->pos++]];
    }
}

const Curve* BinaryMeanTree::get_mean() const {
    return &head->curve;
}
