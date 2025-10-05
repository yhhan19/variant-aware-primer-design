#include "bst.hpp"

BSTNode::BSTNode() {

}

BSTNode::~BSTNode() {

}

BST::BST(std::vector<std::pair<index_t, risk_t>> &input) {
    if (input.size() == 0) 
        root = NULL;
    else 
        root = build_bst(input, 0, input.size() - 1);
}

BSTNode *BST::build_bst(std::vector<std::pair<index_t, risk_t>> &input, index_t left, index_t right) {
    index_t mid = ((std::size_t)left + right) / 2;
    BSTNode *root = new_node(input[mid].first, input[mid].second);
    if (left <= mid - 1) root->left = build_bst(input, left, mid - 1);
    if (right >= mid + 1) root->right = build_bst(input, mid + 1, right);
    update(root);
    return root;
}

BST::~BST() {
    free_bst(root);
}

void BST::free_bst(BSTNode *root) {
    if (root == NULL) return ;
    free_bst(root->left);
    free_bst(root->right);
    delete root;
}

void BST::update(BSTNode *root) {
    root->min_key = root->key;
    root->min_index = root->lower = root->upper = root->index;
    if (root->left != NULL) {
        root->lower = root->left->lower;
        if (root->left->min_key < root->min_key) {
            root->min_key = root->left->key;
            root->min_index = root->left->index;
        }
    }
    if (root->right != NULL) {
        root->upper = root->right->upper;
        if (root->right->min_key < root->min_key) {
            root->min_key = root->right->key;
            root->min_index = root->right->index;
        }
    }
}

BSTNode *BST::new_node(index_t r, risk_t risk) {
    BSTNode *node = new BSTNode();
    node->index = node->lower = node->upper = node->min_index = r;
    node->min_key = node->key = risk;
    node->left = node->right = NULL;
    return node;
}

std::pair<index_t, risk_t> BST::range_query(BSTNode *root, index_t lower, index_t upper) {
    if (lower <= root->lower && upper >= root->upper) 
        return std::make_pair(root->min_index, root->min_key);
    index_t ret_index = 0;
    risk_t ret_key = INF;
    if (root->left != NULL && lower <= root->left->upper) {
        std::pair<index_t, risk_t> temp = range_query(root->left, lower, upper);
        if (temp.second < ret_key) {
            ret_index = temp.first;
            ret_key = temp.second;
        }
    }
    if (lower <= root->index && upper >= root->index) {
        if (root->key < ret_key) {
            ret_index = root->index;
            ret_key = root->key;
        }
    }
    if (root->right != NULL && upper >= root->right->lower) {
        std::pair<index_t, risk_t> temp = range_query(root->right, lower, upper);
        if (temp.second < ret_key) {
            ret_index = temp.first;
            ret_key = temp.second;
        }
    }
    return std::make_pair(ret_index, ret_key);
}

std::pair<index_t, risk_t> BST::range_query(index_t lower, index_t upper) {
    return range_query(root, lower, upper);
}
