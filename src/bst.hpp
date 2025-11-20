#ifndef __BST__
#define __BST__

#include "utility.hpp"

class BSTNode {
    friend class BST;
    public:
        BSTNode();
        ~BSTNode();
    private:
        index_t index, lower, upper, min_index;
        risk_t key, min_key;
        BSTNode *left, *right;
};

class BST {
    public:
        BST(std::vector<std::pair<index_t, risk_t>>& input);
        ~BST();
        std::pair<index_t, risk_t> range_query(index_t lower, index_t upper);
    private:
        BSTNode *root;
        BSTNode *new_node(index_t r, risk_t risk);
        void update(BSTNode *root);
        void free_bst(BSTNode *root);
        BSTNode *build_bst(std::vector<std::pair<index_t, risk_t>> &input, index_t left, index_t right);
        std::pair<index_t, risk_t> range_query(BSTNode *root, index_t lower, index_t upper);
 
};

#endif
