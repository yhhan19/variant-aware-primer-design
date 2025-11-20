#ifndef __RISK_OPTIMIZER__
#define __RISK_OPTIMIZER__

#include "utility.hpp"
#include "bst.hpp"

class RiskOptimizer {
    public:
        RiskOptimizer(std::vector<risk_t> input, index_t len, index_t min, index_t max);
        ~RiskOptimizer();
        std::vector<index_t> random_search(std::size_t limit);
        std::vector<index_t> search(risk_t l, risk_t r, risk_t eps);
        void validate_PDR(std::vector<index_t>);
        risk_t score(std::vector<index_t> PDR, risk_t alpha);
    private:
        std::unordered_map<key_t, risk_t> memo_umap;
        std::unordered_map<key_t, key_t> prev_umap;
        std::vector<std::pair<index_t, risk_t>> *solution_vector;
        BST **solutions;
        risk_t *memo;
        key_t *prev;
        index_t size, len, min, max;
        risk_t *risk, *prefix_sum;
        std::vector<index_t> random_PDR();
        index_t greedy_random_between(index_t cmin, index_t cmax);
        risk_t opt(index_t f, index_t r, index_t pf, index_t pr, risk_t u, risk_t alpha);
        risk_t cost(index_t p, risk_t u, risk_t alpha);
        risk_t top_k_opt(risk_t u, std::vector<index_t> &min_PDR, risk_t alpha);
        risk_t top_k_opt_fast(risk_t u, std::vector<index_t> &min_PDR, risk_t alpha);
        risk_t opt_m(index_t f, index_t r, risk_t u, risk_t alpha);
        risk_t top_k_opt_m(risk_t u, std::vector<index_t> &min_PDR, risk_t alpha);
        bool valid(index_t f, index_t r, index_t f_, index_t r_);
        risk_t top_k_opt_mi(risk_t u, std::vector<index_t> &min_PDR, risk_t alpha);
};

#endif
