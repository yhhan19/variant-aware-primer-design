#include "utility.hpp"

class KPartiteGraph {
    public:
        KPartiteGraph(index_t K, index_t N);
        KPartiteGraph(std::string filename);
        ~KPartiteGraph();
        void random_weights();
        index_t part(index_t i) {return i / N;}
        index_t no(index_t i) {return i % N;}
        index_t index(index_t k, index_t n) {return k * N + n;}
        void display();
        weight_t cost(index_t *solution);
        weight_t cost(std::vector<index_t> solution);
        weight_t solve(std::size_t iter);
        std::vector<index_t> solve_fast(std::size_t iter);
        weight_t solve_sa(std::size_t iter);
    private:
        weight_t *vertices;
        weight_t **graph, **opt_graph, **opt_graph_2;
        index_t K, N;
        weight_t dfs(index_t k, index_t *solution, std::size_t *count);
        void dfs2(index_t k, index_t *solution, weight_t precost, weight_t *max, std::size_t *count);
        weight_t bound(index_t k, index_t *solution);
        void preprocessing();
        weight_t random_search(std::size_t iter);
        std::vector<index_t> random_search_single_fast(std::size_t restarts);
        std::vector<index_t> beam_search_fast(std::size_t beam_width, std::size_t max_iters);
        weight_t simulated_annealing(std::size_t restarts);
};
