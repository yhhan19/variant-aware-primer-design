#include "graph.hpp"

KPartiteGraph::KPartiteGraph(index_t K, index_t N) {
    this->K = K;
    this->N = N;
    graph = new weight_t*[K * N];
    for (index_t i = 0; i < K * N; i ++) 
        graph[i] = new weight_t[K * N];
    opt_graph = new weight_t*[K];
    for (index_t i = 0; i < K; i ++) 
        opt_graph[i] = new weight_t[K];
    opt_graph_2 = new weight_t*[K * N];
    for (index_t i = 0; i < K * N; i ++) 
        opt_graph_2[i] = new weight_t[K];
    vertices = new weight_t[K * N];
}

KPartiteGraph::KPartiteGraph(std::string filename) {
    std::ifstream fin(filename);
    index_t E;
    fin >> K >> N >> E;
    graph = new weight_t*[K * N];
    for (index_t i = 0; i < K * N; i ++) {
        graph[i] = new weight_t[K * N];
        for (index_t j = 0; j < K * N; j ++) 
            graph[i][j] = 0;
    }
    vertices = new weight_t[K * N];
    for (index_t i = 0; i < K * N; i ++) {
        vertices[i] = 0;
    }
    index_t dummy;
    for (index_t i = 0; i < K * N; i ++) 
        fin >> dummy >> vertices[i];
    for (index_t k = 0; k < E; k ++) {
        index_t i, j;
        weight_t w;
        fin >> i >> j >> w;
        graph[i][j] = graph[j][i] = w;
    }
    fin.close();
}

KPartiteGraph::~KPartiteGraph() {
    for (index_t i = 0; i < K * N; i ++) 
        delete [] graph[i];
    delete [] graph;
    for (index_t i = 0; i < K; i ++) 
        delete [] opt_graph[i];
    delete [] opt_graph;
    for (index_t i = 0; i < K; i ++) 
        delete [] opt_graph_2[i];
    delete [] opt_graph_2;
}

void KPartiteGraph::random_weights() {
    for (index_t i = 0; i < K * N; i ++) {
        for (index_t j = i; j < K * N; j ++) {
            if (part(i) == part(j)) 
                graph[i][j] = graph[j][i] = 0;
            else 
                graph[i][j] = graph[j][i] = - (weight_t) rand() / RAND_MAX;
        }
    }
    for (index_t i = 0; i < K * N; i ++) {
        vertices[i] = - (weight_t) rand() / RAND_MAX;
    }
}

void KPartiteGraph::preprocessing() {
    for (index_t i = 0; i < K * N; i ++) {
        for (index_t j = 0; j < K; j ++) {
            weight_t temp = 0;
            for (index_t j_ = 0; j_ < N; j_ ++)
                temp = std::max(temp, graph[i][index(j, j_)]);
            opt_graph_2[i][j] = temp;
        }
    }
    for (index_t i = 0; i < K; i ++) {
        for (index_t j = i + 1; j < K; j ++) {
            weight_t temp = 0;
            for (index_t i_ = 0; i_ < N; i_ ++) 
                for (index_t j_ = 0; j_ < N; j_ ++) 
                    temp = std::max(temp, graph[index(i, i_)][index(j, j_)]);
            opt_graph[i][j] = opt_graph[j][i] = temp;
        }
    }
}

void KPartiteGraph::display() {
    for (index_t i = 0; i < K * N; i ++) {
        for (index_t j = 0; j < K * N; j ++) {
            std::cout << std::setw(10) << graph[i][j];
        }
        std::cout << std::endl;
    }
}

weight_t KPartiteGraph::cost(index_t *solution) {
    weight_t total = 0;
    for (index_t i = 0; i < K; i ++) {
        total += vertices[index(i, solution[i])];
        for (index_t j = i + 1; j < K; j ++) {
            index_t i_ = index(i, solution[i]), j_ = index(j, solution[j]);
            total += graph[i_][j_];
        }
    }
    return total;
}

weight_t KPartiteGraph::cost(std::vector<index_t> solution) {
    weight_t total = 0;
    for (index_t i = 0; i < K; i ++) {
        total += vertices[index(i, solution[i])];
        for (index_t j = i + 1; j < K; j ++) {
            index_t i_ = index(i, solution[i]), j_ = index(j, solution[j]);
            total += graph[i_][j_];
        }
    }
    return total;
}

weight_t KPartiteGraph::dfs(index_t k, index_t *solution, std::size_t *count) {
    (*count) ++;
    if (k == K) {
        return cost(solution);
    }
    else {
        weight_t max = -INF;
        for (index_t i = 0; i < N; i ++) {
            solution[k] = i;
            weight_t temp = dfs(k + 1, solution, count);
            if (temp > max) max = temp;
        }
        return max;
    }
}

std::vector<index_t> KPartiteGraph::solve_fast(std::size_t iter) {
    return random_search_single_fast(iter);
}

weight_t KPartiteGraph::solve(std::size_t iter) {
    index_t *solution = new index_t[K];
    std::size_t count = 0;
    return dfs(0, solution, &count);
}

weight_t KPartiteGraph::solve_sa(std::size_t iter) {
    return simulated_annealing(iter);
}

weight_t KPartiteGraph::random_search(std::size_t iter) {
    index_t *solution = new index_t[K];
    std::size_t i = 0;
    weight_t max = 0;
    while (i < iter) {
        for (index_t k = 0; k < K; k ++) 
            solution[k] = random_between(0, N);
        weight_t current = cost(solution);
        while (true) {
            index_t best_k = -1, best_n;
            weight_t best = current;
            for (index_t j = 0; j < K; j ++) {
                weight_t next = current;
                for (index_t k = 0; k < K; k ++) {
                    if (k == j) continue;
                    next -= graph[index(k, solution[k])][index(j, solution[j])];
                }
                for (index_t n = 0; n < N; n ++) {
                    weight_t delta = 0;
                    for (index_t k = 0; k < K; k ++) {
                        if (k == j) continue;
                        delta += graph[index(k, solution[k])][index(j, n)];
                    }
                    weight_t temp = next + delta;
                    if (temp > best + 1e-6) {
                        best = temp;
                        best_k = j;
                        best_n = n;
                    }
                }
            }
            if (best_k == -1) break;
            current = best;
            solution[best_k] = best_n;
        }
        i += 1;
        if (current > max) {
            max = current;
        }
    }
    delete [] solution;
    return max;
}

std::vector<index_t> KPartiteGraph::random_search_single_fast(std::size_t restarts) {
    // solution[p] stores the within-part index n (0..N-1)
    std::vector<index_t> solution(K), best_solution(K);

    // S[p] = sum_{q!=p} w( index(p,solution[p]), index(q,solution[q]) )
    std::vector<weight_t> S(K, 0);

    // contrib[p][n] = sum_{q!=p} w( index(p,n), index(q,solution[q]) )
    std::vector<std::vector<weight_t>> contrib(K, std::vector<weight_t>(N, 0));

    auto gv = [&](int p) -> index_t { return index(p, solution[p]); };
    auto w  = [&](index_t u, index_t v) -> weight_t { return graph[u][v]; };
    auto v = [&](index_t u) -> weight_t { return vertices[u]; };

    auto init_S_and_contrib = [&]() {
        // S
        for (int p = 0; p < (int)K; ++p) {
            weight_t s = 0; //v(gv(p)) * 2;
            for (int q = 0; q < (int)K; ++q) if (q != p)
                s += w(gv(p), gv(q));
            S[p] = s;
        }
        // contrib
        for (int p = 0; p < (int)K; ++p) {
            for (int n = 0; n < (int)N; ++n) {
                index_t v_ = index(p, n);
                weight_t s = 0; //v(v_) * 2;
                for (int q = 0; q < (int)K; ++q) if (q != p)
                    s += w(v_, gv(q));
                contrib[p][n] = s;
            }
        }
    };

    auto total_from_S = [&]() -> weight_t {
        // each edge counted twice in Σ_p S[p]
        long double sum2 = 0, sum1 = 0;
        for (int p = 0; p < (int)K; ++p) sum2 += (long double)S[p];
        for (int p = 0; p < (int)K; ++p) sum1 += (long double)v(gv(p));
        return (weight_t)(sum2 * 0.5) + sum1;
    };

    weight_t global_best = std::numeric_limits<weight_t>::lowest();

    for (std::size_t r = 0; r < restarts; ++r) {
        // random init: NOTE random_between must be [lo, hi] inclusive.
        for (int p = 0; p < (int)K; ++p) solution[p] = random_between(0, N - 1);

        init_S_and_contrib();
        weight_t current = total_from_S();

        while (true) {
            int best_p = -1, best_n = -1;
            weight_t best_gain = 0; // stop when no positive gain

            // Best-improvement scan with O(1) gain lookups
            for (int p = 0; p < (int)K; ++p) {
                for (int n = 0; n < (int)N; ++n) {
                    if (n == solution[p]) continue;
                    weight_t gain = contrib[p][n] - S[p] + v(index(p, n)) - v(gv(p));
                    if (gain > best_gain) {
                        best_gain = gain;
                        best_p = p;
                        best_n = n;
                    }
                }
            }

            if (best_p == -1) break; // local optimum

            // Apply move (best_p := best_n)
            const index_t old_v = gv(best_p);
            solution[best_p] = best_n;
            const index_t new_v = gv(best_p);
            current += best_gain;

            // 1) Update S for all parts (O(K))
            for (int q = 0; q < (int)K; ++q) if (q != best_p)
                S[q] += w(new_v, gv(q)) - w(old_v, gv(q));
            S[best_p] = contrib[best_p][best_n];

            // 2) Update contrib for all other parts/candidates (O(K·N))
            for (int q = 0; q < (int)K; ++q) if (q != best_p) {
                for (int m = 0; m < (int)N; ++m) {
                    const index_t u = index(q, m);
                    contrib[q][m] += w(u, new_v) - w(u, old_v);
                }
            }

            // 3) Rebuild contrib for the changed part (O(K·N))
            for (int n = 0; n < (int)N; ++n) {
                const index_t v_ = index(best_p, n);
                weight_t s = 0; // v(v_) * 2;
                for (int q = 0; q < (int)K; ++q) if (q != best_p)
                    s += w(v_, gv(q));
                contrib[best_p][n] = s;
            }
        }

        if (current > global_best) {
            global_best = current;
            best_solution = solution;
            // Optional: std::cout << "restart " << r << " best " << global_best << "\n";
        }
    }

    return best_solution;
}

std::vector<index_t>
KPartiteGraph::beam_search_fast(std::size_t beam_width,
                                std::size_t max_iters) {
    if (beam_width == 0) {
        return std::vector<index_t>(K, 0);
    }

    struct BeamNode {
        std::vector<index_t> sol;                 // solution[p] in 0..N-1
        std::vector<weight_t> S;                  // S[p]
        std::vector<std::vector<weight_t>> contrib; // contrib[p][n]
        weight_t score;                           // objective value
    };

    auto gv_idx = [&](const std::vector<index_t>& sol, int p) -> index_t {
        return index(p, sol[p]);
    };
    auto w = [&](index_t u, index_t v_) -> weight_t { return graph[u][v_]; };
    auto v = [&](index_t u) -> weight_t { return vertices[u]; };

    auto total_from_S = [&](const BeamNode& node) -> weight_t {
        long double sum2 = 0, sum1 = 0;
        for (int p = 0; p < (int)K; ++p) sum2 += (long double)node.S[p];
        for (int p = 0; p < (int)K; ++p)
            sum1 += (long double)v(gv_idx(node.sol, p));
        return (weight_t)(sum2 * 0.5L) + sum1;
    };

    auto init_node = [&](BeamNode& node) {
        node.S.assign(K, 0);
        node.contrib.assign(K, std::vector<weight_t>(N, 0));

        // S[p] = sum_{q != p} w(gv(p), gv(q))
        for (int p = 0; p < (int)K; ++p) {
            weight_t s = 0;
            index_t vp = gv_idx(node.sol, p);
            for (int q = 0; q < (int)K; ++q) if (q != p) {
                s += w(vp, gv_idx(node.sol, q));
            }
            node.S[p] = s;
        }
        // contrib[p][n] = sum_{q != p} w(index(p,n), gv(q))
        for (int p = 0; p < (int)K; ++p) {
            for (int n = 0; n < (int)N; ++n) {
                index_t cand = index(p, n);
                weight_t s = 0;
                for (int q = 0; q < (int)K; ++q) if (q != p) {
                    s += w(cand, gv_idx(node.sol, q));
                }
                node.contrib[p][n] = s;
            }
        }
        node.score = total_from_S(node);
    };

    // --- Initialize beam with random solutions ---
    std::vector<BeamNode> beam;
    beam.reserve(beam_width);

    for (std::size_t b = 0; b < beam_width; ++b) {
        BeamNode node;
        node.sol.resize(K);
        for (int p = 0; p < (int)K; ++p) {
            node.sol[p] = random_between(0, N - 1);
        }
        init_node(node);
        beam.push_back(std::move(node));
    }

    // Keep beam sorted by score (descending)
    auto by_score_desc = [](const BeamNode& a, const BeamNode& b) {
        return a.score > b.score;
    };
    std::sort(beam.begin(), beam.end(), by_score_desc);

    struct NeighborCandidate {
        weight_t score;
        std::size_t parent_idx;
        int part;
        int new_n;
    };

    for (std::size_t iter = 0; iter < max_iters; ++iter) {
        std::vector<NeighborCandidate> candidates;
        candidates.reserve(beam.size() * K * (N - 1));

        // --- Expand all neighbors of all beam states ---
        for (std::size_t bi = 0; bi < beam.size(); ++bi) {
            const BeamNode& node = beam[bi];

            for (int p = 0; p < (int)K; ++p) {
                int old_n = node.sol[p];
                index_t old_v = index(p, old_n);

                for (int n = 0; n < (int)N; ++n) {
                    if (n == old_n) continue;

                    // Gain from changing part p to candidate n in this node
                    weight_t gain = node.contrib[p][n] - node.S[p]
                                  + v(index(p, n)) - v(old_v);
                    weight_t new_score = node.score + gain;

                    candidates.push_back(
                        NeighborCandidate{new_score, bi, p, n});
                }
            }
        }

        if (candidates.empty()) break;

        // Sort all neighbors globally by score (descending)
        std::sort(candidates.begin(), candidates.end(),
                  [](const NeighborCandidate& a,
                     const NeighborCandidate& b) {
                      return a.score > b.score;
                  });

        // --- Build next beam from top-k neighbors ---
        std::vector<BeamNode> next_beam;
        next_beam.reserve(beam_width);

        for (std::size_t i = 0;
             i < candidates.size() && next_beam.size() < beam_width;
             ++i) {
            const auto& cand = candidates[i];
            const BeamNode& parent = beam[cand.parent_idx];

            BeamNode child;
            child.sol = parent.sol;
            child.sol[cand.part] = cand.new_n;

            // Compute S, contrib, and score for this child
            init_node(child);
            // Optionally, we could assert that child.score ≈ cand.score

            next_beam.push_back(std::move(child));
        }

        if (next_beam.empty()) break;

        std::sort(next_beam.begin(), next_beam.end(), by_score_desc);
        beam.swap(next_beam);
    }

    // Return the best solution from the final beam
    auto best_it = std::max_element(beam.begin(), beam.end(), by_score_desc);
    if (best_it != beam.end()) {
        return best_it->sol;
    }
    return std::vector<index_t>(K, 0);
}


void KPartiteGraph::dfs2(index_t k, index_t *solution, weight_t precost, weight_t *max, std::size_t *count) {
    (*count) ++;
    if (k == K) {
        if (precost > *max) *max = precost;
    }
    else {
        for (index_t i = 0; i < N; i ++) {
            solution[k] = i;
            weight_t delta = 0;
            for (index_t j = 0; j < k; j ++) 
                delta += graph[index(k, i)][index(j, solution[j])];
            weight_t postcost = bound(k, solution);
            if (postcost + precost + delta > *max)
                dfs2(k + 1, solution, precost + delta, max, count);
        }
    }
}

weight_t KPartiteGraph::simulated_annealing(std::size_t restarts) {
    // --- parameters (tweak as needed) ---
    const std::size_t steps_per_restart   = -1;
    const std::size_t proposals_per_step  = 200;
    const double alpha                    = 0.985;    // geometric cooling
    const double Tmin                     = 1e-6;

    // state
    std::vector<int> solution(K, 0);
    std::vector<weight_t> S(K, 0.0);                         // S[p] = Σ_{q≠p} w(sol[p],sol[q])
    std::vector<std::vector<weight_t>> contrib(              // contrib[p][n] = Σ_{q≠p} w(index(p,n),sol[q])
        K, std::vector<weight_t>(N, 0.0));

    auto w = [&](index_t u, index_t v) -> weight_t {
        return graph[u][v];
    };
    auto gv = [&](int p) -> index_t { return index(p, solution[p]); };

    auto init_random_solution = [&](){
        for (int p = 0; p < K; ++p) solution[p] = random_between(0, N - 1);
    };

    auto init_caches = [&](){
        // Build S and contrib for current solution
        for (int p = 0; p < K; ++p) {
            const index_t vp = gv(p);
            weight_t sp = 0.0;
            for (int q = 0; q < K; ++q) if (q != p) sp += w(vp, gv(q));
            S[p] = sp;

            for (int n = 0; n < N; ++n) {
                const index_t v = index(p, n);
                weight_t c = 0.0;
                for (int q = 0; q < K; ++q) if (q != p) c += w(v, gv(q));
                contrib[p][n] = c;
            }
        }
    };

    auto objective_from_S = [&]() -> weight_t {
        long double sum2 = 0.0L;
        for (int p = 0; p < K; ++p) sum2 += (long double)S[p];
        return (weight_t)(sum2 * 0.5L); // each edge counted twice
    };

    // crude auto-T0: average abs gain of a few random proposals
    auto estimate_T0 = [&]() -> double {
        const size_t samples = std::min<size_t>(1000, static_cast<size_t>(K) * std::max<size_t>(1, N));
        long double s = 0.0L;
        for (int t = 0; t < samples; ++t) {
            const int p = random_between(0, K - 1);
            int n = solution[p];
            if (N > 1) do { n = random_between(0, N - 1); } while (n == solution[p]);
            const weight_t g = contrib[p][n] - S[p];
            s += std::fabs((long double)g);
        }
        double avg = (samples ? (double)(s / samples) : 1.0);
        return (avg > 1e-9 ? avg : 1.0);
    };

    weight_t global_best = std::numeric_limits<weight_t>::lowest();

    for (std::size_t r = 0; r < restarts; ++r) {
        // init
        init_random_solution();
        init_caches();
        weight_t current = objective_from_S();
        weight_t best_this = current;

        double T = estimate_T0();

        for (std::size_t step = 0; step < steps_per_restart; ++step) {
            if (T < Tmin) break;

            for (std::size_t it = 0; it < proposals_per_step; ++it) {
                // propose: flip one part to a new candidate
                const int p = random_between(0, K - 1);
                int n = solution[p];
                if (N > 1) do { n = random_between(0, N - 1); } while (n == solution[p]);

                const weight_t gain = contrib[p][n] - S[p];
                bool accept = (gain >= 0.0);
                if (!accept && T > 0.0) {
                    // Metropolis for maximization: exp(gain / T)
                    const double u = (double)random_between(0, 1000000) / 1000000.0;
                    accept = (std::exp(std::min<weight_t>(gain / T, 40.0)) > u);
                }
                if (!accept) continue;

                // apply move p := n
                const index_t old_v = gv(p);
                solution[p] = n;
                const index_t new_v = gv(p);
                current += gain;

                // update S for all parts (O(K))
                for (int q = 0; q < K; ++q) if (q != p) {
                    const index_t vq = gv(q);
                    S[q] += w(new_v, vq) - w(old_v, vq);
                }
                S[p] = contrib[p][n];

                // update contrib for other parts (O(K·N))
                for (int q = 0; q < K; ++q) if (q != p) {
                    for (int m = 0; m < N; ++m) {
                        const index_t u = index(q, m);
                        contrib[q][m] += w(u, new_v) - w(u, old_v);
                    }
                }
                // rebuild contrib for changed part p (O(K·N))
                for (int m = 0; m < N; ++m) {
                    const index_t v = index(p, m);
                    weight_t c = 0.0;
                    for (int q = 0; q < K; ++q) if (q != p) c += w(v, gv(q));
                    contrib[p][m] = c;
                }

                if (current > best_this) best_this = current;
            }

            // cool
            T *= alpha;
        }

        if (best_this > global_best) global_best = best_this;
    }

    return global_best;
}

weight_t KPartiteGraph::bound(index_t k, index_t *solution) {
    weight_t total = 0;
    for (index_t i = 0; i <= k; i ++) {
        for (index_t j = k + 1; j < K; j ++) {
            total += opt_graph_2[index(i, solution[i])][j];
        }
    }
    for (index_t i = k + 1; i < K; i ++) {
        for (index_t j = i + 1; j < K; j ++) {
            total += opt_graph[i][j];
        }
    }
    return total;
}
