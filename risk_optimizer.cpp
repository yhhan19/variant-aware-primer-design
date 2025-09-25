#include "risk_optimizer.hpp"

RiskOptimizer::RiskOptimizer(std::vector<risk_t> input, index_t len, index_t min, index_t max) {
    size = input.size();
    risk = new risk_t[size];
    prefix_sum = new risk_t[size + 1];
    exclude = new bool[size];
    prefix_sum[0] = 0;
    for (index_t i = 0; i < input.size(); i ++) {
        risk[i] = input[i];
        prefix_sum[i + 1] = prefix_sum[i] + risk[i];
    }
    this->len = len;
    this->min = min;
    this->max = max;
    std::vector<std::pair<index_t, risk_t>> c;
    for (index_t i = 0; i < input.size(); i ++) {
        if (i + len >= size) continue;
        c.push_back(std::make_pair(i, cost(i, 0, 0)));
        exclude[i] = true;
    }
    std::sort(c.begin(), c.end(), [](const std::pair<index_t, risk_t> &a, const std::pair<index_t, risk_t> &b) {
        return a.second < b.second;
    });
    for (index_t i = 0; i < input.size() / 10 * 8; i ++) {
        exclude[c[i].first] = false;
    }
}

risk_t RiskOptimizer::cost(index_t p, risk_t u, risk_t alpha) {
    risk_t sum = prefix_sum[p + len] - prefix_sum[p];
    risk_t sum2 = sum * sum - u;
    if (sum2 < 0) return u * alpha;
    return sum2 + u * alpha;
}

RiskOptimizer::~RiskOptimizer() {
    delete [] risk;
    delete [] prefix_sum;
    delete [] exclude;
}

std::vector<index_t> RiskOptimizer::random_PDR() {
    std::vector<index_t> PDR;
    while (true) {
        index_t k = PDR.size();
        index_t cmin, cmax;
        if (k % 2 == 0) {
            if (k == 0) {
                cmin = 0;
                cmax = 3 * len - 1;
            }
            else if (k == 2) {
                cmin = std::max(PDR[k - 2] + len, PDR[k - 1] + 3 * len - max);
                cmax = PDR[k - 1] - 1;
            }
            else {
                cmin = std::max(PDR[k - 3] + len, PDR[k - 1] + 3 * len - max);
                cmax = PDR[k - 1] - 1;
            }
        }
        else {
            if (k == 1) {
                cmin = PDR[k - 1] + min - len;
                cmax = PDR[k - 1] + max - 1;
            }
            else {
                cmin = std::max(PDR[k - 1] - len + min, PDR[k - 2] + 2 * len);
                cmax = PDR[k - 1] + max - 1;
            }
            if (cmax >= size) {
                PDR.pop_back();
                break;
            }
        }
        PDR.push_back(greedy_random_between(cmin, cmax - len + 2));
    }
    return PDR;
}

std::vector<index_t> RiskOptimizer::random_search(std::size_t limit) {
    risk_t min = -1;
    std::vector<index_t> min_PDR;
    for (std::size_t i = 0; i < limit; i ++) {
        std::vector<index_t> PDR = random_PDR();
        risk_t temp = score(PDR, ALPHA);
        if (min < 0 || temp < min) {
            min = temp;
            min_PDR = PDR;
            std::cout << "iteration: " << i << ", updated: " << min << std::endl;
        }
    }
    return min_PDR;
}

index_t RiskOptimizer::greedy_random_between(index_t cmin, index_t cmax) {
    std::vector<std::pair<index_t, risk_t>> all_sum;
    for (index_t i = cmin; i < cmax; i ++) {
        risk_t sum = 0;
        for (index_t j = i; j < i + len; j ++) 
            sum += risk[j];
        all_sum.push_back(std::make_pair(i, sum * sum));
    }
    std::sort(all_sum.begin(), all_sum.end(), [](const std::pair<index_t, risk_t> &a, const std::pair<index_t, risk_t> &b) {
        return a.second < b.second;
    });
    index_t bound = (cmax - cmin) * 3 / 10;
    if (bound == 0) bound = 1;
    return all_sum[random_between(0, bound)].first;
}

void RiskOptimizer::validate_PDR(std::vector<index_t> PDR) {
    for (index_t i = 0; i < PDR.size(); i ++) {
        if (i % 2 == 0) {
            if (i == 0) {
                continue;
            }
            else if (i == 2) {
                assert(PDR[i] - PDR[i - 2] >= len);
            }
            else {
                assert(PDR[i] - PDR[i - 3] >= len);
            }
            index_t l = PDR[i + 1] - PDR[i];
            assert(l + len <= max);
            assert(l + len >= min);
        }
        else {
            if (i == PDR.size() - 1) {
                assert(PDR[i] - PDR[i - 2] >= len);
            }
            else {
                assert(PDR[i] - PDR[i + 1] >= len);
            }
        }
    }
    std::cout << "uncovered target: " << PDR[0] << " (front), " << size - PDR[PDR.size() - 1] << " (rear)" << std::endl;
    assert(PDR[0] <= len * 2);
    // assert(size - PDR[PDR.size() - 1] <= 3 * len);
    assert(size - PDR[PDR.size() - 1] >= len);
}

risk_t RiskOptimizer::score(std::vector<index_t> PDR, risk_t alpha) {
    std::vector<risk_t> all_sum;
    for (index_t i = 0; i < PDR.size(); i ++) {
        risk_t sum = 0;
        for (index_t j = PDR[i]; j < PDR[i] + len; j ++) 
            sum += risk[j];
        all_sum.push_back(sum * sum);
    }
    std::sort(all_sum.begin(), all_sum.end());
    risk_t total = 0;
    for (index_t i = 0; i < std::ceil(alpha * PDR.size()); i ++) 
        total += all_sum[PDR.size() - i - 1];
    return total;
}

risk_t RiskOptimizer::opt(index_t f, index_t r, index_t f_, index_t r_, risk_t u, risk_t alpha) {
    if (exclude[f] || exclude[r] || exclude[f_] || exclude[r_]) return INF; 
    if (r < r_ + len || r_ < f + len || f < f_ + len) return INF;
    key_t k = to_key(f, r, f_, r_);
    if (memo.find(k) == memo.end()) {
        risk_t min_risk;
        key_t min_key;
        if (f_ <= len * 2) {
            min_risk = cost(f, u, alpha) + cost(r, u, alpha) + cost(f_, u, alpha) + cost(r_, u, alpha);
            min_key = 0;
        }
        else {
            min_risk = INF;
            for (index_t r__ = f_ + len; f >= r__ + len; r__ ++) {
                for (index_t f__ = (r__ >= max ? r__ - max : 0); r__  >= f__ + min; f__ ++) {
                    risk_t temp = opt(f_, r_, f__, r__, u, alpha);
                    if (temp < min_risk) {
                        min_risk = temp;
                        min_key = to_key(f_, r_, f__, r__);
                    }
                }
            }
            min_risk += cost(f, u, alpha) + cost(r, u, alpha);
        }
        memo[k] = min_risk;
        prev[k] = min_key;
    }
    return memo[k];
}

risk_t RiskOptimizer::top_k_opt(risk_t u, std::vector<index_t> &min_PDR, risk_t alpha) {
    min_PDR.clear();
    memo.clear();
    prev.clear();
    risk_t min_risk = INF;
    key_t min_key = 0;
    for (index_t r = size - len * 3; r <= size - len; r ++) {
        for (index_t f = (r >= max ? r - max : 0); r >= f + min; f ++) {
            for (index_t r_ = f + len; r >= r_ + len; r_ ++) {
                for (index_t f_ = (r_ >= max ? r_ - max : 0); r_ >= f_ + min; f_ ++) {
                    risk_t temp = opt(f, r, f_, r_, u, alpha);
                    if (temp < min_risk) {
                        min_risk = temp;
                        min_key = to_key(f, r, f_, r_);
                    }
                }
            }
        }
    }
    while (min_key != 0) {
        index_t a, b, c, d;
        to_index(min_key, a, b, c, d);
        min_PDR.push_back(b);
        min_PDR.push_back(a);
        min_key = prev[min_key];
        if (min_key == 0) {
            min_PDR.push_back(d);
            min_PDR.push_back(c);
        }
    }
    std::reverse(min_PDR.begin(), min_PDR.end());
    index_t i = min_PDR.size();
    return min_risk;
}

std::vector<index_t> RiskOptimizer::search(risk_t lower, risk_t upper, risk_t eps) {
    max -= len;
    min -= len;
    std::vector<index_t> min_PDR;
    risk_t alpha = ALPHA;
    risk_t l = lower, r = upper;
    while (r - l > eps) {
        risk_t u = (l + r) / 2;
        risk_t temp = top_k_opt(u, min_PDR, alpha), grad = alpha * min_PDR.size();
        std::cout << "parameter: " << u << ", relaxed: " << temp << std::endl;
        for (index_t i = 0; i < min_PDR.size(); i ++) 
            if (u < cost(min_PDR[i], 0, 0)) grad --;
        *(&(grad < 0 ? l : r)) = u;
    }
    max += len;
    min += len;
    return min_PDR;
}
