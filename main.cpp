#include "utility.hpp"
#include "risk_optimizer.hpp"

int main(int argc, char **argv) {

    auto seed = time(0);
    srand(seed);
    auto input = random_risk(50000);
    // read_rate("../data/fmd.acc.msa.fasta.rate", input);
    index_t len_amp = 800;
    std::cout << "target length: " << input.size() << " bp" << std::endl;
    std::cout << "max amplicon length: " << len_amp << " bp" << std::endl;
    std::cout << "random seed: " << seed << std::endl;
    std::cout << std::endl;

    std::chrono::high_resolution_clock::time_point t1, t2;

    std::cout << "DP-based relaxed convex optimizer" << std::endl;
    t1 = std::chrono::high_resolution_clock::now();
    RiskOptimizer ro(input, 40, len_amp, len_amp);
    auto PDR = ro.search(0, 20000, 5);
    ro.validate_PDR(PDR);
    risk_t score = ro.score(PDR, ALPHA);
    t2 = std::chrono::high_resolution_clock::now();
    auto time = std::chrono::duration_cast<std::chrono::seconds>(t2 - t1).count();
    std::cout << "loss: " << score << ", ";
    std::cout << "time: " << time << " secs" << std::endl;
    std::cout << std::endl;

    std::cout << "Olivar greedy random optimizer" << std::endl;
    t1 = std::chrono::high_resolution_clock::now();
    RiskOptimizer ro_(input, 40, len_amp, len_amp);
    auto PDR_ = ro_.random_search(ITER_LIMIT * input.size() / len_amp);
    ro_.validate_PDR(PDR_);
    risk_t score_ = ro_.score(PDR_, ALPHA);
    t2 = std::chrono::high_resolution_clock::now();
    auto time_ = std::chrono::duration_cast<std::chrono::seconds>(t2 - t1).count();
    std::cout << "loss: " << score_ << ", ";
    std::cout << "time: " << time_ << " secs" << std::endl;
    std::cout << std::endl;

    std::cout << "improvement: " << (int) ((score_ - score) / score_ * 100) << "%" << std::endl;
    return 0;
}
