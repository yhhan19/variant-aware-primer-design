#include "utility.hpp"
#include "risk_optimizer.hpp"
#include "graph.hpp"

int main(int argc, char **argv) {
    
    time_t seed = 42;
    std::string input_file, output_file, mode;
    index_t len_amp = 420, len_amp_min = 252, len_PDR = 40;
    risk_t u_max = 10000, u_min = 0.1;
    std::size_t iter = 1000;
    for (int i = 0; i < argc; i ++) {
        std::string to_parse(argv[i]);
        if (to_parse == "-i") input_file = std::string(argv[++ i]);
        if (to_parse == "-o") output_file = std::string(argv[++ i]);
        if (to_parse == "-x") mode = std::string(argv[++ i]);
        if (to_parse == "-S") seed = std::stoi(argv[++ i]);
        if (to_parse == "-Ln") len_amp = std::stoi(argv[++ i]);
        if (to_parse == "-Lx") len_amp_min = std::stoi(argv[++ i]);
        if (to_parse == "-Lp") len_PDR = std::stoi(argv[++ i]);
        if (to_parse == "-Ux") u_max = std::stof(argv[++ i]);
        if (to_parse == "-Un") u_min = std::stof(argv[++ i]);
        if (to_parse == "-I") iter = std::stoi(argv[++ i]);
    }
    srand(seed);
    std::string ref;
    std::vector<risk_t> rate;
    if (read_rate_ref(input_file, rate, ref) != 0) return 1;

    std::chrono::high_resolution_clock::time_point t1, t2;
    t1 = std::chrono::high_resolution_clock::now();
    if (mode == "risk-d1") {
        std::cout << "DP-based relaxed convex optimizer (s1)" << std::endl;
        RiskOptimizer ro(rate, len_PDR, len_amp, len_amp);
        auto PDR = ro.search(0, u_max, u_min);
        ro.validate_PDR(PDR);
        risk_t score = ro.score(PDR, ALPHA);
        std::cout << "loss: " << score << std::endl;
        std::ofstream fout(output_file);
        fout << display_PDR(PDR) << std::endl;
        fout.close();
    }
    else if (mode == "risk-d2") {
        std::cout << "DP-based relaxed convex optimizer (s2)" << std::endl;
        RiskOptimizer ro(rate, len_PDR, len_amp_min, len_amp);
        auto PDR = ro.search(0, u_max, u_min);
        ro.validate_PDR(PDR);
        risk_t score = ro.score(PDR, ALPHA);
        std::cout << "loss: " << score << std::endl;
        std::ofstream fout(output_file);
        fout << display_PDR(PDR) << std::endl;
        fout.close();
    }
    else if (mode == "risk-h") {
        std::cout << "Olivar greedy random optimizer" << std::endl;
        RiskOptimizer ro(rate, len_PDR, len_amp_min, len_amp);
        std::cout << "iterations: " << ITER_LIMIT * rate.size() / len_amp << std::endl;
        auto PDR = ro.random_search(ITER_LIMIT * rate.size() / len_amp);
        ro.validate_PDR(PDR);
        risk_t score = ro.score(PDR, ALPHA);
        std::cout << "loss: " << score << std::endl;
        std::ofstream fout(output_file);
        fout << display_PDR(PDR) << std::endl;
        fout.close();
    }
    else if (mode == "dimer-h") {
        KPartiteGraph *g = new KPartiteGraph(input_file);
        auto solution = g->solve_fast(iter);
        std::cout << "loss: " << g->cost(solution) << std::endl;
        // delete g;
        std::ofstream fout(output_file);
        for (index_t i : solution) 
            fout << i << " "; 
        fout << std::endl;
        fout.close();
    }
    else {

    }
    t2 = std::chrono::high_resolution_clock::now();
    long long time = std::chrono::duration_cast<std::chrono::seconds>(t2 - t1).count();
    std::cout << "time: " << time << " secs" << std::endl;
    
    return 0;
}
