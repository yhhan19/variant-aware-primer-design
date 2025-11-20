#include "utility.hpp"

std::vector<risk_t> random_risk(std::size_t size) {
    std::vector<risk_t> input;
    for (std::size_t i = 0; i < size; i ++) 
        input.push_back(4 * (float)(rand()) / (float)(RAND_MAX));
    return input;
}

std::size_t random_between(std::size_t min, std::size_t max) {
    std::size_t r = rand() * RAND_MAX + rand();
    return min + r % (max - min);
}

key_t to_key(index_t f, index_t r, index_t f_, index_t r_) {
    return (key_t) f + ((key_t) r << 16) + ((key_t) f_ << 32) + ((key_t) r_ << 48);
}

key_t to_key_2(index_t r, index_t r_) {
    return (key_t) r + ((key_t) (r - r_) << 16);
}

void to_index(key_t k, index_t &f, index_t &r, index_t &f_, index_t &r_) {
    f = k & 0xffff;
    r = (k >> 16) & 0xffff;
    f_ = (k >> 32) & 0xffff;
    r_ = (k >> 48) & 0xffff;
}

void to_index_2(key_t k, index_t &r, index_t &r_) {
    r = k & 0xffff;
    r_ = r - (k >> 16) & 0xffff;
}

std::size_t read_rate(std::string filename, std::vector<risk_t> &input) {
    std::ifstream fin(filename);
    if (! fin.is_open()) {
        std::cout << "file not open" << std::endl;
        return -1;
    }
    input.clear();
    std::string line;
    while (std::getline(fin, line)) {
        std::stringstream ss(line);
        risk_t rate;
        ss >> rate >> rate;
        input.push_back(rate);
    }
    fin.close();
    return 0;
}

std::size_t read_rate_2(std::string filename, std::vector<risk_t> &input, std::string &ref) {
    std::ifstream fin(filename);
    if (! fin.is_open()) {
        std::cout << "file not open" << std::endl;
        return -1;
    }
    input.clear();
    std::string line;
    std::vector<std::pair<index_t, risk_t>> temp;
    index_t max = 0;
    std::getline(fin, line);
    while (std::getline(fin, line)) {
        std::stringstream ss(line);
        index_t i;
        risk_t rate;
        char comma, nc;
        //    1    ,        A        ,        0.0     ,        0.0     ,        0.0     ,        0.0     ,        0.0
        ss >> i >> comma >> nc >> comma >> rate >> comma >> rate >> comma >> rate >> comma >> rate >> comma >> rate;
        temp.push_back(std::make_pair(i, rate));
        ref += nc;
        max = std::max(max, i);
    }
    for (index_t i = 0; i < max; i ++) 
        input.push_back(0);
    for (auto elem : temp) {
        input[elem.first] = elem.second;
    }
    fin.close();
    return 0;
}

std::size_t read_rate_ref(std::string filename, std::vector<risk_t> &input, std::string &ref) {
    std::ifstream fin(filename);
    if (! fin.is_open()) {
        std::cout << "file not open" << std::endl;
        return 1;
    }
    input.clear();
    std::string line;
    std::vector<std::pair<index_t, risk_t>> temp;
    index_t max = 0;
    std::getline(fin, line);
    while (std::getline(fin, line)) {
        std::stringstream ss(line);
        index_t i;
        risk_t rate;
        char comma, nc;
        ss >> i >> comma >> nc >> comma >> rate;
        temp.push_back(std::make_pair(i, rate));
        ref += nc;
        max = std::max(max, i);
    }
    for (index_t i = 0; i < max; i ++) 
        input.push_back(0);
    for (auto elem : temp) {
        input[elem.first] = elem.second;
    }
    fin.close();
    return 0;
}

std::string display_PDR(std::vector<index_t> PDR) {
    std::string PDRs = "";
    for (index_t i = 0; i < PDR.size(); i ++) 
        PDRs += std::to_string(PDR[i]) + " ";
    return PDRs;
}