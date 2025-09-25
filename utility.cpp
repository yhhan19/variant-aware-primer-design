#include "utility.hpp"

std::vector<risk_t> random_risk(std::size_t size) {
    std::vector<risk_t> input;
    for (std::size_t i = 0; i < size; i ++) 
        input.push_back(5 * (float)(rand()) / (float)(RAND_MAX));
    return input;
}

std::size_t random_between(std::size_t min, std::size_t max) {
    std::size_t r = rand() * RAND_MAX + rand();
    return min + r % (max - min);
}

key_t to_key(index_t f, index_t r, index_t f_, index_t r_) {
    return (key_t) f + ((key_t) r << 16) + ((key_t) f_ << 32) + ((key_t) r_ << 48);
}

void to_index(key_t k, index_t &f, index_t &r, index_t &f_, index_t &r_) {
    f = k & 0xffff;
    r = (k >> 16) & 0xffff;
    f_ = (k >> 32) & 0xffff;
    r_ = (k >> 48) & 0xffff;
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
