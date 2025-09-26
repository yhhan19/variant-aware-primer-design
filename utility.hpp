#ifndef __UTILITY__
#define __UTILITY__

#include <string>
#include <vector>
#include <fstream>
#include <iostream>
#include <sstream>
#include <unordered_map>
#include <cassert>
#include <algorithm>
#include <chrono>
#include <cmath>
#include <ctime>
#include <iomanip>

#define risk_t double
#define ALPHA 0.1
#define key_t uint64_t
#define index_t uint16_t
#define INF 1e18
#define ITER_LIMIT 500
#define KEY_LIMIT (1 << 26)

std::vector<risk_t> random_risk(std::size_t size);
std::size_t random_between(std::size_t min, std::size_t max);
key_t to_key(index_t f, index_t r, index_t f_, index_t r_);
void to_index(key_t k, index_t &f, index_t &r, index_t &f_, index_t &r_);
std::size_t read_rate(std::string filename, std::vector<risk_t> &input);
key_t to_key_2(index_t r, index_t r_);
void to_index_2(key_t k, index_t &r, index_t &r_);

#endif
