#include <iostream>
#include <vector>
#include <random>
#include <chrono>
#include <iomanip>
#include <iomanip>
#include "Timer.h"

const uint32_t seed_val = 1337;

void InitializeRandomNumbers(std::vector<uint32_t>& vec, int interval) {
    std::mt19937 rng;
    rng.seed(seed_val);
    std::uniform_int_distribution<uint32_t> unif(0, std::exp2(interval));
    for (unsigned int i = 0; i < vec.size(); i++) {
        uint32_t a_random_float = unif(rng);
        vec[i] = a_random_float;
    }
}

void InitializeRandomNumbers(std::vector<uint64_t>& vec, int interval) {
    std::mt19937_64 rng;
    rng.seed(seed_val);
    std::uniform_int_distribution<uint64_t> unif(0, std::exp2(interval));
    for (unsigned int i = 0; i < vec.size(); i++) {
        uint64_t a_random_float = unif(rng);
        vec[i] = a_random_float;
    }
}

int main() {
    const int n = 0x1p23;

    std::vector<uint32_t> int32_8(n);
    std::vector<uint32_t> int32_10(n);
    std::vector<uint32_t> int32_12(n);

    std::vector<uint64_t> int64_8(n);
    std::vector<uint64_t> int64_10(n);
    std::vector<uint64_t> int64_12(n);

    InitializeRandomNumbers(int32_8, 5);
    InitializeRandomNumbers(int32_10, 7);
    InitializeRandomNumbers(int32_12, 9);

    InitializeRandomNumbers(int64_8, 5);
    InitializeRandomNumbers(int64_10, 7);
    InitializeRandomNumbers(int64_12, 9);

    // 32

    floatingExp2Integer::Timer timer_int32_8;
    uint32_t sum_int32_8 = 0.0;
    uint32_t sum_int32_8_2 = 0.0;
    for (unsigned int i = 0; i < int32_8.size(); i++) {
        sum_int32_8 += int32_8[i];
        sum_int32_8_2 += sum_int32_8;
    }
    timer_int32_8.stop();

    floatingExp2Integer::Timer timer_int32_10;
    uint32_t sum_int32_10 = 0.0;
    uint32_t sum_int32_10_2 = 0.0;
    for (unsigned int i = 0; i < int32_10.size(); i++) {
        sum_int32_10 += int32_10[i];
        sum_int32_10_2 += sum_int32_10;
    }
    timer_int32_10.stop();

    floatingExp2Integer::Timer timer_int32_12;
    uint32_t sum_int32_12 = 0.0;
    uint32_t sum_int32_12_2 = 0.0;
    for (unsigned int i = 0; i < int32_12.size(); i++) {
        sum_int32_12 += int32_12[i];
        sum_int32_12_2 += sum_int32_12;
    }
    timer_int32_12.stop();

    // 64

    floatingExp2Integer::Timer timer_int64_8;
    uint64_t sum_int64_8 = 0.0;
    uint64_t sum_int64_8_2 = 0.0;
    for (unsigned int i = 0; i < int64_8.size(); i++) {
        sum_int64_8 += int64_8[i];
        sum_int64_8_2 += sum_int64_8;
    }
    timer_int64_8.stop();

    floatingExp2Integer::Timer timer_int64_10;
    uint64_t sum_int64_10 = 0.0;
    uint64_t sum_int64_10_2 = 0.0;
    for (unsigned int i = 0; i < int64_10.size(); i++) {
        sum_int64_10 += int64_10[i];
        sum_int64_10_2 += sum_int64_10;
    }
    timer_int64_10.stop();

    floatingExp2Integer::Timer timer_int64_12;
    uint64_t sum_int64_12 = 0.0;
    uint64_t sum_int64_12_2 = 0.0;
    for (unsigned int i = 0; i < int64_12.size(); i++) {
        sum_int64_12 += int64_12[i];
        sum_int64_12_2 += sum_int64_12;
    }
    timer_int64_12.stop();

    std::cout << "n: " << n << std::endl;

    std::cout << "timer_int32_8:  " << timer_int32_8.time() << std::endl;
    std::cout << "timer_int32_10: " << timer_int32_10.time() << std::endl;
    std::cout << "timer_int32_12: " << timer_int32_12.time() << std::endl;

    std::cout << "timer_int64_8:  " << timer_int64_8.time() << std::endl;
    std::cout << "timer_int64_10: " << timer_int64_10.time() << std::endl;
    std::cout << "timer_int64_12: " << timer_int64_12.time() << std::endl;

    std::cout << std::endl;

    std::cout << "timer_int32_8:  " << sum_int32_8 << std::endl;
    std::cout << "timer_int32_10: " << sum_int32_10 << std::endl;
    std::cout << "timer_int32_12: " << sum_int32_12 << std::endl;

    std::cout << "timer_int64_8:  " << sum_int64_8 << std::endl;
    std::cout << "timer_int64_10: " << sum_int64_10 << std::endl;
    std::cout << "timer_int64_12: " << sum_int64_12 << std::endl;

    std::cout << std::endl;

    std::cout << "timer_int32_8_2:  " << sum_int32_8_2 << std::endl;
    std::cout << "timer_int32_10_2: " << sum_int32_10_2 << std::endl;
    std::cout << "timer_int32_12_2: " << sum_int32_12_2 << std::endl;

    std::cout << "timer_int64_8_2:  " << sum_int64_8_2 << std::endl;
    std::cout << "timer_int64_10_2: " << sum_int64_10_2 << std::endl;
    std::cout << "timer_int64_12_2: " << sum_int64_12_2 << std::endl;
}