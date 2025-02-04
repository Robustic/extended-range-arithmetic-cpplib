#include <iostream>
#include <vector>
#include <random>
#include <chrono>
#include "Float64Exp2Int64.h"

typedef std::mt19937 MyRNG;
uint32_t seed_val = 1337;

void GetRandomNumbers(std::vector<double>& vec);

int main() {
    union {
        unsigned long long num;
        double fp;
    } pun;
    
    // floatingExp2Integer::DoubleExp2Int test0 {0};
    // std::cout << std::dec << test0.mantissa() << " " << test0.exponent() << std::endl;
    // std::cout << std::dec << test0.asDouble() << std::endl;

    // floatingExp2Integer::DoubleExp2Int test1 {1};
    // std::cout << std::dec << test1.mantissa() << " " << test1.exponent() << std::endl;
    // std::cout << std::dec << test1.asDouble() << std::endl;
    
    // floatingExp2Integer::DoubleExp2Int test2 {2};
    // std::cout << std::dec << test2.mantissa() << " " << test2.exponent() << std::endl;
    // std::cout << std::dec << test2.asDouble() << std::endl;

    floatingExp2Integer::Float64Exp2Int64 sum {1};
    double s1 = 0.0;
    double s2 = 0.0;
    double simpleSum = 0.0;
    double test = 0.0;

    std::vector<double> vec(1000000);
    GetRandomNumbers(vec);

    int bad = 0;
    for (int k = 0; k < 10; k++) {
        for (int i = 0; i < vec.size(); i++)
            bad += vec[i];
    }  

    test += bad;

    auto start1 = std::chrono::high_resolution_clock::now();
    for (int i = 0; i < vec.size(); i++) {
        sum *= vec[i];
    }
    auto stop1 = std::chrono::high_resolution_clock::now();

    test += sum.mantissa();

    auto start2 = std::chrono::high_resolution_clock::now();
    for (int i = 0; i < vec.size(); i+=2) {
        s1 += std::log2(vec[i]);
        s2 += std::log2(vec[i+1]);
    }
    auto stop2 = std::chrono::high_resolution_clock::now();

    test += s1 + s2;

    auto start3 = std::chrono::high_resolution_clock::now();
    for (int i = 0; i < vec.size(); i++) {
        simpleSum *= vec[i];//*vec[i];
    }
    auto stop3 = std::chrono::high_resolution_clock::now();

    std::cout << sum.mantissa() << "  " << sum.exponent() << std::endl;
    std::cout << s1+s2 << std::endl;


    auto duration1 = std::chrono::duration_cast<std::chrono::microseconds>(stop1 - start1);
    auto duration2 = std::chrono::duration_cast<std::chrono::microseconds>(stop2 - start2);
    auto duration3 = std::chrono::duration_cast<std::chrono::microseconds>(stop3 - start3);


    std::cout << bad << std::endl;
    std::cout << duration1.count() << std::endl;
    std::cout << duration2.count() << std::endl;
    std::cout << duration3.count() << std::endl;
    std::cout << test << std::endl;
}


void GetRandomNumbers(std::vector<double>& vec) {
    MyRNG rng;
    rng.seed(seed_val);
    std::uniform_real_distribution<double> unif(0, 10);
    for (int i = 0; i < vec.size(); i++) {
        double a_random_double = unif(rng);
        vec[i] = a_random_double;
    }
}



