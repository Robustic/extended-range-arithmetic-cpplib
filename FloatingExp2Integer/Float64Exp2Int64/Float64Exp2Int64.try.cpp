#include <iostream>
#include <vector>
#include <random>
#include <chrono>
#include "Float64Exp2Int64.h"

typedef std::mt19937 MyRNG;
uint32_t seed_val = 1337;

void InitializeRandomNumbers(std::vector<double>& vec);

int main() {
    floatingExp2Integer::Float64Exp2Int64 sum {1};
    double s1 = 0.0;
    double s2 = 0.0;
    double simpleSum = 0.0;
    double test = 0.0;

    std::vector<double> vec(1000000);
    InitializeRandomNumbers(vec);

    int bad = 0;
    for (unsigned int k = 0; k < 10; k++) {
        for (unsigned int i = 0; i < vec.size(); i++)
            bad += vec[i];
    }  

    test += bad;

    auto start1 = std::chrono::high_resolution_clock::now();
    for (unsigned int i = 0; i < vec.size(); i++) {
        sum *= vec[i];
    }
    auto stop1 = std::chrono::high_resolution_clock::now();

    test += sum.mantissa();

    auto start2 = std::chrono::high_resolution_clock::now();
    for (unsigned int i = 0; i < vec.size(); i+=2) {
        s1 += std::log2(vec[i]);
        s2 += std::log2(vec[i+1]);
    }
    auto stop2 = std::chrono::high_resolution_clock::now();

    test += s1 + s2;

    auto start3 = std::chrono::high_resolution_clock::now();
    for (unsigned int i = 0; i < vec.size(); i++) {
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


void InitializeRandomNumbers(std::vector<double>& vec) {
    MyRNG rng;
    rng.seed(seed_val);
    std::uniform_real_distribution<double> unif(0, 10);
    for (unsigned int i = 0; i < vec.size(); i++) {
        double a_random_double = unif(rng);
        vec[i] = a_random_double;
    }
}



