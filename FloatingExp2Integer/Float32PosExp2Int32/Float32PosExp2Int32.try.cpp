#include <iostream>
#include <vector>
#include <random>
#include <chrono>
#include <iomanip>
#include "Flt.h"
#include "Float32PosExp2Int32.h"

const uint32_t seed_val = 1337;

void InitializeRandomNumbers(std::vector<float>& vec);
void FloatToFltValues(std::vector<float>& floatValues, std::vector<floatingExp2Integer::Flt>& fltValues);
void FltToLog2Values(std::vector<floatingExp2Integer::Flt>& fltValues, std::vector<float>& log2Values);
void Log2ToFloat32PosExp2Int32Values(std::vector<float>& log2Values, 
    std::vector<floatingExp2Integer::Float32PosExp2Int32>& float32PosExp2Int32Values);
float LogSumExp2Trick(std::vector<float>& log2Values);

int main() {
    const int n = 10000;

    std::vector<float> floatValues(n);
    std::vector<floatingExp2Integer::Flt> fltValues(n);
    std::vector<float> log2Values(n);
    std::vector<floatingExp2Integer::Float32PosExp2Int32> float32PosExp2Int32Values(n);

    InitializeRandomNumbers(floatValues);
    FloatToFltValues(floatValues, fltValues);
    FltToLog2Values(fltValues, log2Values);
    Log2ToFloat32PosExp2Int32Values(log2Values, float32PosExp2Int32Values);


    auto start0 = std::chrono::high_resolution_clock::now();
    float floatSum = 0.0;
    for (unsigned int i = 0; i < floatValues.size(); i++) {
        floatSum += floatValues[i];
    }
    auto stop0 = std::chrono::high_resolution_clock::now();

    auto start1 = std::chrono::high_resolution_clock::now();
    floatingExp2Integer::Flt fltSum = 0.0;
    for (unsigned int i = 0; i < fltValues.size(); i++) {
        fltSum += fltValues[i];
    }
    auto stop1 = std::chrono::high_resolution_clock::now();

    auto start2 = std::chrono::high_resolution_clock::now();
    float log2sum = 0.0;
    log2sum = LogSumExp2Trick(log2Values);
    auto stop2 = std::chrono::high_resolution_clock::now();

    auto start3 = std::chrono::high_resolution_clock::now();
    floatingExp2Integer::Float32PosExp2Int32 float32PosExp2Int32Sum = float32PosExp2Int32Values[0];
    for (unsigned int i = 1; i < float32PosExp2Int32Values.size(); i++) {
        float32PosExp2Int32Sum += float32PosExp2Int32Values[i];
    }
    auto stop3 = std::chrono::high_resolution_clock::now();


    std::cout << std::setprecision(std::numeric_limits<float>::max_digits10);
    std::cout << "Float sum:                " << floatSum << std::endl;
    std::cout << "Flt sum:                  " << fltSum.asFloat() << std::endl;
    std::cout << "Log2 sum:                 " << std::exp2(log2sum) << std::endl;
    std::cout << "Float32PosExp2Int32 sum:  " << float32PosExp2Int32Sum.asFloat() << std::endl;

    auto duration0 = std::chrono::duration_cast<std::chrono::microseconds>(stop0 - start0);
    auto duration1 = std::chrono::duration_cast<std::chrono::microseconds>(stop1 - start1);
    auto duration2 = std::chrono::duration_cast<std::chrono::microseconds>(stop2 - start2);
    auto duration3 = std::chrono::duration_cast<std::chrono::microseconds>(stop3 - start3);

    std::cout << std::endl;
    std::cout << "Float time:               " << duration0.count() << std::endl;
    std::cout << "Flt time:                 " << duration1.count() << std::endl;
    std::cout << "Log2 time:                " << duration2.count() << std::endl;
    std::cout << "Float32PosExp2Int32 time: " << duration3.count() << std::endl;
    std::cout << std::endl;

    std::cout << "Float32PosExp2Int32 time / Log2 time " << (float)duration3.count() / (float)duration2.count() << std::endl;
    std::cout << "Log2 time / Float32PosExp2Int32 time " << (float)duration2.count() / (float)duration3.count() << std::endl;
}

void InitializeRandomNumbers(std::vector<float>& vec) {
    std::mt19937 rng;
    rng.seed(seed_val);
    std::uniform_real_distribution<float> unif(1e-30, 1e-1);
    for (unsigned int i = 0; i < vec.size(); i++) {
        float a_random_float = unif(rng);
        vec[i] = a_random_float;
    }
}

void FloatToFltValues(std::vector<float>& floatValues, std::vector<floatingExp2Integer::Flt>& fltValues) {
    for (unsigned int i = 0; i < fltValues.size(); i++) {
        fltValues[i] += floatValues[i];
    }
}

void FltToLog2Values(std::vector<floatingExp2Integer::Flt>& fltValues, std::vector<float>& log2Values) {
    for (unsigned int i = 0; i < log2Values.size(); i++) {
        log2Values[i] = std::log2(fltValues[i].asFloat());
    }
}

void Log2ToFloat32PosExp2Int32Values(std::vector<float>& log2Values, 
    std::vector<floatingExp2Integer::Float32PosExp2Int32>& float32PosExp2Int32Values) {
    for (unsigned int i = 0; i < float32PosExp2Int32Values.size(); i++) {
        float32PosExp2Int32Values[i].log2ToFloat32Exp2Int32(log2Values[i]);
    }
}

float LogSumExp2Trick(std::vector<float>& log2Values) {
    float max = *std::max_element(log2Values.begin(), log2Values.end());
    float sum = 0;

    for (float val : log2Values) {
        sum += std::exp2(val - max);
    }

    return max + std::log2(sum);
}


