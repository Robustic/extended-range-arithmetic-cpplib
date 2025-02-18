#include <iostream>
#include <vector>
#include <random>
#include <chrono>
#include <iomanip>
#include "Flt.h"
#include "Flt2.h"
#include "Int32PosExp2Int32.h"
#include "../Float32PosExp2Int32/Float32PosExp2Int32.h"

const uint32_t seed_val = 1337;

void InitializeRandomNumbers(std::vector<float>& vec);
void FloatToFltValues(std::vector<float>& floatValues, std::vector<floatingExp2Integer::Flt>& fltValues);
void FloatToFlt2Values(std::vector<float>& floatValues, std::vector<floatingExp2Integer::Flt2>& flt2Values);
void FltToLog2Values(std::vector<floatingExp2Integer::Flt>& fltValues, std::vector<float>& log2Values);
void Log2ToInt32PosExp2Int32Values(std::vector<float>& log2Values, 
    std::vector<floatingExp2Integer::Int32PosExp2Int32>& int32PosExp2Int32Values);
void Log2ToFloat32PosExp2Int32Values(std::vector<float>& log2Values, 
    std::vector<floatingExp2Integer::Float32PosExp2Int32>& float32PosExp2Int32Values);
float LogSumExp2Trick(std::vector<float>& log2Values);
float Log2Multiply(std::vector<float>& log2Values);

void loop(int n, double results[], std::int64_t time[]);

int main() {

    const int n[] = { 1000, 3000, 10000, 30000, 100000, 300000, 1000000, 3000000, 10000000, 30000000, 100000000 };
    const unsigned int numberOfCases = sizeof(n) / sizeof(n[0]);
    double results[numberOfCases][9];
    std::int64_t time[numberOfCases][9];

    for (unsigned int i = 0; i < numberOfCases; i++) {
        loop(n[i], results[i], time[i]);
    }

    std::cout << std::endl;
    for (unsigned int i = 0; i < numberOfCases; i++) {
        std::cout << n[i] << " ";
        for (unsigned int k = 0; k < 9; k++) {
            std::cout << results[i][k] << " ";
        }
        std::cout << std::endl;
    }

    std::cout << std::endl;
    for (unsigned int i = 0; i < numberOfCases; i++) {
        std::cout << n[i] << " ";
        for (unsigned int k = 0; k < 9; k++) {
            std::cout << time[i][k] << " ";
        }
        std::cout << std::endl;
    }
}

void loop(int n, double results[], std::int64_t time[]) {

    std::vector<float> floatValues(n);
    std::vector<floatingExp2Integer::Flt> fltValues(n);
    std::vector<floatingExp2Integer::Flt2> flt2Values(n);
    std::vector<float> log2Values(n);
    std::vector<floatingExp2Integer::Int32PosExp2Int32> int32PosExp2Int32Values(n);
    std::vector<floatingExp2Integer::Float32PosExp2Int32> float32PosExp2Int32Values(n);

    InitializeRandomNumbers(floatValues);
    FloatToFltValues(floatValues, fltValues);
    FloatToFlt2Values(floatValues, flt2Values);
    FltToLog2Values(fltValues, log2Values);
    Log2ToInt32PosExp2Int32Values(log2Values, int32PosExp2Int32Values);
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

    auto start1_2 = std::chrono::high_resolution_clock::now();
    floatingExp2Integer::Flt2 flt2Sum = 0.0;
    for (unsigned int i = 0; i < flt2Values.size(); i++) {
        flt2Sum += flt2Values[i];
    }
    auto stop1_2 = std::chrono::high_resolution_clock::now();

    auto start2 = std::chrono::high_resolution_clock::now();
    float log2sum = 0.0;
    log2sum = LogSumExp2Trick(log2Values);
    auto stop2 = std::chrono::high_resolution_clock::now();

    auto start3 = std::chrono::high_resolution_clock::now();
    floatingExp2Integer::Int32PosExp2Int32 int32PosExp2Int32Sum = int32PosExp2Int32Values[0];
    for (unsigned int i = 1; i < int32PosExp2Int32Values.size(); i++) {
        int32PosExp2Int32Sum += int32PosExp2Int32Values[i];
    }
    auto stop3 = std::chrono::high_resolution_clock::now();

    auto start4 = std::chrono::high_resolution_clock::now();
    floatingExp2Integer::Float32PosExp2Int32 float32PosExp2Int32Sum = float32PosExp2Int32Values[0];
    for (unsigned int i = 1; i < float32PosExp2Int32Values.size(); i++) {
        float32PosExp2Int32Sum += float32PosExp2Int32Values[i];
    }
    auto stop4 = std::chrono::high_resolution_clock::now();

    auto start2_2 = std::chrono::high_resolution_clock::now();
    float log2sum_2 = 0.0;
    log2sum_2 = Log2Multiply(log2Values);
    auto stop2_2 = std::chrono::high_resolution_clock::now();

    auto start3_2 = std::chrono::high_resolution_clock::now();
    floatingExp2Integer::Int32PosExp2Int32 int32PosExp2Int32Sum_2 = int32PosExp2Int32Values[0];
    for (unsigned int i = 1; i < int32PosExp2Int32Values.size(); i++) {
        int32PosExp2Int32Sum_2 *= int32PosExp2Int32Values[i];
    }
    auto stop3_2 = std::chrono::high_resolution_clock::now();

    auto start4_2 = std::chrono::high_resolution_clock::now();
    floatingExp2Integer::Float32PosExp2Int32 float32PosExp2Int32Sum_2 = float32PosExp2Int32Values[0];
    for (unsigned int i = 1; i < float32PosExp2Int32Values.size(); i++) {
        float32PosExp2Int32Sum_2 *= float32PosExp2Int32Values[i];
    }
    auto stop4_2 = std::chrono::high_resolution_clock::now();


    std::cout << std::setprecision(std::numeric_limits<float>::max_digits10);
    std::cout << "Float sum:                    " << floatSum << std::endl;
    std::cout << "Flt sum:                      " << fltSum.asFloat() << std::endl;
    std::cout << "Flt2 sum:                     " << flt2Sum.asFloat() << std::endl;
    std::cout << "Log2 sum:                     " << std::exp2(log2sum) << std::endl;
    std::cout << "Int32PosExp2Int32 sum:        " << int32PosExp2Int32Sum.asFloat() << std::endl;
    std::cout << "Float32PosExp2Int32 sum:      " << float32PosExp2Int32Sum.asFloat() << std::endl;
    std::cout << "Log2 multiply:                " << log2sum_2 << std::endl;
    std::cout << "Int32PosExp2Int32 multiply:   " << int32PosExp2Int32Sum_2.int32Exp2Int32ToLog2() << std::endl;
    std::cout << "Float32PosExp2Int32 multiply: " << float32PosExp2Int32Sum_2.float32Exp2Int32ToLog2() << std::endl;
    std::cout << int32PosExp2Int32Sum_2.sicnificand() << "  " << int32PosExp2Int32Sum_2.exponent() << std::endl;

    results[0] = floatSum;
    results[1] = fltSum.asFloat();
    results[2] = flt2Sum.asFloat();
    results[3] = std::exp2(log2sum);
    results[4] = int32PosExp2Int32Sum.asFloat();
    results[5] = float32PosExp2Int32Sum.asFloat();
    results[6] = log2sum_2;
    results[7] = int32PosExp2Int32Sum_2.int32Exp2Int32ToLog2();
    results[8] = float32PosExp2Int32Sum_2.float32Exp2Int32ToLog2();

    auto duration0 = std::chrono::duration_cast<std::chrono::microseconds>(stop0 - start0);
    auto duration1 = std::chrono::duration_cast<std::chrono::microseconds>(stop1 - start1);
    auto duration1_2 = std::chrono::duration_cast<std::chrono::microseconds>(stop1_2 - start1_2);
    auto duration2 = std::chrono::duration_cast<std::chrono::microseconds>(stop2 - start2);
    auto duration3 = std::chrono::duration_cast<std::chrono::microseconds>(stop3 - start3);
    auto duration4 = std::chrono::duration_cast<std::chrono::microseconds>(stop4 - start4);
    auto duration2_2 = std::chrono::duration_cast<std::chrono::microseconds>(stop2_2 - start2_2);
    auto duration3_2 = std::chrono::duration_cast<std::chrono::microseconds>(stop3_2 - start3_2);
    auto duration4_2 = std::chrono::duration_cast<std::chrono::microseconds>(stop4_2 - start4_2);

    time[0] = duration0.count();
    time[1] = duration1.count();
    time[2] = duration1_2.count();
    time[3] = duration2.count();
    time[4] = duration3.count();
    time[5] = duration4.count();
    time[6] = duration2_2.count();
    time[7] = duration3_2.count();
    time[8] = duration4_2.count();

    std::cout << std::endl;
    std::cout << "Float time:                   " << duration0.count() << std::endl;
    std::cout << "Flt time:                     " << duration1.count() << std::endl;
    std::cout << "Flt2 time:                    " << duration1_2.count() << std::endl;
    std::cout << "Log2 time:                    " << duration2.count() << std::endl;
    std::cout << "Int32PosExp2Int32 time:       " << duration3.count() << std::endl;
    std::cout << "Float32PosExp2Int32 time:     " << duration4.count() << std::endl;
    std::cout << "Log2 multiply:                " << duration2_2.count() << std::endl;
    std::cout << "Int32PosExp2Int32 multiply:   " << duration3_2.count() << std::endl;
    std::cout << "Float32PosExp2Int32 multiply: " << duration4_2.count() << std::endl;
    std::cout << std::endl;

    std::cout << "Int32PosExp2Int32 time / Log2 time " << (float)duration3.count() / (float)duration2.count() << std::endl;
    std::cout << "Log2 time / Int32PosExp2Int32 time " << (float)duration2.count() / (float)duration3.count() << std::endl;
}

void InitializeRandomNumbers(std::vector<float>& vec) {
    std::mt19937 rng;
    rng.seed(seed_val);
    std::uniform_real_distribution<float> unif(1e-30, 1e-10);
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

void FloatToFlt2Values(std::vector<float>& floatValues, std::vector<floatingExp2Integer::Flt2>& flt2Values) {
    for (unsigned int i = 0; i < flt2Values.size(); i++) {
        flt2Values[i] += floatValues[i];
    }
}

void FltToLog2Values(std::vector<floatingExp2Integer::Flt>& fltValues, std::vector<float>& log2Values) {
    for (unsigned int i = 0; i < log2Values.size(); i++) {
        log2Values[i] = std::log2(fltValues[i].asFloat());
    }
}

void Log2ToInt32PosExp2Int32Values(std::vector<float>& log2Values, 
    std::vector<floatingExp2Integer::Int32PosExp2Int32>& int32PosExp2Int32Values) {
    for (unsigned int i = 0; i < int32PosExp2Int32Values.size(); i++) {
        int32PosExp2Int32Values[i].log2ToInt32Exp2Int32(log2Values[i]);
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

float Log2Multiply(std::vector<float>& log2Values) {
    float sum = 0;
    for (float val : log2Values) {
        sum += val;
    }

    return sum;
}


