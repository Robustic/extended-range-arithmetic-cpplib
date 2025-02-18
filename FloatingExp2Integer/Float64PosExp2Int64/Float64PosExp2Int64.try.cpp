#include <iostream>
#include <vector>
#include <random>
#include <chrono>
#include <iomanip>
#include "Dbl.h"
#include "Float64PosExp2Int64.h"

const uint32_t seed_val = 1337;

void InitializeRandomNumbers(std::vector<double>& vec);
void DoubleToDblValues(std::vector<double>& doubleValues, std::vector<floatingExp2Integer::Dbl>& dblValues);
void DblToLog2Values(std::vector<floatingExp2Integer::Dbl>& dblValues, std::vector<double>& log2Values);
void Log2ToFloat64PosExp2Int64Values(std::vector<double>& log2Values, 
    std::vector<floatingExp2Integer::Float64PosExp2Int64>& float64PosExp2Int64Values);
double LogSumExp2Trick(std::vector<double>& log2Values);

int main() {
    const int n = 100000000;

    std::vector<double> doubleValues(n);
    std::vector<floatingExp2Integer::Dbl> dblValues(n);
    std::vector<double> log2Values(n);
    std::vector<floatingExp2Integer::Float64PosExp2Int64> float64PosExp2Int64Values(n);

    InitializeRandomNumbers(doubleValues);
    DoubleToDblValues(doubleValues, dblValues);
    DblToLog2Values(dblValues, log2Values);
    Log2ToFloat64PosExp2Int64Values(log2Values, float64PosExp2Int64Values);


    auto start0 = std::chrono::high_resolution_clock::now();
    double doubleSum = 0.0;
    for (unsigned int i = 0; i < doubleValues.size(); i++) {
        doubleSum += doubleValues[i];
    }
    auto stop0 = std::chrono::high_resolution_clock::now();

    auto start1 = std::chrono::high_resolution_clock::now();
    floatingExp2Integer::Dbl dblSum = 0.0;
    for (unsigned int i = 0; i < dblValues.size(); i++) {
        dblSum += dblValues[i];
    }
    auto stop1 = std::chrono::high_resolution_clock::now();

    auto start2 = std::chrono::high_resolution_clock::now();
    double log2sum = 0.0;
    log2sum = LogSumExp2Trick(log2Values);
    auto stop2 = std::chrono::high_resolution_clock::now();

    auto start3 = std::chrono::high_resolution_clock::now();
    floatingExp2Integer::Float64PosExp2Int64 float64PosExp2Int64Sum = float64PosExp2Int64Values[0];
    for (unsigned int i = 1; i < float64PosExp2Int64Values.size(); i++) {
        float64PosExp2Int64Sum += float64PosExp2Int64Values[i];
    }
    auto stop3 = std::chrono::high_resolution_clock::now();


    std::cout << std::setprecision(std::numeric_limits<double>::max_digits10);
    std::cout << "Double sum:               " << doubleSum << std::endl;
    std::cout << "Dbl sum:                  " << dblSum.asDouble() << std::endl;
    std::cout << "Log2 sum:                 " << std::exp2(log2sum) << std::endl;
    std::cout << "Float64PosExp2Int64 sum:  " << float64PosExp2Int64Sum.asDouble() << std::endl;

    auto duration0 = std::chrono::duration_cast<std::chrono::microseconds>(stop0 - start0);
    auto duration1 = std::chrono::duration_cast<std::chrono::microseconds>(stop1 - start1);
    auto duration2 = std::chrono::duration_cast<std::chrono::microseconds>(stop2 - start2);
    auto duration3 = std::chrono::duration_cast<std::chrono::microseconds>(stop3 - start3);

    std::cout << std::endl;
    std::cout << "Double time:              " << duration0.count() << std::endl;
    std::cout << "Dbl time:                 " << duration1.count() << std::endl;
    std::cout << "Log2 time:                " << duration2.count() << std::endl;
    std::cout << "Float64PosExp2Int64 time: " << duration3.count() << std::endl;
    std::cout << std::endl;

    std::cout << "Float64PosExp2Int64 time / Log2 time " << (double)duration3.count() / (double)duration2.count() << std::endl;
    std::cout << "Log2 time / Float64PosExp2Int64 time " << (double)duration2.count() / (double)duration3.count() << std::endl;
}

void InitializeRandomNumbers(std::vector<double>& vec) {
    std::mt19937 rng;
    rng.seed(seed_val);
    std::uniform_real_distribution<double> unif(1e-300, 1e-10);
    for (unsigned int i = 0; i < vec.size(); i++) {
        double a_random_double = unif(rng);
        vec[i] = a_random_double;
    }
}

void DoubleToDblValues(std::vector<double>& doubleValues, std::vector<floatingExp2Integer::Dbl>& dblValues) {
    for (unsigned int i = 0; i < dblValues.size(); i++) {
        dblValues[i] += doubleValues[i];
    }
}

void DblToLog2Values(std::vector<floatingExp2Integer::Dbl>& dblValues, std::vector<double>& log2Values) {
    for (unsigned int i = 0; i < log2Values.size(); i++) {
        log2Values[i] = std::log2(dblValues[i].asDouble());
    }
}

void Log2ToFloat64PosExp2Int64Values(std::vector<double>& log2Values, 
    std::vector<floatingExp2Integer::Float64PosExp2Int64>& float64PosExp2Int64Values) {
    for (unsigned int i = 0; i < float64PosExp2Int64Values.size(); i++) {
        float64PosExp2Int64Values[i].Log2ToFloat64Exp2Int64(log2Values[i]);
    }
}

double LogSumExp2Trick(std::vector<double>& log2Values) {
    double max = *std::max_element(log2Values.begin(), log2Values.end());
    double sum = 0;

    for (double val : log2Values) {
        sum += std::exp2(val - max);
    }

    return max + std::log2(sum);
}


