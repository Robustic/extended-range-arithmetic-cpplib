#include <iostream>
#include <vector>
#include <random>
#include <chrono>
#include <iomanip>
#include "Timer.h"
#include "Flt.h"
#include "Flt2.h"
#include "./Int32PosExp2Int32/Int32PosExp2Int32.h"
#include "./Float32PosExp2Int32/Float32PosExp2Int32.h"
#include "./Float32Exp2Int32/Float32Exp2Int32.h"

const uint32_t seed_val = 1337;

void InitializeRandomNumbers(std::vector<float>& vec);
void FloatToFltValues(std::vector<float>& floatValues, std::vector<floatingExp2Integer::Flt>& fltValues);
void FloatToFlt2Values(std::vector<float>& floatValues, std::vector<floatingExp2Integer::Flt2>& flt2Values);
void FloatToLog2Values(std::vector<float>& fltValues, std::vector<float>& log2Values);
void FloatToInt32PosExp2Int32Values(std::vector<float>& floatValues, 
    std::vector<floatingExp2Integer::Int32PosExp2Int32>& int32PosExp2Int32Values);
void FloatToFloat32PosExp2Int32Values(std::vector<float>& floatValues, 
    std::vector<floatingExp2Integer::Float32PosExp2Int32>& float32PosExp2Int32Values);
void FloatToFloat32Exp2Int32Values(std::vector<float>& floatValues, 
    std::vector<floatingExp2Integer::Float32Exp2Int32>& float32Exp2Int32Values);
float LogSumExp2Trick(std::vector<float>& log2Values);
float Log2Multiply(std::vector<float>& log2Values);

void loop(int n, double results[], std::int64_t time[]);

int main() {
    const int number_of_types = 11;
    const int n[] = { 1000, 3000, 10000, 30000, 100000, 300000, 1000000, 3000000, 10000000, 30000000, 100000000 };
    const unsigned int number_of_cases = sizeof(n) / sizeof(n[0]);
    double results[number_of_cases][number_of_types];
    std::int64_t time[number_of_cases][number_of_types];

    for (unsigned int i = 0; i < number_of_cases; i++) {
        loop(n[i], results[i], time[i]);
    }

    std::cout << std::endl;
    for (unsigned int i = 0; i < number_of_cases; i++) {
        std::cout << n[i] << " ";
        for (unsigned int k = 0; k < number_of_types; k++) {
            std::cout << results[i][k] << " ";
        }
        std::cout << std::endl;
    }

    std::cout << std::endl;
    for (unsigned int i = 0; i < number_of_cases; i++) {
        std::cout << n[i] << " ";
        for (unsigned int k = 0; k < number_of_types; k++) {
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
    std::vector<floatingExp2Integer::Float32Exp2Int32> float32Exp2Int32Values(n);

    InitializeRandomNumbers(floatValues);
    FloatToFltValues(floatValues, fltValues);
    FloatToFlt2Values(floatValues, flt2Values);
    FloatToLog2Values(floatValues, log2Values);
    FloatToInt32PosExp2Int32Values(floatValues, int32PosExp2Int32Values);
    FloatToFloat32PosExp2Int32Values(floatValues, float32PosExp2Int32Values);
    FloatToFloat32Exp2Int32Values(floatValues, float32Exp2Int32Values);

    floatingExp2Integer::Timer timer_float;
    float floatSum = 0.0;
    for (unsigned int i = 0; i < floatValues.size(); i++) {
        floatSum += floatValues[i];
    }
    timer_float.stop();

    floatingExp2Integer::Timer timer_Flt;
    floatingExp2Integer::Flt fltSum = 0.0;
    for (unsigned int i = 0; i < fltValues.size(); i++) {
        fltSum += fltValues[i];
    }
    timer_Flt.stop();

    floatingExp2Integer::Timer timer_Flt2;
    floatingExp2Integer::Flt2 flt2Sum = 0.0;
    for (unsigned int i = 0; i < flt2Values.size(); i++) {
        flt2Sum += flt2Values[i];
    }
    timer_Flt2.stop();

    floatingExp2Integer::Timer timer_log2;
    float log2sum = 0.0;
    log2sum = LogSumExp2Trick(log2Values);
    timer_log2.stop();

    floatingExp2Integer::Timer timer_Int32PosExp2Int32;
    floatingExp2Integer::Int32PosExp2Int32 int32PosExp2Int32Sum = int32PosExp2Int32Values[0];
    for (unsigned int i = 1; i < int32PosExp2Int32Values.size(); i++) {
        int32PosExp2Int32Sum += int32PosExp2Int32Values[i];
    }
    timer_Int32PosExp2Int32.stop();

    floatingExp2Integer::Timer timer_Float32PosExp2Int32;
    floatingExp2Integer::Float32PosExp2Int32 float32PosExp2Int32Sum = float32PosExp2Int32Values[0];
    for (unsigned int i = 1; i < float32PosExp2Int32Values.size(); i++) {
        float32PosExp2Int32Sum += float32PosExp2Int32Values[i];
    }
    timer_Float32PosExp2Int32.stop();

    floatingExp2Integer::Timer timer_Float32Exp2Int32;
    floatingExp2Integer::Float32Exp2Int32 float32Exp2Int32Sum = float32Exp2Int32Values[0];
    for (unsigned int i = 1; i < float32Exp2Int32Values.size(); i++) {
        float32Exp2Int32Sum += float32Exp2Int32Values[i];
    }
    timer_Float32Exp2Int32.stop();

    floatingExp2Integer::Timer timer_multiply_log2;
    float log2sum_2 = 0.0;
    log2sum_2 = Log2Multiply(log2Values);
    timer_multiply_log2.stop();

    floatingExp2Integer::Timer timer_multiply_Int32PosExp2Int32;
    floatingExp2Integer::Int32PosExp2Int32 int32PosExp2Int32Sum_2 = int32PosExp2Int32Values[0];
    for (unsigned int i = 1; i < int32PosExp2Int32Values.size(); i++) {
        int32PosExp2Int32Sum_2 *= int32PosExp2Int32Values[i];
    }
    timer_multiply_Int32PosExp2Int32.stop();

    floatingExp2Integer::Timer timer_multiply_Float32PosExp2Int32;
    floatingExp2Integer::Float32PosExp2Int32 float32PosExp2Int32Sum_2 = float32PosExp2Int32Values[0];
    for (unsigned int i = 1; i < float32PosExp2Int32Values.size(); i++) {
        float32PosExp2Int32Sum_2 *= float32PosExp2Int32Values[i];
    }
    timer_multiply_Float32PosExp2Int32.stop();

    floatingExp2Integer::Timer timer_multiply_Float32Exp2Int32;
    floatingExp2Integer::Float32Exp2Int32 float32Exp2Int32Sum_2 = float32Exp2Int32Values[0];
    for (unsigned int i = 1; i < float32Exp2Int32Values.size(); i++) {
        float32Exp2Int32Sum_2 *= float32Exp2Int32Values[i];
    }
    timer_multiply_Float32Exp2Int32.stop();


    std::cout << std::setprecision(std::numeric_limits<float>::max_digits10);
    std::cout << "Float sum:                    " << floatSum << std::endl;
    std::cout << "Flt sum:                      " << fltSum.asFloat() << std::endl;
    std::cout << "Flt2 sum:                     " << flt2Sum.asFloat() << std::endl;
    std::cout << "Log2 sum:                     " << std::exp2(log2sum) << std::endl;
    std::cout << "Int32PosExp2Int32 sum:        " << int32PosExp2Int32Sum.asFloat() << std::endl;
    std::cout << "Float32PosExp2Int32 sum:      " << float32PosExp2Int32Sum.asFloat() << std::endl;
    std::cout << "Float32Exp2Int32 sum:         " << float32Exp2Int32Sum.asFloat() << std::endl;
    std::cout << "Log2 multiply:                " << log2sum_2 << std::endl;
    std::cout << "Int32PosExp2Int32 multiply:   " << int32PosExp2Int32Sum_2.int32PosExp2Int32ToLog2() << std::endl;
    std::cout << "Float32PosExp2Int32 multiply: " << float32PosExp2Int32Sum_2.float32PosExp2Int32ToLog2() << std::endl;
    std::cout << "Float32Exp2Int32 multiply:    " << float32Exp2Int32Sum_2.float32Exp2Int32ToLog2() << std::endl;

    results[0] = floatSum;
    results[1] = fltSum.asFloat();
    results[2] = flt2Sum.asFloat();
    results[3] = std::exp2(log2sum);
    results[4] = int32PosExp2Int32Sum.asFloat();
    results[5] = float32PosExp2Int32Sum.asFloat();
    results[6] = float32Exp2Int32Sum.asFloat();
    results[7] = log2sum_2;
    results[8] = int32PosExp2Int32Sum_2.int32PosExp2Int32ToLog2();
    results[9] = float32PosExp2Int32Sum_2.float32PosExp2Int32ToLog2();
    results[10] = float32Exp2Int32Sum_2.float32Exp2Int32ToLog2();

    time[0] = timer_float.time();
    time[1] = timer_Flt.time();
    time[2] = timer_Flt2.time();
    time[3] = timer_log2.time();
    time[4] = timer_Int32PosExp2Int32.time();
    time[5] = timer_Float32PosExp2Int32.time();
    time[6] = timer_Float32Exp2Int32.time();
    time[7] = timer_multiply_log2.time();
    time[8] = timer_multiply_Int32PosExp2Int32.time();
    time[9] = timer_multiply_Float32PosExp2Int32.time();
    time[10] = timer_multiply_Float32Exp2Int32.time();

    std::cout << std::endl;
    std::cout << "Float time:                   " << timer_float.time() << std::endl;
    std::cout << "Flt time:                     " << timer_Flt.time() << std::endl;
    std::cout << "Flt2 time:                    " << timer_Flt2.time() << std::endl;
    std::cout << "Log2 time:                    " << timer_log2.time() << std::endl;
    std::cout << "Int32PosExp2Int32 time:       " << timer_Int32PosExp2Int32.time() << std::endl;
    std::cout << "Float32PosExp2Int32 time:     " << timer_Float32PosExp2Int32.time() << std::endl;
    std::cout << "Float32Exp2Int32 time:        " << timer_Float32Exp2Int32.time() << std::endl;
    std::cout << "Log2 multiply:                " << timer_multiply_log2.time() << std::endl;
    std::cout << "Int32PosExp2Int32 multiply:   " << timer_multiply_Int32PosExp2Int32.time() << std::endl;
    std::cout << "Float32PosExp2Int32 multiply: " << timer_multiply_Float32PosExp2Int32.time() << std::endl;
    std::cout << "Float32Exp2Int32 multiply:    " << timer_multiply_Float32Exp2Int32.time() << std::endl;
    std::cout << std::endl;
}

void InitializeRandomNumbers(std::vector<float>& vec) {
    std::mt19937 rng;
    rng.seed(seed_val);
    std::uniform_real_distribution<float> unif(std::log2(1e-11), std::log2(1e-10));
    for (unsigned int i = 0; i < vec.size(); i++) {
        float a_random_float = unif(rng);
        vec[i] = std::exp2f(a_random_float);
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

void FloatToLog2Values(std::vector<float>& floatValues, std::vector<float>& log2Values) {
    for (unsigned int i = 0; i < log2Values.size(); i++) {
        log2Values[i] = std::log2(floatValues[i]);
    }
}

void FloatToInt32PosExp2Int32Values(std::vector<float>& floatValues, 
    std::vector<floatingExp2Integer::Int32PosExp2Int32>& int32PosExp2Int32Values) {
    for (unsigned int i = 0; i < int32PosExp2Int32Values.size(); i++) {
        int32PosExp2Int32Values[i].floatToInt32PosExp2Int32(floatValues[i]);
    }
}

void FloatToFloat32PosExp2Int32Values(std::vector<float>& floatValues, 
    std::vector<floatingExp2Integer::Float32PosExp2Int32>& float32PosExp2Int32Values) {
    for (unsigned int i = 0; i < float32PosExp2Int32Values.size(); i++) {
        float32PosExp2Int32Values[i].floatToFloat32PosExp2Int32(floatValues[i]);
    }
}

void FloatToFloat32Exp2Int32Values(std::vector<float>& floatValues, 
    std::vector<floatingExp2Integer::Float32Exp2Int32>& float32Exp2Int32Values) {
    for (unsigned int i = 0; i < float32Exp2Int32Values.size(); i++) {
        float32Exp2Int32Values[i].floatToFloat32Exp2Int32(floatValues[i]);
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


