#include <iostream>
#include <vector>
#include <random>
#include <chrono>
#include <iomanip>
#include "Timer.h"
#include "Dbl.h"
#include "Dbl2.h"
#include "./Int64PosExp2Int64/Int64PosExp2Int64.h"
#include "./Float64PosExp2Int64/Float64PosExp2Int64.h"
#include "./Float64Exp2Int64/Float64Exp2Int64.h"

const uint32_t seed_val = 1337;

void InitializeRandomNumbers(std::vector<double>& vec);
void DoubleToDblValues(std::vector<double>& doubleValues, std::vector<floatingExp2Integer::Dbl>& dblValues);
void DoubleToDbl2Values(std::vector<double>& doubleValues, std::vector<floatingExp2Integer::Dbl2>& dbl2Values);
void DoubleToLog2Values(std::vector<double>& dblValues, std::vector<double>& log2Values);
void DoubleToInt64PosExp2Int64Values(std::vector<double>& doubleValues, 
    std::vector<floatingExp2Integer::Int64PosExp2Int64>& int64PosExp2Int64Values);
void DoubleToFloat64PosExp2Int64Values(std::vector<double>& doubleValues, 
    std::vector<floatingExp2Integer::Float64PosExp2Int64>& float64PosExp2Int64Values);
void DoubleToFloat64Exp2Int64Values(std::vector<double>& doubleValues, 
    std::vector<floatingExp2Integer::Float64Exp2Int64>& float64Exp2Int64Values);
double LogSumExp2Trick(std::vector<double>& log2Values);
double Log2Multiply(std::vector<double>& log2Values);

void loop(int n, double results[], std::int64_t time[]);

int main() {
    const int number_of_types = 12;
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

double sum_double(std::vector<double>& doubleValues) {
    double doubleSum = 0.0;
    for (unsigned int i = 0; i < doubleValues.size(); i++) {
        doubleSum += doubleValues[i];
    }
    return doubleSum;
}

double sum_Dbl(std::vector<floatingExp2Integer::Dbl>& dblValues) {
    floatingExp2Integer::Dbl dblSum = 0.0;
    for (unsigned int i = 0; i < dblValues.size(); i++) {
        dblSum += dblValues[i];
    }
    return dblSum.asDouble();
}

double sum_Dbl2(std::vector<floatingExp2Integer::Dbl2>& dbl2Values) {
    floatingExp2Integer::Dbl2 dbl2Sum = 0.0;
    for (unsigned int i = 0; i < dbl2Values.size(); i++) {
        dbl2Sum += dbl2Values[i];
    }
    return dbl2Sum.asDouble();
}

double sum_Int64PosExp2Int64(std::vector<floatingExp2Integer::Int64PosExp2Int64>& int64PosExp2Int64Values) {
    floatingExp2Integer::Int64PosExp2Int64 int64PosExp2Int64Sum = int64PosExp2Int64Values[0];
    for (unsigned int i = 1; i < int64PosExp2Int64Values.size(); i++) {
        int64PosExp2Int64Sum += int64PosExp2Int64Values[i];
    }
    return int64PosExp2Int64Sum.asDouble();
}

double sum_Float64PosExp2Int64(std::vector<floatingExp2Integer::Float64PosExp2Int64>& float64PosExp2Int64Values) {
    floatingExp2Integer::Float64PosExp2Int64 float64PosExp2Int64Sum = float64PosExp2Int64Values[0];
    for (unsigned int i = 1; i < float64PosExp2Int64Values.size(); i++) {
        float64PosExp2Int64Sum += float64PosExp2Int64Values[i];
    }
    return float64PosExp2Int64Sum.asDouble();
}

double sum_Float64Exp2Int64(std::vector<floatingExp2Integer::Float64Exp2Int64>& float64Exp2Int64Values) {
    floatingExp2Integer::Float64Exp2Int64 float64Exp2Int64Sum = float64Exp2Int64Values[0];
    for (unsigned int i = 1; i < float64Exp2Int64Values.size(); i++) {
        float64Exp2Int64Sum += float64Exp2Int64Values[i];
    }
    return float64Exp2Int64Sum.asDouble();
}

double multiply_Int64PosExp2Int64(std::vector<floatingExp2Integer::Int64PosExp2Int64>& int64PosExp2Int64Values) {
    floatingExp2Integer::Int64PosExp2Int64 int64PosExp2Int64Multiply = int64PosExp2Int64Values[0];
    for (unsigned int i = 1; i < int64PosExp2Int64Values.size(); i++) {
        int64PosExp2Int64Multiply *= int64PosExp2Int64Values[i];
    }
    return int64PosExp2Int64Multiply.int64PosExp2Int64ToLog2();
}

double multiply_Float64PosExp2Int64(std::vector<floatingExp2Integer::Float64PosExp2Int64>& float64PosExp2Int64Values) {
    floatingExp2Integer::Float64PosExp2Int64 float64PosExp2Int64Multiply = float64PosExp2Int64Values[0];
    for (unsigned int i = 1; i < float64PosExp2Int64Values.size(); i++) {
        float64PosExp2Int64Multiply *= float64PosExp2Int64Values[i];
    }
    return float64PosExp2Int64Multiply.float64PosExp2Int64ToLog2();
}

double multiply_Float64Exp2Int64(std::vector<floatingExp2Integer::Float64Exp2Int64>& float64Exp2Int64Values) {
    floatingExp2Integer::Float64Exp2Int64 float64Exp2Int64Multiply = float64Exp2Int64Values[0];
    for (unsigned int i = 1; i < float64Exp2Int64Values.size(); i++) {
        float64Exp2Int64Multiply *= float64Exp2Int64Values[i];
    }
    return float64Exp2Int64Multiply.float64Exp2Int64ToLog2();
}

double sum_vectorization_Float64PosExp2Int64(std::vector<floatingExp2Integer::Float64PosExp2Int64>& float64PosExp2Int64Values) {
    floatingExp2Integer::Float64PosExp2Int64 float64PosExp2Int64ValuesSum(float64PosExp2Int64Values);
    return float64PosExp2Int64ValuesSum.asDouble();
}

void loop(int n, double results[], std::int64_t time[]) {

    std::vector<double> doubleValues(n);
    std::vector<floatingExp2Integer::Dbl> dblValues(n);
    std::vector<floatingExp2Integer::Dbl2> dbl2Values(n);
    std::vector<double> log2Values(n);
    std::vector<floatingExp2Integer::Int64PosExp2Int64> int64PosExp2Int64Values(n);
    std::vector<floatingExp2Integer::Float64PosExp2Int64> float64PosExp2Int64Values(n);
    std::vector<floatingExp2Integer::Float64Exp2Int64> float64Exp2Int64Values(n);

    InitializeRandomNumbers(doubleValues);
    DoubleToDblValues(doubleValues, dblValues);
    DoubleToDbl2Values(doubleValues, dbl2Values);
    DoubleToLog2Values(doubleValues, log2Values);
    DoubleToInt64PosExp2Int64Values(doubleValues, int64PosExp2Int64Values);
    DoubleToFloat64PosExp2Int64Values(doubleValues, float64PosExp2Int64Values);
    DoubleToFloat64Exp2Int64Values(doubleValues, float64Exp2Int64Values);

    floatingExp2Integer::Timer timer_double;
    double doubleSum = sum_double(doubleValues);
    timer_double.stop();

    floatingExp2Integer::Timer timer_Dbl;
    double dblSum = sum_Dbl(dblValues);    
    timer_Dbl.stop();

    floatingExp2Integer::Timer timer_Dbl2;
    double dbl2Sum = sum_Dbl2(dbl2Values);
    timer_Dbl2.stop();

    floatingExp2Integer::Timer timer_log2;
    double log2sum = LogSumExp2Trick(log2Values);
    timer_log2.stop();

    floatingExp2Integer::Timer timer_Int64PosExp2Int64;
    double int64PosExp2Int64Sum = sum_Int64PosExp2Int64(int64PosExp2Int64Values);
    timer_Int64PosExp2Int64.stop();

    floatingExp2Integer::Timer timer_Float64PosExp2Int64;
    double float64PosExp2Int64Sum = sum_Float64PosExp2Int64(float64PosExp2Int64Values);
    timer_Float64PosExp2Int64.stop();

    floatingExp2Integer::Timer timer_Float64Exp2Int64;
    double float64Exp2Int64Sum = sum_Float64Exp2Int64(float64Exp2Int64Values);
    timer_Float64Exp2Int64.stop();

    floatingExp2Integer::Timer timer_multiply_log2;
    double log2sum_2 = Log2Multiply(log2Values);
    timer_multiply_log2.stop();

    floatingExp2Integer::Timer timer_multiply_Int64PosExp2Int64;
    double int64PosExp2Int64Multiply = multiply_Int64PosExp2Int64(int64PosExp2Int64Values);
    timer_multiply_Int64PosExp2Int64.stop();

    floatingExp2Integer::Timer timer_multiply_Float64PosExp2Int64;
    double float64PosExp2Int64Multiply = multiply_Float64PosExp2Int64(float64PosExp2Int64Values);
    timer_multiply_Float64PosExp2Int64.stop();

    floatingExp2Integer::Timer timer_multiply_Float64Exp2Int64;
    double float64Exp2Int64Multiply = multiply_Float64Exp2Int64(float64Exp2Int64Values);
    timer_multiply_Float64Exp2Int64.stop();

    floatingExp2Integer::Timer timer_Float64PosExp2Int64_vec;
    double float64PosExp2Int64VectorizationSum = sum_vectorization_Float64PosExp2Int64(float64PosExp2Int64Values);
    timer_Float64PosExp2Int64_vec.stop();


    std::cout << std::setprecision(std::numeric_limits<double>::max_digits10);
    std::cout << "Double sum:                             " << doubleSum << std::endl;
    std::cout << "Dbl sum:                                " << dblSum << std::endl;
    std::cout << "Dbl2 sum:                               " << dbl2Sum << std::endl;
    std::cout << "Log2 sum:                               " << std::exp2(log2sum) << std::endl;
    std::cout << "Int64PosExp2Int64 sum:                  " << int64PosExp2Int64Sum << std::endl;
    std::cout << "Float64PosExp2Int64 sum:                " << float64PosExp2Int64Sum << std::endl;
    std::cout << "Float64Exp2Int64 sum:                   " << float64Exp2Int64Sum << std::endl;
    std::cout << "Log2 multiply:                          " << log2sum_2 << std::endl;
    std::cout << "Int64PosExp2Int64 multiply:             " << int64PosExp2Int64Multiply << std::endl;
    std::cout << "Float64PosExp2Int64 multiply:           " << float64PosExp2Int64Multiply << std::endl;
    std::cout << "Float64Exp2Int64 multiply:              " << float64Exp2Int64Multiply << std::endl;
    std::cout << "Float64PosExp2Int64 sum vectorization:  " << float64PosExp2Int64VectorizationSum << std::endl;

    results[0] = doubleSum;
    results[1] = dblSum;
    results[2] = dbl2Sum;
    results[3] = std::exp2(log2sum);
    results[4] = int64PosExp2Int64Sum;
    results[5] = float64PosExp2Int64Sum;
    results[6] = float64Exp2Int64Sum;
    results[7] = log2sum_2;
    results[8] = int64PosExp2Int64Multiply;
    results[9] = float64PosExp2Int64Multiply;
    results[10] = float64Exp2Int64Multiply;
    results[11] = float64PosExp2Int64VectorizationSum;

    time[0] = timer_double.time();
    time[1] = timer_Dbl.time();
    time[2] = timer_Dbl2.time();
    time[3] = timer_log2.time();
    time[4] = timer_Int64PosExp2Int64.time();
    time[5] = timer_Float64PosExp2Int64.time();
    time[6] = timer_Float64Exp2Int64.time();
    time[7] = timer_multiply_log2.time();
    time[8] = timer_multiply_Int64PosExp2Int64.time();
    time[9] = timer_multiply_Float64PosExp2Int64.time();
    time[10] = timer_multiply_Float64Exp2Int64.time();
    time[11] = timer_Float64PosExp2Int64_vec.time();

    std::cout << std::endl;
    std::cout << "Double time:                              " << timer_double.time() << std::endl;
    std::cout << "Dbl time:                                 " << timer_Dbl.time() << std::endl;
    std::cout << "Dbl2 time:                                " << timer_Dbl2.time() << std::endl;
    std::cout << "Log2 time:                                " << timer_log2.time() << std::endl;
    std::cout << "Int64PosExp2Int64 time:                   " << timer_Int64PosExp2Int64.time() << std::endl;
    std::cout << "Float64PosExp2Int64 time:                 " << timer_Float64PosExp2Int64.time() << std::endl;
    std::cout << "Float64Exp2Int64 time:                    " << timer_Float64Exp2Int64.time() << std::endl;
    std::cout << "Log2 multiply:                            " << timer_multiply_log2.time() << std::endl;
    std::cout << "Int64PosExp2Int64 multiply:               " << timer_multiply_Int64PosExp2Int64.time() << std::endl;
    std::cout << "Float64PosExp2Int64 multiply:             " << timer_multiply_Float64PosExp2Int64.time() << std::endl;
    std::cout << "Float64Exp2Int64 multiply:                " << timer_multiply_Float64Exp2Int64.time() << std::endl;
    std::cout << "Float64PosExp2Int64 time: vectorization:  " << timer_Float64PosExp2Int64_vec.time() << std::endl;
    std::cout << std::endl;
}

void InitializeRandomNumbers(std::vector<double>& vec) {
    std::mt19937 rng;
    rng.seed(seed_val);
    std::uniform_real_distribution<float> unif(std::log2(1e-11), std::log2(1e-10));
    for (unsigned int i = 0; i < vec.size(); i++) {
        float a_random_float = unif(rng);
        vec[i] = (double)std::exp2f(a_random_float);
    }
}

void DoubleToDblValues(std::vector<double>& doubleValues, std::vector<floatingExp2Integer::Dbl>& dblValues) {
    for (unsigned int i = 0; i < dblValues.size(); i++) {
        dblValues[i] += doubleValues[i];
    }
}

void DoubleToDbl2Values(std::vector<double>& doubleValues, std::vector<floatingExp2Integer::Dbl2>& dbl2Values) {
    for (unsigned int i = 0; i < dbl2Values.size(); i++) {
        dbl2Values[i] += doubleValues[i];
    }
}

void DoubleToLog2Values(std::vector<double>& doubleValues, std::vector<double>& log2Values) {
    for (unsigned int i = 0; i < log2Values.size(); i++) {
        log2Values[i] = std::log2(doubleValues[i]);
    }
}

void DoubleToInt64PosExp2Int64Values(std::vector<double>& doubleValues, 
    std::vector<floatingExp2Integer::Int64PosExp2Int64>& int64PosExp2Int64Values) {
    for (unsigned int i = 0; i < int64PosExp2Int64Values.size(); i++) {
        int64PosExp2Int64Values[i].doubleToInt64PosExp2Int64(doubleValues[i]);
    }
}

void DoubleToFloat64PosExp2Int64Values(std::vector<double>& doubleValues, 
    std::vector<floatingExp2Integer::Float64PosExp2Int64>& float64PosExp2Int64Values) {
    for (unsigned int i = 0; i < float64PosExp2Int64Values.size(); i++) {
        float64PosExp2Int64Values[i].doubleToFloat64PosExp2Int64(doubleValues[i]);
    }
}

void DoubleToFloat64Exp2Int64Values(std::vector<double>& doubleValues, 
    std::vector<floatingExp2Integer::Float64Exp2Int64>& float64Exp2Int64Values) {
    for (unsigned int i = 0; i < float64Exp2Int64Values.size(); i++) {
        float64Exp2Int64Values[i].doubleToFloat64Exp2Int64(doubleValues[i]);
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

double Log2Multiply(std::vector<double>& log2Values) {
    double sum = 0;
    for (double val : log2Values) {
        sum += val;
    }

    return sum;
}


