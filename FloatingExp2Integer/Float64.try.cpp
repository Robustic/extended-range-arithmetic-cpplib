#include <iostream>
#include <vector>
#include <immintrin.h>
#define _GNU_SOURCE      /* These two lines are optional. They are only needed to use the scalar */
#include <math.h> 
#include <random>
#include <chrono>
#include <iomanip>
#include <Eigen/Dense>
#include "Timer.h"
#include "Dbl.h"
#include "Dbl2.h"
#include "Dbl3.h"
#include "./Int64PosExp2Int64/Int64PosExp2Int64.h"
#include "./Float64PosExp2Int64/Float64PosExp2Int64.h"
#include "./Float64Exp2Int64/Float64Exp2Int64.h"

const uint32_t seed_val = 1337;

void InitializeRandomNumbers(std::vector<double>& vec);
void DoubleToEigenValues(std::vector<double>& doubleValues, Eigen::VectorXd& eigen_values);
void DoubleToDblValues(std::vector<double>& doubleValues, std::vector<floatingExp2Integer::Dbl>& dblValues);
void DoubleToDbl2Values(std::vector<double>& doubleValues, std::vector<floatingExp2Integer::Dbl2>& dbl2Values);
void DoubleToDbl3Values(std::vector<double>& doubleValues, std::vector<floatingExp2Integer::Dbl3>& dbl3Values);
void DoubleToLog2Values(std::vector<double>& dblValues, std::vector<double>& log2Values);
void DoubleToInt64PosExp2Int64Values(std::vector<double>& doubleValues, 
    std::vector<floatingExp2Integer::Int64PosExp2Int64>& int64PosExp2Int64Values);
void DoubleToFloat64PosExp2Int64Values(std::vector<double>& doubleValues, 
    std::vector<floatingExp2Integer::Float64PosExp2Int64>& float64PosExp2Int64Values);
void DoubleToFloat64Exp2Int64Values(std::vector<double>& doubleValues, 
    std::vector<floatingExp2Integer::Float64Exp2Int64>& float64Exp2Int64Values);
double LogSumExp2Trick(std::vector<double>& log2Values, int64_t& time);
double Log2Multiply(std::vector<double>& log2Values, int64_t& time);

void loop(int n, double results[], std::int64_t time[]);

int main() {
    const int number_of_types = 15;
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

double vectorizedSum(const std::vector<double>& doubleValues) {
    size_t size = doubleValues.size();
    size_t i = 0;
    __m256d sumVec1 = _mm256_setzero_pd();  // Initialize sum vector to 0
    __m256d sumVec2 = _mm256_setzero_pd();  // Initialize sum vector to 0

    // Process 4 doubles at a time using AVX
    for (; i + 7 < size; i += 8) {
        __m256d dataVec1 = _mm256_loadu_pd(&doubleValues[i]); // Load 4 values
        sumVec1 = _mm256_add_pd(sumVec1, dataVec1);  // Sum them into sumVec
        __m256d dataVec2 = _mm256_loadu_pd(&doubleValues[i+4]); // Load 4 values
        sumVec2 = _mm256_add_pd(sumVec2, dataVec2);  // Sum them into sumVec
    }

    // Reduce sumVec into a scalar
    double sumArray1[4];
    double sumArray2[4];
    _mm256_storeu_pd(sumArray1, sumVec1); // Store back to memory
    _mm256_storeu_pd(sumArray2, sumVec2); // Store back to memory
    double sum1 = sumArray1[0] + sumArray1[1] + sumArray1[2] + sumArray1[3];
    double sum2 = sumArray2[0] + sumArray2[1] + sumArray2[2] + sumArray2[3];

    // Handle remaining elements (if size isn't a multiple of 4)
    for (; i < size; i++) {
        sum1 += doubleValues[i];
    }

    return sum1 + sum2;
}

double sum_double(std::vector<double>& doubleValues, int64_t& time) {
    floatingExp2Integer::Timer timer;
    // double doubleSum = 0.0;
    // for (unsigned int i = 0; i < doubleValues.size(); i++) {
    //     doubleSum += doubleValues[i];
    // }
    // SECOND OPTION
    double doubleSum1 = 0.0;
    double doubleSum2 = 0.0;
    double doubleSum3 = 0.0;
    unsigned int i;
    for (i = 0; i + 2 < doubleValues.size(); i += 3) {
        doubleSum1 += doubleValues[i];
        doubleSum2 += doubleValues[i+1];
        doubleSum3 += doubleValues[i+2];
    }
    for (i = i; i < doubleValues.size(); i++) {
        doubleSum1 += doubleValues[i];
    }
    doubleSum1 += doubleSum2 + doubleSum3;
    // THIRD OPTION
    //double doubleSum1 = vectorizedSum(doubleValues);
    timer.stop();
    time = timer.time();
    return doubleSum1;
}

double sum_Dbl(std::vector<floatingExp2Integer::Dbl>& dblValues, int64_t& time) {
    floatingExp2Integer::Timer timer;
    floatingExp2Integer::Dbl dblSum = 0.0;
    for (unsigned int i = 0; i < dblValues.size(); i++) {
        dblSum += dblValues[i];
    }
    timer.stop();
    time = timer.time();
    return dblSum.asDouble();
}

double sum_Dbl2(std::vector<floatingExp2Integer::Dbl2>& dbl2Values, int64_t& time) {
    floatingExp2Integer::Timer timer;
    floatingExp2Integer::Dbl2 dbl2Sum = 0.0;
    for (unsigned int i = 0; i < dbl2Values.size(); i++) {
        dbl2Sum += dbl2Values[i];
    }
    timer.stop();
    time = timer.time();
    return dbl2Sum.asDouble();
}

double sum_Int64PosExp2Int64(std::vector<floatingExp2Integer::Int64PosExp2Int64>& int64PosExp2Int64Values, int64_t& time) {
    floatingExp2Integer::Timer timer;
    floatingExp2Integer::Int64PosExp2Int64 int64PosExp2Int64Sum = int64PosExp2Int64Values[0];
    for (unsigned int i = 1; i < int64PosExp2Int64Values.size(); i++) {
        int64PosExp2Int64Sum += int64PosExp2Int64Values[i];
    }
    timer.stop();
    time = timer.time();
    return int64PosExp2Int64Sum.asDouble();
}

double sum_Float64PosExp2Int64(std::vector<floatingExp2Integer::Float64PosExp2Int64>& float64PosExp2Int64Values, int64_t& time) {
    floatingExp2Integer::Timer timer;
    floatingExp2Integer::Float64PosExp2Int64 float64PosExp2Int64Sum = float64PosExp2Int64Values[0];
    for (unsigned int i = 1; i < float64PosExp2Int64Values.size(); i++) {
        float64PosExp2Int64Sum += float64PosExp2Int64Values[i];
    }
    timer.stop();
    time = timer.time();
    return float64PosExp2Int64Sum.asDouble();
}

double sum_Float64Exp2Int64(std::vector<floatingExp2Integer::Float64Exp2Int64>& float64Exp2Int64Values, int64_t& time) {
    floatingExp2Integer::Timer timer;
    floatingExp2Integer::Float64Exp2Int64 float64Exp2Int64Sum = float64Exp2Int64Values[0];
    for (unsigned int i = 1; i < float64Exp2Int64Values.size(); i++) {
        float64Exp2Int64Sum += float64Exp2Int64Values[i];
    }
    timer.stop();
    time = timer.time();
    return float64Exp2Int64Sum.asDouble();
}

double multiply_Int64PosExp2Int64(std::vector<floatingExp2Integer::Int64PosExp2Int64>& int64PosExp2Int64Values, int64_t& time) {
    floatingExp2Integer::Timer timer;
    floatingExp2Integer::Int64PosExp2Int64 int64PosExp2Int64Multiply = int64PosExp2Int64Values[0];
    for (unsigned int i = 1; i < int64PosExp2Int64Values.size(); i++) {
        int64PosExp2Int64Multiply *= int64PosExp2Int64Values[i];
    }
    timer.stop();
    time = timer.time();
    return int64PosExp2Int64Multiply.int64PosExp2Int64ToLog2();
}

double multiply_Float64PosExp2Int64(std::vector<floatingExp2Integer::Float64PosExp2Int64>& float64PosExp2Int64Values, int64_t& time) {
    floatingExp2Integer::Timer timer;
    floatingExp2Integer::Float64PosExp2Int64 float64PosExp2Int64Multiply = float64PosExp2Int64Values[0];
    for (unsigned int i = 1; i < float64PosExp2Int64Values.size(); i++) {
        float64PosExp2Int64Multiply *= float64PosExp2Int64Values[i];
    }
    timer.stop();
    time = timer.time();
    return float64PosExp2Int64Multiply.float64PosExp2Int64ToLog2();
}

double multiply_Float64Exp2Int64(std::vector<floatingExp2Integer::Float64Exp2Int64>& float64Exp2Int64Values, int64_t& time) {
    floatingExp2Integer::Timer timer;
    floatingExp2Integer::Float64Exp2Int64 float64Exp2Int64Multiply = float64Exp2Int64Values[0];
    for (unsigned int i = 1; i < float64Exp2Int64Values.size(); i++) {
        float64Exp2Int64Multiply *= float64Exp2Int64Values[i];
    }
    timer.stop();
    time = timer.time();
    return float64Exp2Int64Multiply.float64Exp2Int64ToLog2();
}

double sum_vectorization_Float64PosExp2Int64(std::vector<floatingExp2Integer::Float64PosExp2Int64>& float64PosExp2Int64Values, int64_t& time) {
    floatingExp2Integer::Timer timer;
    floatingExp2Integer::Float64PosExp2Int64 float64PosExp2Int64ValuesSum(float64PosExp2Int64Values);
    timer.stop();
    time = timer.time();
    return float64PosExp2Int64ValuesSum.asDouble();
}

double sum_parallelization_Dbl(std::vector<floatingExp2Integer::Dbl>& dblValues, int64_t& time) {
    floatingExp2Integer::Timer timer;
    floatingExp2Integer::Dbl dblParallelizationSum(dblValues);
    timer.stop();
    time = timer.time();
    return dblParallelizationSum.asDouble();
}

double sum_parallelization_Dbl2(std::vector<floatingExp2Integer::Dbl2>& dbl2Values, int64_t& time) {
    floatingExp2Integer::Timer timer;
    floatingExp2Integer::Dbl2 dbl2ParallelizationSum(dbl2Values);
    timer.stop();
    time = timer.time();
    return dbl2ParallelizationSum.asDouble();
}

double sum_parallelization_Dbl3(std::vector<floatingExp2Integer::Dbl3>& dbl3Values, int64_t& time) {
    floatingExp2Integer::Timer timer;
    floatingExp2Integer::Dbl3 dbl3ParallelizationSum(dbl3Values);
    timer.stop();
    time = timer.time();
    return dbl3ParallelizationSum.asDouble();
}

// double logsumexp_avx2(const std::vector<double>& x) {
//     size_t size = x.size();
//     __m256d maxVec = _mm256_set1_pd(-INFINITY);
    
//     // Compute max_x using AVX2
//     for (size_t i = 0; i < size; i += 4) {
//         __m256d data = _mm256_loadu_pd(&x[i]);
//         maxVec = _mm256_max_pd(maxVec, data);
//     }
//     double max_x[4];
//     _mm256_storeu_pd(max_x, maxVec);
//     double max_val = std::max({max_x[0], max_x[1], max_x[2], max_x[3]});

//     // Compute sum(exp(x - max_x))
//     __m256d sumVec = _mm256_setzero_pd();
//     for (size_t i = 0; i < size; i += 4) {
//         __m256d data = _mm256_loadu_pd(&x[i]);
//         __m256d expData = _mm256_exp2_pd(_mm256_sub_pd(data, _mm256_set1_pd(max_val)));
//         sumVec = _mm256_add_pd(sumVec, expData);
//     }
    
//     double sumExp[4];
//     _mm256_storeu_pd(sumExp, sumVec);
//     double sum = sumExp[0] + sumExp[1] + sumExp[2] + sumExp[3];

//     return max_val + std::log(sum);
// }

// double eigen_logsumexp(Eigen::VectorXd& eigen_values, int64_t& time) {
//     floatingExp2Integer::Timer timer;
//     double eigen_sum = eigen_values.array().exp().sum().log();;
//     timer.stop();
//     time = timer.time();
//     return eigen_sum;
// }

void loop(int n, double results[], std::int64_t time[]) {

    std::vector<double> doubleValues(n);
    Eigen::VectorXd eigen_values(n);
    std::vector<floatingExp2Integer::Dbl> dblValues(n);
    std::vector<floatingExp2Integer::Dbl2> dbl2Values(n);
    std::vector<floatingExp2Integer::Dbl3> dbl3Values(n);
    std::vector<double> log2Values(n);
    std::vector<floatingExp2Integer::Int64PosExp2Int64> int64PosExp2Int64Values(n);
    std::vector<floatingExp2Integer::Float64PosExp2Int64> float64PosExp2Int64Values(n);
    std::vector<floatingExp2Integer::Float64Exp2Int64> float64Exp2Int64Values(n);

    InitializeRandomNumbers(doubleValues);
    DoubleToEigenValues(doubleValues, eigen_values);
    DoubleToDblValues(doubleValues, dblValues);
    DoubleToDbl2Values(doubleValues, dbl2Values);
    DoubleToDbl3Values(doubleValues, dbl3Values);
    DoubleToLog2Values(doubleValues, log2Values);
    DoubleToInt64PosExp2Int64Values(doubleValues, int64PosExp2Int64Values);
    DoubleToFloat64PosExp2Int64Values(doubleValues, float64PosExp2Int64Values);
    DoubleToFloat64Exp2Int64Values(doubleValues, float64Exp2Int64Values);

    int64_t timer_double;
    double doubleSum = sum_double(doubleValues, timer_double);

    int64_t timer_Dbl;
    double dblSum = sum_Dbl(dblValues, timer_Dbl);

    int64_t timer_Dbl2;
    double dbl2Sum = sum_Dbl2(dbl2Values, timer_Dbl2);

    int64_t timer_log2;
    double log2sum = LogSumExp2Trick(log2Values, timer_log2);

    int64_t timer_Int64PosExp2Int64;
    double int64PosExp2Int64Sum = sum_Int64PosExp2Int64(int64PosExp2Int64Values, timer_Int64PosExp2Int64);

    int64_t timer_Float64PosExp2Int64;
    double float64PosExp2Int64Sum = sum_Float64PosExp2Int64(float64PosExp2Int64Values, timer_Float64PosExp2Int64);

    int64_t timer_Float64Exp2Int64;
    double float64Exp2Int64Sum = sum_Float64Exp2Int64(float64Exp2Int64Values, timer_Float64Exp2Int64);

    int64_t timer_multiply_log2;
    double log2sum_2 = Log2Multiply(log2Values, timer_multiply_log2);

    int64_t timer_multiply_Int64PosExp2Int64;
    double int64PosExp2Int64Multiply = multiply_Int64PosExp2Int64(int64PosExp2Int64Values, timer_multiply_Int64PosExp2Int64);

    int64_t timer_multiply_Float64PosExp2Int64;
    double float64PosExp2Int64Multiply = multiply_Float64PosExp2Int64(float64PosExp2Int64Values, timer_multiply_Float64PosExp2Int64);

    int64_t timer_multiply_Float64Exp2Int64;
    double float64Exp2Int64Multiply = multiply_Float64Exp2Int64(float64Exp2Int64Values, timer_multiply_Float64Exp2Int64);

    int64_t timer_Float64PosExp2Int64_vec;
    double float64PosExp2Int64VectorizationSum = sum_vectorization_Float64PosExp2Int64(float64PosExp2Int64Values, timer_Float64PosExp2Int64_vec);

    int64_t timer_Dbl_paral;
    double dblParallelizationSum = sum_parallelization_Dbl(dblValues, timer_Dbl_paral);

    int64_t timer_Dbl2_paral;
    double dbl2ParallelizationSum = sum_parallelization_Dbl2(dbl2Values, timer_Dbl2_paral);

    int64_t timer_Dbl3_paral;
    double dbl3ParallelizationSum = sum_parallelization_Dbl3(dbl3Values, timer_Dbl3_paral);


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
    std::cout << "Dbl sum parallelization:                " << dblParallelizationSum << std::endl;
    std::cout << "Dbl2 sum parallelization:               " << dbl2ParallelizationSum << std::endl;
    std::cout << "Dbl3 sum parallelization:               " << dbl3ParallelizationSum << std::endl;

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
    results[12] = dblParallelizationSum;
    results[13] = dbl2ParallelizationSum;
    results[14] = dbl3ParallelizationSum;

    time[0] = timer_double;
    time[1] = timer_Dbl;
    time[2] = timer_Dbl2;
    time[3] = timer_log2;
    time[4] = timer_Int64PosExp2Int64;
    time[5] = timer_Float64PosExp2Int64;
    time[6] = timer_Float64Exp2Int64;
    time[7] = timer_multiply_log2;
    time[8] = timer_multiply_Int64PosExp2Int64;
    time[9] = timer_multiply_Float64PosExp2Int64;
    time[10] = timer_multiply_Float64Exp2Int64;
    time[11] = timer_Float64PosExp2Int64_vec;
    time[12] = timer_Dbl_paral;
    time[13] = timer_Dbl2_paral;
    time[14] = timer_Dbl3_paral;

    std::cout << std::endl;
    std::cout << "Double time:                              " << timer_double << std::endl;
    std::cout << "Dbl time:                                 " << timer_Dbl << std::endl;
    std::cout << "Dbl2 time:                                " << timer_Dbl2 << std::endl;
    std::cout << "Log2 time:                                " << timer_log2 << std::endl;
    std::cout << "Int64PosExp2Int64 time:                   " << timer_Int64PosExp2Int64 << std::endl;
    std::cout << "Float64PosExp2Int64 time:                 " << timer_Float64PosExp2Int64 << std::endl;
    std::cout << "Float64Exp2Int64 time:                    " << timer_Float64Exp2Int64 << std::endl;
    std::cout << "Log2 multiply:                            " << timer_multiply_log2 << std::endl;
    std::cout << "Int64PosExp2Int64 multiply:               " << timer_multiply_Int64PosExp2Int64 << std::endl;
    std::cout << "Float64PosExp2Int64 multiply:             " << timer_multiply_Float64PosExp2Int64 << std::endl;
    std::cout << "Float64Exp2Int64 multiply:                " << timer_multiply_Float64Exp2Int64 << std::endl;
    std::cout << "Float64PosExp2Int64 time: vectorization:  " << timer_Float64PosExp2Int64_vec << std::endl;
    std::cout << "Dbl time parallelization:                 " << timer_Dbl_paral << std::endl;
    std::cout << "Dbl2 time parallelization:                " << timer_Dbl2_paral << std::endl;
    std::cout << "Dbl3 time parallelization:                " << timer_Dbl3_paral << std::endl;
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

void DoubleToEigenValues(std::vector<double>& doubleValues, Eigen::VectorXd& eigen_values) {
    for (unsigned int i = 0; i < eigen_values.size(); i++) {
        eigen_values[i] += doubleValues[i];
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

void DoubleToDbl3Values(std::vector<double>& doubleValues, std::vector<floatingExp2Integer::Dbl3>& dbl3Values) {
    for (unsigned int i = 0; i < dbl3Values.size(); i++) {
        dbl3Values[i] += doubleValues[i];
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

double LogSumExp2Trick(std::vector<double>& log2Values, int64_t& time) {
    floatingExp2Integer::Timer timer;
    double max = *std::max_element(log2Values.begin(), log2Values.end());
    double sum = 0;

    for (double val : log2Values) {
        sum += std::exp2(val - max);
    }
    timer.stop();
    time = timer.time();
    return max + std::log2(sum);
}

double Log2Multiply(std::vector<double>& log2Values, int64_t& time) {
    floatingExp2Integer::Timer timer;
    double sum = 0;
    for (double val : log2Values) {
        sum += val;
    }
    timer.stop();
    time = timer.time();
    return sum;
}


