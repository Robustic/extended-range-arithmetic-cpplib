#include <iostream>
#include <vector>
#include <immintrin.h>
#include <random>
#include <chrono>
#include <iomanip>
#include <cstring>
#include "exp_table.h"
#include "Timer.h"
#include "Dbl.h"
#include "Dbl2.h"
#include "Dbl3.h"
#include "./Int64PosExp2Int64/Int64PosExp2Int64.h"
#include "./Float64PosExp2Int64/Float64PosExp2Int64.h"
#include "./Float64Exp2Int64/Float64Exp2Int64.h"

// extern "C"{
//     double logsumexp_avx2(const double* x, size_t size);
// }

const uint32_t seed_val = 1337;


void InitializeRandomNumbers(std::vector<double>& vec);
void DoubleToDblValues(std::vector<double>& doubleValues, std::vector<floatingExp2Integer::Dbl>& dblValues);
void DoubleToDbl2Values(std::vector<double>& doubleValues, std::vector<floatingExp2Integer::Dbl2>& dbl2Values);
void DoubleToDbl3Values(std::vector<double>& doubleValues, std::vector<floatingExp2Integer::Dbl3>& dbl3Values);
void DoubleToLog2Values(std::vector<double>& dblValues, std::vector<double>& log2Values);
void DoubleToLogValues(std::vector<double>& dblValues, std::vector<double>& logValues);
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
    const int number_of_types = 16;
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

inline double logsumexp_avx2_cpp(const double* x, size_t size, __m256d* data1, __m256d* data2, __m256d* data3, __m256d* data4) {
    __m256d maxVec1 = _mm256_set1_pd(-INFINITY);
    __m256d maxVec2 = _mm256_set1_pd(-INFINITY);
    __m256d maxVec3 = _mm256_set1_pd(-INFINITY);
    __m256d maxVec4 = _mm256_set1_pd(-INFINITY);

    for (size_t i = 0; i < 128; i++) {
        data1[i] = _mm256_loadu_pd(&x[16*i]);
        data2[i] = _mm256_loadu_pd(&x[16*i+4]);
        data3[i] = _mm256_loadu_pd(&x[16*i+8]);
        data4[i] = _mm256_loadu_pd(&x[16*i+12]);
        maxVec1 = _mm256_max_pd(maxVec1, data1[i]);
        maxVec2 = _mm256_max_pd(maxVec2, data2[i]);
        maxVec3 = _mm256_max_pd(maxVec3, data3[i]);
        maxVec4 = _mm256_max_pd(maxVec4, data4[i]);
    }

    double max_x1[4];
    double max_x2[4];
    double max_x3[4];
    double max_x4[4];
    _mm256_storeu_pd(max_x1, maxVec1);
    _mm256_storeu_pd(max_x2, maxVec2);
    _mm256_storeu_pd(max_x3, maxVec3);
    _mm256_storeu_pd(max_x4, maxVec4);
    double max_val = fmax(fmax(fmax(fmax(max_x1[0], max_x1[1]), fmax(max_x1[2], max_x1[3])),
                               fmax(fmax(max_x2[0], max_x2[1]), fmax(max_x2[2], max_x2[3]))),
                          fmax(fmax(fmax(max_x3[0], max_x3[1]), fmax(max_x3[2], max_x3[3])),
                               fmax(fmax(max_x4[0], max_x4[1]), fmax(max_x4[2], max_x4[3]))));

    __m256d max = _mm256_set1_pd(max_val);

    __m256d sumVec1 = _mm256_setzero_pd();
    __m256d sumVec2 = _mm256_setzero_pd();
    __m256d sumVec3 = _mm256_setzero_pd();
    __m256d sumVec4 = _mm256_setzero_pd();

    __m256d one = _mm256_set1_pd(1.0);
    __m256d half = _mm256_set1_pd(0.5);
    __m256d sixth = _mm256_set1_pd(1.0 / 6.0);
    __m256d th24 = _mm256_set1_pd(1.0 / 24.0);
    __m256d th120 = _mm256_set1_pd(1.0 / 120.0);
    __m256d th720 = _mm256_set1_pd(1.0 / 720.0);

    for (size_t i = 0; i < 128; i++) {
        __m256d x1 = _mm256_sub_pd(data1[i], max);
        __m256d x2 = _mm256_sub_pd(data2[i], max);
        __m256d x3 = _mm256_sub_pd(data3[i], max);
        __m256d x4 = _mm256_sub_pd(data4[i], max);

        __m256d x1_2 = _mm256_mul_pd(x1, x1);
        __m256d x2_2 = _mm256_mul_pd(x2, x2);
        __m256d x3_2 = _mm256_mul_pd(x3, x3);
        __m256d x4_2 = _mm256_mul_pd(x4, x4);
    
        __m256d x1_3 = _mm256_mul_pd(x1, x1_2);
        __m256d x2_3 = _mm256_mul_pd(x2, x2_2);
        __m256d x3_3 = _mm256_mul_pd(x3, x3_2);
        __m256d x4_3 = _mm256_mul_pd(x4, x4_2);
    
        __m256d x1_4 = _mm256_mul_pd(x1, x1_3);
        __m256d x2_4 = _mm256_mul_pd(x2, x2_3);
        __m256d x3_4 = _mm256_mul_pd(x3, x3_3);
        __m256d x4_4 = _mm256_mul_pd(x4, x4_3);
    
        __m256d x1_5 = _mm256_mul_pd(x1, x1_4);
        __m256d x2_5 = _mm256_mul_pd(x2, x2_4);
        __m256d x3_5 = _mm256_mul_pd(x3, x3_4);
        __m256d x4_5 = _mm256_mul_pd(x4, x4_4);
    
        __m256d x1_6 = _mm256_mul_pd(x1, x1_5);
        __m256d x2_6 = _mm256_mul_pd(x2, x2_5);
        __m256d x3_6 = _mm256_mul_pd(x3, x3_5);
        __m256d x4_6 = _mm256_mul_pd(x4, x4_5);
    
        x1 = _mm256_add_pd(x1, one);
        x2 = _mm256_add_pd(x2, one);
        x3 = _mm256_add_pd(x3, one);
        x4 = _mm256_add_pd(x4, one);

        x1 = _mm256_fmadd_pd(x1_2, half, x1);
        x2 = _mm256_fmadd_pd(x2_2, half, x2);
        x3 = _mm256_fmadd_pd(x3_2, half, x3);
        x4 = _mm256_fmadd_pd(x4_2, half, x4);

        x1 = _mm256_fmadd_pd(x1_3, sixth, x1);
        x2 = _mm256_fmadd_pd(x2_3, sixth, x2);
        x3 = _mm256_fmadd_pd(x3_3, sixth, x3);
        x4 = _mm256_fmadd_pd(x4_3, sixth, x4);

        x1 = _mm256_fmadd_pd(x1_4, th24, x1);
        x2 = _mm256_fmadd_pd(x2_4, th24, x2);
        x3 = _mm256_fmadd_pd(x3_4, th24, x3);
        x4 = _mm256_fmadd_pd(x4_4, th24, x4);

        x1 = _mm256_fmadd_pd(x1_5, th120, x1);
        x2 = _mm256_fmadd_pd(x2_5, th120, x2);
        x3 = _mm256_fmadd_pd(x3_5, th120, x3);
        x4 = _mm256_fmadd_pd(x4_5, th120, x4);

        x1 = _mm256_fmadd_pd(x1_6, th720, x1);
        x2 = _mm256_fmadd_pd(x2_6, th720, x2);
        x3 = _mm256_fmadd_pd(x3_6, th720, x3);
        x4 = _mm256_fmadd_pd(x4_6, th720, x4);

        sumVec1 = _mm256_add_pd(sumVec1, x1);
        sumVec2 = _mm256_add_pd(sumVec2, x2);
        sumVec3 = _mm256_add_pd(sumVec3, x3);
        sumVec4 = _mm256_add_pd(sumVec4, x4);
    }

    double sumExp[4];
    _mm256_storeu_pd(sumExp, sumVec1);
    double sum = sumExp[0] + sumExp[1] + sumExp[2] + sumExp[3];
    _mm256_storeu_pd(sumExp, sumVec2);
    sum += sumExp[0] + sumExp[1] + sumExp[2] + sumExp[3];
    _mm256_storeu_pd(sumExp, sumVec3);
    sum += sumExp[0] + sumExp[1] + sumExp[2] + sumExp[3];
    _mm256_storeu_pd(sumExp, sumVec4);
    sum += sumExp[0] + sumExp[1] + sumExp[2] + sumExp[3];

    return max_val + std::log(sum);
}

double LogSumExpTrick(std::vector<double>& logValues) {
    double max = *std::max_element(logValues.begin(), logValues.end());
    double sum = 0;

    for (double val : logValues) {
        sum += std::exp(val - max);
    }
    return max + std::log(sum);
}

double sum_vectorization_log(std::vector<double>& logValues, int64_t& time) {
    double* ptr = logValues.data();
    floatingExp2Integer::Timer timer;

    unsigned int small_size = 4 * 4 * 128;
    unsigned int n_full = logValues.size() / small_size;
    unsigned int n_rest = logValues.size() % small_size;
    std::vector<double> temp(n_full + n_rest, 0.0);

    static __m256d data1[128];
    static __m256d data2[128];
    static __m256d data3[128];
    static __m256d data4[128];

    for (unsigned int i = 0; i < n_full; i++) {
        temp[i] = logsumexp_avx2_cpp(ptr + i * small_size, small_size, data1, data2, data3, data4);
    }
    for (unsigned int i = 0; i < n_rest; i++) {
        temp[n_full + i] = logValues[n_full * small_size + i];
    }

    double logVectorizationSum = LogSumExpTrick(temp);
    timer.stop();
    time = timer.time();

    return logVectorizationSum;
}

void loop(int n, double results[], std::int64_t time[]) {

    std::vector<double> doubleValues(n);
    std::vector<floatingExp2Integer::Dbl> dblValues(n);
    std::vector<floatingExp2Integer::Dbl2> dbl2Values(n);
    std::vector<floatingExp2Integer::Dbl3> dbl3Values(n);
    std::vector<double> log2Values(n);
    std::vector<double> logValues(n);
    std::vector<floatingExp2Integer::Int64PosExp2Int64> int64PosExp2Int64Values(n);
    std::vector<floatingExp2Integer::Float64PosExp2Int64> float64PosExp2Int64Values(n);
    std::vector<floatingExp2Integer::Float64Exp2Int64> float64Exp2Int64Values(n);

    InitializeRandomNumbers(doubleValues);
    DoubleToDblValues(doubleValues, dblValues);
    DoubleToDbl2Values(doubleValues, dbl2Values);
    DoubleToDbl3Values(doubleValues, dbl3Values);
    DoubleToLog2Values(doubleValues, log2Values);
    DoubleToLogValues(doubleValues, logValues);
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

    int64_t timer_log_vector;
    double logVectorizationSum = sum_vectorization_log(logValues, timer_log_vector);


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
    std::cout << "Double sum vectorization:               " << std::exp(logVectorizationSum) << std::endl;


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
    results[15] = std::exp(logVectorizationSum);

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
    time[15] = timer_log_vector;

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
    std::cout << "Float64PosExp2Int64 time vectorization:   " << timer_Float64PosExp2Int64_vec << std::endl;
    std::cout << "Dbl time parallelization:                 " << timer_Dbl_paral << std::endl;
    std::cout << "Dbl2 time parallelization:                " << timer_Dbl2_paral << std::endl;
    std::cout << "Dbl3 time parallelization:                " << timer_Dbl3_paral << std::endl;
    std::cout << "Double time vectorization:                " << timer_log_vector << std::endl;
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

void DoubleToLogValues(std::vector<double>& doubleValues, std::vector<double>& logValues) {
    for (unsigned int i = 0; i < logValues.size(); i++) {
        logValues[i] = std::log(doubleValues[i]);
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


