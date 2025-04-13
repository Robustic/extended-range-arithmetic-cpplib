#include <iostream>
#include <fstream>
#include <limits>
#include <functional>
#include <vector>
#include <immintrin.h>
#include <random>
#include <chrono>
#include <iomanip>
#include <cstring>
#include <cstdlib>
#include "Timer.h"
#include "Dbl.h"
#include "Dbl2.h"
#include "Log2Scale.h"
#include "./Int64PosExp2Int64/Int64PosExp2Int64.h"
#include "./Float64PosExp2Int64/Float64PosExp2Int64.h"
#include "./Float64LargeRangeNumber/Float64LargeRangeNumber.h"
#include "./Fukushima/Fukushima.h"

// *******  PARAMETERS  *******

// CHANGE THESE VALUES TO CALCULATE WITH DIFFERENT VALUE RANGE AND DIFFERENT COUNT OF NUMBERS

//constexpr double min_log2 = -33;//-10000;
//constexpr double max_log2 = -30;//-10;

constexpr size_t n[] = { 1000, 3000, 10000, 30000, 100000, 300000, 1000000, 3000000, 10000000, 30000000, 100000000 };
constexpr size_t n_rounds[] = { 10000, 10000, 10000, 3000, 1000, 300, 100, 30, 10, 3, 1 };
//constexpr size_t n_rounds[] = { 10000, 10000, 10000, 3000, 1000, 1000, 1000, 300, 100, 30, 10 };
//constexpr size_t n[] = { 1, 2, 3, 27, 33, 2000, 3000, 10000000 };
//constexpr size_t n_rounds[] = { 1, 1, 1, 1, 1, 1, 1, 1 };

// *******  COMMON FUNCTIONS  *******

constexpr size_t n_count = sizeof(n) / sizeof(n[0]);

void save_vector_as_binary(const std::vector<double>& vec, double min_log2, double max_log2) {
    std::string filename = "log2range_" + std::to_string(static_cast<int>(min_log2)) + "_to_" + std::to_string(static_cast<int>(max_log2)) + ".bin";
    std::ofstream outFile(filename, std::ios::binary);
    if (!outFile) {
        std::cerr << "Error opening file for writing!" << std::endl;
        return;
    }

    for (const double& value : vec) {
        outFile.write(reinterpret_cast<const char*>(&value), sizeof(double));
    }

    outFile.close();
}

const uint32_t seed_val = 1337;

void InitializeRandomNumbers(std::vector<double>& vec, double min_log2, double max_log2) {
    std::mt19937 rng;
    rng.seed(seed_val);
    std::uniform_real_distribution<double> unif(min_log2, max_log2);
    for (size_t i = 0; i < vec.size(); i++) {
        double a_random_double = unif(rng);
        vec[i] = a_random_double;
    }

    //save_vector_as_binary(vec, min_log2, max_log2);
}

// *******  TEMPLATE  *******

struct ResultCollection {
    std::string case_name;
    std::array<double, n_count> results_as_log2;
    std::array<int64_t, n_count> times;
};

template<typename T>
void calc_perf(const std::vector<double>& values_as_log2, std::vector<ResultCollection>& resultCollections,
        std::function<int64_t(const std::vector<T>&, double&, std::string&)> function) {

    ResultCollection rc;

    for (size_t i = 0; i < n_count; i++) {
        size_t n_current = n[i];
        size_t n_rounds_current = n_rounds[i];

        std::vector<T> values_as_T(n_current);
        T::log2s_to(values_as_log2, values_as_T);

        int64_t time_sum = 0;
        std::string case_name;
        double result_as_log2;

        for (size_t i = 0; i < n_rounds_current; i++) {
            time_sum += function(values_as_T, result_as_log2, case_name);
        }

        rc.case_name = case_name;
        rc.results_as_log2[i] = result_as_log2;
        rc.times[i] = time_sum / n_rounds_current;
    }
    resultCollections.push_back(rc);
}

void calc_perf_double(const std::vector<double>& values_as_log2, std::vector<ResultCollection>& resultCollections,
    std::function<int64_t(const std::vector<double>&, double&, std::string&)> function) {

    ResultCollection rc;

    for (size_t i = 0; i < n_count; i++) {
        size_t n_current = n[i];
        size_t n_rounds_current = n_rounds[i];

        std::vector<double> values_as_double(n_current);
        for (size_t i = 0; i < values_as_double.size(); i++) {
            values_as_double[i] = std::exp2(values_as_log2[i]);
        }

        int64_t time_sum = 0;
        std::string case_name;
        double result_as_log2;

        for (size_t i = 0; i < n_rounds_current; i++) {
            time_sum += function(values_as_double, result_as_log2, case_name);
        }

        rc.case_name = case_name;
        rc.results_as_log2[i] = result_as_log2;
        rc.times[i] = time_sum / n_rounds_current;
    }
    resultCollections.push_back(rc);
}

void calc_perf_log2scale(const std::vector<double>& values_as_log2, std::vector<ResultCollection>& resultCollections,
    std::function<int64_t(const std::vector<double>&, double&, std::string&)> function) {

    ResultCollection rc;

    for (size_t i = 0; i < n_count; i++) {
        size_t n_current = n[i];
        size_t n_rounds_current = n_rounds[i];

        std::vector<double> values_as_log2_limitted(n_current);
        for (size_t i = 0; i < values_as_log2_limitted.size(); i++) {
            values_as_log2_limitted[i] = values_as_log2[i];
        }

        int64_t time_sum = 0;
        std::string case_name;
        double result_as_log2;

        for (size_t i = 0; i < n_rounds_current; i++) {
            time_sum += function(values_as_log2_limitted, result_as_log2, case_name);
        }

        rc.case_name = case_name;
        rc.results_as_log2[i] = result_as_log2;
        rc.times[i] = time_sum / n_rounds_current;
    }
    resultCollections.push_back(rc);
}

void calc_perf_Float64LargeRangeNumber(const std::vector<double>& values_as_log2, std::vector<ResultCollection>& resultCollections,
    std::function<int64_t(const std::vector<double>&, double&, std::string&)> function) {

    ResultCollection rc;

    for (size_t i = 0; i < n_count; i++) {
        size_t n_current = n[i];
        size_t n_rounds_current = n_rounds[i];

        std::vector<double> values_as_double(n_current);
        floatingExp2Integer::Float64LargeRangeNumber::log2s_to(values_as_log2, values_as_double);

        int64_t time_sum = 0;
        std::string case_name;
        double result_as_log2;

        for (size_t i = 0; i < n_rounds_current; i++) {
            time_sum += function(values_as_double, result_as_log2, case_name);
        }

        rc.case_name = case_name;
        rc.results_as_log2[i] = result_as_log2;
        rc.times[i] = time_sum / n_rounds_current;
    }
    resultCollections.push_back(rc);
}

// *******  CASES  *******

// *******  double  *******

int64_t sum_sequential_dbl(const std::vector<double>& values_as_double, double& result_as_log2, std::string& case_name) {
    case_name = __func__;

    floatingExp2Integer::Timer timer;
    double sum = 0.0;
    for (size_t i = 0; i < values_as_double.size(); i++) {
        sum += values_as_double[i];
    }
    timer.stop();

    result_as_log2 = std::log2(sum);
    return timer.time();
}

int64_t sum_parallel_dbl(const std::vector<double>& values_as_double, double& result_as_log2, std::string& case_name) {
    case_name = __func__;

    floatingExp2Integer::Timer timer;    
    double res = 0.0;
    constexpr size_t n_parallel = 4;

    if (values_as_double.size() < 2 * 8 * n_parallel) {
        for (size_t k = 0; k < values_as_double.size(); k++) {
            res += values_as_double[k];
        }
    }
    else {
        __m512d vres[n_parallel];
        __m512d vb[n_parallel];

        for (size_t m = 0; m < n_parallel; m++) {
            vres[m] = _mm512_loadu_pd(&values_as_double[8 * m]);
        }

        size_t i = 8 * n_parallel;
        for (i = 8 * n_parallel; i + (8 * n_parallel - 1) < values_as_double.size(); i += 8 * n_parallel) {
            for (size_t m = 0; m < n_parallel; m++) {
                vb[m] = _mm512_loadu_pd(&values_as_double[i + 8 * m]);
                vres[m] = vres[m] + vb[m];
            }
        }

        for (size_t m = 0; m < n_parallel; m++) {
            for (size_t k = 0; k < 8; k++) {
                res += vres[m][k];
            }
        }

        for (i = i; i < values_as_double.size(); i++) {
            res += values_as_double[i];
        }
    }
    timer.stop();

    result_as_log2 = std::log2(res);
    return timer.time();
}

int64_t multiply_sequential_dbl(const std::vector<double>& values_as_double, double& result_as_log2, std::string& case_name) {
    case_name = __func__;

    floatingExp2Integer::Timer timer;
    double multiplied = 1.0;
    for (size_t i = 0; i < values_as_double.size(); i++) {
        multiplied *= values_as_double[i];
    }
    timer.stop();

    result_as_log2 = std::log2(multiplied);
    return timer.time();
}

int64_t multiply_parallel_dbl(const std::vector<double>& values_as_double, double& result_as_log2, std::string& case_name) {
    case_name = __func__;

    floatingExp2Integer::Timer timer;
    double res = 1.0;
    constexpr size_t n_parallel = 4;

    if (values_as_double.size() < 2 * 8 * n_parallel) {
        for (size_t k = 0; k < values_as_double.size(); k++) {
            res *= values_as_double[k];
        }
    }
    else {
        __m512d vres[n_parallel];
        __m512d vb[n_parallel];

        for (size_t m = 0; m < n_parallel; m++) {
            vres[m] = _mm512_loadu_pd(&values_as_double[8 * m]);
        }

        size_t i = 8 * n_parallel;
        for (i = 8 * n_parallel; i + (8 * n_parallel - 1) < values_as_double.size(); i += 8 * n_parallel) {
            for (size_t m = 0; m < n_parallel; m++) {
                vb[m] = _mm512_loadu_pd(&values_as_double[i + 8 * m]);
                vres[m] = vres[m] * vb[m];
            }
        }

        for (size_t m = 0; m < n_parallel; m++) {
            for (size_t k = 0; k < 8; k++) {            
                res *= vres[m][k];
            }
        }

        for (i = i; i < values_as_double.size(); i++) {
            res *= values_as_double[i];
        }
    }
    timer.stop();

    result_as_log2 = std::log2(res);
    return timer.time();
}

// *******  log2scale  *******

int64_t sum_sequential_log2scale(const std::vector<double>& values_as_log2, double& result_as_log2, std::string& case_name) {
    case_name = __func__;

    floatingExp2Integer::Timer timer;
    double sum = values_as_log2[0];
    for (size_t i = 1; i < values_as_log2.size(); i++) {
        double value = values_as_log2[i];
        if (sum - value > 54) {
            sum = sum;
        }
        else if (value - sum > 54) {
            sum = value;
        }
        else if (sum > value) {
            sum = sum + std::log2(1 + std::exp2(value - sum));
        }
        else {
            sum = value + std::log2(std::exp2(sum - value) + 1);
        }
    }
    timer.stop();

    result_as_log2 = sum;
    return timer.time();
}

// More efficient function compiled with Intel Compiler to utilize SIMD
int64_t sum_parallel_log2scale(const std::vector<double>& values_as_log2, double& result_as_log2, std::string& case_name) {
    case_name = __func__;

    floatingExp2Integer::Timer timer;
    double sum = 0.0;
    double max = *std::max_element(values_as_log2.begin(), values_as_log2.end());
    for (size_t i = 0; i < values_as_log2.size(); i++) {
        sum += std::exp2(values_as_log2[i] - max);
    }
    timer.stop();

    result_as_log2 = max + std::log2(sum);
    return timer.time();
}

//// ********************************************
// 
//// @echo off
////set OMP_NUM_THREADS = 8
////icx /Qiopenmp -O3 -Qstd=c++23 /Qpar-num-threads:1 -march=skylake-avx512 -ffast-math Float64.try.cpp Float64ExtendedExp.cpp -o Float64.exe
////echo Float64.exe
//
//int64_t sum_parallel_log2scale_OPTION_2_SIMD(const std::vector<double>& values_as_log2, double& result_as_log2, std::string& case_name) {
//    case_name = __func__;
//
//    floatingExp2Integer::Timer timer;
//    size_t i = 0;
//    double sum = 0.0;
//    double max = *std::max_element(values_as_log2.begin(), values_as_log2.end());
//    if (8 <= values_as_log2.size()) {
//        __m256d vmax = _mm256_set1_pd(max);
//        __m256d va = _mm256_loadu_pd(&values_as_log2[0]);
//        va = va - vmax;
//        __m256d vsum = _mm256_exp2_pd(va);
//        for (i = 4; i + 3 < values_as_log2.size(); i += 4) {
//            __m256d vb = _mm256_loadu_pd(&values_as_log2[i]);
//            vb = vb - vmax;
//            __m256d vexp = _mm256_exp2_pd(vb);
//            vsum = vsum + vexp;
//        }
//
//        for (size_t k = 0; k < 4; k++) {
//            sum += vsum[k];
//        }
//    }
//
//    for (i = i; i < values_as_log2.size(); i++) {
//        sum += std::exp2(values_as_log2[i] - max);
//    }
//    timer.stop();
//
//    result_as_log2 = max + std::log2(sum);
//    return timer.time();
//}
//
//static const size_t SIZE = 4096;
//
//inline double calculate_array_sum_log2_split(const std::vector<double>& values, size_t start) {
//    double sum = 0.0;
//    double max = *std::max_element(values.begin() + start, values.begin() + start + SIZE);
//
//    __m256d vmax = _mm256_set1_pd(max);
//    __m256d va = _mm256_loadu_pd(&values[start]);
//    va = va - vmax;
//    __m256d vsum = _mm256_exp2_pd(va);
//
//    for (size_t i = start + 4; i + 3 < start + SIZE; i += 4) {
//        __m256d vb = _mm256_loadu_pd(&values[i]);
//        vb = vb - vmax;
//        __m256d vexp = _mm256_exp2_pd(vb);
//        vsum = vsum + vexp;
//    }
//
//    for (size_t k = 0; k < 4; k++) {
//        sum += vsum[k];
//    }
//    double result = max + std::log2(sum);
//
//    return result;
//}
//
//inline double calculate_array_sum_log2_rest(const std::vector<double>& values, size_t start) {
//    size_t i = 0;
//    double sum = 0.0;
//    double max = *std::max_element(values.begin() + start, values.end());
//    if (8 <= values.size() - start) {
//        __m256d vmax = _mm256_set1_pd(max);
//        __m256d va = _mm256_loadu_pd(&values[0]);
//        va = va - vmax;
//        __m256d vsum = _mm256_exp2_pd(va);
//        for (i = 4 + start; i + 3 < values.size(); i += 4) {
//            __m256d vb = _mm256_loadu_pd(&values[i]);
//            vb = vb - vmax;
//            __m256d vexp = _mm256_exp2_pd(vb);
//            vsum = vsum + vexp;
//        }
//
//        for (size_t k = 0; k < 4; k++) {
//            sum += vsum[k];
//        }
//    }
//
//    for (i = i; i < values.size(); i++) {
//        sum += std::exp2(values[i] - max);
//    }
//    double result = max + std::log2(sum);
//
//    return result;
//}
//
//int64_t sum_parallel_log2scale_OPTION_3_SIMD_betterCaching(const std::vector<double>& values_as_log2, double& result_as_log2, std::string& case_name) {
//    case_name = __func__;
//
//    floatingExp2Integer::Timer timer;
//    size_t n_count = ((values_as_log2.size() - 1) / SIZE) + 1;
//
//    std::vector<double> collected_values(n_count);
//
//    for (size_t i = 0; i < n_count - 1; i++) {
//        collected_values[i] = calculate_array_sum_log2_split(values_as_log2, i * SIZE);
//    }
//    collected_values[n_count - 1] = calculate_array_sum_log2_rest(values_as_log2, (n_count - 1) * SIZE);
//
//    double result = calculate_array_sum_log2_rest(collected_values, 0);
//    timer.stop();
//
//    result_as_log2 = result;
//    return timer.time();
//}
//
//int64_t sum_parallel_log2scale_OPTION_4_SIMD_betterCaching_threading(const std::vector<double>& values_as_log2, double& result_as_log2, std::string& case_name) {
//    case_name = __func__;
//
//    floatingExp2Integer::Timer timer;
//    size_t n_count = ((values_as_log2.size() - 1) / SIZE) + 1;
//
//    std::vector<double> collected_values(n_count);
//
//#pragma omp parallel for schedule(static, 1)
//    for (size_t i = 0; i < n_count - 1; i++) {
//        collected_values[i] = calculate_array_sum_log2_split(values_as_log2, i * SIZE);
//    }
//    collected_values[n_count - 1] = calculate_array_sum_log2_rest(values_as_log2, (n_count - 1) * SIZE);
//
//    double result = calculate_array_sum_log2_rest(collected_values, 0);
//    timer.stop();
//
//    result_as_log2 = result;
//    return timer.time();
//}
//
//// ********************************************

int64_t multiply_sequential_log2scale(const std::vector<double>& values_as_log2, double& result_as_log2, std::string& case_name) {
    case_name = __func__;

    floatingExp2Integer::Timer timer;
    double sum = 0.0;
    for (size_t i = 0; i < values_as_log2.size(); i++) {
        // log(A*B) = log(A) + log(B) 
        sum += values_as_log2[i];
    }
    timer.stop();

    result_as_log2 = sum;
    return timer.time();
}

int64_t multiply_parallel_log2scale(const std::vector<double>& values_as_log2, double& result_as_log2, std::string& case_name) {
    case_name = __func__;

    floatingExp2Integer::Timer timer;
    double res = 0.0;
    constexpr size_t n_parallel = 4;

    if (values_as_log2.size() < 2 * 8 * n_parallel) {
        for (size_t k = 0; k < values_as_log2.size(); k++) {
            res += values_as_log2[k];
        }
    }
    else {
        __m512d vres[n_parallel];
        __m512d vb[n_parallel];

        for (size_t m = 0; m < n_parallel; m++) {
            vres[m] = _mm512_loadu_pd(&values_as_log2[8 * m]);
        }

        size_t i = 8 * n_parallel;
        for (i = 8 * n_parallel; i + (8 * n_parallel - 1) < values_as_log2.size(); i += 8 * n_parallel) {
            for (size_t m = 0; m < n_parallel; m++) {
                vb[m] = _mm512_loadu_pd(&values_as_log2[i + 8 * m]);
                vres[m] = vres[m] + vb[m];
            }
        }

        for (size_t m = 0; m < n_parallel; m++) {
            for (size_t k = 0; k < 8; k++) {
                res += vres[m][k];
            }
        }

        for (i = i; i < values_as_log2.size(); i++) {
            res += values_as_log2[i];
        }
    }
    timer.stop();

    result_as_log2 = res;
    return timer.time();
}

// *******  Dbl1  *******

int64_t sum_sequential_Dbl1(const std::vector<floatingExp2Integer::Dbl>& values, double& result_as_log2, std::string& case_name) {
    case_name = __func__;

    floatingExp2Integer::Timer timer;
    floatingExp2Integer::Dbl sum(0.0);
    for (size_t i = 0; i < values.size(); i++) {
        sum += values[i];
    }
    timer.stop();

    result_as_log2 = sum.as_log2();
    return timer.time();
}

int64_t sum_parallel_Dbl1(const std::vector<floatingExp2Integer::Dbl>& values, double& result_as_log2, std::string& case_name) {
    case_name = __func__;

    floatingExp2Integer::Timer timer;
    floatingExp2Integer::Dbl sum;
    sum.sum(values);
    timer.stop();

    result_as_log2 = sum.as_log2();
    return timer.time();
}

int64_t  multiply_sequential_Dbl1(const std::vector<floatingExp2Integer::Dbl>& values, double& result_as_log2, std::string& case_name) {
    case_name = __func__;

    floatingExp2Integer::Timer timer;
    floatingExp2Integer::Dbl multiplied(1.0);
    for (size_t i = 0; i < values.size(); i++) {
        multiplied *= values[i];
    }
    timer.stop();

    result_as_log2 = multiplied.as_log2();
    return timer.time();
}

int64_t  multiply_parallel_Dbl1(const std::vector<floatingExp2Integer::Dbl>& values, double& result_as_log2, std::string& case_name) {
    case_name = __func__;

    floatingExp2Integer::Timer timer;
    floatingExp2Integer::Dbl multiplied;
    multiplied.multiply(values);
    timer.stop();

    result_as_log2 = multiplied.as_log2();
    return timer.time();
}

// *******  Dbl2  *******

int64_t sum_sequential_Dbl2(const std::vector<floatingExp2Integer::Dbl2>& values, double& result_as_log2, std::string& case_name) {
    case_name = __func__;

    floatingExp2Integer::Timer timer;
    floatingExp2Integer::Dbl2 sum(0.0);
    for (size_t i = 0; i < values.size(); i++) {
        sum += values[i];
    }
    timer.stop();

    result_as_log2 = sum.as_log2();
    return timer.time();
}

int64_t  sum_parallel_Dbl2(const std::vector<floatingExp2Integer::Dbl2>& values, double& result_as_log2, std::string& case_name) {
    case_name = __func__;

    floatingExp2Integer::Timer timer;
    floatingExp2Integer::Dbl2 sum;
    sum.sum(values);
    timer.stop();

    result_as_log2 = sum.as_log2();
    return timer.time();
}

int64_t  multiply_sequential_Dbl2(const std::vector<floatingExp2Integer::Dbl2>& values, double& result_as_log2, std::string& case_name) {
    case_name = __func__;

    floatingExp2Integer::Timer timer;
    floatingExp2Integer::Dbl2 multiplied(1.0);
    for (size_t i = 0; i < values.size(); i++) {
        multiplied *= values[i];
    }
    timer.stop();

    result_as_log2 = multiplied.as_log2();
    return timer.time();
}

int64_t  multiply_parallel_Dbl2(const std::vector<floatingExp2Integer::Dbl2>& values, double& result_as_log2, std::string& case_name) {
    case_name = __func__;

    floatingExp2Integer::Timer timer;
    floatingExp2Integer::Dbl2 multiplied;
    multiplied.multiply(values);
    timer.stop();

    result_as_log2 = multiplied.as_log2();
    return timer.time();
}

// *******  Fukushima  *******

int64_t  sum_sequential_Fukushima(const std::vector<floatingExp2Integer::Fukushima>& values, double& result_as_log2, std::string& case_name) {
    case_name = __func__;

    floatingExp2Integer::Timer timer;
    floatingExp2Integer::Fukushima sum = values[0];
    for (size_t i = 1; i < values.size(); i++) {
        sum += values[i];
    }
    timer.stop();

    result_as_log2 = sum.as_log2();
    return timer.time();
}

int64_t  sum_parallel_Fukushima(const std::vector<floatingExp2Integer::Fukushima>& values, double& result_as_log2, std::string& case_name) {
    case_name = __func__;

    floatingExp2Integer::Timer timer;
    floatingExp2Integer::Fukushima collector;
    collector.sum(values);
    timer.stop();

    result_as_log2 = collector.as_log2();
    return timer.time();
}

int64_t  multiply_sequential_Fukushima(const std::vector<floatingExp2Integer::Fukushima>& values, double& result_as_log2, std::string& case_name) {
    case_name = __func__;

    floatingExp2Integer::Timer timer;
    floatingExp2Integer::Fukushima res(1.0);
    for (size_t i = 0; i < values.size(); i++) {
        res *= values[i];
    }
    timer.stop();

    result_as_log2 = res.as_log2();
    return timer.time();
}

int64_t  multiply_parallel_Fukushima(const std::vector<floatingExp2Integer::Fukushima>& values, double& result_as_log2, std::string& case_name) {
    case_name = __func__;

    floatingExp2Integer::Timer timer;
    floatingExp2Integer::Fukushima collector;
    collector.multiply(values);
    timer.stop();

    result_as_log2 = collector.as_log2();
    return timer.time();
}

// *******  Float64LargeRangeNumber  *******

int64_t  sum_sequential_Float64LargeRangeNumber(const std::vector<double>& values_as_lrn, double& result_as_log2, std::string& case_name) {
    case_name = __func__;

    floatingExp2Integer::Timer timer;
    double sum = values_as_lrn[0];
    for (size_t i = 1; i < values_as_lrn.size(); i++) {
        double value = values_as_lrn[i];
        sum = floatingExp2Integer::Float64LargeRangeNumber::sum(sum, value);
    }
    timer.stop();

    result_as_log2 = floatingExp2Integer::Float64LargeRangeNumber::as_log2(sum);
    return timer.time();
}

int64_t  sum_parallel_Float64LargeRangeNumber(const std::vector<double>& values_as_lrn, double& result_as_log2, std::string& case_name) {
    case_name = __func__;

    floatingExp2Integer::Timer timer;
    double sum = floatingExp2Integer::Float64LargeRangeNumber::sum(values_as_lrn);
    timer.stop();

    result_as_log2 = floatingExp2Integer::Float64LargeRangeNumber::as_log2(sum);
    return timer.time();
}

int64_t  multiply_sequential_Float64LargeRangeNumber(const std::vector<double>& values_as_lrn, double& result_as_log2, std::string& case_name) {
    case_name = __func__;

    floatingExp2Integer::Timer timer;
    double res = values_as_lrn[0];
    for (size_t i = 1; i < values_as_lrn.size(); i++) {
        double value = values_as_lrn[i];
        res = floatingExp2Integer::Float64LargeRangeNumber::multiply(res, value);
    }
    timer.stop();

    result_as_log2 = floatingExp2Integer::Float64LargeRangeNumber::as_log2(res);
    return timer.time();
}

int64_t  multiply_parallel_Float64LargeRangeNumber(const std::vector<double>& values_as_lrn, double& result_as_log2, std::string& case_name) {
    case_name = __func__;

    floatingExp2Integer::Timer timer;
    double res = floatingExp2Integer::Float64LargeRangeNumber::multiply(values_as_lrn);
    timer.stop();

    result_as_log2 = floatingExp2Integer::Float64LargeRangeNumber::as_log2(res);
    return timer.time();
}

// *******  Int64PosExp2Int64  *******

int64_t  sum_sequential_Int64PosExp2Int64(const std::vector<floatingExp2Integer::Int64PosExp2Int64>& values, double& result_as_log2, std::string& case_name) {
    case_name = __func__;

    floatingExp2Integer::Timer timer;
    floatingExp2Integer::Int64PosExp2Int64 res = values[0];
    for (size_t i = 1; i < values.size(); i++) {
        res += values[i];
    }
    timer.stop();

    result_as_log2 = res.as_log2();
    return timer.time();
}

int64_t  sum_parallel_Int64PosExp2Int64(const std::vector<floatingExp2Integer::Int64PosExp2Int64>& values, double& result_as_log2, std::string& case_name) {
    case_name = __func__;

    floatingExp2Integer::Timer timer;
    floatingExp2Integer::Int64PosExp2Int64 res;
    res.sum(values);
    timer.stop();

    result_as_log2 = res.as_log2();
    return timer.time();
}

int64_t  multiply_sequential_Int64PosExp2Int64(const std::vector<floatingExp2Integer::Int64PosExp2Int64>& values, double& result_as_log2, std::string& case_name) {
    case_name = __func__;

    floatingExp2Integer::Timer timer;
    floatingExp2Integer::Int64PosExp2Int64 res = values[0];
    for (size_t i = 1; i < values.size(); i++) {
        res *= values[i];
    }
    timer.stop();

    result_as_log2 = res.as_log2();
    return timer.time();
}

int64_t  multiply_parallel_Int64PosExp2Int64(const std::vector<floatingExp2Integer::Int64PosExp2Int64>& values, double& result_as_log2, std::string& case_name) {
    case_name = __func__;

    floatingExp2Integer::Timer timer;
    floatingExp2Integer::Int64PosExp2Int64 res;
    res.multiply(values);
    timer.stop();

    result_as_log2 = res.as_log2();
    return timer.time();
}

// *******  Float64PosExp2Int64  *******

int64_t  sum_sequential_Float64PosExp2Int64(const std::vector<floatingExp2Integer::Float64PosExp2Int64>& values, double& result_as_log2, std::string& case_name) {
    case_name = __func__;

    floatingExp2Integer::Timer timer;
    floatingExp2Integer::Float64PosExp2Int64 res = values[0];
    for (size_t i = 1; i < values.size(); i++) {
        res += values[i];
    }
    timer.stop();

    result_as_log2 = res.as_log2();
    return timer.time();
}

int64_t  sum_parallel_Float64PosExp2Int64(const std::vector<floatingExp2Integer::Float64PosExp2Int64>& values, double& result_as_log2, std::string& case_name) {
    case_name = __func__;

    floatingExp2Integer::Timer timer;
    floatingExp2Integer::Float64PosExp2Int64 res;
    res.sum(values);
    timer.stop();

    result_as_log2 = res.as_log2();
    return timer.time();
}

int64_t  multiply_sequential_Float64PosExp2Int64(const std::vector<floatingExp2Integer::Float64PosExp2Int64>& values, double& result_as_log2, std::string& case_name) {
    case_name = __func__;

    floatingExp2Integer::Timer timer;
    floatingExp2Integer::Float64PosExp2Int64 res = values[0];
    for (size_t i = 1; i < values.size(); i++) {
        res *= values[i];
    }
    timer.stop();

    result_as_log2 = res.as_log2();
    return timer.time();
}

int64_t  multiply_parallel_Float64PosExp2Int64(const std::vector<floatingExp2Integer::Float64PosExp2Int64>& values, double& result_as_log2, std::string& case_name) {
    case_name = __func__;

    floatingExp2Integer::Timer timer;
    floatingExp2Integer::Float64PosExp2Int64 res;
    res.multiply(values);
    timer.stop();

    result_as_log2 = res.as_log2();
    return timer.time();
}

// *******  MAIN  *******

int main(int argc, char* argv[]) {
    std::cout << "START  " << floatingExp2Integer::Timer::current_time() << std::endl;

    double min_log2;
    double max_log2;

    if (argc != 3) {
        min_log2 = -33.0;//-10000;
        max_log2 = -30.0;//-10;
    }
    else {
        min_log2 = std::atof(argv[1]);
        max_log2 = std::atof(argv[2]);
    }

    std::cout << "MIN_limit_log2_of_random_numbers: " << min_log2 << std::endl;
    std::cout << "MAX_limit_log2_of_random_numbers: " << max_log2 << std::endl;

    std::vector<ResultCollection> resultCollections;

    constexpr unsigned int n_max = *std::max_element(n, n + n_count);
    std::vector<double> double_values(n_max);
    InitializeRandomNumbers(double_values, min_log2, max_log2);

    std::cout << "RANDOM NUMBERS DEFINED" << std::endl;


    calc_perf_double(double_values, resultCollections, sum_sequential_dbl);
    calc_perf_log2scale(double_values, resultCollections, sum_sequential_log2scale);
    calc_perf<floatingExp2Integer::Dbl>(double_values, resultCollections, sum_sequential_Dbl1);
    calc_perf<floatingExp2Integer::Dbl2>(double_values, resultCollections, sum_sequential_Dbl2);
    calc_perf<floatingExp2Integer::Fukushima>(double_values, resultCollections, sum_sequential_Fukushima);
    calc_perf_Float64LargeRangeNumber(double_values, resultCollections, sum_sequential_Float64LargeRangeNumber);
    calc_perf<floatingExp2Integer::Int64PosExp2Int64>(double_values, resultCollections, sum_sequential_Int64PosExp2Int64);
    calc_perf<floatingExp2Integer::Float64PosExp2Int64>(double_values, resultCollections, sum_sequential_Float64PosExp2Int64);

    calc_perf_double(double_values, resultCollections, sum_parallel_dbl);
    calc_perf_log2scale(double_values, resultCollections, sum_parallel_log2scale);
    calc_perf<floatingExp2Integer::Dbl>(double_values, resultCollections, sum_parallel_Dbl1);
    calc_perf<floatingExp2Integer::Dbl2>(double_values, resultCollections, sum_parallel_Dbl2);
    calc_perf<floatingExp2Integer::Fukushima>(double_values, resultCollections, sum_parallel_Fukushima);
    calc_perf_Float64LargeRangeNumber(double_values, resultCollections, sum_parallel_Float64LargeRangeNumber);
    calc_perf<floatingExp2Integer::Int64PosExp2Int64>(double_values, resultCollections, sum_parallel_Int64PosExp2Int64);
    calc_perf<floatingExp2Integer::Float64PosExp2Int64>(double_values, resultCollections, sum_parallel_Float64PosExp2Int64);

    calc_perf_double(double_values, resultCollections, multiply_sequential_dbl);
    calc_perf_log2scale(double_values, resultCollections, multiply_sequential_log2scale);
    calc_perf<floatingExp2Integer::Dbl>(double_values, resultCollections, multiply_sequential_Dbl1);
    calc_perf<floatingExp2Integer::Dbl2>(double_values, resultCollections, multiply_sequential_Dbl2);
    calc_perf<floatingExp2Integer::Fukushima>(double_values, resultCollections, multiply_sequential_Fukushima);
    calc_perf_Float64LargeRangeNumber(double_values, resultCollections, multiply_sequential_Float64LargeRangeNumber);
    calc_perf<floatingExp2Integer::Int64PosExp2Int64>(double_values, resultCollections, multiply_sequential_Int64PosExp2Int64);
    calc_perf<floatingExp2Integer::Float64PosExp2Int64>(double_values, resultCollections, multiply_sequential_Float64PosExp2Int64);

    calc_perf_double(double_values, resultCollections, multiply_parallel_dbl);
    calc_perf_log2scale(double_values, resultCollections, multiply_parallel_log2scale);
    calc_perf<floatingExp2Integer::Dbl>(double_values, resultCollections, multiply_parallel_Dbl1);
    calc_perf<floatingExp2Integer::Dbl2>(double_values, resultCollections, multiply_parallel_Dbl2);
    calc_perf<floatingExp2Integer::Fukushima>(double_values, resultCollections, multiply_parallel_Fukushima);
    calc_perf_Float64LargeRangeNumber(double_values, resultCollections, multiply_parallel_Float64LargeRangeNumber);
    calc_perf<floatingExp2Integer::Int64PosExp2Int64>(double_values, resultCollections, multiply_parallel_Int64PosExp2Int64);
    calc_perf<floatingExp2Integer::Float64PosExp2Int64>(double_values, resultCollections, multiply_parallel_Float64PosExp2Int64);


    std::cout << std::setprecision(std::numeric_limits<double>::digits10 + 1);

    std::cout << std::endl << "RESULTS" << std::endl << "n ";
    for (size_t i = 0; i < n_count; i++) {
        std::cout << n[i] << " ";
    }

    std::cout << std::endl;
    for (size_t f = 0; f < resultCollections.size(); f++) {
        ResultCollection& current = resultCollections[f];
        std::cout << current.case_name << " ";
        for (size_t i = 0; i < n_count; i++) {
            std::cout << current.results_as_log2[i] << " ";
        }
        std::cout << std::endl;
    }

    std::cout << std::endl << "CALCULATION_TIME" << std::endl << "Rounds ";
    for (size_t i = 0; i < n_count; i++) {
        std::cout << n_rounds[i] << " ";
    }
    std::cout << std::endl << "n ";
    for (size_t i = 0; i < n_count; i++) {
        std::cout << n[i] << " ";
    }

    std::cout << std::endl;
    for (size_t f = 0; f < resultCollections.size(); f++) {
        ResultCollection& current = resultCollections[f];
        std::cout << current.case_name << " ";
        for (size_t i = 0; i < n_count; i++) {
            std::cout << current.times[i] << " ";
        }
        std::cout << std::endl;
    }

    //std::cout << "Press any key to exit... ";
    //std::cin.get();
    //std::cout << std::endl;

    resultCollections.clear();

    std::cout << std::endl << "END  " << floatingExp2Integer::Timer::current_time() << std::endl;;

    return 0;
}
