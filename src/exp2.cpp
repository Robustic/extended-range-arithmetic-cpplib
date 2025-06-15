#include <iostream>
#include <limits>
#include <functional>
#include <vector>
#include <immintrin.h>
#include <random>
#include <chrono>
#include <iomanip>
#include <cstring>
#include <dvec.h>
#include <immintrin.h>
#include "Timer.h"
#include "Float64ExtendedExp.h"

const uint32_t seed_val = 1337;

void InitializeRandomNumbers(std::vector<double>& vec) {
    std::mt19937 rng;
    rng.seed(seed_val);
    std::uniform_real_distribution<double> unif(std::log2(1e-11), std::log2(1e-10));
    for (unsigned int i = 0; i < vec.size(); i++) {
        float a_random_doable = unif(rng);
        vec[i] = std::exp2(a_random_doable);
    }
}

void DoubleToLogValues(const std::vector<double>& doubleValues, std::vector<double>& logValues) {
    for (unsigned int i = 0; i < logValues.size(); i++) {
        logValues[i] = std::log2(doubleValues[i]);
    }
}

void DoubleToFloat64ExtendedExpValues(const std::vector<double>& doubleValues,
    std::vector<double>& extended_doubles) {
    floatingExp2Integer::Float64ExtendedExp encoderi;
    encoderi.doubleToFloat64ExtendedExp((unsigned int)doubleValues.size(), doubleValues.data(), extended_doubles.data());
}

/////////////////////////////////


/// Log2 sum trick and (sequential and array)

int64_t calculate_sequential_sum_log2(const std::vector<double>& values, double& result) {
    floatingExp2Integer::Timer timer;
    double sum = values[0];
    for (unsigned int i = 1; i < values.size(); i++) {
        double value = values[i];
        double max = std::max(sum, value);
        sum = max + std::log2(std::exp2(sum - max) + std::exp2(value - max));
    }
    timer.stop();
    result = std::exp2(sum);
    return timer.time();
}

const unsigned int SIZE = 1024;

inline double calculate_array_sum_log2_split(const std::vector<double>& values, unsigned int start) {
    double sum = 0.0;
    double max = *std::max_element(values.begin() + start, values.begin() + start + SIZE);

    __m256d vmax = _mm256_set1_pd(max);
    __m256d va = _mm256_loadu_pd(&values[start]);
    va = va - vmax;
    __m256d vsum = _mm256_exp2_pd(va);

    for (unsigned int i = start + 4; i + 3 < start + SIZE; i += 4) {
        __m256d vb = _mm256_loadu_pd(&values[i]);
        vb = vb - vmax;
        __m256d vexp = _mm256_exp2_pd(vb);
        vsum = vsum + vexp;
    }

    for (unsigned int k = 0; k < 4; k++) {
        sum += vsum[k];
    }
    double result = max + std::log2(sum);

    return result;
}

inline double calculate_array_sum_log2_rest(const std::vector<double>& values, unsigned int start) {
    unsigned int i = 0;
    double sum = 0.0;
    double max = *std::max_element(values.begin() + start, values.end());
    if (8 <= values.size() - start) {
        __m256d vmax = _mm256_set1_pd(max);
        __m256d va = _mm256_loadu_pd(&values[0]);
        va = va - vmax;
        __m256d vsum = _mm256_exp2_pd(va);
        for (i = 4 + start; i + 3 < values.size(); i += 4) {
            __m256d vb = _mm256_loadu_pd(&values[i]);
            vb = vb - vmax;
            __m256d vexp = _mm256_exp2_pd(vb);
            vsum = vsum + vexp;
        }

        for (unsigned int k = 0; k < 4; k++) {
            sum += vsum[k];
        }
    }

    for (i = i; i < values.size(); i++) {
        sum += std::exp2(values[i] - max);
    }
    double result = max + std::log2(sum);

    return result;
}

int64_t calculate_array_sum_log2_smart(const std::vector<double>& values, double& result) {
    floatingExp2Integer::Timer timer;

    unsigned int n_count = ((values.size() - 1) / SIZE) + 1;

    std::vector<double> collected_values(n_count);

    for (unsigned int i = 0; i < n_count - 1; i++) {
        collected_values[i] = calculate_array_sum_log2_split(values, i * SIZE);
    }
    collected_values[n_count - 1] = calculate_array_sum_log2_rest(values, (n_count - 1) * SIZE);

    result = std::exp2(calculate_array_sum_log2_rest(collected_values, 0));
    timer.stop();
    return timer.time();
}

int64_t calculate_array_sum_log2(const std::vector<double>& values, double& result) {
    floatingExp2Integer::Timer timer;
    unsigned int i = 0;
    double sum = 0.0;
    double max = *std::max_element(values.begin(), values.end());
    if (8 <= values.size()) {
        __m256d vmax = _mm256_set1_pd(max);
        __m256d va = _mm256_loadu_pd(&values[0]);
        va = va - vmax;
        __m256d vsum = _mm256_exp2_pd(va);
        for (i = 4; i + 3 < values.size(); i += 4) {
            __m256d vb = _mm256_loadu_pd(&values[i]);
            vb = vb - vmax;
            __m256d vexp = _mm256_exp2_pd(vb);
            vsum = vsum + vexp;
        }

        for (unsigned int k = 0; k < 4; k++) {
            sum += vsum[k];
        }
    }

    for (i = i; i < values.size(); i++) {
        sum += std::exp2(values[i] - max);
    }
    result = std::exp2(max + std::log2(sum));
    timer.stop();
    return timer.time();
}

////////////

int64_t calculate_array_sum_Float64ExtendedExp(unsigned int n, double* values, double& result) {
    floatingExp2Integer::Timer timer;
    floatingExp2Integer::Float64ExtendedExp collector(1.0);
    double encoded = collector.sumFloat64ExtendedExp(n, values);
    double sum = collector.decodeFloat64ExtendedExp(encoded);
    timer.stop();
    result = sum;
    return timer.time();
}

int64_t calculate_sequential_sum_Float64ExtendedExp(unsigned int n, double* values, double& result) {
    floatingExp2Integer::Timer timer;
    floatingExp2Integer::Float64ExtendedExp collector(1.0);
    collector.encodedToFloat64ExtendedExp(values[0]);
    for (unsigned int i = 1; i < n; i++) {
        unsigned k = collector.encoded > 100.001 ? i - 1 : i;
        volatile double value = values[k];
        collector += value;
    }
    double sum = collector.asDouble();
    timer.stop();
    result = sum;
    return timer.time();
}

//////////////

int main() {
    unsigned int n_current = 100000000;
    std::vector<double> double_values(n_current);
    InitializeRandomNumbers(double_values);

    std::vector<double> values_converted(double_values.size());
    DoubleToLogValues(double_values, values_converted);

    std::vector<double> values_extended_exp(double_values.size());
    DoubleToFloat64ExtendedExpValues(double_values, values_extended_exp);

    std::cout << std::setprecision(std::numeric_limits<double>::digits10 + 1);

    double result = 0.0;
    int64_t time1 = calculate_sequential_sum_log2(values_converted, result);
    std::cout << "calculate_sequential_sum_log2:                 " << time1 << "  " << result << std::endl;

    int64_t time2 = calculate_array_sum_log2(values_converted, result);
    std::cout << "calculate_avg_array_sum_log2:                   " << time2 << "  " << result << std::endl;

    int64_t time3 = calculate_array_sum_log2_smart(values_converted, result);
    std::cout << "calculate_array_sum_log2_smart:                 " << time3 << "  " << result << std::endl;

    int64_t time4 = calculate_sequential_sum_Float64ExtendedExp(values_extended_exp.size(), values_extended_exp.data(), result);
    std::cout << "calculate_sequential_sum_Float64ExtendedExp:   " << time4 << "  " << result << std::endl;

    int64_t time5 = calculate_array_sum_Float64ExtendedExp(values_extended_exp.size(), values_extended_exp.data(), result);
    std::cout << "calculate_array_sum_Float64ExtendedExp:          " << time5 << "  " << result << std::endl;

    return 0;
}

