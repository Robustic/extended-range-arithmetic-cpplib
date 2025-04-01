#include <iostream>
#include <limits>
#include <functional>
#include <vector>
#include <immintrin.h>
#include <random>
#include <chrono>
#include <iomanip>
#include <cstring>
#include "Timer.h"
#include "Dbl.h"
#include "Dbl2.h"
#include "Dbl3.h"
#include "./Int64PosExp2Int64/Int64PosExp2Int64.h"
#include "./Float64PosExp2Int64/Float64PosExp2Int64.h"
#include "./Float64Exp2Int64/Float64Exp2Int64.h"
#include "./Float64ExtendedExp/Float64ExtendedExp.h"
#include "./Fukushima/Fukushima.h"

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

void DoubleToDblValues(const std::vector<double>& doubleValues, std::vector<floatingExp2Integer::Dbl>& dblValues) {
    for (unsigned int i = 0; i < dblValues.size(); i++) {
        dblValues[i] += doubleValues[i];
    }
}

void DoubleToDbl2Values(const std::vector<double>& doubleValues, std::vector<floatingExp2Integer::Dbl2>& dbl2Values) {
    for (unsigned int i = 0; i < dbl2Values.size(); i++) {
        dbl2Values[i] += doubleValues[i];
    }
}

void DoubleToLogValues(const std::vector<double>& doubleValues, std::vector<double>& logValues) {
    for (unsigned int i = 0; i < logValues.size(); i++) {
        logValues[i] = std::log2(doubleValues[i]);
    }
}

void DoubleToFukushimaValues(const std::vector<double>& doubleValues,
    std::vector<floatingExp2Integer::Fukushima>& fukushimaValues) {
    for (unsigned int i = 0; i < fukushimaValues.size(); i++) {
        fukushimaValues[i].doubleToFukushima(doubleValues[i]);
    }
}

/////////////////////////////////

/// double sum and multiply

std::int64_t calculate_array_sum_dbl(const std::vector<double>& values, double& result) {
    floatingExp2Integer::Timer timer;
    unsigned int i = 0;
    double sum = 0.0; 
    if (16 <= values.size()) {
        __m512d vsum = _mm512_loadu_pd(&values[0]);
        for (i = 8; i + 7 < values.size(); i += 8) {
            __m512d vb = _mm512_loadu_pd(&values[i]);
            vsum = vsum + vb;
        }

        for (unsigned int k = 0; k < 8; k++) {
            sum += vsum[k];
        }
    }

    for (i = i; i < values.size(); i++) {
        sum += values[i];
    }
    timer.stop();
    result = sum;
    return timer.time();
}

std::int64_t calculate_avg_array_sum_dbl(std::string& header, unsigned int n_rounds, const std::vector<double>& values, double& result)
{
    header = "double_array_sum:";

    std::int64_t time_sum = 0.0;
    for (unsigned int i = 0; i < n_rounds; i++) {
        time_sum += calculate_array_sum_dbl(values, result);
    }
    return time_sum / n_rounds;
}

std::int64_t calculate_array_multiply_dbl(const std::vector<double>& values, double& result) {
    floatingExp2Integer::Timer timer;
    unsigned int i = 0;
    double multiplied = 1.0;
    if (16 <= values.size()) {
        __m512d vmul = _mm512_set1_pd(1.0);
        for (i = 0; i + 7 < values.size(); i += 8) {
            __m512d vb = _mm512_loadu_pd(&values[i]);
            vmul = vmul * vb;
        }

        for (unsigned int k = 0; k < 8; k++) {
            multiplied *= vmul[k];
        }
    }

    for (i = i; i < values.size(); i++) {
        multiplied *= values[i];
    }
    timer.stop();
    result = multiplied;
    return timer.time();
}

std::int64_t calculate_avg_array_multiply_dbl(std::string& header, unsigned int n_rounds, const std::vector<double>& values, double& result)
{
    header = "double_array_multiply:";

    std::int64_t time_sum = 0.0;
    for (unsigned int i = 0; i < n_rounds; i++) {
        time_sum += calculate_array_multiply_dbl(values, result);
    }
    return time_sum / n_rounds;
}

/// Log2 sum trick and (sequential and array)

std::int64_t calculate_sequential_sum_log2(const std::vector<double>& values, double& result) {
    floatingExp2Integer::Timer timer;
    double sum = values[0];
    for (unsigned int i = 1; i < values.size(); i++) {
        double value = values[i];
        if (sum > value) {
            sum = sum + std::log2(1 + std::exp2(value - sum));
        }
        else {
            sum = value + std::log2(std::exp2(sum - value) + 1);
        }
    }
    timer.stop();
    result = std::exp2(sum);
    return timer.time();
}

std::int64_t calculate_avg_sequential_sum_log2(std::string& header, unsigned int n_rounds, const std::vector<double>& values, double& result)
{
    header = "log2_sequential_sum:";

    std::vector<double> values_converted(values.size());
    DoubleToLogValues(values, values_converted);

    std::int64_t time_sum = 0.0;
    for (unsigned int i = 0; i < n_rounds; i++) {
        time_sum += calculate_sequential_sum_log2(values_converted, result);
    }
    return time_sum / n_rounds;
}

std::int64_t calculate_array_sum_log2(const std::vector<double>& values, double& result) {
    floatingExp2Integer::Timer timer;
    unsigned int i = 0;
    double sum = 0.0;
    double max = *std::max_element(values.begin(), values.end());
    for (i = 0; i < values.size(); i++) {
        sum += std::exp2(values[i] - max);
    }
    timer.stop();
    result = max + std::log2(sum);
    return timer.time();
}

std::int64_t calculate_avg_array_sum_log2(std::string& header, unsigned int n_rounds, const std::vector<double>& values, double& result)
{
    header = "log2_array_sum:";

    std::vector<double> values_converted(values.size());
    DoubleToLogValues(values, values_converted);

    std::int64_t time_sum = 0.0;
    for (unsigned int i = 0; i < n_rounds; i++) {
        time_sum += calculate_array_sum_log2(values_converted, result);
    }
    return time_sum / n_rounds;
}

std::int64_t calculate_array_multiply_log2(const std::vector<double>& values, double& result) {
    floatingExp2Integer::Timer timer;
    unsigned int i = 0;
    double sum = 0.0;
    for (i = 0; i < values.size(); i++) {
        sum += values[i];
    }
    timer.stop();
    result = sum;
    return timer.time();
}

std::int64_t calculate_avg_array_multiply_log2(std::string& header, unsigned int n_rounds, const std::vector<double>& values, double& result)
{
    header = "log2_array_multiply:";

    std::vector<double> values_converted(values.size());
    DoubleToLogValues(values, values_converted);

    std::int64_t time_sum = 0.0;
    for (unsigned int i = 0; i < n_rounds; i++) {
        time_sum += calculate_array_multiply_log2(values_converted, result);
    }
    return time_sum / n_rounds;
}

/// Dbl1 and Dbl2 sum

std::int64_t calculate_sequential_sum_Dbl1(const std::vector<floatingExp2Integer::Dbl>& values, double& result) {
    floatingExp2Integer::Timer timer;
    floatingExp2Integer::Dbl sum = 0.0;
    for (unsigned int i = 0; i < values.size(); i++) {
        sum += values[i];
    }
    timer.stop();
    result = sum.asDouble();
    return timer.time();
}

std::int64_t calculate_avg_sequential_sum_Dbl1(std::string& header, unsigned int n_rounds, const std::vector<double>& values, double& result)
{
    header = "Dbl1_sequential_sum:";

    std::vector<floatingExp2Integer::Dbl> values_converted(values.size());
    DoubleToDblValues(values, values_converted);

    std::int64_t time_sum = 0.0;
    for (unsigned int i = 0; i < n_rounds; i++) {
        time_sum += calculate_sequential_sum_Dbl1(values_converted, result);
    }
    return time_sum / n_rounds;
}

std::int64_t calculate_sequential_sum_Dbl2(const std::vector<floatingExp2Integer::Dbl2>& values, double& result) {
    floatingExp2Integer::Timer timer;
    floatingExp2Integer::Dbl2 sum = 0.0;
    for (unsigned int i = 0; i < values.size(); i++) {
        sum += values[i];
    }
    timer.stop();
    result = sum.asDouble();
    return timer.time();
}

std::int64_t calculate_avg_sequential_sum_Dbl2(std::string& header, unsigned int n_rounds, const std::vector<double>& values, double& result)
{
    header = "Dbl2_sequential_sum:";

    std::vector<floatingExp2Integer::Dbl2> values_converted(values.size());
    DoubleToDbl2Values(values, values_converted);

    std::int64_t time_sum = 0.0;
    for (unsigned int i = 0; i < n_rounds; i++) {
        time_sum += calculate_sequential_sum_Dbl2(values_converted, result);
    }
    return time_sum / n_rounds;
}

std::int64_t calculate_array_sum_Dbl1(const std::vector<floatingExp2Integer::Dbl>& values, double& result) {
    floatingExp2Integer::Timer timer;
    floatingExp2Integer::Dbl sum = 0.0;
    sum.sumDbl(values);
    timer.stop();
    result = sum.asDouble();
    return timer.time();
}

std::int64_t calculate_avg_array_sum_Dbl1(std::string& header, unsigned int n_rounds, const std::vector<double>& values, double& result)
{
    header = "Dbl1_array_sum:";

    std::vector<floatingExp2Integer::Dbl> values_converted(values.size());
    DoubleToDblValues(values, values_converted);

    std::int64_t time_sum = 0.0;
    for (unsigned int i = 0; i < n_rounds; i++) {
        time_sum += calculate_array_sum_Dbl1(values_converted, result);
    }
    return time_sum / n_rounds;
}

std::int64_t calculate_array_sum_Dbl2(const std::vector<floatingExp2Integer::Dbl2>& values, double& result) {
    floatingExp2Integer::Timer timer;
    floatingExp2Integer::Dbl2 sum = 0.0;
    sum.sumDbl2(values);
    timer.stop();
    result = sum.asDouble();
    return timer.time();
}

std::int64_t calculate_avg_array_sum_Dbl2(std::string& header, unsigned int n_rounds, const std::vector<double>& values, double& result)
{
    header = "Dbl2_array_sum:";

    std::vector<floatingExp2Integer::Dbl2> values_converted(values.size());
    DoubleToDbl2Values(values, values_converted);

    std::int64_t time_sum = 0.0;
    for (unsigned int i = 0; i < n_rounds; i++) {
        time_sum += calculate_array_sum_Dbl2(values_converted, result);
    }
    return time_sum / n_rounds;
}

/// Dbl1 and Dbl2 multiply

std::int64_t calculate_sequential_multiply_Dbl1(const std::vector<floatingExp2Integer::Dbl>& values, double& result) {
    floatingExp2Integer::Timer timer;
    floatingExp2Integer::Dbl multiplied = 1.0;
    for (unsigned int i = 0; i < values.size(); i++) {
        multiplied *= values[i];
    }
    timer.stop();
    result = multiplied.asDouble();
    return timer.time();
}

std::int64_t calculate_avg_sequential_multiply_Dbl1(std::string& header, unsigned int n_rounds, const std::vector<double>& values, double& result)
{
    header = "Dbl1_sequential_multiply:";

    std::vector<floatingExp2Integer::Dbl> values_converted(values.size());
    DoubleToDblValues(values, values_converted);

    std::int64_t time_sum = 0.0;
    for (unsigned int i = 0; i < n_rounds; i++) {
        time_sum += calculate_sequential_multiply_Dbl1(values_converted, result);
    }
    return time_sum / n_rounds;
}

std::int64_t calculate_sequential_multiply_Dbl2(const std::vector<floatingExp2Integer::Dbl2>& values, double& result) {
    floatingExp2Integer::Timer timer;
    floatingExp2Integer::Dbl2 multiplied = 1.0;
    for (unsigned int i = 0; i < values.size(); i++) {
        multiplied *= values[i];
    }
    timer.stop();
    result = multiplied.asDouble();
    return timer.time();
}

std::int64_t calculate_avg_sequential_multiply_Dbl2(std::string& header, unsigned int n_rounds, const std::vector<double>& values, double& result)
{
    header = "Dbl2_sequential_multiply:";

    std::vector<floatingExp2Integer::Dbl2> values_converted(values.size());
    DoubleToDbl2Values(values, values_converted);

    std::int64_t time_sum = 0.0;
    for (unsigned int i = 0; i < n_rounds; i++) {
        time_sum += calculate_sequential_multiply_Dbl2(values_converted, result);
    }
    return time_sum / n_rounds;
}

std::int64_t calculate_array_multiply_Dbl1(const std::vector<floatingExp2Integer::Dbl>& values, double& result) {
    floatingExp2Integer::Timer timer;
    floatingExp2Integer::Dbl multiplied = 1.0;
    multiplied.multiplyDbl(values);
    timer.stop();
    result = multiplied.asDouble();
    return timer.time();
}

std::int64_t calculate_avg_array_multiply_Dbl1(std::string& header, unsigned int n_rounds, const std::vector<double>& values, double& result)
{
    header = "Dbl1_array_multiply:";

    std::vector<floatingExp2Integer::Dbl> values_converted(values.size());
    DoubleToDblValues(values, values_converted);

    std::int64_t time_sum = 0.0;
    for (unsigned int i = 0; i < n_rounds; i++) {
        time_sum += calculate_array_multiply_Dbl1(values_converted, result);
    }
    return time_sum / n_rounds;
}

std::int64_t calculate_array_multiply_Dbl2(const std::vector<floatingExp2Integer::Dbl2>& values, double& result) {
    floatingExp2Integer::Timer timer;
    floatingExp2Integer::Dbl2 multiplied = 1.0;
    multiplied.multiplyDbl2(values);
    timer.stop();
    result = multiplied.asDouble();
    return timer.time();
}

std::int64_t calculate_avg_array_multiply_Dbl2(std::string& header, unsigned int n_rounds, const std::vector<double>& values, double& result)
{
    header = "Dbl2_array_multiply:";

    std::vector<floatingExp2Integer::Dbl2> values_converted(values.size());
    DoubleToDbl2Values(values, values_converted);

    std::int64_t time_sum = 0.0;
    for (unsigned int i = 0; i < n_rounds; i++) {
        time_sum += calculate_array_multiply_Dbl2(values_converted, result);
    }
    return time_sum / n_rounds;
}

/// 

std::int64_t calculate_array_sum_Float64PosExp2Int64(const std::vector<floatingExp2Integer::Float64PosExp2Int64>& values, double& result) {
    floatingExp2Integer::Timer timer;
    floatingExp2Integer::Float64PosExp2Int64 sum(values);
    timer.stop();
    result = sum.asDouble();
    return timer.time();
}

void DoubleToFloat64PosExp2Int64Values(const std::vector<double>& doubleValues, 
    std::vector<floatingExp2Integer::Float64PosExp2Int64>& float64PosExp2Int64Values) {
    for (unsigned int i = 0; i < float64PosExp2Int64Values.size(); i++) {
        float64PosExp2Int64Values[i].doubleToFloat64PosExp2Int64(doubleValues[i]);
    }
}

std::int64_t calculate_avg_array_sum_Float64PosExp2Int64(std::string& header, unsigned int n_rounds, const std::vector<double>& values, double& result)
{
    header = "Float64PosExp2Int64_array_sum:";

    std::vector<floatingExp2Integer::Float64PosExp2Int64> values_converted(values.size());
    DoubleToFloat64PosExp2Int64Values(values, values_converted);

    std::int64_t time_sum = 0.0;
    for (unsigned int i = 0; i < n_rounds; i++) {
        time_sum += calculate_array_sum_Float64PosExp2Int64(values_converted, result);
    }
    return time_sum / n_rounds;
}

std::int64_t calculate_sequential_sum_Fukushima(const std::vector<floatingExp2Integer::Fukushima>& values, double& result) {
    floatingExp2Integer::Timer timer;
    floatingExp2Integer::Fukushima sum = 0.0;
    for (unsigned int i = 0; i < values.size(); i++) {
        unsigned k = sum.exp > 1LL ? (i == 0 ? i : i - 1) : i;
        sum += values[k];
    }
    timer.stop();
    result = sum.asDouble();
    return timer.time();
}

std::int64_t calculate_avg_sequential_sum_Fukushima(std::string& header, unsigned int n_rounds, const std::vector<double>& values, double& result)
{
    header = "Fukushima_sequential_sum:";

    std::vector<floatingExp2Integer::Fukushima> values_converted(values.size());
    DoubleToFukushimaValues(values, values_converted);

    std::int64_t time_sum = 0.0;
    for (unsigned int i = 0; i < n_rounds; i++) {
        time_sum += calculate_sequential_sum_Fukushima(values_converted, result);
    }
    return time_sum / n_rounds;
}

std::int64_t calculate_array_sum_Fukushima(const std::vector<floatingExp2Integer::Fukushima>& values, double& result) {
    floatingExp2Integer::Timer timer;
    floatingExp2Integer::Fukushima collector;
    collector.sum(values);
    timer.stop();
    result = collector.asDouble();
    return timer.time();
}

std::int64_t calculate_avg_array_sum_Fukushima(std::string& header, unsigned int n_rounds, const std::vector<double>& values, double& result)
{
    header = "Fukushima_array_sum:";

    std::vector<floatingExp2Integer::Fukushima> values_converted(values.size());
    DoubleToFukushimaValues(values, values_converted);

    std::int64_t time_sum = 0.0;
    for (unsigned int i = 0; i < n_rounds; i++) {
        time_sum += calculate_array_sum_Fukushima(values_converted, result);
    }
    return time_sum / n_rounds;
}

std::int64_t calculate_sequential_multiply_Fukushima(const std::vector<floatingExp2Integer::Fukushima>& values, double& result) {
    floatingExp2Integer::Timer timer;
    floatingExp2Integer::Fukushima res = 1.0;
    for (unsigned int i = 0; i < values.size(); i++) {
        unsigned k = res.exp > 1LL ? (i == 0 ? i : i - 1) : i;
        res *= values[k];
    }
    timer.stop();
    result = res.fukushimaToLog2();
    return timer.time();
}

std::int64_t calculate_avg_sequential_multiply_Fukushima(std::string& header, unsigned int n_rounds, const std::vector<double>& values, double& result)
{
    header = "Fukushima_sequential_multiply:";

    std::vector<floatingExp2Integer::Fukushima> values_converted(values.size());
    DoubleToFukushimaValues(values, values_converted);

    std::int64_t time_sum = 0.0;
    for (unsigned int i = 0; i < n_rounds; i++) {
        time_sum += calculate_sequential_multiply_Fukushima(values_converted, result);
    }
    return time_sum / n_rounds;
}

std::int64_t calculate_array_multiply_Fukushima(const std::vector<floatingExp2Integer::Fukushima>& values, double& result) {
    floatingExp2Integer::Timer timer;
    floatingExp2Integer::Fukushima collector;
    collector.multiply(values);
    timer.stop();
    result = collector.fukushimaToLog2();
    return timer.time();
}

std::int64_t calculate_avg_array_multiply_Fukushima(std::string& header, unsigned int n_rounds, const std::vector<double>& values, double& result)
{
    header = "Fukushima_array_multiply:";

    std::vector<floatingExp2Integer::Fukushima> values_converted(values.size());
    DoubleToFukushimaValues(values, values_converted);

    std::int64_t time_sum = 0.0;
    for (unsigned int i = 0; i < n_rounds; i++) {
        time_sum += calculate_array_multiply_Fukushima(values_converted, result);
    }
    return time_sum / n_rounds;
}

std::int64_t calculate_array_sum_Float64ExtendedExp(unsigned int n, double* values, double& result) {
    floatingExp2Integer::Timer timer;
    floatingExp2Integer::Float64ExtendedExp collector(1.0);
    double encoded = collector.sumFloat64ExtendedExp(n, values);
    double sum = collector.decodeFloat64ExtendedExp(encoded);
    timer.stop();
    result = sum;
    return timer.time();
}

void DoubleToFloat64ExtendedExpValues(const std::vector<double>& doubleValues, 
    std::vector<floatingExp2Integer::Float64ExtendedExp>& fukushimaValues) {
    for (unsigned int i = 0; i < fukushimaValues.size(); i++) {
        fukushimaValues[i].doubleToFloat64ExtendedExp(doubleValues[i]);
    }
}

std::int64_t calculate_avg_array_sum_Float64ExtendedExp(std::string& header, unsigned int n_rounds, const std::vector<double>& values, double& result)
{
    header = "Float64ExtendedExp_array_sum:";

    floatingExp2Integer::Float64ExtendedExp extendedExp(1.0);

    double* values_converted = new double[values.size()];
    extendedExp.doubleToFloat64ExtendedExp(values.size(), values.data(), values_converted);

    std::int64_t time_sum = 0.0;
    for (unsigned int i = 0; i < n_rounds; i++) {
        time_sum += calculate_array_sum_Float64ExtendedExp(values.size(), values_converted, result);
    }

    delete[] values_converted;
    return time_sum / n_rounds;
}

std::int64_t calculate_sequential_sum_Float64ExtendedExp(unsigned int n, double* values, double& result) {
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

std::int64_t calculate_avg_sequential_sum_Float64ExtendedExp(std::string& header, unsigned int n_rounds, const std::vector<double>& values, double& result)
{
    header = "Float64ExtendedExp_sequential_sum:";

    floatingExp2Integer::Float64ExtendedExp extendedExp(1.0);

    double* values_converted = new double[values.size()];
    extendedExp.doubleToFloat64ExtendedExp(values.size(), values.data(), values_converted);

    std::int64_t time_sum = 0.0;
    for (unsigned int i = 0; i < n_rounds; i++) {
        time_sum += calculate_sequential_sum_Float64ExtendedExp(values.size(), values_converted, result);
    }

    delete[] values_converted;
    return time_sum / n_rounds;
}

int main() {
    constexpr unsigned int n[] = { 1000, 3000, 10000, 30000, 100000, 300000, 1000000, 3000000, 10000000, 30000000, 100000000 };
    constexpr unsigned int n_rounds[] = { 100000, 30000, 10000, 3000, 1000, 300, 100, 30, 10, 3, 1 };
    constexpr unsigned int n_count = sizeof(n) / sizeof(n[0]);

    std::vector<std::function<std::int64_t(std::string&, int, const std::vector<double>&, double&)>> functions;
    //functions.push_back(calculate_avg_array_sum_dbl);
    //functions.push_back(calculate_avg_array_multiply_dbl);
    //functions.push_back(calculate_avg_sequential_sum_Dbl1);
    //functions.push_back(calculate_avg_sequential_sum_Dbl2);
    //functions.push_back(calculate_avg_array_sum_Dbl1);
    //functions.push_back(calculate_avg_array_sum_Dbl2);
    //functions.push_back(calculate_avg_sequential_multiply_Dbl1);
    //functions.push_back(calculate_avg_sequential_multiply_Dbl2);
    //functions.push_back(calculate_avg_array_multiply_Dbl1);
    //functions.push_back(calculate_avg_array_multiply_Dbl2);
    //functions.push_back(calculate_avg_sequential_sum_log2);
    //functions.push_back(calculate_avg_array_sum_log2);
    functions.push_back(calculate_avg_array_multiply_log2);
    functions.push_back(calculate_avg_sequential_sum_Fukushima);
    functions.push_back(calculate_avg_array_sum_Fukushima);
    functions.push_back(calculate_avg_sequential_multiply_Fukushima);
    functions.push_back(calculate_avg_array_multiply_Fukushima);
    //functions.push_back(calculate_avg_array_sum_Float64PosExp2Int64);
    //functions.push_back(calculate_avg_array_sum_Float64ExtendedExp);
    //functions.push_back(calculate_avg_sequential_sum_Float64ExtendedExp);

    int functions_count = functions.size();

    std::vector<std::string> all_headers(functions_count);
    std::vector<std::array<double, n_count>> all_results(functions_count);
    std::vector<std::array<std::int64_t, n_count>> all_times(functions_count);

    for (unsigned int i = 0; i < n_count; i++) {
        int n_current = n[i];
        int n_rounds_current = n_rounds[i];

        std::vector<double> double_values(n_current);
        InitializeRandomNumbers(double_values);

        for (unsigned int f = 0; f < functions.size(); f++) {
            std::function<std::int64_t(std::string&, int, const std::vector<double>&, double&)> function = functions[f];
            std::string header;
            double result;
            std::int64_t time = function(header, n_rounds_current, double_values, result);
            all_headers[f] = header;
            all_results[f][i] = result;
            all_times[f][i] = time;
        }
    }

    std::cout << std::setprecision(std::numeric_limits<double>::digits10 + 1);

    std::cout << std::endl << "Result" << std::endl << "n ";
    for (unsigned int i = 0; i < n_count; i++) {
        std::cout << n[i] << " ";
    }

    std::cout << std::endl;
    for (unsigned int f = 0; f < functions.size(); f++) {
        std::cout << all_headers[f] << " ";
        for (unsigned int i = 0; i < n_count; i++) {
            std::cout << all_results[f][i] << " ";
        }
        std::cout << std::endl;
    }

    std::cout << std::endl << "Time" << std::endl << "n ";
    for (unsigned int i = 0; i < n_count; i++) {
        std::cout << n[i] << " ";
    }

    std::cout << std::endl;
    for (unsigned int f = 0; f < functions.size(); f++) {
        std::cout << all_headers[f] << " ";
        for (unsigned int i = 0; i < n_count; i++) {
            std::cout << all_times[f][i] << " ";
        }
        std::cout << std::endl;
    }

    //std::cout << "Press any key to exit... ";
    //std::cin.get();
    //std::cout << std::endl;

    return 0;
}
