#include <cstdint>
#include <cmath>
#include <vector>
#include <cstring>
#include "Float64PosExp2Int64.h"

typedef double double4_t __attribute__((vector_size(4 * sizeof(double))));
typedef uint64_t uint64_4_t __attribute__((vector_size(4 * sizeof(uint64_t))));
typedef int64_t int64_4_t __attribute__((vector_size(4 * sizeof(int64_t))));

constexpr int64_4_t int64_4_1023 { 1023LL, 1023LL, 1023LL, 1023LL };
constexpr double4_t double4_0 { 0.0, 0.0, 0.0, 0.0 };

namespace floatingExp2Integer
{
    Float64PosExp2Int64::Float64PosExp2Int64() {
        scnfcnd = 1.0;
        exp = 0;
        this->scale();
    }

    Float64PosExp2Int64::Float64PosExp2Int64(double dbl) {
        scnfcnd = dbl;
        exp = 0;
        this->scale();
    }

    Float64PosExp2Int64::Float64PosExp2Int64(double sicnificand, int64_t exponent) {
        scnfcnd = sicnificand;
        exp = exponent;
        this->scale();
    }

    void Float64PosExp2Int64::sum(const std::vector<floatingExp2Integer::Float64PosExp2Int64>& vector) {
        const unsigned int parallel_count = 4;

        double scnfcndSum[parallel_count];
        int64_t expSum[parallel_count];

        for (unsigned int k = 0; k < parallel_count; k++) {
            scnfcndSum[k] = vector[k].scnfcnd;
            expSum[k] = vector[k].exp;
        }

        double scnfcndCurrent[parallel_count];
        int64_t expCurrent[parallel_count];

        for (unsigned int i = parallel_count; i + (parallel_count - 1) < vector.size(); i += parallel_count) {
            for (unsigned int k = 0; k < parallel_count; k++) {
                scnfcndCurrent[k] = vector[i + k].scnfcnd;
                expCurrent[k] = vector[i + k].exp;
            }

            for (unsigned int k = 0; k < parallel_count; k++) {
                int64_t exp_diff = expSum[k] - expCurrent[k];
                uint64_t* sgnfcnd_bits = reinterpret_cast<uint64_t*>(&scnfcndSum[k]);
                uint64_t* sgnfcnd_bits_z = reinterpret_cast<uint64_t*>(&scnfcndCurrent[k]);

                if (exp_diff >= 0) {
                    if (exp_diff > 64) {
                        //
                    }
                    else {
                        *sgnfcnd_bits_z -= exp_diff << 52;
                        scnfcndSum[k] += scnfcndCurrent[k];
                    }
                }
                else {
                    if (exp_diff < -64) {
                        scnfcndSum[k] = scnfcndCurrent[k];
                        expSum[k] = expCurrent[k];
                    }
                    else {
                        expSum[k] = expCurrent[k];
                        *sgnfcnd_bits -= (-exp_diff) << 52;
                        scnfcndSum[k] += scnfcndCurrent[k];
                    }
                }

                if (scnfcndSum[k] >= 2.0) {
                    scnfcndSum[k] *= 0.5;
                    expSum[k]++;
                }
            }
        }

        floatingExp2Integer::Float64PosExp2Int64 sum(scnfcndSum[0], expSum[0]);

        for (unsigned int k = 1; k < parallel_count; k++) {
            floatingExp2Integer::Float64PosExp2Int64 current(scnfcndSum[k], expSum[k]);
            sum += current;
        }

        scnfcnd = sum.scnfcnd;
        exp = sum.exp;

        this->scale();
    }

    void Float64PosExp2Int64::double_to(double dbl) {
        scnfcnd = dbl;
        exp = 0;
        this->scale();
    }

    void Float64PosExp2Int64::log2_to(double log2) {
        exp = (int64_t)log2;
        scnfcnd = std::exp2(log2 - exp);
        this->scale();
    }

    double Float64PosExp2Int64::as_log2() const {
        return std::log2(scnfcnd) + exp;
    }

    double Float64PosExp2Int64::as_double() const {
        return scnfcnd * std::pow(2.0, exp);
    }

    Float64PosExp2Int64& Float64PosExp2Int64::operator+=(Float64PosExp2Int64 z) {
        int64_t exp_diff = exp - z.exp;
        uint64_t* sgnfcnd_bits = reinterpret_cast<uint64_t*>(&scnfcnd);
        uint64_t* sgnfcnd_bits_z = reinterpret_cast<uint64_t*>(&z.scnfcnd);

        if (exp_diff >= 0) {
            if (exp_diff > 64) {
                return *this;
            }
            *sgnfcnd_bits_z -= exp_diff << 52;
            scnfcnd += z.scnfcnd;
        }
        else {
            if (exp_diff < -64) {
                scnfcnd = z.scnfcnd;
                exp = z.exp;
                return *this;
            }
            exp = z.exp;
            *sgnfcnd_bits -= (-exp_diff) << 52;
            scnfcnd += z.scnfcnd;
        }

        if (scnfcnd >= 2.0) {
            scnfcnd *= 0.5;
            exp++;
        }
        return *this;
    }

    void Float64PosExp2Int64::multiply(const std::vector<floatingExp2Integer::Float64PosExp2Int64>& vector) {
        const unsigned int parallel_count = 4;

        double scnfcndSum[parallel_count];
        int64_t expSum[parallel_count];

        for (unsigned int k = 0; k < parallel_count; k++) {
            scnfcndSum[k] = vector[k].scnfcnd;
            expSum[k] = vector[k].exp;
        }

        double scnfcndCurrent[parallel_count];
        int64_t expCurrent[parallel_count];

        for (unsigned int i = parallel_count; i + (parallel_count - 1) < vector.size(); i += parallel_count) {
            for (unsigned int k = 0; k < parallel_count; k++) {
                scnfcndCurrent[k] = vector[i + k].scnfcnd;
                expCurrent[k] = vector[i + k].exp;
            }

            for (unsigned int k = 0; k < parallel_count; k++) {
                scnfcndSum[k] *= scnfcndCurrent[k];
                expSum[k] += expCurrent[k];
                if (scnfcndSum[k] >= 2.0) {
                    scnfcndSum[k] *= 0.5;
                    expSum[k]++;
                }
            }
        }

        floatingExp2Integer::Float64PosExp2Int64 sum(scnfcndSum[0], expSum[0]);

        for (unsigned int k = 1; k < parallel_count; k++) {
            floatingExp2Integer::Float64PosExp2Int64 current(scnfcndSum[k], expSum[k]);
            sum *= current;
        }

        scnfcnd = sum.scnfcnd;
        exp = sum.exp;

        this->scale();
    }

    Float64PosExp2Int64& Float64PosExp2Int64::operator*=(Float64PosExp2Int64 z) {
        scnfcnd *= z.scnfcnd;
        exp += z.exp;
        if (scnfcnd >= 2.0) {
            scnfcnd *= 0.5;
            exp++;
        }
        return *this;
    }

    Float64PosExp2Int64& Float64PosExp2Int64::operator*=(double& dbl) {
        Float64PosExp2Int64 z(dbl);
        scnfcnd *= z.scnfcnd;
        exp += z.exp;
        this->checkRuleForScale();
        return *this;
    }

    inline void Float64PosExp2Int64::checkRuleForScale() {
        if (0x1p12 <= scnfcnd) {
            this->scale();
        }
    }

    inline void Float64PosExp2Int64::scale() {
        uint64_t* sgnfcnd_bits = reinterpret_cast<uint64_t*>(&scnfcnd);
        exp += ((*sgnfcnd_bits & 0x7FF0000000000000ull) >> 52) - 1023;
        *sgnfcnd_bits &= 0x800FFFFFFFFFFFFFull;
        *sgnfcnd_bits |= 0x3FF0000000000000ull;
        scnfcnd = *reinterpret_cast<double*>(sgnfcnd_bits);
    }

    Float64PosExp2Int64 operator+(Float64PosExp2Int64 a, const Float64PosExp2Int64 b) { return a+=b; }
    Float64PosExp2Int64 operator*(Float64PosExp2Int64 a, const Float64PosExp2Int64 b) { return a*=b; }
}

