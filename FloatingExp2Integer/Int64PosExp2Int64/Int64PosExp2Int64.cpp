#include <cstdint>
#include <cmath>
#include <vector>
#include "Int64PosExp2Int64.h"

namespace floatingExp2Integer
{
    Int64PosExp2Int64::Int64PosExp2Int64() {
        this->fromDouble(1.0);
        this->scale();
    }

    Int64PosExp2Int64::Int64PosExp2Int64(double dbl) {
        this->fromDouble(dbl);
        this->scale();
    }

    Int64PosExp2Int64::Int64PosExp2Int64(double significand, uint64_t exponent) {
        scnfcnd = significand;
        exp = exponent;
        this->scale();
    }

    void Int64PosExp2Int64::double_to(double dbl) {
        this->fromDouble(dbl);
        this->scale();
    }

    void Int64PosExp2Int64::log2_to(double log2) {
        int64_t exponent = (int64_t)log2;
        double dbl = std::exp2(log2 - exponent);
        this->fromDouble(dbl);
        exp += exponent;
        this->scale();
    }

    double Int64PosExp2Int64::as_log2() const {
        return std::log2(scnfcnd) + exp;
    }

    double Int64PosExp2Int64::as_double() const { 
        return (double)scnfcnd * std::pow(2.0, exp);
    }

    void Int64PosExp2Int64::sum(const std::vector<floatingExp2Integer::Int64PosExp2Int64>& vector) {
        const unsigned int parallel_count = 4;

        uint64_t scnfcndSum[parallel_count];
        int64_t expSum[parallel_count];

        for (unsigned int k = 0; k < parallel_count; k++) {
            scnfcndSum[k] = vector[k].scnfcnd;
            expSum[k] = vector[k].exp;
        }

        uint64_t scnfcndCurrent[parallel_count];
        int64_t expCurrent[parallel_count];

        for (unsigned int i = parallel_count; i + (parallel_count - 1) < vector.size(); i += parallel_count) {
            for (unsigned int k = 0; k < parallel_count; k++) {
                scnfcndCurrent[k] = vector[i + k].scnfcnd;
                expCurrent[k] = vector[i + k].exp;
            }

            for (unsigned int k = 0; k < parallel_count; k++) {
                int64_t exp_diff = expSum[k] - expCurrent[k];

                if (exp_diff >= 0) {
                    if (exp_diff > 63) {
                        //
                    }
                    else {
                        scnfcndSum[k] += scnfcndCurrent[k] >> exp_diff;
                    }
                }
                else {
                    if (exp_diff < -63) {
                        scnfcndSum[k] = scnfcndCurrent[k];
                        expSum[k] = expCurrent[k];
                    }
                    else {
                        expSum[k] = expCurrent[k];
                        scnfcndSum[k] = scnfcndSum[k] >> (-exp_diff);
                        scnfcndSum[k] += scnfcndCurrent[k];
                    }
                }

                if (scnfcndSum[k] >= 0x0020000000000000ull) {
                    scnfcndSum[k] = scnfcndSum[k] >> 1;
                    expSum[k]++;
                }
            }
        }

        floatingExp2Integer::Int64PosExp2Int64 sum(scnfcndSum[0], expSum[0]);

        for (unsigned int k = 1; k < parallel_count; k++) {
            floatingExp2Integer::Int64PosExp2Int64 current(scnfcndSum[k], expSum[k]);
            sum += current;
        }

        scnfcnd = sum.scnfcnd;
        exp = sum.exp;

        this->scale();
    }

    Int64PosExp2Int64& Int64PosExp2Int64::operator+=(Int64PosExp2Int64 z) {
        int64_t exp_diff = exp - z.exp;

        if (exp_diff >= 0) {
            if (exp_diff > 63) {
                return *this;
            }
            scnfcnd += z.scnfcnd >> exp_diff;
        }
        else {
            if (exp_diff < -63) {
                scnfcnd = z.scnfcnd;
                exp = z.exp;
                return *this;
            }
            exp = z.exp;
            scnfcnd = scnfcnd >> (-exp_diff);
            scnfcnd += z.scnfcnd;
        }        
        
        if (scnfcnd >= 0x0020000000000000ull) {
            scnfcnd = scnfcnd >> 1;
            exp++;
        }
        return *this;
    }

    void Int64PosExp2Int64::multiply(const std::vector<floatingExp2Integer::Int64PosExp2Int64>& vector) {
        const unsigned int parallel_count = 4;

        uint64_t scnfcndSum[parallel_count];
        int64_t expSum[parallel_count];

        for (unsigned int k = 0; k < parallel_count; k++) {
            scnfcndSum[k] = vector[k].scnfcnd;
            expSum[k] = vector[k].exp;
        }

        uint64_t scnfcndCurrent[parallel_count];
        int64_t expCurrent[parallel_count];

        for (unsigned int i = parallel_count; i + (parallel_count - 1) < vector.size(); i += parallel_count) {
            for (unsigned int k = 0; k < parallel_count; k++) {
                scnfcndCurrent[k] = vector[i + k].scnfcnd;
                expCurrent[k] = vector[i + k].exp;
            }

            for (unsigned int k = 0; k < parallel_count; k++) {
                uint32_t offset = 32 - __builtin_clzll(scnfcndSum[k]);
                uint32_t offset_z = 32 - __builtin_clzll(scnfcndCurrent[k]);
                scnfcndSum[k] = (scnfcndSum[k] >> offset) * (scnfcndCurrent[k] >> offset_z);
                expSum[k] += expCurrent[k] + offset + offset_z;

                if (scnfcndSum[k] >= 0x8000000000000000ull) {
                    scnfcndSum[k] = scnfcndSum[k] >> 11;
                    expSum[k] += 11;
                }
                else {
                    scnfcndSum[k] = scnfcndSum[k] >> 10;
                    expSum[k] += 10;
                }
            }
        }

        floatingExp2Integer::Int64PosExp2Int64 sum(scnfcndSum[0], expSum[0]);
        
        for (unsigned int k = 1; k < parallel_count; k++) {
            floatingExp2Integer::Int64PosExp2Int64 current(scnfcndSum[k], expSum[k]);
            sum *= current;
        }

        scnfcnd = sum.scnfcnd;
        exp = sum.exp;

        this->scale();
    }

    Int64PosExp2Int64& Int64PosExp2Int64::operator*=(Int64PosExp2Int64 z) {
        uint32_t offset = 32 - __builtin_clzll(scnfcnd);
        uint32_t offset_z = 32 - __builtin_clzll(z.scnfcnd);
        scnfcnd = (scnfcnd >> offset) * (z.scnfcnd >> offset_z);
        exp += z.exp + offset + offset_z;

        if (scnfcnd >= 0x8000000000000000ull) {
            scnfcnd = scnfcnd >> 11;
            exp += 11;
        }
        else {
            scnfcnd = scnfcnd >> 10;
            exp += 10;
        }
        return *this;
    }

    inline void Int64PosExp2Int64::fromDouble(double dbl) {
        uint64_t* sgnfcnd_bits = reinterpret_cast<uint64_t*>(&dbl);
        scnfcnd = *sgnfcnd_bits & 0x000FFFFFFFFFFFFFull;
        scnfcnd |= 0x0010000000000000ull;
        exp = ((*sgnfcnd_bits & 0x7FF0000000000000ull) >> 52) - (1023 + 52);
    }

    inline void Int64PosExp2Int64::checkRuleForScale() {
        if (0x8000000000000000ull <= scnfcnd) {
            this->scale();
        }
    }

    inline void Int64PosExp2Int64::scale() {
        uint32_t offset = 11 - __builtin_clzll(scnfcnd);
        exp += offset;
        scnfcnd = scnfcnd >> offset;
    }

    Int64PosExp2Int64 operator+(Int64PosExp2Int64 a, const Int64PosExp2Int64 b) { return a+=b; }
    Int64PosExp2Int64 operator*(Int64PosExp2Int64 a, const Int64PosExp2Int64 b) { return a*=b; }
}
