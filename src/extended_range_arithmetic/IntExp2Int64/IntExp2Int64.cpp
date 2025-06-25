#include <cstdint>
#include <cmath>
#include <vector>
#include "IntExp2Int64.h"

namespace extended_range_arithmetic
{
    IntExp2Int64::IntExp2Int64() {
        this->fromDouble(1.0);
        this->scale();
    }

    IntExp2Int64::IntExp2Int64(double dbl) {
        this->fromDouble(dbl);
        this->scale();
    }

    IntExp2Int64::IntExp2Int64(uint64_t principal_part, int64_t auxiliary_index) {
        principal = principal_part;
        aux = auxiliary_index;
        this->scale();
    }

    void IntExp2Int64::double_to(double dbl) {
        this->fromDouble(dbl);
        this->scale();
    }

    void IntExp2Int64::log2_to(double log2) {
        int64_t exponent = (int64_t)log2;
        double dbl = std::exp2(log2 - exponent);
        this->fromDouble(dbl);
        aux += exponent;
        this->scale();
    }

    double IntExp2Int64::as_log2() const {
        return std::log2(principal) + aux;
    }

    double IntExp2Int64::as_double() const { 
        return (double)principal * std::pow(2.0, aux);
    }

    void IntExp2Int64::sum(const std::vector<extended_range_arithmetic::IntExp2Int64>& vector) {
        const size_t parallel_count = 4;

        if (2 * parallel_count > vector.size()) {
            extended_range_arithmetic::IntExp2Int64 sum = vector[0];
            for (size_t i = 1; i < vector.size(); i++) {
                sum += vector[i];
            }
            principal = sum.principal;
            aux = sum.aux;
            return;
        }

        uint64_t principalSum[parallel_count];
        int64_t auxSum[parallel_count];

        for (size_t k = 0; k < parallel_count; k++) {
            principalSum[k] = vector[k].principal;
            auxSum[k] = vector[k].aux;
        }

        uint64_t principalCurrent[parallel_count];
        int64_t auxCurrent[parallel_count];

        size_t i = parallel_count;

        for (i = parallel_count; i + (parallel_count - 1) < vector.size(); i += parallel_count) {
            for (size_t k = 0; k < parallel_count; k++) {
                principalCurrent[k] = vector[i + k].principal;
                auxCurrent[k] = vector[i + k].aux;

                int64_t aux_diff = auxSum[k] - auxCurrent[k];

                if (aux_diff >= 0) {
                    if (aux_diff > 62) {
                        //
                    }
                    else {
                        principalSum[k] += principalCurrent[k] >> aux_diff;
                    }
                }
                else {
                    if (aux_diff < -62) {
                        principalSum[k] = principalCurrent[k];
                        auxSum[k] = auxCurrent[k];
                    }
                    else {
                        auxSum[k] = auxCurrent[k];
                        principalSum[k] = principalSum[k] >> (-aux_diff);
                        principalSum[k] += principalCurrent[k];
                    }
                }

                if (principalSum[k] >= 0x8000000000000000ull) {
                    principalSum[k] = principalSum[k] >> 11;
                    auxSum[k] += 11LL;
                }
            }
        }

        for (size_t k = 0; k < parallel_count; k++) {
            uint32_t offset = 11 - __builtin_clzll(principalSum[k]);
            principalSum[k] = principalSum[k] >> offset;
            auxSum[k] += offset;
        }

        extended_range_arithmetic::IntExp2Int64 sum(principalSum[0], auxSum[0]);

        for (size_t k = 1; k < parallel_count; k++) {
            extended_range_arithmetic::IntExp2Int64 current(principalSum[k], auxSum[k]);
            sum += current;
        }

        for (i = i; i < vector.size(); i++) {
            sum += vector[i];
        }

        principal = sum.principal;
        aux = sum.aux;

        // already scaled ok
    }

    IntExp2Int64& IntExp2Int64::operator+=(IntExp2Int64 z) {
        int64_t aux_diff = aux - z.aux;

        if (aux_diff >= 0) {
            if (aux_diff > 52) {
                return *this;
            }
            principal += z.principal >> aux_diff;
        }
        else {
            if (aux_diff < -52) {
                principal = z.principal;
                aux = z.aux;
                return *this;
            }
            aux = z.aux;
            principal = principal >> (-aux_diff);
            principal += z.principal;
        }        
        
        if (principal >= 0x0020000000000000ull) {
            principal = principal >> 1;
            aux++;
        }
        return *this;
    }

    void IntExp2Int64::multiply(const std::vector<extended_range_arithmetic::IntExp2Int64>& vector) {
        const size_t parallel_count = 4;

        if (2 * parallel_count > vector.size()) {
            extended_range_arithmetic::IntExp2Int64 sum = vector[0];
            for (size_t i = 1; i < vector.size(); i++) {
                sum *= vector[i];
            }
            principal = sum.principal;
            aux = sum.aux;
            return;
        }

        uint64_t principalSum[parallel_count];
        int64_t auxSum[parallel_count];

        for (size_t k = 0; k < parallel_count; k++) {
            principalSum[k] = vector[k].principal;
            auxSum[k] = vector[k].aux;
        }

        uint64_t principalCurrent[parallel_count];
        int64_t auxCurrent[parallel_count];

        size_t i = parallel_count;

        for (i = parallel_count; i + (parallel_count - 1) < vector.size(); i += parallel_count) {
            for (size_t k = 0; k < parallel_count; k++) {
                principalCurrent[k] = vector[i + k].principal;
                auxCurrent[k] = vector[i + k].aux;

                uint32_t offset = 32 - __builtin_clzll(principalSum[k]);
                uint32_t offset_z = 32 - __builtin_clzll(principalCurrent[k]);
                principalSum[k] = (principalSum[k] >> offset) * (principalCurrent[k] >> offset_z);
                auxSum[k] += auxCurrent[k] + offset + offset_z;

                if (principalSum[k] >= 0x8000000000000000ull) {
                    principalSum[k] = principalSum[k] >> 11;
                    auxSum[k] += 11;
                }
                else {
                    principalSum[k] = principalSum[k] >> 10;
                    auxSum[k] += 10;
                }
            }
        }

        extended_range_arithmetic::IntExp2Int64 sum(principalSum[0], auxSum[0]);
        
        for (size_t k = 1; k < parallel_count; k++) {
            extended_range_arithmetic::IntExp2Int64 current(principalSum[k], auxSum[k]);
            sum *= current;
        }

        for (i = i; i < vector.size(); i++) {
            sum *= vector[i];
        }

        principal = sum.principal;
        aux = sum.aux;

        // already scaled ok
    }

    IntExp2Int64& IntExp2Int64::operator*=(IntExp2Int64 z) {
        uint32_t offset = 32 - __builtin_clzll(principal);
        uint32_t offset_z = 32 - __builtin_clzll(z.principal);
        principal = (principal >> offset) * (z.principal >> offset_z);
        aux += z.aux + offset + offset_z;

        if (principal >= 0x8000000000000000ull) {
            principal = principal >> 11;
            aux += 11;
        }
        else {
            principal = principal >> 10;
            aux += 10;
        }
        return *this;
    }

    inline void IntExp2Int64::fromDouble(double dbl) {
        uint64_t* principal_bits = reinterpret_cast<uint64_t*>(&dbl);
        principal = *principal_bits & 0x000FFFFFFFFFFFFFull;
        principal |= 0x0010000000000000ull;
        aux = ((*principal_bits & 0x7FF0000000000000ull) >> 52) - (1023 + 52);
    }

    inline void IntExp2Int64::checkRuleForScale() {
        if (0x8000000000000000ull <= principal) {
            this->scale();
        }
    }

    inline void IntExp2Int64::scale() {
        uint32_t offset = 11 - __builtin_clzll(principal);
        aux += offset;
        principal = principal >> offset;
    }

    IntExp2Int64 operator+(IntExp2Int64 a, const IntExp2Int64 b) { return a+=b; }
    IntExp2Int64 operator*(IntExp2Int64 a, const IntExp2Int64 b) { return a*=b; }
}
