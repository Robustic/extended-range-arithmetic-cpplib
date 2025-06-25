#include <cstdint>
#include <cmath>
#include <vector>
#include <cstring>
#include "FloatExp2Int64.h"

typedef double double4_t __attribute__((vector_size(4 * sizeof(double))));
typedef uint64_t uint64_4_t __attribute__((vector_size(4 * sizeof(uint64_t))));
typedef int64_t int64_4_t __attribute__((vector_size(4 * sizeof(int64_t))));

constexpr int64_4_t int64_4_1023 { 1023LL, 1023LL, 1023LL, 1023LL };
constexpr double4_t double4_0 { 0.0, 0.0, 0.0, 0.0 };

namespace extended_range_arithmetic
{
    FloatExp2Int64::FloatExp2Int64() {
        principal = 1.0;
        aux = 0;
        this->scale();
    }

    FloatExp2Int64::FloatExp2Int64(double dbl) {
        principal = dbl;
        aux = 0;
        this->scale();
    }

    FloatExp2Int64::FloatExp2Int64(double sicnificand, int64_t exponent) {
        principal = sicnificand;
        aux = exponent;
        this->scale();
    }

    void FloatExp2Int64::sum(const std::vector<extended_range_arithmetic::FloatExp2Int64>& vector) {
        const size_t parallel_count = 4;

        if (2 * parallel_count > vector.size()) {
            extended_range_arithmetic::FloatExp2Int64 sum = vector[0];
            for (size_t i = 1; i < vector.size(); i++) {
                sum += vector[i];
            }
            principal = sum.principal;
            aux = sum.aux;
            return;
        }

        double principalSum[parallel_count];
        int64_t auxSum[parallel_count];

        for (size_t k = 0; k < parallel_count; k++) {
            principalSum[k] = vector[k].principal;
            auxSum[k] = vector[k].aux;
        }

        double principalCurrent[parallel_count];
        int64_t auxCurrent[parallel_count];

        size_t i = parallel_count;

        for (i = parallel_count; i + (parallel_count - 1) < vector.size(); i += parallel_count) {
            for (size_t k = 0; k < parallel_count; k++) {
                principalCurrent[k] = vector[i + k].principal;
                auxCurrent[k] = vector[i + k].aux;

                int64_t aux_diff = auxSum[k] - auxCurrent[k];
                uint64_t* principal_bits = reinterpret_cast<uint64_t*>(&principalSum[k]);
                uint64_t* principal_bits_z = reinterpret_cast<uint64_t*>(&principalCurrent[k]);

                if (aux_diff >= 0) {
                    if (aux_diff > 64) {
                        //
                    }
                    else {
                        *principal_bits_z -= aux_diff << 52;
                        principalSum[k] += principalCurrent[k];
                    }
                }
                else {
                    if (aux_diff < -64) {
                        principalSum[k] = principalCurrent[k];
                        auxSum[k] = auxCurrent[k];
                    }
                    else {
                        auxSum[k] = auxCurrent[k];
                        *principal_bits -= (-aux_diff) << 52;
                        principalSum[k] += principalCurrent[k];
                    }
                }

                if (principalSum[k] >= 0x1p12) {
                    uint64_t* principal_bits = reinterpret_cast<uint64_t*>(&principalSum[k]);
                    auxSum[k] += ((*principal_bits & 0x7FF0000000000000ull) >> 52) - 1023;
                    *principal_bits &= 0x800FFFFFFFFFFFFFull;
                    *principal_bits |= 0x3FF0000000000000ull;
                    principalSum[k] = *reinterpret_cast<double*>(principal_bits);
                }
            }
        }

        for (size_t k = 0; k < parallel_count; k++) {
            uint64_t* principal_bits = reinterpret_cast<uint64_t*>(&principalSum[k]);
            auxSum[k] += ((*principal_bits & 0x7FF0000000000000ull) >> 52) - 1023;
            *principal_bits &= 0x800FFFFFFFFFFFFFull;
            *principal_bits |= 0x3FF0000000000000ull;
            principalSum[k] = *reinterpret_cast<double*>(principal_bits);
        }

        extended_range_arithmetic::FloatExp2Int64 sum(principalSum[0], auxSum[0]);

        for (size_t k = 1; k < parallel_count; k++) {
            extended_range_arithmetic::FloatExp2Int64 current(principalSum[k], auxSum[k]);
            sum += current;
        }

        for (i = i; i < vector.size(); i++) {
            sum += vector[i];
        }

        principal = sum.principal;
        aux = sum.aux;

        // scaling already done
    }

    void FloatExp2Int64::double_to(double dbl) {
        principal = dbl;
        aux = 0;
        this->scale();
    }

    void FloatExp2Int64::log2_to(double log2) {
        aux = (int64_t)log2;
        principal = std::exp2(log2 - aux);
        this->scale();
    }

    double FloatExp2Int64::as_log2() const {
        return std::log2(principal) + aux;
    }

    double FloatExp2Int64::as_double() const {
        return principal * std::pow(2.0, aux);
    }

    FloatExp2Int64& FloatExp2Int64::operator+=(FloatExp2Int64 z) {
        int64_t aux_diff = aux - z.aux;
        uint64_t* principal_bits = reinterpret_cast<uint64_t*>(&principal);
        uint64_t* principal_bits_z = reinterpret_cast<uint64_t*>(&z.principal);

        if (aux_diff >= 0) {
            if (aux_diff > 52) {
                return *this;
            }
            *principal_bits_z -= aux_diff << 52;
            principal += z.principal;
        }
        else {
            if (aux_diff < -52) {
                principal = z.principal;
                aux = z.aux;
                return *this;
            }
            aux = z.aux;
            *principal_bits -= (-aux_diff) << 52;
            principal += z.principal;
        }

        if (principal >= 2.0) {
            principal *= 0.5;
            aux++;
        }
        return *this;
    }

    void FloatExp2Int64::multiply(const std::vector<extended_range_arithmetic::FloatExp2Int64>& vector) {
        const size_t parallel_count = 4;

        if (2 * parallel_count > vector.size()) {
            extended_range_arithmetic::FloatExp2Int64 sum = vector[0];
            for (size_t i = 1; i < vector.size(); i++) {
                sum *= vector[i];
            }
            principal = sum.principal;
            aux = sum.aux;
            return;
        }

        double principalSum[parallel_count];
        int64_t auxSum[parallel_count];

        for (size_t k = 0; k < parallel_count; k++) {
            principalSum[k] = vector[k].principal;
            auxSum[k] = vector[k].aux;
        }

        double principalCurrent[parallel_count];
        int64_t auxCurrent[parallel_count];

        size_t i = parallel_count;

        for (i = parallel_count; i + (parallel_count - 1) < vector.size(); i += parallel_count) {
            for (size_t k = 0; k < parallel_count; k++) {
                principalCurrent[k] = vector[i + k].principal;
                auxCurrent[k] = vector[i + k].aux;

                principalSum[k] *= principalCurrent[k];
                auxSum[k] += auxCurrent[k];

                if (principalSum[k] >= 0x1p500) {
                    uint64_t* principal_bits = reinterpret_cast<uint64_t*>(&principalSum[k]);
                    auxSum[k] += ((*principal_bits & 0x7FF0000000000000ull) >> 52) - 1023;
                    *principal_bits &= 0x800FFFFFFFFFFFFFull;
                    *principal_bits |= 0x3FF0000000000000ull;
                    principalSum[k] = *reinterpret_cast<double*>(principal_bits);
                }
            }
        }

        for (size_t k = 0; k < parallel_count; k++) {
            uint64_t* principal_bits = reinterpret_cast<uint64_t*>(&principalSum[k]);
            auxSum[k] += ((*principal_bits & 0x7FF0000000000000ull) >> 52) - 1023;
            *principal_bits &= 0x800FFFFFFFFFFFFFull;
            *principal_bits |= 0x3FF0000000000000ull;
            principalSum[k] = *reinterpret_cast<double*>(principal_bits);
        }

        extended_range_arithmetic::FloatExp2Int64 sum(principalSum[0], auxSum[0]);

        for (size_t k = 1; k < parallel_count; k++) {
            extended_range_arithmetic::FloatExp2Int64 current(principalSum[k], auxSum[k]);
            sum *= current;
        }

        for (i = i; i < vector.size(); i++) {
            sum *= vector[i];
        }

        principal = sum.principal;
        aux = sum.aux;

        // scaling already done
    }

    FloatExp2Int64& FloatExp2Int64::operator*=(FloatExp2Int64 z) {
        principal *= z.principal;
        aux += z.aux;
        if (principal >= 2.0) {
            principal *= 0.5;
            aux++;
        }
        return *this;
    }

    FloatExp2Int64& FloatExp2Int64::operator*=(double& dbl) {
        FloatExp2Int64 z(dbl);
        principal *= z.principal;
        aux += z.aux;
        this->checkRuleForScale();
        return *this;
    }

    inline void FloatExp2Int64::checkRuleForScale() {
        if (0x1p12 <= principal) {
            this->scale();
        }
    }

    inline void FloatExp2Int64::scale() {
        uint64_t* principal_bits = reinterpret_cast<uint64_t*>(&principal);
        aux += ((*principal_bits & 0x7FF0000000000000ull) >> 52) - 1023;
        *principal_bits &= 0x800FFFFFFFFFFFFFull;
        *principal_bits |= 0x3FF0000000000000ull;
        principal = *reinterpret_cast<double*>(principal_bits);
    }

    FloatExp2Int64 operator+(FloatExp2Int64 a, const FloatExp2Int64 b) { return a+=b; }
    FloatExp2Int64 operator*(FloatExp2Int64 a, const FloatExp2Int64 b) { return a*=b; }
}

