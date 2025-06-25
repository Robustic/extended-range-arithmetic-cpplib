#include <cmath>
#include <cstdint>
#include <bit>
#include <bitset>
#include <iostream>

#include "Xnumber64.h"

namespace extended_range_arithmetic
{
    Xnumber64::Xnumber64() {
        principal = 1.0;
        aux = 0LL;
    }

    Xnumber64::Xnumber64(double dbl) {
        principal = dbl;
        aux = 0LL;
    }

    Xnumber64::Xnumber64(double principal_part, int64_t auxiliary_index) {
        principal = principal_part;
        aux = auxiliary_index;
    }

    void Xnumber64::log2_to(const double from) {
        double from_floored = std::floor(from);
        int64_t exponent = (int64_t)from_floored;
        double multiplier = std::exp2(from - from_floored);

        if (exponent >= 0) {
            aux = (exponent + EXP_MULTIPLIER / 2) / EXP_MULTIPLIER;
            principal = multiplier * std::exp2(exponent - EXP_MULTIPLIER * aux);
        }
        else {
            aux = (exponent - EXP_MULTIPLIER / 2) / EXP_MULTIPLIER;
            if ((exponent - EXP_MULTIPLIER / 2) - EXP_MULTIPLIER * aux == 0LL && multiplier > 1.0) {
                aux++;
            }
            principal = multiplier * std::exp2(exponent - EXP_MULTIPLIER * aux);
        }
    }

    double Xnumber64::as_double() const {
        return principal * std::exp2(EXP_MULTIPLIER * aux);
    }

    double Xnumber64::as_log2() const {
        return std::log2(principal) + EXP_MULTIPLIER * aux;
    }

    void Xnumber64::sum(const std::vector<extended_range_arithmetic::Xnumber64>& vector) {
        const size_t parallel_count = 4;

        if (2 * parallel_count > vector.size()) {
            extended_range_arithmetic::Xnumber64 sum = vector[0];
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

                if (auxSum[k] == auxCurrent[k]) {
                    principalSum[k] = principalSum[k] + principalCurrent[k];
                }
                else {
                    int64_t id = auxSum[k] - auxCurrent[k];
                    if (id == 1) {
                        principalSum[k] = principalSum[k] + (principalCurrent[k] * BIGI);
                    }
                    else if (id > 1) {
                        principalSum[k] = principalSum[k];
                    }
                    else if (id == -1) {
                        principalSum[k] = principalCurrent[k] + (principalSum[k] * BIGI);
                        auxSum[k] = auxCurrent[k];
                    }
                    else {
                        principalSum[k] = principalCurrent[k];
                        auxSum[k] = auxCurrent[k];
                    }
                }

                if (principalSum[k] >= BIGS) {
                    principalSum[k] *= BIGI;
                    auxSum[k] += 1;
                }
            }
        }

        extended_range_arithmetic::Xnumber64 sum(principalSum[0], auxSum[0]);

        for (size_t k = 1; k < parallel_count; k++) {
            extended_range_arithmetic::Xnumber64 current(principalSum[k], auxSum[k]);
            sum += current;
        }

        for (i = i; i < vector.size(); i++) {
            sum += vector[i];
        }

        principal = sum.principal;
        aux = sum.aux;
    }

    void Xnumber64::multiply(const std::vector<extended_range_arithmetic::Xnumber64>& vector) {
        const size_t parallel_count = 6;

        if (2 * parallel_count > vector.size()) {
            extended_range_arithmetic::Xnumber64 res = vector[0];
            for (size_t i = 1; i < vector.size(); i++) {
                res *= vector[i];
            }
            principal = res.principal;
            aux = res.aux;
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

                principalSum[k] = principalSum[k] * principalCurrent[k];
                auxSum[k] = auxSum[k] + auxCurrent[k];

                if (principalSum[k] >= BIGS) {
                    principalSum[k] *= BIGI;
                    auxSum[k] += 1;
                }
                else if (principalSum[k] < BIGSI) {
                    principalSum[k] *= BIG;
                    auxSum[k] -= 1;
                }
            }
        }

        extended_range_arithmetic::Xnumber64 sum(principalSum[0], auxSum[0]);

        for (size_t k = 1; k < parallel_count; k++) {
            extended_range_arithmetic::Xnumber64 current(principalSum[k], auxSum[k]);
            sum *= current;
        }

        for (i = i; i < vector.size(); i++) {
            sum *= vector[i];
        }

        principal = sum.principal;
        aux = sum.aux;
    }

    Xnumber64& Xnumber64::operator+=(Xnumber64 z) {
        if (aux == z.aux) {
            principal = principal + z.principal;
        } 
        else {
            int64_t id = aux - z.aux;
            if (id == 1) {
                principal = principal + (z.principal * BIGI);
            } else if (id == -1) {
                principal = z.principal + (principal * BIGI);
                aux = z.aux;
            } else if (id < -1) {
                principal = z.principal;
                aux = z.aux;
            }
        }

        if (principal >= BIGS) {
            principal *= BIGI;
            aux += 1;
        }

        return *this;
    }

    Xnumber64& Xnumber64::operator*=(Xnumber64 z) {
        principal = principal * z.principal;
        aux = aux + z.aux;

        xnorm();

        return *this;
    }

    inline void Xnumber64::xnorm() {
        if (principal >= BIGS) {
            principal *= BIGI;
            aux += 1;
        } else if (principal < BIGSI) {
            principal *= BIG;
            aux -= 1;
        }
    }

    // Xnumber64 operator+(Xnumber64 a, const Xnumber64 b) { return a+=b; }
    // Xnumber64 operator*(Xnumber64 a, const Xnumber64 b) { return a*=b; }
}

