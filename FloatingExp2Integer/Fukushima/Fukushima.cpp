#include <cmath>
#include <cstdint>
#include <bit>
#include <bitset>
#include <iostream>

#include "Fukushima.h"

namespace floatingExp2Integer
{
    Fukushima::Fukushima() {
        scnfcnd = 1.0;
        exp = 0LL;
    }

    Fukushima::Fukushima(double dbl) {
        scnfcnd = dbl;
        exp = 0LL;
    }

    Fukushima::Fukushima(double dbl, int64_t exponent) {
        scnfcnd = dbl;
        exp = exponent;
    }

    void Fukushima::log2_to(const double from) {
        double from_floored = std::floor(from);
        int64_t exponent = (int64_t)from_floored;
        double multiplier = std::exp2(from - from_floored);

        if (exponent >= 0) {
            exp = (exponent + EXP_MULTIPLIER / 2) / EXP_MULTIPLIER;
            scnfcnd = multiplier * std::exp2(exponent - EXP_MULTIPLIER * exp);
        }
        else {
            exp = (exponent - EXP_MULTIPLIER / 2) / EXP_MULTIPLIER;
            if ((exponent - EXP_MULTIPLIER / 2) - EXP_MULTIPLIER * exp == 0LL && multiplier > 1.0) {
                exp++;
            }
            scnfcnd = multiplier * std::exp2(exponent - EXP_MULTIPLIER * exp);
        }
    }

    double Fukushima::as_double() const {
        return scnfcnd * std::exp2(EXP_MULTIPLIER * exp);
    }

    double Fukushima::as_log2() const {
        return std::log2(scnfcnd) + EXP_MULTIPLIER * exp;
    }

    void Fukushima::sum(const std::vector<floatingExp2Integer::Fukushima>& vector) {
        const size_t parallel_count = 4;

        if (2 * parallel_count > vector.size()) {
            floatingExp2Integer::Fukushima sum = vector[0];
            for (size_t i = 1; i < vector.size(); i++) {
                sum += vector[i];
            }
            scnfcnd = sum.scnfcnd;
            exp = sum.exp;
            return;
        }

        double scnfcndSum[parallel_count];
        int64_t expSum[parallel_count];

        for (size_t k = 0; k < parallel_count; k++) {
            scnfcndSum[k] = vector[k].scnfcnd;
            expSum[k] = vector[k].exp;
        }

        double scnfcndCurrent[parallel_count];
        int64_t expCurrent[parallel_count];

        size_t i = parallel_count;

        for (i = parallel_count; i + (parallel_count - 1) < vector.size(); i += parallel_count) {
            for (size_t k = 0; k < parallel_count; k++) {
                scnfcndCurrent[k] = vector[i + k].scnfcnd;
                expCurrent[k] = vector[i + k].exp;

                if (expSum[k] == expCurrent[k]) {
                    scnfcndSum[k] = scnfcndSum[k] + scnfcndCurrent[k];
                }
                else {
                    int64_t id = expSum[k] - expCurrent[k];
                    if (id == 1) {
                        scnfcndSum[k] = scnfcndSum[k] + (scnfcndCurrent[k] * BIGI);
                    }
                    else if (id > 1) {
                        scnfcndSum[k] = scnfcndSum[k];
                    }
                    else if (id == -1) {
                        scnfcndSum[k] = scnfcndCurrent[k] + (scnfcndSum[k] * BIGI);
                        expSum[k] = expCurrent[k];
                    }
                    else {
                        scnfcndSum[k] = scnfcndCurrent[k];
                        expSum[k] = expCurrent[k];
                    }
                }

                if (scnfcndSum[k] >= BIGS) {
                    scnfcndSum[k] *= BIGI;
                    expSum[k] += 1;
                }
            }
        }

        floatingExp2Integer::Fukushima sum(scnfcndSum[0], expSum[0]);

        for (size_t k = 1; k < parallel_count; k++) {
            floatingExp2Integer::Fukushima current(scnfcndSum[k], expSum[k]);
            sum += current;
        }

        for (i = i; i < vector.size(); i++) {
            sum += vector[i];
        }

        scnfcnd = sum.scnfcnd;
        exp = sum.exp;
    }

    void Fukushima::multiply(const std::vector<floatingExp2Integer::Fukushima>& vector) {
        const size_t parallel_count = 6;

        if (2 * parallel_count > vector.size()) {
            floatingExp2Integer::Fukushima res = vector[0];
            for (size_t i = 1; i < vector.size(); i++) {
                res *= vector[i];
            }
            scnfcnd = res.scnfcnd;
            exp = res.exp;
            return;
        }

        double scnfcndSum[parallel_count];
        int64_t expSum[parallel_count];

        for (size_t k = 0; k < parallel_count; k++) {
            scnfcndSum[k] = vector[k].scnfcnd;
            expSum[k] = vector[k].exp;
        }

        double scnfcndCurrent[parallel_count];
        int64_t expCurrent[parallel_count];

        size_t i = parallel_count;

        for (i = parallel_count; i + (parallel_count - 1) < vector.size(); i += parallel_count) {
            for (size_t k = 0; k < parallel_count; k++) {
                scnfcndCurrent[k] = vector[i + k].scnfcnd;
                expCurrent[k] = vector[i + k].exp;

                scnfcndSum[k] = scnfcndSum[k] * scnfcndCurrent[k];
                expSum[k] = expSum[k] + expCurrent[k];

                if (scnfcndSum[k] >= BIGS) {
                    scnfcndSum[k] *= BIGI;
                    expSum[k] += 1;
                }
                else if (scnfcndSum[k] < BIGSI) {
                    scnfcndSum[k] *= BIG;
                    expSum[k] -= 1;
                }
            }
        }

        floatingExp2Integer::Fukushima sum(scnfcndSum[0], expSum[0]);

        for (size_t k = 1; k < parallel_count; k++) {
            floatingExp2Integer::Fukushima current(scnfcndSum[k], expSum[k]);
            sum *= current;
        }

        for (i = i; i < vector.size(); i++) {
            sum *= vector[i];
        }

        scnfcnd = sum.scnfcnd;
        exp = sum.exp;
    }

    Fukushima& Fukushima::operator+=(Fukushima z) {
        if (exp == z.exp) {
            scnfcnd = scnfcnd + z.scnfcnd;
        } 
        else {
            int64_t id = exp - z.exp;
            if (id == 1) {
                scnfcnd = scnfcnd + (z.scnfcnd * BIGI);
            } else if (id > 1) {
                scnfcnd = scnfcnd;
            } else if (id == -1) {
                scnfcnd = z.scnfcnd + (scnfcnd * BIGI);
                exp = z.exp;
            } else {
                scnfcnd = z.scnfcnd;
                exp = z.exp;
            }
        }

        if (scnfcnd >= BIGS) {
            scnfcnd *= BIGI;
            exp += 1;
        }

        return *this;
    }

    Fukushima& Fukushima::operator*=(Fukushima z) {
        scnfcnd = scnfcnd * z.scnfcnd;
        exp = exp + z.exp;

        xnorm();

        return *this;
    }

    inline void Fukushima::xnorm() {
        if (scnfcnd >= BIGS) {
            scnfcnd *= BIGI;
            exp += 1;
        } else if (scnfcnd < BIGSI) {
            scnfcnd *= BIG;
            exp -= 1;
        }
    }

    // Fukushima operator+(Fukushima a, const Fukushima b) { return a+=b; }
    // Fukushima operator*(Fukushima a, const Fukushima b) { return a*=b; }
}

