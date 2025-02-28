#include <cstdint>
#include <cmath>
#include <vector>
#include <cstring>

#include <cstdint>
#include <bit>
#include <iostream>
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

    Float64PosExp2Int64::Float64PosExp2Int64(double sicnificand, std::int64_t exponent) {
        scnfcnd = sicnificand;
        exp = exponent;
        this->scale();
    }

    Float64PosExp2Int64::Float64PosExp2Int64(std::vector<floatingExp2Integer::Float64PosExp2Int64>& vector) {
        const int pn = 3;

        if (vector.size() < 4 * pn) {
            Float64PosExp2Int64 a(0, std::numeric_limits<int64_t>::min());
            for (unsigned int i = i; i < vector.size(); i++) {
                a += vector[i];
            }    
            scnfcnd = a.scnfcnd;
            exp = a.exp;
        }

        double4_t sa[pn];
        int64_4_t ea[pn];
    
        for (int p = 0; p < pn; p++) {
            for (int k = 0; k < 4; k++) {
                sa[p][k] = vector[p*4 + k].scnfcnd;
                ea[p][k] = vector[p*4 + k].exp;
            }
        }

        double4_t sb[pn];
        int64_4_t eb[pn];

        int64_4_t exp_diff[pn];
        int64_4_t sign_exp_diff[pn];

        uint64_4_t ib1[pn];
        uint64_4_t add[pn];

        uint64_4_t ia1[pn];
    
        unsigned int index;
        int counter = 0;
        for (index = 4 * pn; index + (4 * pn - 1) < vector.size(); index += 4 * pn) {    
            for (int p = 0; p < pn; p++) {
                for (int k = 0; k < 4; k++) {
                    int vector_index = index + p*4 + k;
                    sb[p][k] = vector[vector_index].scnfcnd;
                    eb[p][k] = vector[vector_index].exp;
                }
            }
            for (int p = 0; p < pn; p++) {
                exp_diff[p] = ea[p] - eb[p];
                sign_exp_diff[p] = exp_diff[p] >= 0LL ? 1LL : -1LL;
            // }
            // //uint64_4_t& ib1 = reinterpret_cast<uint64_4_t&>(sb[p]);
            // std::memcpy(&ib1, &sb, sizeof ib1); 
            // for (int p = 0; p < pn; p++) {
            
                ib1[p] = std::bit_cast<uint64_4_t>(sb[p]);

                add[p] = (uint64_4_t)((exp_diff[p] * sign_exp_diff[p]) << 52) & 0x7FF0000000000000ull;
                ib1[p] -= add[p] * sign_exp_diff[p];

                sb[p] = std::bit_cast<double4_t>(ib1[p]);
            // }
            // std::memcpy(&sb, &ib1, sizeof sb);
            // for (int p = 0; p < pn; p++) {
        
                sa[p] = exp_diff[p] <= -64LL ? double4_0 : sa[p];
                ea[p] = exp_diff[p] <= -64LL ? eb[p] : ea[p];
                sb[p] = exp_diff[p] >= 64LL ? double4_0 : sb[p];
        
                sa[p] += sb[p];
            }
            counter++;
            if (counter == 10) {
                //std::memcpy(&ia1, &sa, sizeof ia1);
                for (int p = 0; p < pn; p++) {
                    //uint64_4_t& ia1 = reinterpret_cast<uint64_4_t&>(sa[p]);
                    ia1[p] = std::bit_cast<uint64_4_t>(sa[p]);
                    ea[p] = (ea[p] - int64_4_1023) + (int64_4_t)((ia1[p] & 0x7FF0000000000000ull) >> 52);
                    ia1[p] &= 0x800FFFFFFFFFFFFFull;
                    ia1[p] |= 0x3FF0000000000000ull;
                    sa[p] = std::bit_cast<double4_t>(ia1[p]);
                }
                //std::memcpy(&sa, &ia1, sizeof sa);
                counter = 0;
            }
        }        
        
        Float64PosExp2Int64 a(1, std::numeric_limits<int64_t>::min());
        for (int p = 0; p < pn; p++) {
            for (int k = 0; k < 4; k++) {
                Float64PosExp2Int64 z(sa[p][k], ea[p][k]);
                a += z;
            }
            std::cout << std::endl;
        }

        for (index = index; index < vector.size(); index++) {
            a += vector[index];
        }

        scnfcnd = a.scnfcnd;
        exp = a.exp;
    }

    void Float64PosExp2Int64::doubleToFloat64PosExp2Int64(double dbl) {
        scnfcnd = dbl;
        exp = 0;
        this->scale();
    }

    void Float64PosExp2Int64::log2ToFloat64PosExp2Int64(double log2) {
        exp = (std::int64_t)log2;
        scnfcnd = std::exp2(log2 - exp);
        this->scale();
    }

    double Float64PosExp2Int64::float64PosExp2Int64ToLog2() const {
        return std::log2(scnfcnd) + exp;
    }

    double Float64PosExp2Int64::sicnificand() {
        this->scale();
        return scnfcnd; 
    }

    std::int64_t Float64PosExp2Int64::exponent() {
        this->scale();
        return exp;
    }

    double Float64PosExp2Int64::asDouble() const {
        return scnfcnd * std::pow(2.0, exp);
    }

    Float64PosExp2Int64& Float64PosExp2Int64::operator+=(Float64PosExp2Int64 z) {
        std::int64_t exp_diff = exp - z.exp;
        std::uint64_t* sgnfcnd_bits_z = reinterpret_cast<std::uint64_t*>(&z.scnfcnd);

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
            *sgnfcnd_bits_z += (-exp_diff) << 52;
            scnfcnd += z.scnfcnd;
        }

        this->checkRuleForScale();    
        return *this;
    }

    Float64PosExp2Int64& Float64PosExp2Int64::operator*=(Float64PosExp2Int64 z) {
        scnfcnd *= z.scnfcnd;
        exp += z.exp;
        this->checkRuleForScale();
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

