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

    Float64PosExp2Int64::Float64PosExp2Int64(double sicnificand, std::int64_t exponent) {
        scnfcnd = sicnificand;
        exp = exponent;
        this->scale();
    }

    Float64PosExp2Int64::Float64PosExp2Int64(std::vector<floatingExp2Integer::Float64PosExp2Int64>& vector) {
        double4_t sa1;
        int64_4_t ea1;
    
        for (int k = 0; k < 4; k++) {
            sa1[k] = vector[k].scnfcnd;
            ea1[k] = vector[k].exp;
        }
    
        unsigned int i;
        for (i = 4; i + 3 < vector.size(); i += 4) {
            double4_t sb1;
            int64_4_t eb1;
    
            for (int k = 0; k < 4; k++) {
                sb1[k] = vector[i + k].scnfcnd;
                eb1[k] = vector[i + k].exp;
            }
    
            int64_4_t exp_diff = ea1 - eb1;
            int64_4_t sign_exp_diff = exp_diff >= 0LL ? 1LL : -1LL;
            
            uint64_4_t ib1;
            std::memcpy(&ib1, &sb1, sizeof ib1);
    
            uint64_4_t add = (uint64_4_t)((exp_diff * sign_exp_diff) << 52) & 0x7FF0000000000000ull;
            ib1 -= add * sign_exp_diff;
    
            std::memcpy(&sb1, &ib1, sizeof sb1);
    
            sa1 = exp_diff <= -64LL ? double4_0 : sa1;
            sb1 = exp_diff >= 64LL ? double4_0 : sb1;
    
            sa1 += sb1;
    
            uint64_4_t ia1;
            std::memcpy(&ia1, &sa1, sizeof ia1);
            ea1 = (ea1 - int64_4_1023) + (int64_4_t)((ia1 & 0x7FF0000000000000ull) >> 52);
            ia1 &= 0x800FFFFFFFFFFFFFull;
            ia1 |= 0x3FF0000000000000ull;
            std::memcpy(&sa1, &ia1, sizeof sa1);
        }        
        
        Float64PosExp2Int64 a(sa1[0], ea1[0]);
        for (int k = 1; k < 4; k++) {
            Float64PosExp2Int64 z(sa1[k], ea1[k]);
            a += z;
        }

        for (i = i; i < vector.size(); i++) {
            a += vector[i];
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

