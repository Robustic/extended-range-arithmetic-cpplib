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

    Float64PosExp2Int64::Float64PosExp2Int64(const std::vector<floatingExp2Integer::Float64PosExp2Int64>& vector) {
        if (vector.size() < 4 * 2) {
            Float64PosExp2Int64 a(vector[0].scnfcnd, vector[0].exp);
            for (unsigned int i = 0; i < vector.size(); i++) {
                a += vector[i];
            }    
            scnfcnd = a.scnfcnd;
            exp = a.exp;
        }

        double4_t sa1;
        int64_4_t ea1;
        double4_t sa2;
        int64_4_t ea2;
        double4_t sa3;
        int64_4_t ea3;
        double4_t sa4;
        int64_4_t ea4;

        double4_t sb1;
        int64_4_t eb1;
        double4_t sb2;
        int64_4_t eb2;
        double4_t sb3;
        int64_4_t eb3;
        double4_t sb4;
        int64_4_t eb4;
    
        for (int k = 0; k < 4; k++) {
            sa1[k] = vector[k].scnfcnd;
            sa2[k] = vector[k + 4].scnfcnd;
            ea1[k] = vector[k].exp;
            ea2[k] = vector[k + 4].exp;
            sa3[k] = vector[k + 8].scnfcnd;
            sa4[k] = vector[k + 12].scnfcnd;
            ea3[k] = vector[k + 8].exp;
            ea4[k] = vector[k + 12].exp;
        }

        const auto* vec_ptr = vector.data();
    
        int counter = 0;
        unsigned int i;
        for (i = 16; i + 15 < vector.size(); i += 16) {
    
            for (int k = 0; k < 4; k++) {
                int ik = i + k;
                int ik4 = ik + 4;
                int ik8 = ik + 8;
                int ik12 = ik + 12;
                sb1[k] = vector[ik].scnfcnd;
                sb2[k] = vector[ik4].scnfcnd;
                sb3[k] = vector[ik8].scnfcnd;
                sb4[k] = vector[ik12].scnfcnd;
                eb1[k] = vector[ik].exp;
                eb2[k] = vector[ik4].exp;
                eb3[k] = vector[ik8].exp;
                eb4[k] = vector[ik12].exp;
            }
    
            int64_4_t exp_diff1 = ea1 - eb1;
            int64_4_t exp_diff2 = ea2 - eb2;
            int64_4_t exp_diff3 = ea3 - eb3;
            int64_4_t exp_diff4 = ea4 - eb4;

            uint64_4_t& ib1 = reinterpret_cast<uint64_4_t&>(sb1); 
            uint64_4_t& ib2 = reinterpret_cast<uint64_4_t&>(sb2);
            uint64_4_t& ib3 = reinterpret_cast<uint64_4_t&>(sb3); 
            uint64_4_t& ib4 = reinterpret_cast<uint64_4_t&>(sb4);
            // ib1 = ((((ib1 & 0x7FF0000000000000ull) >> 52) - exp_diff1) << 52) | (ib1 & 0x800FFFFFFFFFFFFFull);
            // ib2 = ((((ib2 & 0x7FF0000000000000ull) >> 52) - exp_diff2) << 52) | (ib2 & 0x800FFFFFFFFFFFFFull);
            // ib3 = ((((ib3 & 0x7FF0000000000000ull) >> 52) - exp_diff3) << 52) | (ib3 & 0x800FFFFFFFFFFFFFull);
            // ib4 = ((((ib4 & 0x7FF0000000000000ull) >> 52) - exp_diff4) << 52) | (ib4 & 0x800FFFFFFFFFFFFFull);
            ib1 = ib1 - (exp_diff1 << 52);
            ib2 = ib2 - (exp_diff2 << 52);
            ib3 = ib3 - (exp_diff3 << 52);
            ib4 = ib4 - (exp_diff4 << 52);
    
            sa1 = exp_diff1 <= -64LL ? double4_0 : sa1;
            sb1 = exp_diff1 >= 64LL ? double4_0 : sb1;
            sa2 = exp_diff2 <= -64LL ? double4_0 : sa2;
            sb2 = exp_diff2 >= 64LL ? double4_0 : sb2;
            sa3 = exp_diff3 <= -64LL ? double4_0 : sa3;
            sb3 = exp_diff3 >= 64LL ? double4_0 : sb3;
            sa4 = exp_diff4 <= -64LL ? double4_0 : sa4;
            sb4 = exp_diff4 >= 64LL ? double4_0 : sb4;
    
            sa1 += sb1;
            sa2 += sb2;
            sa3 += sb3;
            sa4 += sb4;

            counter++;
            if (counter > 30) {
                counter = 0;
                uint64_4_t& ia1 = reinterpret_cast<uint64_4_t&>(sa1); 
                uint64_4_t& ia2 = reinterpret_cast<uint64_4_t&>(sa2);
                uint64_4_t& ia3 = reinterpret_cast<uint64_4_t&>(sa3); 
                uint64_4_t& ia4 = reinterpret_cast<uint64_4_t&>(sa4);
                ea1 = (ea1 - int64_4_1023) + (int64_4_t)((ia1 & 0x7FF0000000000000ull) >> 52);
                ea2 = (ea2 - int64_4_1023) + (int64_4_t)((ia2 & 0x7FF0000000000000ull) >> 52);
                ea3 = (ea3 - int64_4_1023) + (int64_4_t)((ia3 & 0x7FF0000000000000ull) >> 52);
                ea4 = (ea4 - int64_4_1023) + (int64_4_t)((ia4 & 0x7FF0000000000000ull) >> 52);
                ia1 &= 0x800FFFFFFFFFFFFFull;
                ia1 |= 0x3FF0000000000000ull;
                ia2 &= 0x800FFFFFFFFFFFFFull;
                ia2 |= 0x3FF0000000000000ull;
                ia3 &= 0x800FFFFFFFFFFFFFull;
                ia3 |= 0x3FF0000000000000ull;
                ia4 &= 0x800FFFFFFFFFFFFFull;
                ia4 |= 0x3FF0000000000000ull;
            }
        }        
        
        Float64PosExp2Int64 a1(sa1[0], ea1[0]);
        Float64PosExp2Int64 a2(sa2[0], ea2[0]);
        Float64PosExp2Int64 a3(sa3[0], ea3[0]);
        Float64PosExp2Int64 a4(sa4[0], ea4[0]);
        for (int k = 1; k < 4; k++) {
            Float64PosExp2Int64 z1(sa1[k], ea1[k]);
            Float64PosExp2Int64 z2(sa2[k], ea2[k]);
            Float64PosExp2Int64 z3(sa3[k], ea3[k]);
            Float64PosExp2Int64 z4(sa4[k], ea4[k]);
            a1 += z1;
            a2 += z2;
            a3 += z3;
            a4 += z4;
        }
        a1 += a2 + a3 + a4;

        for (i = i; i < vector.size(); i++) {
            a1 += vector[i];
        }

        scnfcnd = a1.scnfcnd;
        exp = a1.exp;
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

