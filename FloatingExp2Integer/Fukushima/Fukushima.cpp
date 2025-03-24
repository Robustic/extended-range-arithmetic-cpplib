#include <cmath>
#include <cstdint>
#include <immintrin.h>

#include "Fukushima.h"

namespace floatingExp2Integer
{
    Fukushima::Fukushima() {
        this->doubleToFukushima(1.0);
    }

    Fukushima::Fukushima(double dbl) {
        this->doubleToFukushima(dbl);
    }

    Fukushima::Fukushima(double dbl, std::int64_t exponent) {
        scnfcnd = dbl;
        exp = exponent;
    }

    void Fukushima::doubleToFukushima(double dbl) {
        scnfcnd = dbl;
        exp = 0LL;
    }

    double Fukushima::asDouble() const {
        return scnfcnd * std::exp2(920 * exp);
    }

    Fukushima::Fukushima(const std::vector<floatingExp2Integer::Fukushima>& vector) {
        const double BIGI = 0x1p-960;
        const double BIGS = 0x1p480;
        const __m256d BIGI_4 = _mm256_set1_pd(BIGI);
        const __m256d BIGS_4 = _mm256_set1_pd(BIGS);

        const __m256d zerosd_4 = _mm256_set1_pd(0.0);
        const __m256i zerosi_4 = _mm256_set1_epi64x(0LL);
        const __m256i minus1i_4 = _mm256_set1_epi64x(-1LL);
        const __m256i plus1i_4 = _mm256_set1_epi64x(1LL);

        __m256d as1;
        __m256i ae1;
        __m256d bs1;
        __m256i be1;

        __m256d as2;
        __m256i ae2;
        __m256d bs2;
        __m256i be2;

        for (int k = 0; k < 4; k++) {
            as1[k] = vector[k].scnfcnd;
            ae1[k] = vector[k].exp;
            as2[k] = vector[k + 4].scnfcnd;
            ae2[k] = vector[k + 4].exp;
        }

        int i = 8;
        for (i = 8; i + 7 < vector.size(); i+=8) {
            for (int k = 0; k < 4; k++) {
                bs1[k] = vector[i + k].scnfcnd;
                // be1[k] = vector[i + k].exp;
                bs2[k] = vector[i + k + 4].scnfcnd;
                // be2[k] = vector[i + k + 4].exp;
            }

            // __m256i exp_diff1 = ae1 - be1;
            // __m256i exp_diff2 = ae2 - be2;

            // ae1 = exp_diff1 < zerosi_4 ? be1 : ae1;
            // ae2 = exp_diff2 < zerosi_4 ? be2 : ae2;

            // bs1 = exp_diff1 == plus1i_4 ? bs1 * BIGI_4 : bs1;
            // as1 = exp_diff1 == minus1i_4 ? as1 * BIGI_4 : as1;
            // bs1 = exp_diff1 > plus1i_4 ? zerosd_4 : bs1;
            // as1 = exp_diff1 < minus1i_4 ? zerosd_4 : as1;

            // bs2 = exp_diff2 == plus1i_4 ? bs2 * BIGI_4 : bs2;
            // as2 = exp_diff2 == minus1i_4 ? as2 * BIGI_4 : as2;
            // bs2 = exp_diff2 > plus1i_4 ? zerosd_4 : bs2;
            // as2 = exp_diff2 < minus1i_4 ? zerosd_4 : as2;

    // // Integer comparisons (producing masks for blending)
    // __m256i mask_eq_plus1i_1 = _mm256_cmpeq_epi64(exp_diff1, plus1i_4);
    // __m256i mask_eq_minus1i_1 = _mm256_cmpeq_epi64(exp_diff1, minus1i_4);
    // __m256i mask_gt_plus1i_1 = _mm256_cmpgt_epi64(exp_diff1, plus1i_4);
    // __m256i mask_lt_minus1i_1 = _mm256_cmpgt_epi64(minus1i_4, exp_diff1); // (A < B) == (B > A)

    // __m256i mask_eq_plus1i_2 = _mm256_cmpeq_epi64(exp_diff2, plus1i_4);
    // __m256i mask_eq_minus1i_2 = _mm256_cmpeq_epi64(exp_diff2, minus1i_4);
    // __m256i mask_gt_plus1i_2 = _mm256_cmpgt_epi64(exp_diff2, plus1i_4);
    // __m256i mask_lt_minus1i_2 = _mm256_cmpgt_epi64(minus1i_4, exp_diff2); // (A < B) == (B > A)

    // // Convert integer masks to double masks for blending
    // __m256d mask_eq_plus1i_1_d = _mm256_castsi256_pd(mask_eq_plus1i_1);
    // __m256d mask_eq_minus1i_1_d = _mm256_castsi256_pd(mask_eq_minus1i_1);
    // __m256d mask_gt_plus1i_1_d = _mm256_castsi256_pd(mask_gt_plus1i_1);
    // __m256d mask_lt_minus1i_1_d = _mm256_castsi256_pd(mask_lt_minus1i_1);

    // __m256d mask_eq_plus1i_2_d = _mm256_castsi256_pd(mask_eq_plus1i_2);
    // __m256d mask_eq_minus1i_2_d = _mm256_castsi256_pd(mask_eq_minus1i_2);
    // __m256d mask_gt_plus1i_2_d = _mm256_castsi256_pd(mask_gt_plus1i_2);
    // __m256d mask_lt_minus1i_2_d = _mm256_castsi256_pd(mask_lt_minus1i_2);

    // // Use masked blend to avoid branches
    // bs1 = _mm256_blendv_pd(bs1, _mm256_mul_pd(bs1, BIGI_4), mask_eq_plus1i_1_d);
    // as1 = _mm256_blendv_pd(as1, _mm256_mul_pd(as1, BIGI_4), mask_eq_minus1i_1_d);
    // bs1 = _mm256_blendv_pd(bs1, zerosd_4, mask_gt_plus1i_1_d);
    // as1 = _mm256_blendv_pd(as1, zerosd_4, mask_lt_minus1i_1_d);

    // bs2 = _mm256_blendv_pd(bs2, _mm256_mul_pd(bs2, BIGI_4), mask_eq_plus1i_2_d);
    // as2 = _mm256_blendv_pd(as2, _mm256_mul_pd(as2, BIGI_4), mask_eq_minus1i_2_d);
    // bs2 = _mm256_blendv_pd(bs2, zerosd_4, mask_gt_plus1i_2_d);
    // as2 = _mm256_blendv_pd(as2, zerosd_4, mask_lt_minus1i_2_d);

            as1 += bs1;
            as2 += bs2;
        }

        as1 += as2;
        
        floatingExp2Integer::Fukushima f0(as1[0], ae1[0]);
        floatingExp2Integer::Fukushima f1(as1[1], ae1[1]);
        floatingExp2Integer::Fukushima f2(as1[2], ae1[2]);
        floatingExp2Integer::Fukushima f3(as1[3], ae1[3]);

        f0 += f1;
        f0 += f2;
        f0 += f3;

        for (i = i; i < vector.size(); i++) {
            f0 += vector[i];
        }

        scnfcnd = f0.scnfcnd;
        exp = f0.exp;
    }

    Fukushima& Fukushima::operator+=(Fukushima z) {
        //const int IND = 960;
        const double BIGI = 0x1p-960; // pow(2.0, -IND);
        const double BIGS = 0x1p480; // std::pow(2.0, IND / 2);

        if (exp == z.exp) {
            scnfcnd = scnfcnd + z.scnfcnd;
        } 
        else {
            int id = exp - z.exp;
            if (id == 1) {
                scnfcnd = scnfcnd + (z.scnfcnd * BIGI);
            } else if (id == -1) {
                scnfcnd = z.scnfcnd + (scnfcnd * BIGI);
                exp = z.exp;
            } else if (id > 1) {
                scnfcnd = scnfcnd;
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

    inline void Fukushima::xnorm() {
        //constexpr int IND = 960;
        const double BIG = 0x1p960; // std::pow(2.0, IND);
        const double BIGI = 0x1p-960; // std::pow(2.0, -IND);
        const double BIGS = 0x1p480; // std::pow(2.0, IND / 2);
        const double BIGSI = 0x1p-480; // std::pow(2.0, -(IND / 2));
    
        //double w = std::abs(scnfcnd);
    
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

