/* When main: */
/* gcc -Wall -O3 -m64 -march=native libmvec_ex.c -lm -lmvec -o libmvec_ex */
/* When not main: */
/* gcc -shared -fPIC -Wall -O3 -m64 -march=native libmvec_ex.c -o libmvec_ex.so -lm -lmvec */

#include <immintrin.h>
#include <stdio.h>
#include <stdint.h>
#define _GNU_SOURCE
#include <math.h>
                  
__m256d _ZGVdN4v_exp(__m256d x);        
__m256  _ZGVdN8v_logf(__m256 x);

double logsumexp_avx2(const double* x, size_t size) {
  __m256d maxVec1 = _mm256_set1_pd(-INFINITY);
  __m256d maxVec2 = _mm256_set1_pd(-INFINITY);

  // Compute max_x using AVX2
  for (size_t i = 0; i < size; i += 8) {
      __m256d data1 = _mm256_loadu_pd(&x[i]);
      __m256d data2 = _mm256_loadu_pd(&x[i+4]);
      maxVec1 = _mm256_max_pd(maxVec1, data1);
      maxVec2 = _mm256_max_pd(maxVec2, data2);
  }

  double max_x1[4];
  double max_x2[4];
  _mm256_storeu_pd(max_x1, maxVec1);
  _mm256_storeu_pd(max_x2, maxVec2);
  double max_val = fmax(fmax(fmax(max_x1[0], max_x1[1]), fmax(max_x1[2], max_x1[3])),
                        fmax(fmax(max_x2[0], max_x2[1]), fmax(max_x2[2], max_x2[3])));

  __m256d max = _mm256_set1_pd(max_val);

  // Compute sum(exp(x - max_x))
  __m256d sumVec1 = _mm256_setzero_pd();
  __m256d sumVec2 = _mm256_setzero_pd();
  for (size_t i = 0; i < size; i += 8) {
      __m256d data1 = _mm256_loadu_pd(&x[i]);
      __m256d data2 = _mm256_loadu_pd(&x[i+4]);
      __m256d expData1 = _ZGVdN4v_exp(_mm256_sub_pd(data1, max));
      __m256d expData2 = _ZGVdN4v_exp(_mm256_sub_pd(data2, max));
      sumVec1 = _mm256_add_pd(sumVec1, expData1);
      sumVec2 = _mm256_add_pd(sumVec2, expData2);
  }

  double sumExp[4];
  _mm256_storeu_pd(sumExp, sumVec1);
  double sum = sumExp[0] + sumExp[1] + sumExp[2] + sumExp[3];
  _mm256_storeu_pd(sumExp, sumVec2);
  sum += sumExp[0] + sumExp[1] + sumExp[2] + sumExp[3];

  return max_val + log(sum);
}

double logsumexp_Taylor(const double* x, size_t size, __m256d* data1, __m256d* data2, __m256d* data3, __m256d* data4) {
  __m256d maxVec1 = _mm256_set1_pd(-INFINITY);
  __m256d maxVec2 = _mm256_set1_pd(-INFINITY);
  __m256d maxVec3 = _mm256_set1_pd(-INFINITY);
  __m256d maxVec4 = _mm256_set1_pd(-INFINITY);

  for (size_t i = 0; i < 128; i++) {
      data1[i] = _mm256_loadu_pd(&x[16*i]);
      data2[i] = _mm256_loadu_pd(&x[16*i+4]);
      data3[i] = _mm256_loadu_pd(&x[16*i+8]);
      data4[i] = _mm256_loadu_pd(&x[16*i+12]);
      maxVec1 = _mm256_max_pd(maxVec1, data1[i]);
      maxVec2 = _mm256_max_pd(maxVec2, data2[i]);
      maxVec3 = _mm256_max_pd(maxVec3, data3[i]);
      maxVec4 = _mm256_max_pd(maxVec4, data4[i]);
  }

  double max_x1[4];
  double max_x2[4];
  double max_x3[4];
  double max_x4[4];
  _mm256_storeu_pd(max_x1, maxVec1);
  _mm256_storeu_pd(max_x2, maxVec2);
  _mm256_storeu_pd(max_x3, maxVec3);
  _mm256_storeu_pd(max_x4, maxVec4);
  double max_val = fmax(fmax(fmax(fmax(max_x1[0], max_x1[1]), fmax(max_x1[2], max_x1[3])),
                              fmax(fmax(max_x2[0], max_x2[1]), fmax(max_x2[2], max_x2[3]))),
                      fmax(fmax(fmax(max_x3[0], max_x3[1]), fmax(max_x3[2], max_x3[3])),
                              fmax(fmax(max_x4[0], max_x4[1]), fmax(max_x4[2], max_x4[3]))));

  __m256d max = _mm256_set1_pd(max_val);

  __m256d sumVec1 = _mm256_setzero_pd();
  __m256d sumVec2 = _mm256_setzero_pd();
  __m256d sumVec3 = _mm256_setzero_pd();
  __m256d sumVec4 = _mm256_setzero_pd();

  __m256d one = _mm256_set1_pd(1.0);
  __m256d half = _mm256_set1_pd(0.5);
  __m256d sixth = _mm256_set1_pd(1.0 / 6.0);
  __m256d th24 = _mm256_set1_pd(1.0 / 24.0);
  __m256d th120 = _mm256_set1_pd(1.0 / 120.0);
  __m256d th720 = _mm256_set1_pd(1.0 / 720.0);

  for (size_t i = 0; i < 128; i++) {
      __m256d x1 = _mm256_sub_pd(data1[i], max);
      __m256d x2 = _mm256_sub_pd(data2[i], max);
      __m256d x3 = _mm256_sub_pd(data3[i], max);
      __m256d x4 = _mm256_sub_pd(data4[i], max);

      __m256d x1_2 = _mm256_mul_pd(x1, x1);
      __m256d x2_2 = _mm256_mul_pd(x2, x2);
      __m256d x3_2 = _mm256_mul_pd(x3, x3);
      __m256d x4_2 = _mm256_mul_pd(x4, x4);

      __m256d x1_3 = _mm256_mul_pd(x1, x1_2);
      __m256d x2_3 = _mm256_mul_pd(x2, x2_2);
      __m256d x3_3 = _mm256_mul_pd(x3, x3_2);
      __m256d x4_3 = _mm256_mul_pd(x4, x4_2);

      __m256d x1_4 = _mm256_mul_pd(x1, x1_3);
      __m256d x2_4 = _mm256_mul_pd(x2, x2_3);
      __m256d x3_4 = _mm256_mul_pd(x3, x3_3);
      __m256d x4_4 = _mm256_mul_pd(x4, x4_3);

      __m256d x1_5 = _mm256_mul_pd(x1, x1_4);
      __m256d x2_5 = _mm256_mul_pd(x2, x2_4);
      __m256d x3_5 = _mm256_mul_pd(x3, x3_4);
      __m256d x4_5 = _mm256_mul_pd(x4, x4_4);

      __m256d x1_6 = _mm256_mul_pd(x1, x1_5);
      __m256d x2_6 = _mm256_mul_pd(x2, x2_5);
      __m256d x3_6 = _mm256_mul_pd(x3, x3_5);
      __m256d x4_6 = _mm256_mul_pd(x4, x4_5);

      __m256d res = _mm256_add_pd(x1, one);
      __m256d res = _mm256_add_pd(x2, one);
      __m256d res = _mm256_add_pd(x3, one);
      __m256d res = _mm256_add_pd(x4, one);

      x1 = _mm256_fmadd_pd(x1_2, half, x1);
      x2 = _mm256_fmadd_pd(x2_2, half, x2);
      x3 = _mm256_fmadd_pd(x3_2, half, x3);
      x4 = _mm256_fmadd_pd(x4_2, half, x4);

      x1 = _mm256_fmadd_pd(x1_3, sixth, x1);
      x2 = _mm256_fmadd_pd(x2_3, sixth, x2);
      x3 = _mm256_fmadd_pd(x3_3, sixth, x3);
      x4 = _mm256_fmadd_pd(x4_3, sixth, x4);

      x1 = _mm256_fmadd_pd(x1_4, th24, x1);
      x2 = _mm256_fmadd_pd(x2_4, th24, x2);
      x3 = _mm256_fmadd_pd(x3_4, th24, x3);
      x4 = _mm256_fmadd_pd(x4_4, th24, x4);

      x1 = _mm256_fmadd_pd(x1_5, th120, x1);
      x2 = _mm256_fmadd_pd(x2_5, th120, x2);
      x3 = _mm256_fmadd_pd(x3_5, th120, x3);
      x4 = _mm256_fmadd_pd(x4_5, th120, x4);

      x1 = _mm256_fmadd_pd(x1_6, th720, x1);
      x2 = _mm256_fmadd_pd(x2_6, th720, x2);
      x3 = _mm256_fmadd_pd(x3_6, th720, x3);
      x4 = _mm256_fmadd_pd(x4_6, th720, x4);

      sumVec1 = _mm256_add_pd(sumVec1, x1);
      sumVec2 = _mm256_add_pd(sumVec2, x2);
      sumVec3 = _mm256_add_pd(sumVec3, x3);
      sumVec4 = _mm256_add_pd(sumVec4, x4);
  }

  double sumExp[4];
  _mm256_storeu_pd(sumExp, sumVec1);
  double sum = sumExp[0] + sumExp[1] + sumExp[2] + sumExp[3];
  _mm256_storeu_pd(sumExp, sumVec2);
  sum += sumExp[0] + sumExp[1] + sumExp[2] + sumExp[3];
  _mm256_storeu_pd(sumExp, sumVec3);
  sum += sumExp[0] + sumExp[1] + sumExp[2] + sumExp[3];
  _mm256_storeu_pd(sumExp, sumVec4);
  sum += sumExp[0] + sumExp[1] + sumExp[2] + sumExp[3];

  return max_val + log(sum);
}
