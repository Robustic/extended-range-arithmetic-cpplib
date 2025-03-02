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
