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
  __m256d maxVec = _mm256_set1_pd(-INFINITY);

  // Compute max_x using AVX2
  for (size_t i = 0; i < size; i += 4) {
      __m256d data = _mm256_loadu_pd(&x[i]);
      maxVec = _mm256_max_pd(maxVec, data);
  }

  double max_x[4];
  _mm256_storeu_pd(max_x, maxVec);
  double max_val = fmax(fmax(max_x[0], max_x[1]), fmax(max_x[2], max_x[3]));

  __m256d max = _mm256_set1_pd(max_val);

  // Compute sum(exp(x - max_x))
  __m256d sumVec = _mm256_setzero_pd();
  for (size_t i = 0; i < size; i += 4) {
      __m256d data = _mm256_loadu_pd(&x[i]);
      __m256d expData = _ZGVdN4v_exp(_mm256_sub_pd(data, max));
      sumVec = _mm256_add_pd(sumVec, expData);
  }

  double sumExp[4];
  _mm256_storeu_pd(sumExp, sumVec);
  double sum = sumExp[0] + sumExp[1] + sumExp[2] + sumExp[3];

  return max_val + log(sum);
}
