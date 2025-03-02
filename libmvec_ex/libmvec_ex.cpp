/* gcc -Wall -O3 -m64 -march=skylake libmvec_ex.c -lm */ 
/* gcc -Wall -O3 -m64 -march=native libmvec_ex.c -lm -lmvec -o libmvec_ex */
#include <immintrin.h>
#include <stdio.h>
#include <stdint.h>
#define _GNU_SOURCE      /* These two lines are optional. They are only needed to use the scalar */
#include <math.h>        /* functions sin, exp, sincos,...                                       */

__m256d _ZGVdN4v_cos(__m256d x);                    
__m256d _ZGVdN4v_exp(__m256d x);                    
__m256d _ZGVdN4v_log(__m256d x);                    
__m256d _ZGVdN4v_sin(__m256d x);                    
__m256d _ZGVdN4vv_pow(__m256d x, __m256d y); 
void    _ZGVdN4vvv_sincos(__m256d x, __m256i ptrs, __m256i ptrc);

__m256  _ZGVdN8v_cosf(__m256 x);            
__m256  _ZGVdN8v_expf(__m256 x);            
__m256  _ZGVdN8v_logf(__m256 x);            
__m256  _ZGVdN8v_sinf(__m256 x);            
__m256  _ZGVdN8vv_powf(__m256 x, __m256 y); 
void    _ZGVdN8vvv_sincosf(__m256 x, __m256i ptrs_lo, __m256i ptrs_hi, __m256i ptrc_lo, __m256i ptrc_hi);

int mm256_print_pd(__m256d x);
int mm256_print_ps(__m256 x);

int main()
{
  __m256d x, y, z;
  __m256 xf, yf, zf;
  double z_s[4], z_c[4];          /* sincos output                              */
  float zf_s[8], zf_c[8];         /* sincosf output                             */
  __m256i ptrs, ptrc;             /* Pointers to the elements of z_s and z_c    */ 
  __m256i ptrs_lo, ptrs_hi, ptrc_lo, ptrc_hi;

  x = _mm256_set_pd (0.04, 0.03, 0.02, 0.01);
  y = _mm256_set_pd (2.0, 1.5, 1.0, 0.5);
  xf = _mm256_set_ps (0.08f, 0.07f, 0.06f, 0.05f, 0.04f, 0.03f, 0.02f, 0.01f);
  yf = _mm256_set_ps (4.0f, 3.5f, 3.0f, 2.5f, 2.0f, 1.5f, 1.0f, 0.5f);

  printf("AVX2 Double precision examples\n");
                             printf("x             "); mm256_print_pd(x);
                             printf("y             "); mm256_print_pd(y);
  z =_ZGVdN4v_cos(x);        printf("cos(x)        "); mm256_print_pd(z);
  z =_ZGVdN4v_exp(x);        printf("exp(x)        "); mm256_print_pd(z);
  z =_ZGVdN4v_log(x);        printf("log(x)        "); mm256_print_pd(z);
  z =_ZGVdN4v_sin(x);        printf("sin(x)        "); mm256_print_pd(z);
  z =_ZGVdN4vv_pow(x, y);    printf("pow(x,y)      "); mm256_print_pd(z);

  ptrs = _mm256_set_epi64x((uint64_t)&z_s[3],(uint64_t)&z_s[2],(uint64_t)&z_s[1],(uint64_t)&z_s[0]);
  ptrc = _mm256_set_epi64x((uint64_t)&z_c[3],(uint64_t)&z_c[2],(uint64_t)&z_c[1],(uint64_t)&z_c[0]);
/* Alternative: ptrs = _mm256_add_epi64(_mm256_set1_epi64x((uint64_t)&z_s[0]),_mm256_set_epi64x(24,16,8,0));  */
/* This might be more efficient if the destination addresses are contiguous in memory.                        */
  _ZGVdN4vvv_sincos(x, ptrs, ptrc);                                 /* The results of _ZGVdN4vvv_sincos are scattered into the adresses in ptrs and ptrc */ 
  printf("sincos cos(x) %12.8f %12.8f %12.8f %12.8f  \n", z_c[3], z_c[2], z_c[1], z_c[0]);
  printf("sincos sin(x) %12.8f %12.8f %12.8f %12.8f\n\n", z_s[3], z_s[2], z_s[1], z_s[0]);

  printf("AVX2 Single precision examples\n");
                             printf("x             "); mm256_print_ps(xf);
                             printf("y             "); mm256_print_ps(yf);
  zf =_ZGVdN8v_cosf(xf);     printf("cosf(x)       "); mm256_print_ps(zf);
  zf =_ZGVdN8v_expf(xf);     printf("expf(x)       "); mm256_print_ps(zf);
  zf =_ZGVdN8v_logf(xf);     printf("logf(x)       "); mm256_print_ps(zf);
  zf =_ZGVdN8v_sinf(xf);     printf("sinf(x)       "); mm256_print_ps(zf);
  zf =_ZGVdN8vv_powf(xf, yf);printf("powf(x,y)     "); mm256_print_ps(zf);

  ptrs_lo = _mm256_set_epi64x((uint64_t)&zf_s[3],(uint64_t)&zf_s[2],(uint64_t)&zf_s[1],(uint64_t)&zf_s[0]);
  ptrs_hi = _mm256_set_epi64x((uint64_t)&zf_s[7],(uint64_t)&zf_s[6],(uint64_t)&zf_s[5],(uint64_t)&zf_s[4]);
  ptrc_lo = _mm256_set_epi64x((uint64_t)&zf_c[3],(uint64_t)&zf_c[2],(uint64_t)&zf_c[1],(uint64_t)&zf_c[0]);
  ptrc_hi = _mm256_set_epi64x((uint64_t)&zf_c[7],(uint64_t)&zf_c[6],(uint64_t)&zf_c[5],(uint64_t)&zf_c[4]);
  _ZGVdN8vvv_sincosf(xf, ptrs_lo, ptrs_hi, ptrc_lo, ptrc_hi);       /* The results of _ZGVdN8vvv_sincosf are scattered to the adresses in ptrs_lo, ptrs_hi, ptrc_lo, and ptrc_hi */ 
  printf("sincosf cos(x)%12.8f %12.8f %12.8f %12.8f %12.8f %12.8f %12.8f %12.8f \n", 
           zf_c[7], zf_c[6], zf_c[5], zf_c[4], zf_c[3], zf_c[2], zf_c[1], zf_c[0]);
  printf("sincosf sin(x)%12.8f %12.8f %12.8f %12.8f %12.8f %12.8f %12.8f %12.8f  \n",
           zf_s[7], zf_s[6], zf_s[5], zf_s[4], zf_s[3], zf_s[2], zf_s[1], zf_s[0]);

  return 0;
}


__attribute__ ((noinline)) int mm256_print_pd(__m256d x){
    double vec_x[4];
    _mm256_storeu_pd(vec_x,x);
    printf("%12.8f %12.8f %12.8f %12.8f  \n", vec_x[3], vec_x[2], vec_x[1], vec_x[0]);
    return 0;
}


__attribute__ ((noinline)) int mm256_print_ps(__m256 x){
    float vec_x[8];
    _mm256_storeu_ps(vec_x,x);
    printf("%12.8f %12.8f %12.8f %12.8f %12.8f %12.8f %12.8f %12.8f \n", vec_x[7], vec_x[6], vec_x[5], vec_x[4],
                     vec_x[3], vec_x[2], vec_x[1], vec_x[0]);
    return 0;
}
