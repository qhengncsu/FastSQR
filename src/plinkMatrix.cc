#include <math.h>
#include "glmnetMatrix.h"
#include "omp.h"
#include "R_ext/Print.h"
// Computation heavily adapted from
// https://github.com/chrchang/plink-ng/blob/master/2.0/pgenlib_ffi_support.cc
inline uintptr_t DivUp(uintptr_t val, uint32_t divisor)
{
    return (val + divisor - 1) / divisor;
}
static constexpr uintptr_t kMask5555 = (~static_cast<uintptr_t>(0)) / 3;

static constexpr int32_t kBitsPerWord = 64;
static constexpr int32_t kBitsPerWordD2 = kBitsPerWord / 2;
static constexpr int32_t kCacheline = 64;
static constexpr int32_t kBytesPerWord = kBitsPerWord / 8;
static constexpr int32_t kWordsPerCacheline = kCacheline / kBytesPerWord;
static constexpr int32_t kNypsPerCacheline = kCacheline * 4;

#ifdef _WIN64
inline uint32_t ctzw(unsigned long long ullii)
{
    return __builtin_ctzll(ullii);
}
#else
inline uint32_t ctzw(unsigned long ulii)
{
    return __builtin_ctzl(ulii);
}
#endif

static constexpr uintptr_t k1LU = static_cast<uintptr_t>(1);

inline uintptr_t bzhi(uintptr_t ww, uint32_t idx) {
  return ww & ((k1LU << idx) - k1LU);
}

inline void ZeroTrailingBits(uintptr_t bit_ct, uintptr_t* bitarr) {
  const uint32_t trail_ct = bit_ct % kBitsPerWord;
  if (trail_ct) {
    bitarr[bit_ct / kBitsPerWord] = bzhi(bitarr[bit_ct / kBitsPerWord], trail_ct);
  }
}


static inline void ZeroTrailingNyps(uintptr_t nyp_ct, uintptr_t* bitarr) {
  ZeroTrailingBits(nyp_ct * 2, bitarr);
}

void eigen_get_eta(double *eta, const double *cov, const double *weights, int no, int ncov);


void GetWeightsByValueNoDosage(const double* weights, const uintptr_t* genoarr,
                               uint32_t sample_ct, double* buf) {
    const uint32_t word_ct = DivUp(sample_ct, kBitsPerWordD2);
    double result = 0.0;
    double result2 = 0.0;
    double miss_weight = 0.0;
#pragma omp parallel
    {
        int total_threads = omp_get_num_threads();
        int threadid = omp_get_thread_num();
        int size = (word_ct + total_threads * 8 - 1) / (total_threads * 8);
        size *= 8;
        uint32_t start = threadid * size;
        uint32_t end = (1 + threadid) * size;
        double result_local = 0.0;
        double result2_local = 0.0;
        double miss_weight_local = 0.0;
        if (end > word_ct)
        {
            end = word_ct;
        }
        for (uint32_t widx = start; widx < end; ++widx)
        {
            const uintptr_t geno_word = genoarr[widx];
            if (!geno_word){
                continue;
            }
            const double *cur_weights = &(weights[widx * kBitsPerWordD2]);
            uintptr_t geno_word1 = geno_word & kMask5555;
            uintptr_t geno_word2 = (geno_word >> 1) & kMask5555;
            uintptr_t geno_missing_word = geno_word1 & geno_word2;
            geno_word1 ^= geno_missing_word;
            while (geno_word1)
            {
                const uint32_t sample_idx_lowbits = ctzw(geno_word1) / 2;
                result_local += cur_weights[sample_idx_lowbits];
                geno_word1 &= geno_word1 - 1;
            }
            geno_word2 ^= geno_missing_word;
            while (geno_word2)
            {
                const uint32_t sample_idx_lowbits = ctzw(geno_word2) / 2;
                result2_local += cur_weights[sample_idx_lowbits];
                geno_word2 &= geno_word2 - 1;
            }
            while (geno_missing_word)
            {
                const uint32_t sample_idx_lowbits = ctzw(geno_missing_word) / 2;
                miss_weight_local += cur_weights[sample_idx_lowbits];
                geno_missing_word &= geno_missing_word - 1;
            }
        }
        #pragma omp critical
        {
            result += result_local;
            result2 += result2_local;
            miss_weight += miss_weight_local;
        }
    }

    buf[0] = result;
    buf[1] = result2;
    buf[2] = miss_weight;
}

PlinkMatrix::PlinkMatrix(int no, int ni, const uintptr_t* x,
                         const double* xim, int intr, int ncov, const double *cov) {
    this->no = no;
    this->ni = ni;
    data = x;
    this->xim = xim;
    this->center = (bool)intr;
    const uint32_t cache_line_ct = DivUp(no, kNypsPerCacheline);
    word_ct = kWordsPerCacheline * cache_line_ct;

    this->ncov = ncov;
    this->cov = cov;

}
PlinkMatrix::~PlinkMatrix()
{
    data = nullptr;
    xim = nullptr;
    cov = nullptr;
}

double PlinkMatrix::dot_product(int j, const double *v, double vsum)
{
    // assert((j > 0) && (j < ni));
    if(j >= ncov){
        double buf[3];
        GetWeightsByValueNoDosage(v, &(data[((uintptr_t)(j-ncov)) * word_ct]), no, buf);
        double result = buf[0] + 2 * buf[1] + buf[2] * xim[j - ncov];
        if(center) {
            result -= xim[j - ncov] * vsum;
        }
        return result;
    }

    double result = 0;
    for(int i = 0; i < no; ++i) {
        result += v[i] * cov[j*no + i];
    }
    return result;
}

double PlinkMatrix::vx2(int j, const double* v, double vsum, double *xm) {
    if(j >= ncov){
        double buf[3];
        GetWeightsByValueNoDosage(v, &(data[((uintptr_t)(j-ncov)) * word_ct]), no, buf);
        double result = buf[0] + 4 * buf[1] + buf[2] * xim[j - ncov] * xim[j - ncov];
        if(center) {
            double l1 = (buf[0] + 2 * buf[1] + buf[2] * xim[j - ncov]) * 2 * xim[j - ncov];
            double l2 = xim[j - ncov] * xim[j - ncov] *vsum;
            result += (l2 - l1);
            (*xm) = buf[0] + 2 * buf[1] + buf[2] * xim[j - ncov];
        }
        return result;
    }

    double result = 0;
    for(int i = 0; i < no; ++i) {
        result += v[i] * cov[j*no + i] * cov[j*no + i];
    }
    (*xm) = 0;
    return result;

}

void PlinkMatrix::update_res(int j, double d, const double *weights,
                             double *r, double *rsum, double vsum, double vx)
{
    if(j >= ncov){
        const uintptr_t *genoarr = &(data[((uintptr_t)(j-ncov)) * word_ct]);

        const uint32_t word_ct_local = DivUp(no, kBitsPerWordD2);
#pragma omp parallel
    {
        int total_threads = omp_get_num_threads();
        int threadid = omp_get_thread_num();
        int size = (word_ct_local + total_threads * 8 - 1) / (total_threads * 8);
        size *= 8;
        uint32_t start = threadid * size;
        uint32_t end = (1 + threadid) * size;
        if (end > word_ct_local)
        {
            end = word_ct_local;
        }
        for (uint32_t widx = start; widx < end; ++widx)
        {
            const uintptr_t geno_word = genoarr[widx];
            if (!geno_word)
            {
                continue;
            }
            const double *cur_weights = &(weights[widx * kBitsPerWordD2]);
            uintptr_t geno_word1 = geno_word & kMask5555;
            uintptr_t geno_word2 = (geno_word >> 1) & kMask5555;
            uintptr_t geno_missing_word = geno_word1 & geno_word2;
            geno_word1 ^= geno_missing_word;
            while (geno_word1)
            {
                const uint32_t sample_idx_lowbits = ctzw(geno_word1) / 2;
                r[widx * kBitsPerWordD2 + sample_idx_lowbits] -=
                    d * cur_weights[sample_idx_lowbits];
                geno_word1 &= geno_word1 - 1;
            }
            geno_word2 ^= geno_missing_word;
            while (geno_word2)
            {
                const uint32_t sample_idx_lowbits = ctzw(geno_word2) / 2;
                r[widx * kBitsPerWordD2 + sample_idx_lowbits] -=
                    2 * d * cur_weights[sample_idx_lowbits];
                geno_word2 &= geno_word2 - 1;
            }
            while (geno_missing_word)
            {
                const uint32_t sample_idx_lowbits = ctzw(geno_missing_word) / 2;
                r[widx * kBitsPerWordD2 + sample_idx_lowbits] -=
                    xim[j - ncov] * d * cur_weights[sample_idx_lowbits];
                geno_missing_word &= geno_missing_word - 1;
            }
        }
    }
        if(center) {
            MatrixGlmnet::update_res_eigen(r, weights, d*xim[j - ncov], no);
            (*rsum) += d * (vsum * xim[j - ncov] - vx);
        }
        return;
    }

    for (int i = 0; i < no; ++i) {
        r[i] -= d * weights[i] * cov[j * no + i];
        (*rsum) -= d * weights[i] * cov[j * no + i];
    }
    return;
}

void PlinkMatrix::multv(double *eta, const double *weights, double aint) {
    for (int i = 0; i < no; ++i)
    {
        eta[i] = aint;
        //Rprintf("Index: %d, eta: %f\n", i, aint);
    }

    // This is only useful for small no, maybe I should just not use it at all?
    const uint32_t word_ct_local = DivUp(no, kBitsPerWordD2);
#pragma omp parallel
    {
        int total_threads = omp_get_num_threads();
        int threadid = omp_get_thread_num();
        int size = (word_ct_local + total_threads * 8 - 1) / (total_threads * 8);
        size *= 8;
        uint32_t start = threadid * size;
        uint32_t end = (1 + threadid) * size;
        if (end > word_ct_local)
        {
            end = word_ct_local;
        }
        for (int j = ncov; j < ni; ++j)
        {
            const uintptr_t *genoarr = &(data[((uintptr_t)(j-ncov)) * word_ct]);
            // it's easy to forget that xim should be extended when
            // we add covariates to it
            double ximpute = xim[j - ncov];
            double wj = weights[j];

            for (uint32_t widx = start; widx < end; ++widx)
            {
                const uintptr_t geno_word = genoarr[widx];
                if (!geno_word)
                {
                    continue;
                }
                uintptr_t geno_word1 = geno_word & kMask5555;
                uintptr_t geno_word2 = (geno_word >> 1) & kMask5555;
                uintptr_t geno_missing_word = geno_word1 & geno_word2;
                geno_word1 ^= geno_missing_word;
                while (geno_word1)
                {
                    const uint32_t sample_idx_lowbits = ctzw(geno_word1) / 2;
                    eta[widx * kBitsPerWordD2 + sample_idx_lowbits] += wj;
                    geno_word1 &= geno_word1 - 1;
                }
                geno_word2 ^= geno_missing_word;
                while (geno_word2)
                {
                    const uint32_t sample_idx_lowbits = ctzw(geno_word2) / 2;
                    eta[widx * kBitsPerWordD2 + sample_idx_lowbits] += 2 * wj;
                    geno_word2 &= geno_word2 - 1;
                }
                while (geno_missing_word)
                {
                    const uint32_t sample_idx_lowbits = ctzw(geno_missing_word) / 2;
                    eta[widx * kBitsPerWordD2 + sample_idx_lowbits] += ximpute * wj;
                    geno_missing_word &= geno_missing_word - 1;
                }
            }
        }
    }

    if(ncov > 0) {
        eigen_get_eta(eta, cov, weights, no, ncov);
    }
}

void PlinkMatrix::compute_eta(double *eta, const double *weights, double aint,
                              bool has_offset, const double *offset)
{
    multv(eta, weights, aint);
    if (has_offset)
    {
        for (int i = 0; i < no; ++i)
        {
            eta[i] += offset[i];
        }
    }
    if(center){
        double inner = 0;
        for(int j = ncov; j < ni; ++j) {
            inner += xim[j - ncov] * weights[j];
        }
        for(int i = 0; i < no; ++i) {
            eta[i] -= inner;
        }
    }
}
