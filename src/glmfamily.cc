#include "glmfamily.h"
#include <vector>
#include <algorithm>
#include <stdexcept>

#include <math.h>
#define _USE_MATH_DEFINES
#include <cmath>

GlmFamily::~GlmFamily() {}
Gaussian::Gaussian() {}

// Copy these vectors, though maybe only need to do it once?
void Gaussian::get_workingset(const double *eta, const double *y,
                              const double *v, double *w, double *z, int len, double* sumbuf) {
    sumbuf[0] = 0;
    sumbuf[1] = 0;
    for (int i = 0; i < len; ++i) {
        w[i] = v[i];
        z[i] = v[i] * (y[i] - eta[i]);
        sumbuf[0] += z[i];
        sumbuf[1] += w[i];
    }
}

double Gaussian::get_deviance(const double *y, const double *eta,
                              const double *v, int len) {
    double result = 0.0;
    for (int i = 0; i < len; ++i) {
        result += (y[i] - eta[i]) * (y[i] - eta[i]) * v[i];
    }
    return result;
}

void Gaussian::get_residual(const double *y, const double *eta, const double *v,
                            double *r, int len) {
    for (int i = 0; i < len; ++i) {
        r[i] = v[i] * (y[i] - eta[i]);
    }
}

double Gaussian::null_deviance(const double *y, const double *v, double *r,
                               int intr, double *eta, bool has_offset,
                               const double *offset, double *aint, int len) {
    if (!has_offset) {
        double weightysquare = 0;
        double weighty = 0;
        double weightsum = 0;
        for (int i = 0; i < len; ++i) {
            weightysquare += v[i] * y[i] * y[i];
            weighty += v[i] * y[i];
            weightsum += v[i];
        }
        if (!intr) {
            for (int i = 0; i < len; ++i) {
                eta[i] = 0;
                r[i] = v[i] * y[i];
            }
            return weightysquare;
        }

        double beta0 = weighty / weightsum;
        for (int i = 0; i < len; ++i) {
            eta[i] = beta0;
            r[i] = v[i] * (y[i] - beta0);
        }
        *aint = beta0;
        return (weightysquare - 2 * weighty * beta0 +
                weightsum * beta0 * beta0);
    }

    double weightysquare = 0;
    double weighty = 0;
    double weightsum = 0;
    for (int i = 0; i < len; ++i) {
        double y0 = y[i] - offset[i];
        weightysquare += v[i] * y0 * y0;
        weighty += v[i] * y0;
        weightsum += v[i];
    }
    if (!intr) {
        for (int i = 0; i < len; ++i) {
            eta[i] = offset[i];
            r[i] = v[i] * (y[i] - offset[i]);
        }
        return weightysquare;
    }
    double beta0 = weighty / weightsum;
    for (int i = 0; i < len; ++i) {
        eta[i] = offset[i] + beta0;
        r[i] = v[i] * (y[i] - beta0 - offset[i]);
    }
    return (weightysquare - 2 * weighty * beta0 + weightsum * beta0 * beta0);
}

double Gaussian::get_obj(const double *y, const double *v, const double *eta, const double *a,
        const int *ju, const double *vp, int no, int ni, double lambda){
  return 0.0;
};

Logistic::Logistic() {}
void Logistic::get_workingset(const double *eta, const double *y,
                              const double *v, double *w, double *z, int len, double* sumbuf) {
    sumbuf[0] = 0;
    sumbuf[1] = 0;
    for (int i = 0; i < len; ++i) {
        double p = 1 / (1 + exp(-eta[i]));
        w[i] = p * (1 - p) * v[i];
        z[i] = (y[i] - p) * v[i];
        sumbuf[0] += z[i];
        sumbuf[1] += w[i];
    }
}

double Logistic::get_deviance(const double *y, const double *eta,
                              const double *v, int len) {
    double result = 0;
    for (int i = 0; i < len; ++i) {
        result += v[i] * (log(1 + exp(eta[i])) - y[i] * eta[i]);
    }
    result *= 2;
    return result;
}

void Logistic::get_residual(const double *y, const double *eta, const double *v,
                            double *r, int len) {
    for (int i = 0; i < len; ++i) {
        double eeta = exp(eta[i]);
        r[i] = v[i] * (y[i] - eeta / (1 + eeta));
    }
}

double Logistic::null_deviance(const double *y, const double *v, double *r,
                               int intr, double *eta, bool has_offset,
                               const double *offset, double *aint, int len) {
    if (!has_offset) {
        double count0 = 0;
        double count1 = 0;
        for (int i = 0; i < len; ++i) {
            count0 += v[i] * y[i];
            count1 += v[i] * (1 - y[i]);
        }
        if (!intr) {
            for (int i = 0; i < len; ++i) {
                eta[i] = 0;
                r[i] = v[i] * (y[i] - 0.5);
            }
            return (2 * log(2) * (count0 + count1));
        }
        double p0inv = 1 + (count1 / count0);
        double p0 = 1/p0inv;
        double aint_val = -log(p0inv - 1);
        *aint = aint_val;
        for (int i = 0; i < len; ++i) {
            eta[i] = aint_val;
            r[i] = v[i] * (y[i] - p0);
        }
        return (2 * log(p0inv) * (count0 + count1));
    }

    throw "Offset has not been implemented for logstic regression\n";
}

double Logistic::get_obj(const double *y, const double *v, const double *eta, const double *a,
                         const int *ju, const double *vp, int no, int ni, double lambda){
  return 0.0;
};

double normal_cdf(double x, double h, double mu = 0.0) {
  // Handle degenerate case (h = 0)
  if (h <= 0.0) {
    throw "bandwidth could not be negative\n";
    }

  // Standardize x: adjust for mean and scale by standard deviation
  double z = (x - mu) / (h * std::sqrt(2.0));

  // Compute CDF using the error function (erf)
  return 0.5 * (1.0 + std::erf(z));
}

void Quantile::get_workingset(const double *eta, const double *y,
                              const double *v, double *w, double *z, int len, double* sumbuf) {
  sumbuf[0] = 0;
  sumbuf[1] = 0;
  double ch = 1 / (std::sqrt(2 * M_PI) * this->h);
  for (int i = 0; i < len; ++i) {
    double res = y[i] - eta[i];
    w[i] = ch * v[i];
    z[i] = (this->tau - normal_cdf(-res, h))  * v[i];
    sumbuf[0] += z[i];
    sumbuf[1] += w[i];
  }
}

double Quantile::get_deviance(const double *y, const double *eta,
                              const double *v, int len) {
  double result = 0;
  double shift = this->tau - 0.5;
  for (int i = 0; i < len; ++i) {
    double res = y[i] - eta[i];
    result += (0.5*fabs(res)+shift*res)*v[i];
  }
  return result;
}

void Quantile::get_residual(const double *y, const double *eta, const double *v,
                            double *r, int len) {
  for (int i = 0; i < len; ++i) {
    r[i] = v[i] * (y[i] - eta[i]);
  }
}

double quantile(const double* y, std::size_t len, double tau) {
  if (len == 0) throw std::invalid_argument("Array length must be > 0");
  if (tau < 0.0 || tau > 1.0) throw std::invalid_argument("tau must be between 0 and 1");

  // Copy to vector to sort
  std::vector<double> data(y, y + len);
  std::sort(data.begin(), data.end());

  // Position in sorted array
  double idx = tau * (len - 1);
  std::size_t lower = static_cast<std::size_t>(idx);
  double frac = idx - lower;

  // Linear interpolation if needed
  if (lower + 1 < len) {
    return data[lower] * (1 - frac) + data[lower + 1] * frac;
  } else {
    return data[lower];  // tau == 1 case
  }
}

double Quantile::null_deviance(const double *y, const double *v, double *r,
                               int intr, double *eta, bool has_offset,
                               const double *offset, double *aint, int len) {
  if (!has_offset) {
    if (!intr) {
      for (int i = 0; i < len; ++i) {
        eta[i] = 0;
        r[i] = v[i] * (normal_cdf(-y[i],this->h)- this->tau);
      }
    }else{
      double eta1 = quantile(y,len,this->tau);
      *aint = eta1;
      for (int i = 0; i < len; ++i) {
        eta[i] = eta1;
        r[i] = v[i] * (normal_cdf(eta1-y[i],this->h)- this->tau);
      }
    }
    return this->get_deviance(y, eta, v, len);
  }
  throw "Offset has not been implemented for quantile regression\n";
}

double lGu(double u){
  return std::sqrt(2/M_PI)*std::exp(-0.5*u*u) + u*(1-2*normal_cdf(-u,1.0));
}

double Quantile::get_obj(const double *y, const double *v, const double *eta, const double *a,
                         const int *ju, const double *vp, int no, int ni, double lambda){
  double h = this->h;
  double tau = this->tau;
  double result = 0;
  double res;
  for(int i=0;i<no;i++){
    res = y[i] - eta[i];
    result += v[i]* (0.5*h*lGu(res/h)+(tau-0.5)*res);
  }
  for(int j=0;j<ni;j++){
    if(!ju[j]){
      continue;
    }
    result += lambda*fabs(a[j])*vp[j];
  }
  return result;
}
