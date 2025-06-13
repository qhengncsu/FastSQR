#include "glmnetMatrix.h"

#include <math.h>

// This will probably be fixed in a new version of Eigen, for now
#ifndef __clang__
#if __GNUC__ < 8
#define _mm256_set_m128d(vh, vl) \
    _mm256_insertf128_pd(_mm256_castpd128_pd256(vl), (vh), 1)
#endif
#endif

#include <Eigen/Core>

MatrixGlmnet::~MatrixGlmnet() {}

int MatrixGlmnet::get_ni() { return this->ni; }
int MatrixGlmnet::get_no() { return this->no; }

double MatrixGlmnet::max_grad(const double *r, const int *ju,
                              const double *vp) {
    double result = 0.0;
    for (int i = 0; i < ni; ++i) {
        if ((!ju[i]) || (vp[i] <= 0.0)) {
            continue;
        }
        result = fmax(result, fabs(this->dot_product(i, r)) / vp[i]);
    }
    return result;
}

double MatrixGlmnet::sumv(const double *v, int len) {
    Eigen::Map<const Eigen::ArrayXd> vmap(v, len);
    return vmap.sum();
}

void MatrixGlmnet::update_res_eigen(double* r, const double *v, double d, int len) {

    Eigen::Map<Eigen::ArrayXd> rmap(r, len);
    Eigen::Map<const Eigen::ArrayXd> vmap(v, len);
    rmap += d * vmap;
}

DenseM::DenseM(int no, int ni, const double *x) {
    this->no = no;
    this->ni = ni;
    data = x;
}
DenseM::~DenseM() { data = nullptr; }

double DenseM::dot_product(int j, const double *v, double vsum) {
    Eigen::Map<const Eigen::VectorXd> x(data + j * no, no);
    Eigen::Map<const Eigen::VectorXd> y(v, no);
    return x.dot(y);
}

double DenseM::vx2(int j, const double *v, double vsum, double *xm) {
    Eigen::Map<const Eigen::ArrayXd> x(data + j * no, no);
    Eigen::Map<const Eigen::ArrayXd> y(v, no);
    return (y * x.square()).sum();
}

void DenseM::update_res(int j, double d, const double *v,
                        double *__restrict r, double *rsum, double vsum, double vx) {
    for (int i = 0; i < no; ++i) {
        r[i] -= d * v[i] * data[j * no + i];
        (*rsum) -= d * v[i] * data[j * no + i];
    }
}

void DenseM::compute_eta(double *__restrict eta, const double *a, double aint,
                         bool has_offset, const double *offset) {
    Eigen::Map<Eigen::VectorXd> eta_map(eta, no);
    Eigen::Map<const Eigen::MatrixXd> X(data, no, ni);
    Eigen::Map<const Eigen::VectorXd> a_map(a, ni);

    eta_map = (X * a_map).array() + aint;
    if (has_offset) {
        Eigen::Map<const Eigen::VectorXd> offset_map(offset, no);
        eta_map += offset_map;
    }
}

void eigen_get_eta(double *eta, const double *cov, const double *weights, int no, int ncov)
{
    Eigen::Map<Eigen::VectorXd> eta_map(eta, no);
    Eigen::Map<const Eigen::MatrixXd> X(cov, no, ncov);
    Eigen::Map<const Eigen::VectorXd> y(weights, ncov);
    eta_map += X * y;
    return;
}