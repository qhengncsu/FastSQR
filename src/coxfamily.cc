#include "glmfamily.h"
#include <math.h>
# include<stdlib.h>

// This will probably be fixed in a new version of Eigen, for now
#ifndef __clang__
#if __GNUC__ < 8
#define _mm256_set_m128d(vh, vl) \
    _mm256_insertf128_pd(_mm256_castpd128_pd256(vl), (vh), 1)
#endif
#endif


#include <Eigen/Core>

typedef const Eigen::PermutationMatrix<Eigen::Dynamic,Eigen::Dynamic> PermMat;

Cox::Cox(const int *order, const int *rankmin, const int *rankmax, int len) : order(order), rankmin(rankmin), rankmax(rankmax) {
    rskden = (double*)malloc(sizeof(double)*len);
}

Cox::~Cox(){
    free(rskden);
    rskden=nullptr;
}

void Cox::get_workingset(const double *eta, const double *y, const double *v,
                         double *w, double *z, int len, double *sumbuf)
{
    Eigen::Map<const Eigen::VectorXi> ordermap(order, len);
    Eigen::Map<const Eigen::VectorXd> etamap(eta, len);
    Eigen::Map<const Eigen::VectorXd> vmap(v, len);
    Eigen::Map<Eigen::VectorXd> wmap(w, len);
    Eigen::Map<Eigen::VectorXd> zmap(z, len);
    // borrow rskden to get an ordered v
    Eigen::Map<Eigen::VectorXd> orderedv(rskden, len);
    PermMat p(ordermap);

    orderedv = p.transpose() * vmap;

    // exp(eta) in increasing time order
    zmap.noalias() = p.transpose() * etamap.array().exp().matrix();


    // reverse cumulative sum
    double current = 0;
    for(int i = 0; i < len; ++i){
        current += zmap[len - 1 - i];
        zmap[len - 1 - i] = current;
    }

    // adjust ties, won't aliase
    for(int i = 0; i < len; ++i){
        zmap[i] = zmap[rankmin[i]];
    }

    // Now zmap is riskdenom, make wmap cumsum v/riskdenom^2
    current = 0;
    for(int i = 0; i < len; ++i){
        current += orderedv[i] /(zmap[i] * zmap[i]);
        wmap[i] = current;
    }

    // set zmap to be cumsum v/riskdenom
    current = 0;
    for(int i = 0; i < len; ++i){
        current += orderedv[i]/zmap[i];
        zmap[i] = current;
    }

    // adjust the ties
    for(int i = 0; i < len; ++i){
        wmap[i] = wmap[rankmax[i]];
        zmap[i] = zmap[rankmax[i]];
    }


    // Now wmap is cumsum v/riskdnom^2 zmap is cumsum v/riskdenom
    // We are ready to get the working weight and working response
    double wsum = 0;
    double zsum = 0;
    for(int i = 0; i < len; ++i){
        double local_eta = eta[order[i]];
        double eeta = exp(local_eta);
        wmap[i] = eeta * (zmap[i]-wmap[i] * eeta);
        // double grad = eeta * zmap[i] - orderedv[i];
        zmap[i] = orderedv[i] - eeta * zmap[i];
        wsum += wmap[i];
        zsum += zmap[i];
    }
    sumbuf[0] = zsum;
    sumbuf[1] = wsum;

    // Get the order back
    wmap = p * wmap;
    zmap = p * zmap;

}

double Cox::get_deviance(const double *y, const double *eta, const double *v,
                         int len)
{
    Eigen::Map<const Eigen::VectorXi> ordermap(order, len);
    PermMat p(ordermap);
    Eigen::Map<const Eigen::VectorXd> etamap(eta, len);
    Eigen::Map<Eigen::VectorXd> risk(rskden, len);
    risk = p.transpose() * etamap.array().exp().matrix();
    double current = 0;
    for(int i = 0;i < len; ++i){
        current += risk[len- 1 -i];
        risk[len - 1 - i] = current;
    }

    for(int i = 0; i < len; ++i){
        risk[i] = risk[rankmin[i]];
    }

    double result = 0;
    for(int i = 0; i < len; ++i){
        double local_v = v[order[i]];
        result += local_v * log(risk[i]) - v[i] * eta[i];
    }
    result *= 2;
    return result;

}

void Cox::get_residual(const double *y, const double *eta, const double *v,
                       double *r, int len)

{
    Eigen::Map<const Eigen::VectorXi> ordermap(order, len);
    Eigen::Map<const Eigen::VectorXd> etamap(eta, len);
    Eigen::Map<const Eigen::VectorXd> vmap(v, len);
    Eigen::Map<Eigen::VectorXd> rmap(r, len);

    // borrow rskden to get an ordered v
    Eigen::Map<Eigen::VectorXd> orderedv(rskden, len);
    PermMat p(ordermap);

    orderedv = p.transpose() * vmap;

    // exp(eta) in increasing time order
    rmap.noalias() = p.transpose() * etamap.array().exp().matrix();
    double current = 0;
    for(int i = 0;i < len; ++i){
        current += rmap[len- 1 -i];
        rmap[len - 1 - i] = current;
    }

    for(int i = 0; i < len; ++i){
        rmap[i] = rmap[rankmin[i]];
    }

    current = 0;
    for(int i = 0; i < len; ++i){
        current += orderedv[i]/rmap[i];
        orderedv[i] = current;
    }

    for(int i = 0; i < len; ++i){
        orderedv[i] = orderedv[rankmax[i]];
    }

    orderedv = p * orderedv;

    for(int i = 0; i < len; ++i){
        r[i] = v[i] - exp(eta[i]) * orderedv[i];
    }

}

double Cox::null_deviance(const double *y, const double *v, double *r, int intr,
                          double *eta, bool has_offset, const double *offset,
                          double *aint, int len)
{
    if (!has_offset)
    {
        for(int i = 0; i < len; ++i){
            r[i] = v[i];
            eta[i] = 0;
        }
        double result = 0;
        double current = 0;
        for (int i = 0; i < len; ++i)
        {
            double local_v = v[order[i]];
            result += local_v * log(len - rankmin[i]);
            // borrow rskden to store residual terms
            current += local_v/(len - rankmin[i]);
            rskden[i] = current;
        }

        for(int i = 0; i < len; ++i){
            rskden[i] = rskden[rankmax[i]];
            r[order[i]] -= rskden[i];
        }

        return result * 2;
    }
    throw "Offset has not been implemented for Cox regression\n";
}

double Cox::get_obj(const double *y, const double *v, const double *eta, const double *a,
                         const int *ju, const double *vp, int no, int ni, double lambda){
  return 0.0;
};
