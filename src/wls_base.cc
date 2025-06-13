#include <math.h>
#include "glmnetMatrix.h"
void wls_base(double alm0, double almc, double alpha, int m, int no, int ni,
              MatrixGlmnet *X, double *__restrict r, const double *v, int intr,
              const int *ju, const double *vp, const double *cl, int nx,
              double thr, int maxit, double *__restrict a, double *aint,
              double *__restrict g, int *__restrict ia, int *__restrict iy,
              int *iz, int *__restrict mm, int *nino, double *rsqc, int *nlp,
              double *__restrict xv, double * xm,int *jerr, int irls_iter, double rsum, double vsum) {

    const double xmz = vsum;
    if(!intr) {
        vsum = 0;
        rsum = 0;
    }
    double ab = almc * alpha;
    double dem = almc * (1.0 - alpha);
    double tlam = alpha * (2.0 * almc - alm0);

    for (int j = 0; j < ni; ++j) {
        if(!ju[j]) {
            continue;
        }

        if ((irls_iter == 0) && (iy[j] == 0)) {
            g[j] = fabs(X->dot_product(j, r, rsum));
            if(g[j] > tlam * vp[j]) {
                iy[j] = 1;
            }
        } 

        if(iy[j]) {
            xv[j] = X->vx2(j, v, vsum, (xm+j));
        }
    }

    bool jz = true;

    while (true) {
        if (!((*iz) && jz)) {
            (*nlp)++;
            double dlx = 0.0;
            for (int j = 0; j < ni; ++j) {
                if (!iy[j]) {
                    continue;
                }

                double gj = X->dot_product(j, r, rsum);
                double aj = a[j];
                double u = gj + aj * xv[j];
                double au = fabs(u) - vp[j] * ab;
                if (au < 0.0) {
                    a[j] = 0.0;
                } else {
                    a[j] = fmax(cl[2 * j],
                                fmin(cl[2 * j + 1],
                                     copysign(au, u) / (xv[j] + vp[j] * dem)));
                }

                if (a[j] == aj) {
                    continue;
                }

                if (mm[j] == 0) {
                    (*nino)++;
                    mm[j] = (*nino);
                    if ((*nino) > nx) {
                        break;
                    }
                    ia[(*nino) - 1] = j;
                }
                double d = a[j] - aj;
                (*rsqc) += d * (2.0 * gj - d * xv[j]);
                X->update_res(j, d, v, r, &rsum, vsum, xm[j]);
                dlx = fmax(xv[j] * d * d, dlx);
            }
            if ((*nino) > nx) {
                break;
            }
            if (intr) {
                double d = rsum / xmz;
                (*aint) += d;
                (*rsqc) += d * (2.0 * rsum - d * xmz);

                dlx = fmax(dlx, xmz * d * d);

                for (int i = 0; i < no; ++i) {
                    r[i] -= d * v[i];
                }
                rsum = 0.0;
            }

            // KKT checking here
            if (dlx < thr) {
                bool ixx = false;
                for (int j = 0; j < ni; ++j) {
                    if (iy[j] || (!ju[j])) {
                        continue;
                    }
                    g[j] = fabs(X->dot_product(j, r, rsum));
                    if (g[j] > ab * vp[j]) {
                        iy[j] = 1;
                        xv[j] = X->vx2(j, v, vsum, xm+j);
                        ixx = true;
                    }
                }

                if (ixx) {
                    continue;
                }
                break;
            }

            if ((*nlp) > maxit) {
                *jerr = -m;
                return;
            }
        }
        (*iz) = 1;

        while (true) {
            (*nlp)++;
            double dlx = 0.0;
            for (int l = 0; l < (*nino); ++l) {
                int k = ia[l];
                double gk = X->dot_product(k, r, rsum);
                double ak = a[k];
                double u = gk + ak * xv[k];
                double au = fabs(u) - vp[k] * ab;
                if (au < 0.0) {
                    a[k] = 0.0;
                } else {
                    a[k] = fmax(cl[2 * k],
                                fmin(cl[2 * k + 1],
                                     copysign(au, u) / (xv[k] + vp[k] * dem)));
                }

                if (ak == a[k]) {
                    continue;
                }
                double d = a[k] - ak;
                (*rsqc) += d * (2.0 * gk - d * xv[k]);
                X->update_res(k, d, v, r, &rsum, vsum, xm[k]);
                
                dlx = fmax(xv[k] * d * d, dlx);
            }

            if (intr) {
                double d = rsum / xmz;
                (*aint) += d;
                (*rsqc) += d * (2.0 * rsum - d * xmz);

                dlx = fmax(dlx, xmz * d * d);
                for (int i = 0; i < no; ++i) {
                    r[i] -= d * v[i];
                }
                rsum = 0.0;

            }

            if (dlx < thr) {
                break;
            }

            if ((*nlp) > maxit) {
                *jerr = -m;
                return;
            }
        }
        jz = false;
    }
    // free(xv);
    return;
}