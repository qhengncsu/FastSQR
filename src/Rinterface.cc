#include "R.h"
#include "R_ext/Print.h"
#include "Rinternals.h"
#include "glmfamily.h"
#include "glmnetMatrix.h"
# include<stdlib.h>
// #include "gperftools/profiler.h"

void glmnetPath(double alpha, MatrixGlmnet *X, const double *y, const double *v,
                int intr, const int *ju, const double *vp, const double *cl,
                int nx, double thr, int maxit, GlmFamily *fam, bool has_offset,
                double *offset, const double *lambdas, int nlambda, int mxitnr,
                const double flmin, int *lmu, double *a0, double *ca, int *ia,
                int *nin, double *devratio_vec, double *alm, int *nlp,
                double *nulldev, int *jerr, double *a, int *iy, int *mm, int nino, int warm, double *residuals);

double get_max_lambda(double alpha,  MatrixGlmnet *X, const double *y, const double *v,
                      int intr, const int *ju, const double *vp,
                      GlmFamily *fam, bool has_offset, double *offset);

GlmFamily *get_family(const char *family, const int *order=nullptr, const int *rankmin=nullptr, const int *rankmax=nullptr, int len=0, double h=0.25, double tau=0.5) {
    if (strcmp(family, "gaussian") == 0) {
        return new Gaussian();
    }
    if (strcmp(family, "logistic") == 0) {
        return new Logistic();
    }
    if(strcmp(family, "cox") == 0) {
        return new Cox(order, rankmin, rankmax, len);
    }
    if(strcmp(family, "quantile") == 0) {
      return new Quantile(h,tau);
    }
    error("invalid family name\n");
}

GlmFamily *get_family2(const char *family, double h=0.25, double tau=0.5) {
  if(strcmp(family, "quantile") == 0) {
    return new Quantile(h,tau);
  }
  error("invalid family name\n");
}

void get_xmxs_dense(const double *x, const double *v, int *ju, double *xm,
                    double *xs, int no, int ni) {
    for (int j = 0; j < ni; ++j) {
        if (!ju[j]) {
            continue;
        }

        double ex = 0.0;
        double ex2 = 0.0;
        for (int i = 0; i < no; ++i) {
            ex += v[i] * x[j * no + i];
            ex2 += v[i] * x[j * no + i] * x[j * no + i];
        }
        // Assumes sum_{i} v_i = 1
        double variance = ex2 - ex * ex;
        if (variance <= 0) {
            ju[j] = 0;
            continue;
        }
        xm[j] = ex;
        xs[j] = sqrt(variance);
    }
}

void standardize(double *x, const int *ju, double *xm, double *xs, int intr,
                 int isd, int no, int ni) {
  if (intr) {
    for (int j = 0; j < ni; ++j) {
      if (!ju[j]) {
        continue;
      }
      for (int i = 0; i < no; ++i) {
        x[j * no + i] -= xm[j];
      }
    }

    if (isd) {
      for (int j = 0; j < ni; ++j) {
        if (!ju[j]) {
          continue;
        }
        for (int i = 0; i < no; ++i) {
          x[j * no + i] /= xs[j];
        }
      }
    }
    return;
  }

  if (isd) {
    for (int j = 0; j < ni; ++j) {
      if (!ju[j]) {
        continue;
      }
      for (int i = 0; i < no; ++i) {
        x[j * no + i] /= xs[j];
      }
    }
  }
}

MatrixGlmnet *get_matrix(SEXP xptr, const char *mattype, int no, int ni,
                         int isd, int intr, int *ju, double *xm, double *xs,
                         const double *v, const double *xim, int ncov, const double *cov) {
    if (strcmp(mattype, "Dense") == 0) {
        double *x = REAL(xptr);
        get_xmxs_dense(x, v, ju, xm, xs, no, ni);
        standardize(x, ju, xm, xs, intr, isd, no, ni);
        MatrixGlmnet *result = new DenseM(no, ni, x);
        return result;
    }

    if (strcmp(mattype, "Plink") == 0) {
        uintptr_t *x = (uintptr_t *)R_ExternalPtrAddr(xptr);
        MatrixGlmnet *result = new PlinkMatrix(no, ni, x, xim, intr, ncov, cov);
        return result;
    }

    error("invalid matrix type\n");
}

#ifdef __cplusplus
extern "C" {
#endif

SEXP testplink(SEXP x2, SEXP no2, SEXP ni2, SEXP xim2, SEXP v2, SEXP eta2) {
    uintptr_t *xptr = (uintptr_t *)R_ExternalPtrAddr(x2);
    MatrixGlmnet *x =
        new PlinkMatrix(asInteger(no2), asInteger(ni2), xptr, REAL(xim2), 0, 0, nullptr);
    double *eta = REAL(eta2);
    double *v = REAL(v2);
    x->compute_eta(eta, v, 0.0, false, nullptr);
    return R_NilValue;
}

SEXP testwz(SEXP eta2, SEXP order2, SEXP rankmax2, SEXP rankmin2, SEXP v2){
    const double *v = REAL(v2);
    int len = length(eta2);
    Cox cox(INTEGER(order2), INTEGER(rankmin2), INTEGER(rankmax2),len);
    double *z = (double*)malloc(sizeof(double)*len);
    double *w = (double*)malloc(sizeof(double)*len);
    double sumbuf[2];
    cox.get_workingset(REAL(eta2), nullptr, v, w, z, len, sumbuf);

    for(int i = 0; i < len; ++i){
        Rprintf("z[%d] is %f and w[%d] is %f \n", i, z[i], i, w[i]);
    }
    Rprintf("zsum is %f and wsum is %f\n", sumbuf[0], sumbuf[1]);
    double dev = cox.get_deviance(nullptr, REAL(eta2), v, len);
    Rprintf("deviance is %f\n", dev);

    double *r = (double*)malloc(sizeof(double)*len);
    double *eta3 = (double*)malloc(sizeof(double)*len);

    double nulldev = cox.null_deviance(nullptr, v, r,0,eta3,false,nullptr, nullptr, len);
    Rprintf("nulldeviance is %f \n", nulldev);
    for(int i = 0; i < len; ++i){
        Rprintf("residual[%d] is %f\n", i, r[i]);
    }

    cox.get_residual(nullptr, REAL(eta2), v, r, len);
    for(int i = 0; i < len; ++i){
        Rprintf("get residual residual[%d] is %f\n", i, r[i]);
    }


    free(r);
    free(eta3);
    free(z);
    free(w);
    return R_NilValue;
}

  SEXP testquantile(SEXP eta2, SEXP y2, SEXP h2, SEXP tau2, SEXP v2){
    const double *v = REAL(v2);
    int len = length(eta2);
    Quantile quantile(asReal(h2),asReal(tau2));
    double *z = (double*)malloc(sizeof(double)*len);
    double *w = (double*)malloc(sizeof(double)*len);
    double sumbuf[2];

    // 计算工作集
    quantile.get_workingset(REAL(eta2), REAL(y2), v, w, z, len, sumbuf);

    // 打印 w 和 z 数组
    Rprintf("\nw values: ");
    for(int i = 0; i < len; i++) Rprintf("%f ", w[i]);

    Rprintf("\nz values: ");
    for(int i = 0; i < len; i++) Rprintf("%f ", z[i]);

    double dev = quantile.get_deviance(REAL(y2), REAL(eta2), v, len);
    Rprintf("\ntau is %f\n", quantile.tau);
    Rprintf("deviance is %f\n", dev);

    dev = quantile.get_deviance(REAL(y2), REAL(eta2), v, len);
    Rprintf("deviance is %f\n", dev);

    double *r = (double*)malloc(sizeof(double)*len);
    double *eta3 = (double*)malloc(sizeof(double)*len);
    double aint = 0;

    // 计算零偏差
    double nulldev = quantile.null_deviance(REAL(y2), v, r, 1, eta3, false, nullptr, &aint, len);

    // 打印 r 和 eta3 数组
    Rprintf("\nr values: ");
    for(int i = 0; i < len; i++) Rprintf("%f ", r[i]);

    Rprintf("\neta3 values: ");
    for(int i = 0; i < len; i++) Rprintf("%f ", eta3[i]);

    Rprintf("\nnulldeviance is %f \n", nulldev);
    Rprintf("aint is %f \n", aint);

    // 清理内存
    free(r);
    free(eta3);
    free(z);
    free(w);
    return R_NilValue;
  }

SEXP max_lambda(SEXP alpha2, SEXP x2, SEXP y2, SEXP weights2,
                SEXP intr2, SEXP isd2, SEXP ju2, SEXP vp2,
                SEXP family2, SEXP h2, SEXP tau2,
                SEXP offset2, SEXP has_offset2, SEXP max_lambda){
  const char *mattype = CHAR(STRING_ELT(VECTOR_ELT(x2, 0), 0));
  int no = asInteger(VECTOR_ELT(x2, 1));
  int ni = asInteger(VECTOR_ELT(x2, 2));
  SEXP xptr = VECTOR_ELT(x2, 3);

  double *xim = nullptr;  // only for plink matrix
  double *cov = nullptr; // only for plink matrix
  int ncov = 0;
  if (strcmp(mattype, "Plink") == 0) {
    xim = REAL(VECTOR_ELT(x2, 4));
    ncov = asInteger(VECTOR_ELT(x2, 5));
    if(ncov > 0) {
      cov = REAL(VECTOR_ELT(x2, 6));
    }
  }
  int intr = asInteger(intr2);
  int isd = asInteger(isd2);
  int *ju = INTEGER(ju2);
  double *v = REAL(weights2);
  bool dup_x = (bool)(intr + isd);
  dup_x = dup_x && (strcmp(mattype, "Dense") == 0);
  if (dup_x) {
    xptr = PROTECT(duplicate(xptr));
  }

  // always compute xm and xs, maybe won't be used
  double *xm = (double *)malloc(sizeof(double) * ni);
  double *xs = (double *)malloc(sizeof(double) * ni);

  MatrixGlmnet *X =
    get_matrix(xptr, mattype, no, ni, isd, intr, ju, xm, xs, v, xim, ncov, cov);

  // Create family object
  const char *family = CHAR(STRING_ELT(family2, 0));
  double *y = nullptr;
  GlmFamily *fam = nullptr;
  if(strcmp(family, "cox") == 0){
    fam = get_family(family, INTEGER(VECTOR_ELT(y2, 0)), INTEGER(VECTOR_ELT(y2, 1)),INTEGER(VECTOR_ELT(y2, 2)), no);
  }else if(strcmp(family, "quantile") == 0){
    double h = asReal(h2);
    double tau = asReal(tau2);
    fam = get_family2(family, h, tau);
    y = REAL(y2);
  }else {
    fam = get_family(family);
    y = REAL(y2);
  }
  double alpha = asReal(alpha2);
  bool has_offset = (bool)asInteger(has_offset2);
  double *offset = nullptr;
  if (has_offset) {
    offset = REAL(offset2);
  }
  double *vp = REAL(vp2);
  double lambda_max = get_max_lambda(alpha, X, y, v, intr, ju, vp, fam, has_offset, offset);
  delete fam;
  delete X;
  Rprintf("Maximum lambda = %f. \n", lambda_max);
  double *lambdamax = REAL(max_lambda);
  free(xm);
  free(xs);
  *lambdamax = lambda_max;
  if (dup_x) UNPROTECT(1);
  return R_NilValue;
}

SEXP solve(SEXP alpha2, SEXP x2, SEXP y2, SEXP weights2, SEXP ju2, SEXP vp2,
           SEXP cl2, SEXP nx2, SEXP nlam2, SEXP flmin2, SEXP ulam2,
           SEXP thresh2, SEXP isd2, SEXP intr2, SEXP maxit2, SEXP lmu2,
           SEXP a02, SEXP ca2, SEXP ia2, SEXP nin2, SEXP devratio2, SEXP alm2,
           SEXP nlp2, SEXP family2, SEXP offset2, SEXP has_offset2,
           SEXP mxitnr2, SEXP nulldev2, SEXP jerr2, SEXP beta02, SEXP iy2,
           SEXP mm2, SEXP nino2, SEXP warm2, SEXP residuals2, SEXP h2, SEXP tau2){
    // Create matrix object
    const char *mattype = CHAR(STRING_ELT(VECTOR_ELT(x2, 0), 0));
    int no = asInteger(VECTOR_ELT(x2, 1));
    int ni = asInteger(VECTOR_ELT(x2, 2));
    SEXP xptr = VECTOR_ELT(x2, 3);

    double *xim = nullptr;  // only for plink matrix
    double *cov = nullptr; // only for plink matrix
    int ncov = 0;
    if (strcmp(mattype, "Plink") == 0) {
        xim = REAL(VECTOR_ELT(x2, 4));
        ncov = asInteger(VECTOR_ELT(x2, 5));
        if(ncov > 0) {
            cov = REAL(VECTOR_ELT(x2, 6));
        }
    }

    int intr = asInteger(intr2);
    int isd = asInteger(isd2);
    int *ju = INTEGER(ju2);
    double *v = REAL(weights2);
    bool dup_x = (bool)(intr + isd);
    dup_x = dup_x && (strcmp(mattype, "Dense") == 0);
    if (dup_x) {
        xptr = PROTECT(duplicate(xptr));
    }

    // always compute xm and xs, maybe won't be used
    double *xm = (double *)malloc(sizeof(double) * ni);
    double *xs = (double *)malloc(sizeof(double) * ni);

    MatrixGlmnet *X =
        get_matrix(xptr, mattype, no, ni, isd, intr, ju, xm, xs, v, xim, ncov, cov);

    // Create family object
    const char *family = CHAR(STRING_ELT(family2, 0));
    double *y = nullptr;
    GlmFamily *fam = nullptr;
    if(strcmp(family, "cox") == 0){
        fam = get_family(family, INTEGER(VECTOR_ELT(y2, 0)), INTEGER(VECTOR_ELT(y2, 1)),INTEGER(VECTOR_ELT(y2, 2)), no);
    }else if(strcmp(family, "quantile") == 0){
        double h = asReal(h2);
        double tau = asReal(tau2);
        fam = get_family2(family, h, tau);
        y = REAL(y2);
    }else {
        fam = get_family(family);
        y = REAL(y2);
    }
    double flmin = asReal(flmin2);
    int nlambda = asInteger(nlam2);
    double *lambdas = REAL(ulam2);

    double alpha = asReal(alpha2);

    double *vp = REAL(vp2);
    double *cl = REAL(cl2);
    int nx = asInteger(nx2);
    double thr = asReal(thresh2);
    int maxit = asInteger(maxit2);
    bool has_offset = (bool)asInteger(has_offset2);
    double *offset = nullptr;
    if (has_offset) {
        offset = REAL(offset2);
    }
    int mxitnr = asInteger(mxitnr2);

    // Things that needs to be modified
    int *lmu = INTEGER(lmu2);
    double *a0 = REAL(a02);
    double *ca = REAL(ca2);
    int *ia = INTEGER(ia2);
    int *nin = INTEGER(nin2);
    double *devratio = REAL(devratio2);
    double *alm = REAL(alm2);
    int *nlp = INTEGER(nlp2);
    int *jerr = INTEGER(jerr2);
    double *nulldev = REAL(nulldev2);
    double *residuals = REAL(residuals2);

    double *beta0 = REAL(PROTECT(duplicate(beta02)));
    int *iy = INTEGER(iy2);
    int *mm = INTEGER(mm2);
    int nino = asInteger(nino2);
    int warm = asInteger(warm2);

    // Forgetting to exclude plink matrix here gave me headache
    if(warm && isd && (strcmp(mattype, "Plink") != 0)) {
        Rprintf("standardizing variables here!\n");
        for(int j = 0; j < ni; ++j) {
            beta0[j] *= xs[j];
        }
    }
    glmnetPath(alpha, X, y, v, intr, ju, vp, cl, nx, thr, maxit, fam,
               has_offset, offset, lambdas, nlambda, mxitnr, flmin, lmu, a0, ca,
               ia, nin, devratio, alm, nlp, nulldev, jerr, beta0, iy, mm, nino, warm, residuals);

    // scale parameters back if standardization happend
    if (isd && (strcmp(mattype, "Plink") != 0)) {
      for (int m = 0; m < (*lmu); ++m) {
        for (int k = 0; k < nin[m]; ++k) {
          ca[m * nx + k] /= xs[ia[k]];
        }
      }
    }

    if (intr) {
      if (strcmp(mattype, "Plink") == 0) {
        for (int m = 0; m < (*lmu); ++m) {
          for (int k = 0; k < nin[m]; ++k) {
            if (ia[k] >= ncov) {
              a0[m] -= ca[m * nx + k] * xim[ia[k] - ncov];
            }
          }
        }
      } else {
        for (int m = 0; m < (*lmu); ++m) {
          for (int k = 0; k < nin[m]; ++k) {
            a0[m] -= ca[m * nx + k] * xm[ia[k]];
          }
        }
      }
    }

    delete fam;
    free(xm);
    free(xs);
    UNPROTECT(1);
    if (dup_x) {
        UNPROTECT(1);
    }
    delete X;
    return R_NilValue;
}


SEXP PlinkMultvC(SEXP x2, SEXP xim2, SEXP cov2, SEXP n2, SEXP p2, SEXP ncov2, SEXP beta2, SEXP a2, SEXP result2){
    int n = asInteger(n2);
    int p = asInteger(p2);
    int ncov = asInteger(ncov2);
    const double *cov = nullptr;
    if(ncov > 0){
        cov = REAL(cov2);
    }
    double *result = REAL(result2);
    const double *beta = REAL(beta2);
    int nlam = ncols(beta2);
    const double *a = REAL(a2);
    const double *xim = REAL(xim2);

    uintptr_t *x = (uintptr_t *)R_ExternalPtrAddr(x2);
    PlinkMatrix PM(n, p, x, xim, 0, ncov, cov);
    //Rprintf("n = %d\n", n);
    for(int j = 0; j < nlam; ++j){
        PM.multv(result+n*j, beta+p*j, a[j]);
    }
    return R_NilValue;
}

#ifdef __cplusplus
}
#endif
