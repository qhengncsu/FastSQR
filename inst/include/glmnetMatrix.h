#ifndef GLMNET_MATRIX
#define GLMNET_MATRIX

#include <stdint.h>
class MatrixGlmnet {
   public:
    // Compute the inner product of X[,j] and v
    virtual double dot_product(int j, const double* v, double vsum=0) = 0;

    // Compute the inner product of X[,j]^2 and v
    // A bit bulky here, also compute X[, j] * v and store to xm
    virtual double vx2(int j, const double* v, double vsum=0, double *xm=nullptr) = 0;

    // Set r = r - d*v*x[,j]
    virtual void update_res(int j, double d, const double* v, double* r, double *rsum=nullptr, double vsum=0, double vx=0) = 0;

    // Set eta = X * a + aint + offset, offset is optional
    virtual void compute_eta(double* eta, const double* a, double aint,
                             bool has_offset, const double* offset) = 0;
    

    // Compute weighted mean and standard deviation of each variable,
    // ignore variables with ju[j] == 0
    // if variable is constant, set ju[j] = 0
    // virtual void get_xmxs(const double* v, const int* ju, double* xm,
    //                       double* xs) = 0;

    static double sumv(const double* v, int len);
    static void update_res_eigen(double* r, const double *v, double d, int len);
    virtual ~MatrixGlmnet();

    int get_no();
    int get_ni();

    double max_grad(const double* r, const int* ju, const double* vp);

   protected:
    int no;  // Number of rows
    int ni;  // Number of variables
};

class DenseM : public MatrixGlmnet {
   public:
    DenseM(int no, int ni, const double* x);
    ~DenseM();

    double dot_product(int j, const double* v,double vsum=0);

    double vx2(int j, const double* v,double vsum=0, double *xm=nullptr);

    void update_res(int j, double d, const double* v, double* r, double *rsum=nullptr, double vsum=0, double vx=0);

    void compute_eta(double* eta, const double* a, double aint, bool has_offset,
                     const double* offset);

   private:
    const double* data;
};

class PlinkMatrix : public MatrixGlmnet {
   public:
    PlinkMatrix(int no, int ni, const uintptr_t* x, const double* xim, int intr, int ncov, const double *cov = nullptr);
    ~PlinkMatrix();

    double dot_product(int j, const double* v,double vsum=0);

    double vx2(int j, const double* v,double vsum=0, double *xm=nullptr);

    void update_res(int j, double d, const double* v, double* r, double *rsum=nullptr, double vsum=0, double vx=0);

    void compute_eta(double* eta, const double* a, double aint, bool has_offset,
                     const double* offset);

    void multv(double *eta, const double *weights, double aint);
    
   private:
     const uintptr_t* data;
     uint32_t word_ct;
     const double *xim; // mean imputation
     bool center;
     int ncov; // number of covariates
     const double *cov;
};

#endif