#ifndef GLM_FAMILY
#define GLM_FAMILY

class GlmFamily {
   public:
     virtual int is_quantile() const { return 0; }
    // Compute the z and w in a IRLS, z here is actually a weighted working
    // residual, also compute sum of z and w
    virtual void get_workingset(const double *eta, const double *y,
                                const double *v, double *w, double *z,
                                int len, double* sumbuf) = 0;
    virtual double get_deviance(const double *y, const double *eta,
                                const double *v, int len) = 0;

    // Residual is vector such that r^TX is the gradient of the log-likelihood
    virtual void get_residual(const double *y, const double *eta,
                              const double *v, double *r, int len) = 0;

    // This function does three things:
    // 1. Compute the null deviance
    // 2. Get parameter estimate for the intercept, if needed
    // 3. Compute weighted residual of a model with intercept only
    virtual double null_deviance(const double *y, const double *v, double *r,
                                 int intr, double *eta, bool has_offset,
                                 const double *offset, double *aint,
                                 int len) = 0;

    virtual double get_obj(const double *y, const double *v, const double *eta, const double *a,
                   const int *ju, const double *vp, int no, int ni, double lambda) = 0;

    virtual ~GlmFamily();
};

class Gaussian : public GlmFamily {
   public:
    Gaussian();
     int is_quantile() const { return 0; }
    // Do nothing
    void get_workingset(const double *eta, const double *y, const double *v,
                        double *w, double *z, int len, double* sumbuf);

    double get_deviance(const double *y, const double *eta, const double *v,
                        int len);

    void get_residual(const double *y, const double *eta, const double *v,
                      double *r, int len);

    double null_deviance(const double *y, const double *v, double *r, int intr,
                         double *eta, bool has_offset, const double *offset,
                         double *aint, int len);
    double get_obj(const double *y, const double *v, const double *eta, const double *a,
                   const int *ju, const double *vp, int no, int ni, double lambda);
};

class Logistic : public GlmFamily {
   public:
    Logistic();
    int is_quantile() const { return 0; }
    void get_workingset(const double *eta, const double *y, const double *v,
                        double *w, double *z, int len, double* sumbuf);

    double get_deviance(const double *y, const double *eta, const double *v,
                        int len);

    void get_residual(const double *y, const double *eta, const double *v,
                      double *r, int len);

    double null_deviance(const double *y, const double *v, double *r, int intr,
                         double *eta, bool has_offset, const double *offset,
                         double *aint, int len);
    double get_obj(const double *y, const double *v, const double *eta, const double *a,
                   const int *ju, const double *vp, int no, int ni, double lambda);
};

class Quantile : public GlmFamily{
  public:
    double h;
    double tau;
    int is_quantile() const { return 1; }
    void get_workingset(const double *eta, const double *y, const double *v,
                        double *w, double *z, int len, double* sumbuf);

    double get_deviance(const double *y, const double *eta, const double *v,
                        int len);

    void get_residual(const double *y, const double *eta, const double *v,
                      double *r, int len);

    double null_deviance(const double *y, const double *v, double *r, int intr,
                         double *eta, bool has_offset, const double *offset,
                         double *aint, int len);

    double get_obj(const double *y, const double *v, const double *eta, const double *a,
                   const int *ju, const double *vp, int no, int ni, double lambda);

    Quantile(double bandwidth, double quantile): h(bandwidth), tau(quantile){}
};

class Cox : public GlmFamily {
   public:
    Cox(const int *order, const int * rankmin, const int * rankmax, int len);
    ~Cox();
    void get_workingset(const double *eta, const double *y, const double *v,
                        double *w, double *z, int len, double* sumbuf);

    double get_deviance(const double *y, const double *eta, const double *v,
                        int len);

    void get_residual(const double *y, const double *eta, const double *v,
                      double *r, int len);

    double null_deviance(const double *y, const double *v, double *r, int intr,
                         double *eta, bool has_offset, const double *offset,
                         double *aint, int len);
    double get_obj(const double *y, const double *v, const double *eta, const double *a,
                   const int *ju, const double *vp, int no, int ni, double lambda);
    private:
      const int * order;
      const int * rankmin;
      const int * rankmax;
      double *rskden;

};

#endif
