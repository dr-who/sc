/* stats.h - Statistical functions
 * C89 portable for DOS, Linux
 * 
 * Data entry: Σ+, Σ-, clrΣ
 * Statistics: mean, sdev, sum, n, median
 * Regression: slope, intercept, r
 * Tests: t-test (one-sample and two-sample)
 */

#ifndef STATS_H
#define STATS_H

#include "config.h"

#ifdef HAVE_STATS

#include "apf.h"
#ifdef HAVE_MATRIX
#include "matrix.h"
#endif

/* Maximum data points for median and t-test */
#ifndef STATS_MAX_DATA
  #define STATS_MAX_DATA 1000
#endif

/* Far memory for DOS - puts large arrays outside 64KB DGROUP */
#ifdef __WATCOMC__
  #define SC_FAR __far
#else
  #define SC_FAR
#endif

/* Statistical registers - summary stats in DGROUP, raw data in far memory */
typedef struct {
    int n;              /* Count of data points */
    apf sum_x;          /* Σx */
    apf sum_x2;         /* Σx² */
    apf sum_y;          /* Σy (for 2-var stats) */
    apf sum_y2;         /* Σy² */
    apf sum_xy;         /* Σxy */
    int has_y;          /* Whether 2-var data was entered */
} stats_summary_t;

extern stats_summary_t stats_summary;

/* Raw data arrays declared in stats.c with SC_FAR */
extern apf SC_FAR stats_data_x[STATS_MAX_DATA];
extern apf SC_FAR stats_data_y[STATS_MAX_DATA];

/* Clear statistical registers */
void stats_clear(void);

/* Add data point (1-var) */
void stats_add(const apf *x);

/* Add data point (2-var) */
void stats_add2(const apf *x, const apf *y);

/* Remove data point (1-var) */
void stats_sub(const apf *x);

/* Get count */
int stats_count(void);

/* Get sum */
void stats_sum(apf *r);

/* Get mean */
void stats_mean(apf *r);

/* Get sample standard deviation (n-1) */
void stats_sdev(apf *r);

/* Get population standard deviation (n) */
void stats_sdevp(apf *r);

/* Get variance (sample) */
void stats_var(apf *r);

/* Get median */
void stats_median(apf *r);

/* Linear regression - slope */
void stats_slope(apf *r);

/* Linear regression - y-intercept */
void stats_intercept(apf *r);

/* Correlation coefficient */
void stats_corr(apf *r);

/* Predict y from x using regression */
void stats_predict_y(apf *r, const apf *x);

/* Predict x from y using regression */
void stats_predict_x(apf *r, const apf *y);

/* One-sample t-test: returns t-statistic comparing mean to mu0 */
void stats_ttest1(apf *t, apf *df, const apf *mu0);

/* Two-sample t-test: x data vs y data, returns t-statistic and df */
void stats_ttest2(apf *t, apf *df);

/* Process stats commands */
void cmd_stats(const char *args);

/* ============== UNIVERSITY STATISTICS FUNCTIONS ============== */

/* Descriptive statistics (double versions) */
double stat_trimean(const double *data, int n);
double stat_sem(const double *data, int n);

/* Probability distributions */
double stat_normcdf(double x);
double stat_tpdf(double x, double df);
double stat_tcdf(double x, double df);
double stat_tinv(double p, double df);
double stat_chi2pdf(double x, double df);
double stat_chi2cdf(double x, double df);
double stat_chi2inv(double p, double df);
double stat_fpdf(double x, double d1, double d2);
double stat_fcdf(double x, double d1, double d2);
double stat_finv(double p, double d1, double d2);
double stat_binopdf(int k, int n, double p);
double stat_binocdf(int k, int n, double p);
double stat_poisspdf(int k, double lambda);
double stat_poisscdf(int k, double lambda);
double stat_exppdf(double x, double lambda);
double stat_expcdf(double x, double lambda);
double stat_expinv(double p, double lambda);
double stat_unifpdf(double x, double a, double b);
double stat_unifcdf(double x, double a, double b);
double stat_unifinv(double p, double a, double b);
double stat_betainc(double x, double a, double b);

/* Hypothesis tests */
void stat_ttest2(const double *x, int nx, const double *y, int ny,
                 double *t_stat, double *p_value, double *df);
void stat_ztest(const double *x, int n, double mu0, double sigma,
                double *z_stat, double *p_value);
void stat_chi2test(const double *observed, const double *expected, int n,
                   double *chi2_stat, double *p_value, int *df);
void stat_anova1(const double *data, const int *groups, int n, int ngroups,
                 double *f_stat, double *p_value);

/* Effect sizes */
double stat_cohend(const double *x, int nx, const double *y, int ny);
double stat_rsquare(const double *y_actual, const double *y_pred, int n);

/* Numerical methods */
double num_simpson(const double *y, int n, double h);

#ifdef HAVE_MATRIX
/* ASCII scatter plot */
void mat_scatter(const matrix_t *x, const matrix_t *y, const matrix_t *groups, const char *chars);
#endif

#endif /* HAVE_STATS */
#endif /* STATS_H */
