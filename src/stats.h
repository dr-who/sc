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

#endif /* HAVE_STATS */
#endif /* STATS_H */
