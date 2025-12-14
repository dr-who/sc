/* stats.h - Statistical functions
 * C89 portable for DOS, Linux
 * 
 * Data entry: Σ+, Σ-, clrΣ
 * Statistics: mean, sdev, sum, n
 * Regression: slope, intercept, r
 */

#ifndef STATS_H
#define STATS_H

#include "config.h"

#ifdef HAVE_STATS

#include "apf.h"

/* Maximum data points (memory limited on VIC-20) */
#ifndef STATS_MAX_DATA
#define STATS_MAX_DATA 100
#endif

/* Statistical registers */
typedef struct {
    int n;              /* Count of data points */
    apf sum_x;          /* Σx */
    apf sum_x2;         /* Σx² */
    apf sum_y;          /* Σy (for 2-var stats) */
    apf sum_y2;         /* Σy² */
    apf sum_xy;         /* Σxy */
} stats_regs_t;

extern stats_regs_t stats_regs;

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

/* Process stats commands */
void cmd_stats(const char *args);

#endif /* HAVE_STATS */
#endif /* STATS_H */
