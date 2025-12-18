/* stats.c - Statistical functions
 * C89 portable for DOS, Linux
 * Uses __far for large arrays to stay under 64KB DGROUP limit on DOS
 */

#include "config.h"

#ifdef HAVE_STATS

#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include <stdlib.h>
#include "stats.h"
#include "apf.h"
#include "apfx.h"
#include "apf_native.h"
#include "matrix.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#ifndef INFINITY
#define INFINITY HUGE_VAL
#endif

/* C89-compatible approximations for special functions */
static double c89_erf(double x)
{
    /* Horner form approximation */
    double t = 1.0 / (1.0 + 0.3275911 * (x < 0 ? -x : x));
    double y = 1.0 - (((((1.061405429 * t - 1.453152027) * t) + 1.421413741)
                        * t - 0.284496736) * t + 0.254829592) * t * exp(-x * x);
    return x < 0 ? -y : y;
}

static double c89_lgamma(double x)
{
    /* Lanczos approximation */
    static double c[7] = {
        1.000000000190015,
        76.18009172947146,
        -86.50532032941677,
        24.01409824083091,
        -1.231739572450155,
        0.1208650973866179e-2,
        -0.5395239384953e-5
    };
    double g = 5.0;
    double sum = c[0];
    int i;
    double t;
    
    if (x <= 0 && x == floor(x)) return 1e308;  /* Pole */
    
    for (i = 1; i < 7; i++) sum += c[i] / (x + i - 1);
    t = x + g - 0.5;
    return 0.5 * log(2 * M_PI) + (x - 0.5) * log(t) - t + log(sum);
}

static double c89_tgamma(double x)
{
    if (x <= 0 && x == floor(x)) return 1e308;  /* Pole */
    return exp(c89_lgamma(x));
}

/* Summary statistics (small, stays in DGROUP) */
stats_summary_t stats_summary;

/* Raw data arrays - use __far on DOS to put in far memory */
apf SC_FAR stats_data_x[STATS_MAX_DATA];
apf SC_FAR stats_data_y[STATS_MAX_DATA];

void stats_clear(void) {
    int i;
    stats_summary.n = 0;
    stats_summary.has_y = 0;
    apf_zero(&stats_summary.sum_x);
    apf_zero(&stats_summary.sum_x2);
    apf_zero(&stats_summary.sum_y);
    apf_zero(&stats_summary.sum_y2);
    apf_zero(&stats_summary.sum_xy);
    for (i = 0; i < STATS_MAX_DATA; i++) {
        apf_zero(&stats_data_x[i]);
        apf_zero(&stats_data_y[i]);
    }
}

void stats_add(const apf *x) {
    apf tmp;
    
    if (stats_summary.n >= STATS_MAX_DATA) {
        printf("Error: max data points (%d) reached\n", STATS_MAX_DATA);
        return;
    }
    
    /* Store raw data in far memory */
    apf_copy(&stats_data_x[stats_summary.n], x);
    
    stats_summary.n++;
    
    /* sum_x += x */
    apf_add(&tmp, &stats_summary.sum_x, x);
    apf_copy(&stats_summary.sum_x, &tmp);
    
    /* sum_x2 += x^2 */
    apf_mul(&tmp, x, x);
    apf_add(&stats_summary.sum_x2, &stats_summary.sum_x2, &tmp);
}

void stats_add2(const apf *x, const apf *y) {
    apf tmp;
    
    if (stats_summary.n >= STATS_MAX_DATA) {
        printf("Error: max data points (%d) reached\n", STATS_MAX_DATA);
        return;
    }
    
    /* Store raw data in far memory */
    apf_copy(&stats_data_x[stats_summary.n], x);
    apf_copy(&stats_data_y[stats_summary.n], y);
    stats_summary.has_y = 1;
    
    stats_summary.n++;
    
    /* sum_x += x */
    apf_add(&tmp, &stats_summary.sum_x, x);
    apf_copy(&stats_summary.sum_x, &tmp);
    
    /* sum_x2 += x^2 */
    apf_mul(&tmp, x, x);
    apf_add(&stats_summary.sum_x2, &stats_summary.sum_x2, &tmp);
    
    /* sum_y += y */
    apf_add(&tmp, &stats_summary.sum_y, y);
    apf_copy(&stats_summary.sum_y, &tmp);
    
    /* sum_y2 += y^2 */
    apf_mul(&tmp, y, y);
    apf_add(&stats_summary.sum_y2, &stats_summary.sum_y2, &tmp);
    
    /* sum_xy += x*y */
    apf_mul(&tmp, x, y);
    apf_add(&stats_summary.sum_xy, &stats_summary.sum_xy, &tmp);
}

void stats_sub(const apf *x) {
    apf tmp;
    int i, found = -1;
    
    if (stats_summary.n > 0) {
        /* Find and remove the data point */
        for (i = 0; i < stats_summary.n; i++) {
            if (apf_eq(&stats_data_x[i], x)) {
                found = i;
                break;
            }
        }
        
        if (found >= 0) {
            /* Shift remaining data */
            for (i = found; i < stats_summary.n - 1; i++) {
                apf_copy(&stats_data_x[i], &stats_data_x[i+1]);
                apf_copy(&stats_data_y[i], &stats_data_y[i+1]);
            }
        }
        
        stats_summary.n--;
        
        /* sum_x -= x */
        apf_sub(&tmp, &stats_summary.sum_x, x);
        apf_copy(&stats_summary.sum_x, &tmp);
        
        /* sum_x2 -= x^2 */
        apf_mul(&tmp, x, x);
        apf_sub(&stats_summary.sum_x2, &stats_summary.sum_x2, &tmp);
    }
}

int stats_count(void) {
    return stats_summary.n;
}

void stats_sum(apf *r) {
    apf_copy(r, &stats_summary.sum_x);
}

void stats_mean(apf *r) {
    apf n;
    if (stats_summary.n == 0) {
        apf_zero(r);
        return;
    }
    apf_from_int(&n, stats_summary.n);
    apf_div(r, &stats_summary.sum_x, &n);
}

/* Sample variance: (Sx2 - (Sx)^2/n) / (n-1) */
void stats_var(apf *r) {
    apf n, n1, sum2_over_n, numer;
    
    if (stats_summary.n < 2) {
        apf_zero(r);
        return;
    }
    
    apf_from_int(&n, stats_summary.n);
    apf_from_int(&n1, stats_summary.n - 1);
    
    /* (Sx)^2 / n */
    apf_mul(&sum2_over_n, &stats_summary.sum_x, &stats_summary.sum_x);
    apf_div(&sum2_over_n, &sum2_over_n, &n);
    
    /* Sx2 - (Sx)^2/n */
    apf_sub(&numer, &stats_summary.sum_x2, &sum2_over_n);
    
    /* Divide by n-1 */
    apf_div(r, &numer, &n1);
}

/* Sample standard deviation: sqrt(variance) */
void stats_sdev(apf *r) {
    stats_var(r);
    apf_sqrt(r, r);
}

/* Population standard deviation: sqrt((Sx2 - (Sx)^2/n) / n) */
void stats_sdevp(apf *r) {
    apf n, sum2_over_n, numer;
    
    if (stats_summary.n < 1) {
        apf_zero(r);
        return;
    }
    
    apf_from_int(&n, stats_summary.n);
    
    /* (Sx)^2 / n */
    apf_mul(&sum2_over_n, &stats_summary.sum_x, &stats_summary.sum_x);
    apf_div(&sum2_over_n, &sum2_over_n, &n);
    
    /* Sx2 - (Sx)^2/n */
    apf_sub(&numer, &stats_summary.sum_x2, &sum2_over_n);
    
    /* Divide by n */
    apf_div(r, &numer, &n);
    
    /* sqrt */
    apf_sqrt(r, r);
}

/* Median - requires sorting a copy of the data */
void stats_median(apf *r) {
    /* Use small local buffer for sorting - limits median to first 100 points */
    apf sorted[100];
    int i, j, n_use;
    apf tmp;
    
    if (stats_summary.n == 0) {
        apf_zero(r);
        return;
    }
    
    n_use = stats_summary.n;
    if (n_use > 100) n_use = 100;  /* Limit for stack safety */
    
    /* Copy data from far memory to local */
    for (i = 0; i < n_use; i++) {
        apf_copy(&sorted[i], &stats_data_x[i]);
    }
    
    /* Simple insertion sort */
    for (i = 1; i < n_use; i++) {
        apf_copy(&tmp, &sorted[i]);
        j = i - 1;
        while (j >= 0 && apf_cmp(&sorted[j], &tmp) > 0) {
            apf_copy(&sorted[j+1], &sorted[j]);
            j--;
        }
        apf_copy(&sorted[j+1], &tmp);
    }
    
    /* Get median */
    if (n_use % 2 == 1) {
        /* Odd: middle element */
        apf_copy(r, &sorted[n_use / 2]);
    } else {
        /* Even: average of two middle elements */
        apf two;
        apf_add(&tmp, &sorted[n_use/2 - 1], &sorted[n_use/2]);
        apf_from_int(&two, 2);
        apf_div(r, &tmp, &two);
    }
}

/* Linear regression slope: (nSxy - SxSy) / (nSx2 - (Sx)^2) */
void stats_slope(apf *r) {
    apf n, numer, denom, tmp1;
    
    if (stats_summary.n < 2) {
        apf_zero(r);
        return;
    }
    
    apf_from_int(&n, stats_summary.n);
    
    /* nSxy */
    apf_mul(&numer, &n, &stats_summary.sum_xy);
    /* SxSy */
    apf_mul(&tmp1, &stats_summary.sum_x, &stats_summary.sum_y);
    /* nSxy - SxSy */
    apf_sub(&numer, &numer, &tmp1);
    
    /* nSx2 */
    apf_mul(&denom, &n, &stats_summary.sum_x2);
    /* (Sx)^2 */
    apf_mul(&tmp1, &stats_summary.sum_x, &stats_summary.sum_x);
    /* nSx2 - (Sx)^2 */
    apf_sub(&denom, &denom, &tmp1);
    
    apf_div(r, &numer, &denom);
}

/* Y-intercept: (Sy - slope*Sx) / n */
void stats_intercept(apf *r) {
    apf slope, n, tmp;
    
    if (stats_summary.n < 2) {
        apf_zero(r);
        return;
    }
    
    stats_slope(&slope);
    apf_from_int(&n, stats_summary.n);
    
    /* slope * Sx */
    apf_mul(&tmp, &slope, &stats_summary.sum_x);
    /* Sy - slope*Sx */
    apf_sub(r, &stats_summary.sum_y, &tmp);
    /* Divide by n */
    apf_div(r, r, &n);
}

/* Correlation coefficient r */
void stats_corr(apf *r) {
    apf n, numer, denom_x, denom_y, denom, tmp1;
    
    if (stats_summary.n < 2) {
        apf_zero(r);
        return;
    }
    
    apf_from_int(&n, stats_summary.n);
    
    /* Numerator: nSxy - SxSy */
    apf_mul(&numer, &n, &stats_summary.sum_xy);
    apf_mul(&tmp1, &stats_summary.sum_x, &stats_summary.sum_y);
    apf_sub(&numer, &numer, &tmp1);
    
    /* Denominator: sqrt[(nSx2 - (Sx)^2)(nSy2 - (Sy)^2)] */
    apf_mul(&denom_x, &n, &stats_summary.sum_x2);
    apf_mul(&tmp1, &stats_summary.sum_x, &stats_summary.sum_x);
    apf_sub(&denom_x, &denom_x, &tmp1);
    
    apf_mul(&denom_y, &n, &stats_summary.sum_y2);
    apf_mul(&tmp1, &stats_summary.sum_y, &stats_summary.sum_y);
    apf_sub(&denom_y, &denom_y, &tmp1);
    
    apf_mul(&denom, &denom_x, &denom_y);
    apf_sqrt(&denom, &denom);
    
    apf_div(r, &numer, &denom);
}

/* Predict y = slope*x + intercept */
void stats_predict_y(apf *r, const apf *x) {
    apf slope, intercept, tmp;
    stats_slope(&slope);
    stats_intercept(&intercept);
    apf_mul(&tmp, &slope, x);
    apf_add(r, &tmp, &intercept);
}

/* Predict x = (y - intercept) / slope */
void stats_predict_x(apf *r, const apf *y) {
    apf slope, intercept, tmp;
    stats_slope(&slope);
    stats_intercept(&intercept);
    apf_sub(&tmp, y, &intercept);
    apf_div(r, &tmp, &slope);
}

/* One-sample t-test: t = (mean - mu0) / (sdev / sqrt(n)) */
void stats_ttest1(apf *t, apf *df, const apf *mu0) {
    apf mean, sdev, n_apf, sqrt_n, se;
    
    if (stats_summary.n < 2) {
        apf_zero(t);
        apf_zero(df);
        return;
    }
    
    stats_mean(&mean);
    stats_sdev(&sdev);
    apf_from_int(&n_apf, stats_summary.n);
    apf_sqrt(&sqrt_n, &n_apf);
    
    /* Standard error = sdev / sqrt(n) */
    apf_div(&se, &sdev, &sqrt_n);
    
    /* t = (mean - mu0) / se */
    apf_sub(t, &mean, mu0);
    apf_div(t, t, &se);
    
    /* df = n - 1 */
    apf_from_int(df, stats_summary.n - 1);
}

/* Two-sample t-test (Welch's t-test for unequal variances) */
void stats_ttest2(apf *t, apf *df) {
    apf mean_x, mean_y, var_x, var_y;
    apf n_x, n_y, tmp1, tmp2, se2, se, numer;
    int i;
    apf sum_y_local, sum_y2_local;
    
    if (!stats_summary.has_y || stats_summary.n < 2) {
        apf_zero(t);
        apf_zero(df);
        printf("Error: need 2-var data for two-sample t-test\n");
        return;
    }
    
    /* Calculate mean_x from x data */
    apf_from_int(&n_x, stats_summary.n);
    apf_div(&mean_x, &stats_summary.sum_x, &n_x);
    
    /* Calculate variance of x */
    stats_var(&var_x);
    
    /* Calculate mean_y and variance_y from y data in far memory */
    apf_zero(&sum_y_local);
    apf_zero(&sum_y2_local);
    for (i = 0; i < stats_summary.n; i++) {
        apf_add(&sum_y_local, &sum_y_local, &stats_data_y[i]);
        apf_mul(&tmp1, &stats_data_y[i], &stats_data_y[i]);
        apf_add(&sum_y2_local, &sum_y2_local, &tmp1);
    }
    apf_from_int(&n_y, stats_summary.n);
    apf_div(&mean_y, &sum_y_local, &n_y);
    
    /* var_y = (sum_y2 - sum_y^2/n) / (n-1) */
    {
        apf n1;
        apf_mul(&tmp1, &sum_y_local, &sum_y_local);
        apf_div(&tmp1, &tmp1, &n_y);
        apf_sub(&tmp2, &sum_y2_local, &tmp1);
        apf_from_int(&n1, stats_summary.n - 1);
        apf_div(&var_y, &tmp2, &n1);
    }
    
    /* se^2 = var_x/n + var_y/n */
    apf_div(&tmp1, &var_x, &n_x);
    apf_div(&tmp2, &var_y, &n_y);
    apf_add(&se2, &tmp1, &tmp2);
    apf_sqrt(&se, &se2);
    
    /* t = (mean_x - mean_y) / se */
    apf_sub(&numer, &mean_x, &mean_y);
    apf_div(t, &numer, &se);
    
    /* Welch-Satterthwaite df approximation */
    {
        apf vx_nx, vy_ny, num, d1, d2, denom, n1;
        
        apf_div(&vx_nx, &var_x, &n_x);
        apf_div(&vy_ny, &var_y, &n_y);
        
        /* Numerator: (vx/nx + vy/ny)^2 */
        apf_add(&tmp1, &vx_nx, &vy_ny);
        apf_mul(&num, &tmp1, &tmp1);
        
        /* Denominator terms */
        apf_from_int(&n1, stats_summary.n - 1);
        apf_mul(&d1, &vx_nx, &vx_nx);
        apf_div(&d1, &d1, &n1);
        apf_mul(&d2, &vy_ny, &vy_ny);
        apf_div(&d2, &d2, &n1);
        apf_add(&denom, &d1, &d2);
        
        apf_div(df, &num, &denom);
    }
}

void cmd_stats(const char *args) {
    char buf[64];
    apf val, val2;
    
    /* Skip whitespace */
    while (*args && isspace((unsigned char)*args)) args++;
    
    if (*args == '\0') {
        /* Show current stats */
        printf("n      = %d\n", stats_summary.n);
        apf_to_str(buf, sizeof(buf), &stats_summary.sum_x, 10);
        printf("Sx     = %s\n", buf);
        stats_mean(&val);
        apf_to_str(buf, sizeof(buf), &val, 10);
        printf("mean   = %s\n", buf);
        if (stats_summary.n >= 1) {
            stats_median(&val);
            apf_to_str(buf, sizeof(buf), &val, 10);
            printf("median = %s\n", buf);
        }
        if (stats_summary.n >= 2) {
            stats_sdev(&val);
            apf_to_str(buf, sizeof(buf), &val, 10);
            printf("sdev   = %s\n", buf);
        }
        if (stats_summary.has_y && stats_summary.n >= 2) {
            stats_corr(&val);
            apf_to_str(buf, sizeof(buf), &val, 10);
            printf("r      = %s\n", buf);
        }
        return;
    }
    
    if (strcmp(args, "clr") == 0 || strcmp(args, "clear") == 0) {
        stats_clear();
        printf("Statistics cleared.\n");
        return;
    }
    
    if (strcmp(args, "mean") == 0) {
        stats_mean(&val);
        apf_to_str(buf, sizeof(buf), &val, 10);
        printf("mean = %s\n", buf);
        return;
    }
    
    if (strcmp(args, "median") == 0 || strcmp(args, "med") == 0) {
        stats_median(&val);
        apf_to_str(buf, sizeof(buf), &val, 10);
        printf("median = %s\n", buf);
        return;
    }
    
    if (strcmp(args, "sdev") == 0 || strcmp(args, "sd") == 0) {
        stats_sdev(&val);
        apf_to_str(buf, sizeof(buf), &val, 10);
        printf("sdev = %s\n", buf);
        return;
    }
    
    if (strcmp(args, "var") == 0 || strcmp(args, "variance") == 0) {
        stats_var(&val);
        apf_to_str(buf, sizeof(buf), &val, 10);
        printf("var = %s\n", buf);
        return;
    }
    
    if (strcmp(args, "sum") == 0) {
        stats_sum(&val);
        apf_to_str(buf, sizeof(buf), &val, 10);
        printf("sum = %s\n", buf);
        return;
    }
    
    if (strcmp(args, "n") == 0 || strcmp(args, "count") == 0) {
        printf("n = %d\n", stats_count());
        return;
    }
    
    if (strcmp(args, "corr") == 0 || strcmp(args, "r") == 0) {
        stats_corr(&val);
        apf_to_str(buf, sizeof(buf), &val, 10);
        printf("r = %s\n", buf);
        return;
    }
    
    if (strcmp(args, "slope") == 0 || strcmp(args, "m") == 0) {
        stats_slope(&val);
        apf_to_str(buf, sizeof(buf), &val, 10);
        printf("slope = %s\n", buf);
        return;
    }
    
    if (strcmp(args, "intercept") == 0 || strcmp(args, "b") == 0) {
        stats_intercept(&val);
        apf_to_str(buf, sizeof(buf), &val, 10);
        printf("intercept = %s\n", buf);
        return;
    }
    
    if (strcmp(args, "reg") == 0 || strcmp(args, "regression") == 0) {
        stats_slope(&val);
        stats_intercept(&val2);
        apf_to_str(buf, sizeof(buf), &val, 10);
        printf("y = %s * x", buf);
        apf_to_str(buf, sizeof(buf), &val2, 10);
        if (val2.sign) {
            printf(" - ");
            val2.sign = 0;
            apf_to_str(buf, sizeof(buf), &val2, 10);
        } else {
            printf(" + ");
        }
        printf("%s\n", buf);
        stats_corr(&val);
        apf_to_str(buf, sizeof(buf), &val, 10);
        printf("r = %s\n", buf);
        return;
    }
    
    if (strncmp(args, "ttest", 5) == 0) {
        const char *p = args + 5;
        while (*p && isspace((unsigned char)*p)) p++;
        
        if (*p == '\0') {
            apf zero;
            apf_zero(&zero);
            stats_ttest1(&val, &val2, &zero);
            apf_to_str(buf, sizeof(buf), &val, 10);
            printf("t = %s\n", buf);
            apf_to_str(buf, sizeof(buf), &val2, 10);
            printf("df = %s\n", buf);
            return;
        }
        
        if (*p == '2' || strncmp(p, "two", 3) == 0) {
            stats_ttest2(&val, &val2);
            apf_to_str(buf, sizeof(buf), &val, 10);
            printf("t = %s\n", buf);
            apf_to_str(buf, sizeof(buf), &val2, 10);
            printf("df = %s (Welch-Satterthwaite)\n", buf);
            return;
        }
        
        {
            apf mu0;
            char mu_str[32];
            if (sscanf(p, "%31s", mu_str) == 1) {
                apf_from_str(&mu0, mu_str);
                stats_ttest1(&val, &val2, &mu0);
                apf_to_str(buf, sizeof(buf), &val, 10);
                printf("t = %s\n", buf);
                apf_to_str(buf, sizeof(buf), &val2, 10);
                printf("df = %s\n", buf);
                return;
            }
        }
    }
    
    printf("Usage: stat [clr|mean|median|sdev|var|sum|n|corr|slope|intercept|reg|ttest]\n");
}

#endif /* HAVE_STATS */

/* ============================================================
 * UNIVERSITY STATISTICS & NUMERICAL METHODS
 * Additional functions for 1st/2nd year courses
 * ============================================================ */

/* Standard normal CDF (double version for internal use) */
double stat_normcdf(double x)
{
    /* Approximation using error function */
    return 0.5 * (1.0 + c89_erf(x / sqrt(2.0)));
}

/* ============== DESCRIPTIVE STATISTICS ============== */

/* Trimean: (Q1 + 2*median + Q3) / 4 */
double stat_trimean(const double *data, int n)
{
    double *sorted, q1, med, q3;
    int i, j;
    
    if (n < 1) return 0;
    sorted = (double*)malloc(n * sizeof(double));
    if (!sorted) return 0;
    
    for (i = 0; i < n; i++) sorted[i] = data[i];
    for (i = 0; i < n - 1; i++) {
        for (j = 0; j < n - i - 1; j++) {
            if (sorted[j] > sorted[j + 1]) {
                double tmp = sorted[j];
                sorted[j] = sorted[j + 1];
                sorted[j + 1] = tmp;
            }
        }
    }
    
    /* Quartiles */
    if (n % 2 == 0) {
        med = (sorted[n/2 - 1] + sorted[n/2]) / 2.0;
    } else {
        med = sorted[n/2];
    }
    
    if ((n/2) % 2 == 0 && n/2 > 0) {
        q1 = (sorted[n/4 - 1] + sorted[n/4]) / 2.0;
        q3 = (sorted[3*n/4 - 1] + sorted[3*n/4]) / 2.0;
    } else {
        q1 = sorted[n/4];
        q3 = sorted[3*n/4];
    }
    
    free(sorted);
    return (q1 + 2*med + q3) / 4.0;
}

/* Standard Error of Mean: std / sqrt(n) */
double stat_sem(const double *data, int n)
{
    double mean = 0, var = 0;
    int i;
    if (n < 2) return 0;
    
    for (i = 0; i < n; i++) mean += data[i];
    mean /= n;
    
    for (i = 0; i < n; i++) {
        double d = data[i] - mean;
        var += d * d;
    }
    var /= (n - 1);  /* Sample variance */
    
    return sqrt(var / n);
}

/* ============== PROBABILITY DISTRIBUTIONS ============== */

/* Helper: Regularized incomplete beta function (for F and t distributions) */
static double betainc_cf(double a, double b, double x)
{
    /* Continued fraction for incomplete beta */
    double qab = a + b;
    double qap = a + 1.0;
    double qam = a - 1.0;
    double c = 1.0;
    double d = 1.0 - qab * x / qap;
    double aa, del, h;
    int m;
    
    if (fabs(d) < 1e-30) d = 1e-30;
    d = 1.0 / d;
    h = d;
    
    for (m = 1; m <= 100; m++) {
        int m2 = 2 * m;
        aa = m * (b - m) * x / ((qam + m2) * (a + m2));
        d = 1.0 + aa * d;
        if (fabs(d) < 1e-30) d = 1e-30;
        c = 1.0 + aa / c;
        if (fabs(c) < 1e-30) c = 1e-30;
        d = 1.0 / d;
        h *= d * c;
        aa = -(a + m) * (qab + m) * x / ((a + m2) * (qap + m2));
        d = 1.0 + aa * d;
        if (fabs(d) < 1e-30) d = 1e-30;
        c = 1.0 + aa / c;
        if (fabs(c) < 1e-30) c = 1e-30;
        d = 1.0 / d;
        del = d * c;
        h *= del;
        if (fabs(del - 1.0) < 1e-10) break;
    }
    return h;
}

double stat_betainc(double x, double a, double b)
{
    double bt;
    if (x < 0 || x > 1) return 0;
    if (x == 0 || x == 1) return x;
    
    bt = exp(c89_lgamma(a + b) - c89_lgamma(a) - c89_lgamma(b) + a * log(x) + b * log(1 - x));
    
    if (x < (a + 1) / (a + b + 2)) {
        return bt * betainc_cf(a, b, x) / a;
    } else {
        return 1.0 - bt * betainc_cf(b, a, 1 - x) / b;
    }
}

/* Student's t distribution PDF */
double stat_tpdf(double x, double df)
{
    double c = exp(c89_lgamma((df + 1) / 2) - c89_lgamma(df / 2)) / sqrt(df * M_PI);
    return c * pow(1 + x * x / df, -(df + 1) / 2);
}

/* Student's t distribution CDF */
double stat_tcdf(double x, double df)
{
    double t = df / (df + x * x);
    double p = stat_betainc(t, df / 2, 0.5) / 2;
    return x >= 0 ? 1 - p : p;
}

/* Student's t inverse CDF (quantile) - Newton's method */
double stat_tinv(double p, double df)
{
    double x = 0, low = -100, high = 100;
    int i;
    
    if (p <= 0) return -1e308;
    if (p >= 1) return 1e308;
    
    /* Bisection to get close */
    for (i = 0; i < 50; i++) {
        x = (low + high) / 2;
        if (stat_tcdf(x, df) < p) low = x;
        else high = x;
    }
    return x;
}

/* Chi-square PDF */
double stat_chi2pdf(double x, double df)
{
    if (x < 0) return 0;
    if (x == 0) return (df == 2) ? 0.5 : 0;
    return pow(x, df/2 - 1) * exp(-x/2) / (pow(2, df/2) * c89_tgamma(df/2));
}

/* Chi-square CDF using incomplete gamma */
double stat_chi2cdf(double x, double df)
{
    double k = df / 2;
    double sum = 0, term;
    int i;
    
    if (x <= 0) return 0;
    
    /* Lower incomplete gamma / gamma(k) */
    term = pow(x/2, k) * exp(-x/2) / c89_tgamma(k + 1);
    sum = term;
    for (i = 1; i < 100; i++) {
        term *= (x / 2) / (k + i);
        sum += term;
        if (fabs(term) < 1e-15 * fabs(sum)) break;
    }
    return sum;
}

/* Chi-square inverse CDF */
double stat_chi2inv(double p, double df)
{
    double x = df, low = 0, high = df * 10;
    int i;
    
    if (p <= 0) return 0;
    if (p >= 1) return 1e308;
    
    /* Bisection */
    for (i = 0; i < 60; i++) {
        x = (low + high) / 2;
        if (stat_chi2cdf(x, df) < p) low = x;
        else high = x;
    }
    return x;
}

/* F distribution PDF */
double stat_fpdf(double x, double d1, double d2)
{
    double num, den;
    if (x < 0) return 0;
    if (x == 0) return (d1 == 2) ? 1 : 0;
    
    num = pow(d1/d2, d1/2) * pow(x, d1/2 - 1);
    den = pow(1 + d1*x/d2, (d1+d2)/2) * exp(c89_lgamma(d1/2) + c89_lgamma(d2/2) - c89_lgamma((d1+d2)/2));
    return num / den;
}

/* F distribution CDF */
double stat_fcdf(double x, double d1, double d2)
{
    if (x <= 0) return 0;
    return stat_betainc(d1 * x / (d1 * x + d2), d1 / 2, d2 / 2);
}

/* F distribution inverse */
double stat_finv(double p, double d1, double d2)
{
    double x = 1, low = 0, high = 100;
    int i;
    
    if (p <= 0) return 0;
    if (p >= 1) return 1e308;
    
    for (i = 0; i < 60; i++) {
        x = (low + high) / 2;
        if (stat_fcdf(x, d1, d2) < p) low = x;
        else high = x;
    }
    return x;
}

/* Binomial coefficient */
double stat_binomcoef(int n, int k)
{
    if (k < 0 || k > n) return 0;
    if (k == 0 || k == n) return 1;
    return exp(c89_lgamma(n + 1) - c89_lgamma(k + 1) - c89_lgamma(n - k + 1));
}

/* Binomial PDF: P(X = k) */
double stat_binopdf(int k, int n, double p)
{
    if (k < 0 || k > n) return 0;
    return stat_binomcoef(n, k) * pow(p, k) * pow(1 - p, n - k);
}

/* Binomial CDF: P(X <= k) */
double stat_binocdf(int k, int n, double p)
{
    double sum = 0;
    int i;
    for (i = 0; i <= k && i <= n; i++) {
        sum += stat_binopdf(i, n, p);
    }
    return sum;
}

/* Poisson PDF */
double stat_poisspdf(int k, double lambda)
{
    if (k < 0 || lambda < 0) return 0;
    return exp(k * log(lambda) - lambda - c89_lgamma(k + 1));
}

/* Poisson CDF */
double stat_poisscdf(int k, double lambda)
{
    double sum = 0;
    int i;
    for (i = 0; i <= k; i++) {
        sum += stat_poisspdf(i, lambda);
    }
    return sum;
}

/* Exponential PDF */
double stat_exppdf(double x, double lambda)
{
    if (x < 0 || lambda <= 0) return 0;
    return lambda * exp(-lambda * x);
}

/* Exponential CDF */
double stat_expcdf(double x, double lambda)
{
    if (x < 0 || lambda <= 0) return 0;
    return 1 - exp(-lambda * x);
}

/* Exponential inverse CDF */
double stat_expinv(double p, double lambda)
{
    if (p <= 0 || lambda <= 0) return 0;
    if (p >= 1) return 1e308;
    return -log(1 - p) / lambda;
}

/* Uniform PDF */
double stat_unifpdf(double x, double a, double b)
{
    if (x < a || x > b || a >= b) return 0;
    return 1.0 / (b - a);
}

/* Uniform CDF */
double stat_unifcdf(double x, double a, double b)
{
    if (x <= a) return 0;
    if (x >= b) return 1;
    return (x - a) / (b - a);
}

/* Uniform inverse CDF */
double stat_unifinv(double p, double a, double b)
{
    if (p <= 0) return a;
    if (p >= 1) return b;
    return a + p * (b - a);
}

/* ============== HYPOTHESIS TESTS ============== */

/* Two-sample t-test (Welch's) */
void stat_ttest2(const double *x, int nx, const double *y, int ny,
                 double *t_stat, double *p_value, double *df)
{
    double mx = 0, my = 0, vx = 0, vy = 0;
    double se, nu;
    int i;
    
    /* Means */
    for (i = 0; i < nx; i++) mx += x[i];
    for (i = 0; i < ny; i++) my += y[i];
    mx /= nx;
    my /= ny;
    
    /* Variances */
    for (i = 0; i < nx; i++) vx += (x[i] - mx) * (x[i] - mx);
    for (i = 0; i < ny; i++) vy += (y[i] - my) * (y[i] - my);
    vx /= (nx - 1);
    vy /= (ny - 1);
    
    /* Welch's t-test */
    se = sqrt(vx / nx + vy / ny);
    *t_stat = (mx - my) / se;
    
    /* Welch-Satterthwaite degrees of freedom */
    nu = (vx / nx + vy / ny) * (vx / nx + vy / ny);
    nu /= (vx * vx / (nx * nx * (nx - 1)) + vy * vy / (ny * ny * (ny - 1)));
    *df = nu;
    
    /* Two-tailed p-value */
    *p_value = 2 * (1 - stat_tcdf(fabs(*t_stat), nu));
}

/* Z-test (known variance or large sample) */
void stat_ztest(const double *x, int n, double mu0, double sigma,
                double *z_stat, double *p_value)
{
    double mx = 0;
    int i;
    
    for (i = 0; i < n; i++) mx += x[i];
    mx /= n;
    
    *z_stat = (mx - mu0) / (sigma / sqrt(n));
    *p_value = 2 * (1 - stat_normcdf(fabs(*z_stat)));
}

/* Chi-square goodness of fit test */
void stat_chi2test(const double *observed, const double *expected, int n,
                   double *chi2_stat, double *p_value, int *df)
{
    double chi2 = 0;
    int i;
    
    for (i = 0; i < n; i++) {
        if (expected[i] > 0) {
            double d = observed[i] - expected[i];
            chi2 += d * d / expected[i];
        }
    }
    
    *chi2_stat = chi2;
    *df = n - 1;
    *p_value = 1 - stat_chi2cdf(chi2, *df);
}

/* One-way ANOVA */
void stat_anova1(const double *data, const int *groups, int n, int ngroups,
                 double *f_stat, double *p_value)
{
    double *means, grand_mean = 0;
    int *counts;
    double ssb = 0, ssw = 0;  /* Sum of squares between/within */
    int i, df1, df2;
    
    means = (double*)calloc(ngroups, sizeof(double));
    counts = (int*)calloc(ngroups, sizeof(int));
    if (!means || !counts) {
        *f_stat = *p_value = 0;
        free(means); free(counts);
        return;
    }
    
    /* Group means */
    for (i = 0; i < n; i++) {
        int g = groups[i];
        if (g >= 0 && g < ngroups) {
            means[g] += data[i];
            counts[g]++;
            grand_mean += data[i];
        }
    }
    grand_mean /= n;
    for (i = 0; i < ngroups; i++) {
        if (counts[i] > 0) means[i] /= counts[i];
    }
    
    /* Sum of squares */
    for (i = 0; i < ngroups; i++) {
        ssb += counts[i] * (means[i] - grand_mean) * (means[i] - grand_mean);
    }
    for (i = 0; i < n; i++) {
        int g = groups[i];
        if (g >= 0 && g < ngroups) {
            double d = data[i] - means[g];
            ssw += d * d;
        }
    }
    
    df1 = ngroups - 1;
    df2 = n - ngroups;
    
    if (df1 > 0 && df2 > 0 && ssw > 0) {
        *f_stat = (ssb / df1) / (ssw / df2);
        *p_value = 1 - stat_fcdf(*f_stat, df1, df2);
    } else {
        *f_stat = *p_value = 0;
    }
    
    free(means);
    free(counts);
}

/* ============== EFFECT SIZES ============== */

/* Cohen's d (pooled standard deviation) */
double stat_cohend(const double *x, int nx, const double *y, int ny)
{
    double mx = 0, my = 0, vx = 0, vy = 0, sp;
    int i;
    
    for (i = 0; i < nx; i++) mx += x[i];
    for (i = 0; i < ny; i++) my += y[i];
    mx /= nx;
    my /= ny;
    
    for (i = 0; i < nx; i++) vx += (x[i] - mx) * (x[i] - mx);
    for (i = 0; i < ny; i++) vy += (y[i] - my) * (y[i] - my);
    
    /* Pooled standard deviation */
    sp = sqrt((vx + vy) / (nx + ny - 2));
    
    return (mx - my) / sp;
}

/* Pearson correlation r (as effect size) - already have corrcoef */

/* R-squared for regression */
double stat_rsquare(const double *y_actual, const double *y_pred, int n)
{
    double ss_res = 0, ss_tot = 0, mean_y = 0;
    int i;
    
    for (i = 0; i < n; i++) mean_y += y_actual[i];
    mean_y /= n;
    
    for (i = 0; i < n; i++) {
        double res = y_actual[i] - y_pred[i];
        double tot = y_actual[i] - mean_y;
        ss_res += res * res;
        ss_tot += tot * tot;
    }
    
    return 1 - ss_res / ss_tot;
}

/* ============== NUMERICAL METHODS ============== */

/* Bisection method for root finding */
double num_bisect(double (*f)(double), double a, double b, double tol, int maxiter)
{
    double fa = f(a), fb = f(b), mid = (a + b) / 2, fmid;
    int i;
    
    if (fa * fb > 0) return (a + b) / 2;  /* No root bracketed */
    
    for (i = 0; i < maxiter; i++) {
        mid = (a + b) / 2;
        fmid = f(mid);
        
        if (fabs(fmid) < tol || (b - a) / 2 < tol) return mid;
        
        if (fa * fmid < 0) {
            b = mid;
            fb = fmid;
        } else {
            a = mid;
            fa = fmid;
        }
    }
    return mid;
}

/* Secant method */
double num_secant(double (*f)(double), double x0, double x1, double tol, int maxiter)
{
    double f0 = f(x0), f1 = f(x1), x2;
    int i;
    
    for (i = 0; i < maxiter; i++) {
        if (fabs(f1 - f0) < 1e-15) break;
        x2 = x1 - f1 * (x1 - x0) / (f1 - f0);
        
        if (fabs(x2 - x1) < tol) return x2;
        
        x0 = x1;
        f0 = f1;
        x1 = x2;
        f1 = f(x2);
    }
    return x1;
}

/* Simpson's rule integration */
double num_simpson(const double *y, int n, double h)
{
    double sum;
    int i;
    
    if (n < 3) {
        /* Fall back to trapezoidal */
        sum = 0;
        for (i = 0; i < n - 1; i++) sum += (y[i] + y[i + 1]) * h / 2;
        return sum;
    }
    
    sum = y[0] + y[n - 1];
    for (i = 1; i < n - 1; i++) {
        sum += (i % 2 == 1) ? 4 * y[i] : 2 * y[i];
    }
    return sum * h / 3;
}

/* Romberg integration (Richardson extrapolation on trapezoidal rule) */
double num_romberg(double (*f)(double), double a, double b, double tol, int maxiter)
{
    double R[10][10];
    double h;
    int i, j, k, n;
    
    h = b - a;
    R[0][0] = 0.5 * h * (f(a) + f(b));
    
    for (i = 1; i < maxiter && i < 10; i++) {
        h /= 2;
        n = 1 << (i - 1);  /* 2^(i-1) */
        
        /* Trapezoidal with new points */
        R[i][0] = 0.5 * R[i-1][0];
        for (k = 1; k <= n; k++) {
            R[i][0] += h * f(a + (2 * k - 1) * h);
        }
        
        /* Richardson extrapolation */
        for (j = 1; j <= i; j++) {
            double factor = 1 << (2 * j);  /* 4^j */
            R[i][j] = (factor * R[i][j-1] - R[i-1][j-1]) / (factor - 1);
        }
        
        if (i > 1 && fabs(R[i][i] - R[i-1][i-1]) < tol) {
            return R[i][i];
        }
    }
    return R[maxiter < 10 ? maxiter - 1 : 9][maxiter < 10 ? maxiter - 1 : 9];
}

/* Euler method for ODE: dy/dt = f(t, y) */
void num_euler(double (*f)(double, double), double t0, double y0, double tf, 
               int nsteps, double *t_out, double *y_out)
{
    double h = (tf - t0) / nsteps;
    double t = t0, y = y0;
    int i;
    
    t_out[0] = t;
    y_out[0] = y;
    
    for (i = 1; i <= nsteps; i++) {
        y = y + h * f(t, y);
        t = t + h;
        t_out[i] = t;
        y_out[i] = y;
    }
}

/* RK4 (4th order Runge-Kutta) for ODE */
void num_rk4(double (*f)(double, double), double t0, double y0, double tf,
             int nsteps, double *t_out, double *y_out)
{
    double h = (tf - t0) / nsteps;
    double t = t0, y = y0;
    double k1, k2, k3, k4;
    int i;
    
    t_out[0] = t;
    y_out[0] = y;
    
    for (i = 1; i <= nsteps; i++) {
        k1 = h * f(t, y);
        k2 = h * f(t + h/2, y + k1/2);
        k3 = h * f(t + h/2, y + k2/2);
        k4 = h * f(t + h, y + k3);
        
        y = y + (k1 + 2*k2 + 2*k3 + k4) / 6;
        t = t + h;
        t_out[i] = t;
        y_out[i] = y;
    }
}

/* Golden section search for minimum */
double num_golden(double (*f)(double), double a, double b, double tol, int maxiter)
{
    double phi = (1 + sqrt(5)) / 2;
    double resphi = 2 - phi;
    double x1 = a + resphi * (b - a);
    double x2 = b - resphi * (b - a);
    double f1 = f(x1), f2 = f(x2);
    int i;
    
    for (i = 0; i < maxiter && (b - a) > tol; i++) {
        if (f1 < f2) {
            b = x2;
            x2 = x1;
            f2 = f1;
            x1 = a + resphi * (b - a);
            f1 = f(x1);
        } else {
            a = x1;
            x1 = x2;
            f1 = f2;
            x2 = b - resphi * (b - a);
            f2 = f(x2);
        }
    }
    return (a + b) / 2;
}

/* Shuffle array (Fisher-Yates) */
void num_shuffle(double *arr, int n)
{
    int i, j;
    double tmp;
    
    for (i = n - 1; i > 0; i--) {
        j = rand() % (i + 1);
        tmp = arr[i];
        arr[i] = arr[j];
        arr[j] = tmp;
    }
}

/* Random sample without replacement */
void num_randsample(const double *pop, int npop, double *sample, int nsample)
{
    int *indices;
    int i, j, tmp;
    
    if (nsample > npop) nsample = npop;
    
    indices = (int*)malloc(npop * sizeof(int));
    if (!indices) return;
    
    for (i = 0; i < npop; i++) indices[i] = i;
    
    /* Partial Fisher-Yates */
    for (i = 0; i < nsample; i++) {
        j = i + rand() % (npop - i);
        tmp = indices[i];
        indices[i] = indices[j];
        indices[j] = tmp;
        sample[i] = pop[indices[i]];
    }
    
    free(indices);
}
/* ============================================================
 * HYPOTHESIS TESTING
 * ============================================================ */

/* ztest(x, mu, sigma) - One-sample z-test
 * Returns z-statistic. Use normcdf to get p-value.
 */
double sf_ztest_stat(double sample_mean, double pop_mean, double pop_std, double n)
{
    if (pop_std <= 0 || n <= 0) return 0.0;
    return (sample_mean - pop_mean) / (pop_std / sqrt(n));
}

/* chi2gof - Chi-square goodness of fit statistic
 * observed and expected are vectors
 */
void mat_chi2stat(matrix_t *r, const matrix_t *observed, const matrix_t *expected)
{
    int n = observed->rows * observed->cols;
    int i;
    double chi2 = 0;
    
    mat_zero(r, 1, 1);
    
    for (i = 0; i < n && i < expected->rows * expected->cols; i++) {
        double o = apf_to_double(&observed->data[i].re);
        double e = apf_to_double(&expected->data[i].re);
        if (e > 0) {
            chi2 += (o - e) * (o - e) / e;
        }
    }
    apf_from_double(&r->data[0].re, chi2);
}

/* ttest2 - Two-sample t-test (independent samples)
 * Returns [t-stat, p-value, df]
 */
void mat_ttest2(matrix_t *r, const matrix_t *x, const matrix_t *y)
{
    int nx = x->rows * x->cols;
    int ny = y->rows * y->cols;
    int i;
    double mean_x = 0, mean_y = 0;
    double var_x = 0, var_y = 0;
    double sp, t, df;
    extern double stat_tcdf(double, double);
    
    mat_zero(r, 1, 3);
    if (nx < 2 || ny < 2) return;
    
    /* Compute means */
    for (i = 0; i < nx; i++) mean_x += apf_to_double(&x->data[i].re);
    for (i = 0; i < ny; i++) mean_y += apf_to_double(&y->data[i].re);
    mean_x /= nx;
    mean_y /= ny;
    
    /* Compute variances */
    for (i = 0; i < nx; i++) {
        double d = apf_to_double(&x->data[i].re) - mean_x;
        var_x += d * d;
    }
    for (i = 0; i < ny; i++) {
        double d = apf_to_double(&y->data[i].re) - mean_y;
        var_y += d * d;
    }
    var_x /= (nx - 1);
    var_y /= (ny - 1);
    
    /* Pooled standard error */
    sp = sqrt(var_x / nx + var_y / ny);
    if (sp < 1e-15) return;
    
    /* Welch's t-test (unequal variances) */
    t = (mean_x - mean_y) / sp;
    
    /* Welch-Satterthwaite degrees of freedom */
    {
        double v1 = var_x / nx;
        double v2 = var_y / ny;
        df = (v1 + v2) * (v1 + v2) / 
             (v1 * v1 / (nx - 1) + v2 * v2 / (ny - 1));
    }
    
    apf_from_double(&r->data[0].re, t);
    apf_from_double(&r->data[1].re, 2 * (1 - stat_tcdf(fabs(t), df)));  /* Two-tailed p-value */
    apf_from_double(&r->data[2].re, df);
}

/* anova1 - One-way ANOVA
 * Input: matrix where each column is a group
 * Returns [F-statistic, p-value, df_between, df_within]
 */
void mat_anova1(matrix_t *r, const matrix_t *data)
{
    int k = data->cols;  /* Number of groups */
    int n = data->rows;  /* Observations per group */
    int i, j;
    double grand_mean = 0, ss_between = 0, ss_within = 0;
    double *group_means;
    double F, df1, df2;
    
    mat_zero(r, 1, 4);
    if (k < 2 || n < 2) return;
    
    group_means = (double*)malloc(k * sizeof(double));
    if (!group_means) return;
    
    /* Compute group means and grand mean */
    for (j = 0; j < k; j++) {
        group_means[j] = 0;
        for (i = 0; i < n; i++) {
            group_means[j] += apf_to_double(&MAT_AT(data, i, j).re);
        }
        group_means[j] /= n;
        grand_mean += group_means[j];
    }
    grand_mean /= k;
    
    /* Sum of squares between groups */
    for (j = 0; j < k; j++) {
        ss_between += n * (group_means[j] - grand_mean) * (group_means[j] - grand_mean);
    }
    
    /* Sum of squares within groups */
    for (j = 0; j < k; j++) {
        for (i = 0; i < n; i++) {
            double d = apf_to_double(&MAT_AT(data, i, j).re) - group_means[j];
            ss_within += d * d;
        }
    }
    
    free(group_means);
    
    df1 = k - 1;
    df2 = k * (n - 1);
    
    if (ss_within < 1e-15) {
        F = INFINITY;
    } else {
        F = (ss_between / df1) / (ss_within / df2);
    }
    
    apf_from_double(&r->data[0].re, F);
    apf_from_double(&r->data[1].re, 1 - stat_fcdf(F, df1, df2));
    apf_from_double(&r->data[2].re, df1);
    apf_from_double(&r->data[3].re, df2);
}

/* ============================================================
 * CORRELATION MEASURES
 * ============================================================ */

/* spearman - Spearman rank correlation */
void mat_spearman(matrix_t *r, const matrix_t *x, const matrix_t *y)
{
    int n = x->rows * x->cols;
    int i, j;
    double *rank_x, *rank_y;
    double mean_rx = 0, mean_ry = 0;
    double sum_xy = 0, sum_xx = 0, sum_yy = 0;
    double rho;
    
    mat_zero(r, 1, 1);
    if (n < 2 || n != y->rows * y->cols) return;
    
    rank_x = (double*)malloc(n * sizeof(double));
    rank_y = (double*)malloc(n * sizeof(double));
    if (!rank_x || !rank_y) {
        if (rank_x) free(rank_x);
        if (rank_y) free(rank_y);
        return;
    }
    
    /* Compute ranks for x */
    for (i = 0; i < n; i++) {
        double xi = apf_to_double(&x->data[i].re);
        double rank = 1;
        for (j = 0; j < n; j++) {
            if (apf_to_double(&x->data[j].re) < xi) rank++;
        }
        /* Handle ties by averaging ranks */
        {
            int ties = 0;
            for (j = 0; j < n; j++) {
                if (fabs(apf_to_double(&x->data[j].re) - xi) < 1e-12) ties++;
            }
            rank_x[i] = rank + (ties - 1) / 2.0;
        }
    }
    
    /* Compute ranks for y */
    for (i = 0; i < n; i++) {
        double yi = apf_to_double(&y->data[i].re);
        double rank = 1;
        for (j = 0; j < n; j++) {
            if (apf_to_double(&y->data[j].re) < yi) rank++;
        }
        {
            int ties = 0;
            for (j = 0; j < n; j++) {
                if (fabs(apf_to_double(&y->data[j].re) - yi) < 1e-12) ties++;
            }
            rank_y[i] = rank + (ties - 1) / 2.0;
        }
    }
    
    /* Pearson correlation of ranks */
    for (i = 0; i < n; i++) {
        mean_rx += rank_x[i];
        mean_ry += rank_y[i];
    }
    mean_rx /= n;
    mean_ry /= n;
    
    for (i = 0; i < n; i++) {
        double dx = rank_x[i] - mean_rx;
        double dy = rank_y[i] - mean_ry;
        sum_xy += dx * dy;
        sum_xx += dx * dx;
        sum_yy += dy * dy;
    }
    
    free(rank_x);
    free(rank_y);
    
    if (sum_xx * sum_yy < 1e-15) {
        rho = 0;
    } else {
        rho = sum_xy / sqrt(sum_xx * sum_yy);
    }
    
    apf_from_double(&r->data[0].re, rho);
}

/* kendall - Kendall tau correlation */
void mat_kendall(matrix_t *r, const matrix_t *x, const matrix_t *y)
{
    int n = x->rows * x->cols;
    int i, j;
    int concordant = 0, discordant = 0;
    double tau;
    
    mat_zero(r, 1, 1);
    if (n < 2 || n != y->rows * y->cols) return;
    
    for (i = 0; i < n - 1; i++) {
        double xi = apf_to_double(&x->data[i].re);
        double yi = apf_to_double(&y->data[i].re);
        for (j = i + 1; j < n; j++) {
            double xj = apf_to_double(&x->data[j].re);
            double yj = apf_to_double(&y->data[j].re);
            double sign_x = (xj > xi) - (xj < xi);
            double sign_y = (yj > yi) - (yj < yi);
            double product = sign_x * sign_y;
            if (product > 0) concordant++;
            else if (product < 0) discordant++;
        }
    }
    
    tau = (double)(concordant - discordant) / (n * (n - 1) / 2);
    apf_from_double(&r->data[0].re, tau);
}

/* ============================================================
 * EFFECT SIZE
 * ============================================================ */

/* cohend - Cohen's d effect size */
void mat_cohend(matrix_t *r, const matrix_t *x, const matrix_t *y)
{
    int nx = x->rows * x->cols;
    int ny = y->rows * y->cols;
    int i;
    double mean_x = 0, mean_y = 0;
    double var_x = 0, var_y = 0;
    double pooled_std, d;
    
    mat_zero(r, 1, 1);
    if (nx < 2 || ny < 2) return;
    
    /* Compute means */
    for (i = 0; i < nx; i++) mean_x += apf_to_double(&x->data[i].re);
    for (i = 0; i < ny; i++) mean_y += apf_to_double(&y->data[i].re);
    mean_x /= nx;
    mean_y /= ny;
    
    /* Compute variances */
    for (i = 0; i < nx; i++) {
        double diff = apf_to_double(&x->data[i].re) - mean_x;
        var_x += diff * diff;
    }
    for (i = 0; i < ny; i++) {
        double diff = apf_to_double(&y->data[i].re) - mean_y;
        var_y += diff * diff;
    }
    var_x /= (nx - 1);
    var_y /= (ny - 1);
    
    /* Pooled standard deviation */
    pooled_std = sqrt(((nx - 1) * var_x + (ny - 1) * var_y) / (nx + ny - 2));
    
    if (pooled_std < 1e-15) {
        d = 0;
    } else {
        d = (mean_x - mean_y) / pooled_std;
    }
    
    apf_from_double(&r->data[0].re, d);
}

/* ============================================================
 * NON-PARAMETRIC TESTS
 * ============================================================ */

/* ranksum - Mann-Whitney U test (Wilcoxon rank-sum)
 * Returns [U-statistic, z-score]
 */
void mat_ranksum(matrix_t *r, const matrix_t *x, const matrix_t *y)
{
    int nx = x->rows * x->cols;
    int ny = y->rows * y->cols;
    int n = nx + ny;
    int i, j;
    double *combined, *ranks;
    int *group;  /* 0 for x, 1 for y */
    double rank_sum_x = 0;
    double U, mean_U, std_U, z;
    
    mat_zero(r, 1, 2);
    if (nx < 1 || ny < 1) return;
    
    combined = (double*)malloc(n * sizeof(double));
    ranks = (double*)malloc(n * sizeof(double));
    group = (int*)malloc(n * sizeof(int));
    if (!combined || !ranks || !group) {
        if (combined) free(combined);
        if (ranks) free(ranks);
        if (group) free(group);
        return;
    }
    
    /* Combine data */
    for (i = 0; i < nx; i++) {
        combined[i] = apf_to_double(&x->data[i].re);
        group[i] = 0;
    }
    for (i = 0; i < ny; i++) {
        combined[nx + i] = apf_to_double(&y->data[i].re);
        group[nx + i] = 1;
    }
    
    /* Compute ranks */
    for (i = 0; i < n; i++) {
        double rank = 1;
        int ties = 0;
        for (j = 0; j < n; j++) {
            if (combined[j] < combined[i]) rank++;
            if (fabs(combined[j] - combined[i]) < 1e-12) ties++;
        }
        ranks[i] = rank + (ties - 1) / 2.0;
    }
    
    /* Sum of ranks for group x */
    for (i = 0; i < n; i++) {
        if (group[i] == 0) rank_sum_x += ranks[i];
    }
    
    /* U statistic */
    U = rank_sum_x - nx * (nx + 1) / 2.0;
    
    /* Normal approximation for large samples */
    mean_U = nx * ny / 2.0;
    std_U = sqrt(nx * ny * (n + 1) / 12.0);
    z = (std_U > 0) ? (U - mean_U) / std_U : 0;
    
    free(combined);
    free(ranks);
    free(group);
    
    apf_from_double(&r->data[0].re, U);
    apf_from_double(&r->data[1].re, z);
}

/* signrank - Wilcoxon signed-rank test
 * For paired samples or one-sample against hypothesized median
 * Returns [W-statistic, z-score]
 */
void mat_signrank(matrix_t *r, const matrix_t *x, const matrix_t *y)
{
    int n = x->rows * x->cols;
    int i, j, count;
    double *diff, *absrank;
    double W_pos = 0, W_neg = 0, W;
    double mean_W, std_W, z;
    
    mat_zero(r, 1, 2);
    if (n < 1) return;
    
    diff = (double*)malloc(n * sizeof(double));
    absrank = (double*)malloc(n * sizeof(double));
    if (!diff || !absrank) {
        if (diff) free(diff);
        if (absrank) free(absrank);
        return;
    }
    
    /* Compute differences */
    count = 0;
    for (i = 0; i < n; i++) {
        double d;
        if (y && y->rows * y->cols > i) {
            d = apf_to_double(&x->data[i].re) - apf_to_double(&y->data[i].re);
        } else {
            d = apf_to_double(&x->data[i].re);  /* One-sample against 0 */
        }
        if (fabs(d) > 1e-12) {  /* Exclude zeros */
            diff[count++] = d;
        }
    }
    
    if (count < 1) {
        free(diff);
        free(absrank);
        return;
    }
    
    /* Rank absolute differences */
    for (i = 0; i < count; i++) {
        double rank = 1;
        int ties = 0;
        double absi = fabs(diff[i]);
        for (j = 0; j < count; j++) {
            double absj = fabs(diff[j]);
            if (absj < absi) rank++;
            if (fabs(absj - absi) < 1e-12) ties++;
        }
        absrank[i] = rank + (ties - 1) / 2.0;
    }
    
    /* Sum positive and negative ranks */
    for (i = 0; i < count; i++) {
        if (diff[i] > 0) W_pos += absrank[i];
        else W_neg += absrank[i];
    }
    
    W = (W_pos < W_neg) ? W_pos : W_neg;
    
    /* Normal approximation */
    mean_W = count * (count + 1) / 4.0;
    std_W = sqrt(count * (count + 1) * (2 * count + 1) / 24.0);
    z = (std_W > 0) ? (W - mean_W) / std_W : 0;
    
    free(diff);
    free(absrank);
    
    apf_from_double(&r->data[0].re, W);
    apf_from_double(&r->data[1].re, z);
}
