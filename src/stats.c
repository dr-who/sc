/* stats.c - Statistical functions
 * C89 portable for DOS, Linux
 * Uses __far for large arrays to stay under 64KB DGROUP limit on DOS
 */

#include "config.h"

#ifdef HAVE_STATS

#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include "stats.h"
#include "apf.h"
#include "apfx.h"

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
