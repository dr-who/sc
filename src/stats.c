/* stats.c - Statistical functions
 * C89 portable for DOS, Linux
 */

#include "config.h"

#ifdef HAVE_STATS

#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include "stats.h"
#include "apf.h"

/* Global statistical registers */
stats_regs_t stats_regs;

void stats_clear(void) {
    stats_regs.n = 0;
    apf_zero(&stats_regs.sum_x);
    apf_zero(&stats_regs.sum_x2);
    apf_zero(&stats_regs.sum_y);
    apf_zero(&stats_regs.sum_y2);
    apf_zero(&stats_regs.sum_xy);
}

void stats_add(const apf *x) {
    apf tmp;
    
    stats_regs.n++;
    
    /* sum_x += x */
    apf_add(&tmp, &stats_regs.sum_x, x);
    apf_copy(&stats_regs.sum_x, &tmp);
    
    /* sum_x2 += x² */
    apf_mul(&tmp, x, x);
    apf_add(&stats_regs.sum_x2, &stats_regs.sum_x2, &tmp);
}

void stats_add2(const apf *x, const apf *y) {
    apf tmp;
    
    stats_regs.n++;
    
    /* sum_x += x */
    apf_add(&tmp, &stats_regs.sum_x, x);
    apf_copy(&stats_regs.sum_x, &tmp);
    
    /* sum_x2 += x² */
    apf_mul(&tmp, x, x);
    apf_add(&stats_regs.sum_x2, &stats_regs.sum_x2, &tmp);
    
    /* sum_y += y */
    apf_add(&tmp, &stats_regs.sum_y, y);
    apf_copy(&stats_regs.sum_y, &tmp);
    
    /* sum_y2 += y² */
    apf_mul(&tmp, y, y);
    apf_add(&stats_regs.sum_y2, &stats_regs.sum_y2, &tmp);
    
    /* sum_xy += x*y */
    apf_mul(&tmp, x, y);
    apf_add(&stats_regs.sum_xy, &stats_regs.sum_xy, &tmp);
}

void stats_sub(const apf *x) {
    apf tmp;
    
    if (stats_regs.n > 0) {
        stats_regs.n--;
        
        /* sum_x -= x */
        apf_sub(&tmp, &stats_regs.sum_x, x);
        apf_copy(&stats_regs.sum_x, &tmp);
        
        /* sum_x2 -= x² */
        apf_mul(&tmp, x, x);
        apf_sub(&stats_regs.sum_x2, &stats_regs.sum_x2, &tmp);
    }
}

int stats_count(void) {
    return stats_regs.n;
}

void stats_sum(apf *r) {
    apf_copy(r, &stats_regs.sum_x);
}

void stats_mean(apf *r) {
    apf n;
    if (stats_regs.n == 0) {
        apf_zero(r);
        return;
    }
    apf_from_int(&n, stats_regs.n);
    apf_div(r, &stats_regs.sum_x, &n);
}

/* Sample variance: (Σx² - (Σx)²/n) / (n-1) */
void stats_var(apf *r) {
    apf n, n1, sum2_over_n, numer;
    
    if (stats_regs.n < 2) {
        apf_zero(r);
        return;
    }
    
    apf_from_int(&n, stats_regs.n);
    apf_from_int(&n1, stats_regs.n - 1);
    
    /* (Σx)² / n */
    apf_mul(&sum2_over_n, &stats_regs.sum_x, &stats_regs.sum_x);
    apf_div(&sum2_over_n, &sum2_over_n, &n);
    
    /* Σx² - (Σx)²/n */
    apf_sub(&numer, &stats_regs.sum_x2, &sum2_over_n);
    
    /* Divide by n-1 */
    apf_div(r, &numer, &n1);
}

/* Sample standard deviation: sqrt(variance) */
void stats_sdev(apf *r) {
    stats_var(r);
    apf_sqrt(r, r);
}

/* Population standard deviation: sqrt((Σx² - (Σx)²/n) / n) */
void stats_sdevp(apf *r) {
    apf n, sum2_over_n, numer;
    
    if (stats_regs.n < 1) {
        apf_zero(r);
        return;
    }
    
    apf_from_int(&n, stats_regs.n);
    
    /* (Σx)² / n */
    apf_mul(&sum2_over_n, &stats_regs.sum_x, &stats_regs.sum_x);
    apf_div(&sum2_over_n, &sum2_over_n, &n);
    
    /* Σx² - (Σx)²/n */
    apf_sub(&numer, &stats_regs.sum_x2, &sum2_over_n);
    
    /* Divide by n */
    apf_div(r, &numer, &n);
    
    /* sqrt */
    apf_sqrt(r, r);
}

/* Linear regression slope: (nΣxy - ΣxΣy) / (nΣx² - (Σx)²) */
void stats_slope(apf *r) {
    apf n, numer, denom, tmp1;
    
    if (stats_regs.n < 2) {
        apf_zero(r);
        return;
    }
    
    apf_from_int(&n, stats_regs.n);
    
    /* nΣxy */
    apf_mul(&numer, &n, &stats_regs.sum_xy);
    /* ΣxΣy */
    apf_mul(&tmp1, &stats_regs.sum_x, &stats_regs.sum_y);
    /* nΣxy - ΣxΣy */
    apf_sub(&numer, &numer, &tmp1);
    
    /* nΣx² */
    apf_mul(&denom, &n, &stats_regs.sum_x2);
    /* (Σx)² */
    apf_mul(&tmp1, &stats_regs.sum_x, &stats_regs.sum_x);
    /* nΣx² - (Σx)² */
    apf_sub(&denom, &denom, &tmp1);
    
    apf_div(r, &numer, &denom);
}

/* Y-intercept: (Σy - slope*Σx) / n */
void stats_intercept(apf *r) {
    apf slope, n, tmp;
    
    if (stats_regs.n < 2) {
        apf_zero(r);
        return;
    }
    
    stats_slope(&slope);
    apf_from_int(&n, stats_regs.n);
    
    /* slope * Σx */
    apf_mul(&tmp, &slope, &stats_regs.sum_x);
    /* Σy - slope*Σx */
    apf_sub(r, &stats_regs.sum_y, &tmp);
    /* Divide by n */
    apf_div(r, r, &n);
}

/* Correlation coefficient r */
void stats_corr(apf *r) {
    apf n, numer, denom_x, denom_y, denom, tmp1;
    
    if (stats_regs.n < 2) {
        apf_zero(r);
        return;
    }
    
    apf_from_int(&n, stats_regs.n);
    
    /* Numerator: nΣxy - ΣxΣy */
    apf_mul(&numer, &n, &stats_regs.sum_xy);
    apf_mul(&tmp1, &stats_regs.sum_x, &stats_regs.sum_y);
    apf_sub(&numer, &numer, &tmp1);
    
    /* Denominator: sqrt((nΣx² - (Σx)²)(nΣy² - (Σy)²)) */
    apf_mul(&denom_x, &n, &stats_regs.sum_x2);
    apf_mul(&tmp1, &stats_regs.sum_x, &stats_regs.sum_x);
    apf_sub(&denom_x, &denom_x, &tmp1);
    
    apf_mul(&denom_y, &n, &stats_regs.sum_y2);
    apf_mul(&tmp1, &stats_regs.sum_y, &stats_regs.sum_y);
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

void cmd_stats(const char *args) {
    char buf[64];
    apf val;
    
    /* Skip whitespace */
    while (*args && isspace((unsigned char)*args)) args++;
    
    if (*args == '\0') {
        /* Show current stats */
        printf("n    = %d\n", stats_regs.n);
        apf_to_str(buf, sizeof(buf), &stats_regs.sum_x, 10);
        printf("Σx   = %s\n", buf);
        stats_mean(&val);
        apf_to_str(buf, sizeof(buf), &val, 10);
        printf("mean = %s\n", buf);
        if (stats_regs.n >= 2) {
            stats_sdev(&val);
            apf_to_str(buf, sizeof(buf), &val, 10);
            printf("sdev = %s\n", buf);
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
    
    if (strcmp(args, "sdev") == 0 || strcmp(args, "sd") == 0) {
        stats_sdev(&val);
        apf_to_str(buf, sizeof(buf), &val, 10);
        printf("sdev = %s\n", buf);
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
    
    printf("Usage: stat [clr|mean|sdev|sum|n]\n");
    printf("  Enter data with: data+ x or data+ x y\n");
}

#endif /* HAVE_STATS */
