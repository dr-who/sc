/* tvm.c - Time Value of Money functions
 * C89 portable for DOS, Linux
 * 
 * TVM equation (end mode):
 *   PV + PMT * ((1-(1+i)^-n)/i) + FV * (1+i)^-n = 0
 * 
 * TVM equation (begin mode):
 *   PV + PMT * (1+i) * ((1-(1+i)^-n)/i) + FV * (1+i)^-n = 0
 */

#include "config.h"

#ifdef HAVE_TVM

#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include "tvm.h"
#include "apf.h"
#include "apfx.h"

/* External display digits */
extern int display_digits;

/* Global TVM registers */
tvm_regs_t tvm_regs;

void tvm_clear(void) {
    apf_zero(&tvm_regs.n);
    apf_zero(&tvm_regs.i);
    apf_zero(&tvm_regs.pv);
    apf_zero(&tvm_regs.pmt);
    apf_zero(&tvm_regs.fv);
    tvm_regs.beg = 0;
}

void tvm_set_n(const apf *val) {
    apf_copy(&tvm_regs.n, val);
}

void tvm_set_i(const apf *val) {
    /* Convert from percent to decimal */
    apf hundred;
    apf_from_int(&hundred, 100);
    apf_div(&tvm_regs.i, val, &hundred);
}

void tvm_set_pv(const apf *val) {
    apf_copy(&tvm_regs.pv, val);
}

void tvm_set_pmt(const apf *val) {
    apf_copy(&tvm_regs.pmt, val);
}

void tvm_set_fv(const apf *val) {
    apf_copy(&tvm_regs.fv, val);
}

void tvm_set_beg(int beg) {
    tvm_regs.beg = beg ? 1 : 0;
}

/* Helper: Calculate (1+i)^n */
static void calc_compound(apf *r, const apf *i, const apf *n) {
    apf one_plus_i;
    apf_from_int(&one_plus_i, 1);
    apf_add(&one_plus_i, &one_plus_i, i);
    apfx_pow(r, &one_plus_i, n);
}

/* Helper: Calculate annuity factor ((1-(1+i)^-n)/i) */
static void calc_annuity(apf *r, const apf *i, const apf *n) {
    apf compound, neg_n, one, numer;
    
    /* Check for i = 0 */
    if (apf_is_zero(i)) {
        apf_copy(r, n);  /* Annuity = n when i = 0 */
        return;
    }
    
    /* (1+i)^-n */
    apf_neg(&neg_n, n);
    calc_compound(&compound, i, &neg_n);
    
    /* 1 - (1+i)^-n */
    apf_from_int(&one, 1);
    apf_sub(&numer, &one, &compound);
    
    /* Divide by i */
    apf_div(r, &numer, i);
}

/* Calculate N given I%, PV, PMT, FV */
int tvm_calc_n(apf *result) {
    apf i, pv, pmt, fv, one, one_plus_i;
    apf numer, denom, ratio, ln_ratio, ln_1pi;
    
    i = tvm_regs.i;
    pv = tvm_regs.pv;
    pmt = tvm_regs.pmt;
    fv = tvm_regs.fv;
    
    /* N = ln((PMT - FV*i) / (PMT + PV*i)) / ln(1+i) */
    apf_from_int(&one, 1);
    apf_add(&one_plus_i, &one, &i);
    
    /* Begin mode adjustment */
    if (tvm_regs.beg) {
        apf_mul(&pmt, &pmt, &one_plus_i);
    }
    
    /* Numerator: PMT - FV*i */
    apf_mul(&numer, &fv, &i);
    apf_sub(&numer, &pmt, &numer);
    
    /* Denominator: PMT + PV*i */
    apf_mul(&denom, &pv, &i);
    apf_add(&denom, &pmt, &denom);
    
    /* ratio */
    apf_div(&ratio, &numer, &denom);
    
    /* ln(ratio) / ln(1+i) */
    apfx_log(&ln_ratio, &ratio);
    apfx_log(&ln_1pi, &one_plus_i);
    apf_div(result, &ln_ratio, &ln_1pi);
    
    return 1;
}

/* Calculate I% given N, PV, PMT, FV (uses Newton-Raphson) */
int tvm_calc_i(apf *result) {
    apf n, pv, pmt, fv;
    apf i, i_new, f, df, delta;
    apf one, compound, neg_n, annuity;
    apf tmp1, tmp2;
    apf hundred, tolerance;
    int iter;
    
    n = tvm_regs.n;
    pv = tvm_regs.pv;
    pmt = tvm_regs.pmt;
    fv = tvm_regs.fv;
    
    apf_from_int(&one, 1);
    apf_from_int(&hundred, 100);
    
    /* Initial guess: 10% */
    apf_from_int(&i, 10);
    apf_div(&i, &i, &hundred);
    
    /* Tolerance */
    apf_from_str(&tolerance, "1e-12");
    
    /* Newton-Raphson iteration */
    for (iter = 0; iter < 100; iter++) {
        apf one_plus_i;
        
        apf_add(&one_plus_i, &one, &i);
        
        /* (1+i)^-n */
        apf_neg(&neg_n, &n);
        apfx_pow(&compound, &one_plus_i, &neg_n);
        
        /* Annuity factor */
        if (!apf_is_zero(&i)) {
            apf_sub(&tmp1, &one, &compound);
            apf_div(&annuity, &tmp1, &i);
        } else {
            apf_copy(&annuity, &n);
        }
        
        /* f(i) = PV + PMT*annuity*(1+beg*i) + FV*compound */
        if (tvm_regs.beg) {
            apf_mul(&tmp1, &pmt, &one_plus_i);
        } else {
            apf_copy(&tmp1, &pmt);
        }
        apf_mul(&tmp1, &tmp1, &annuity);
        apf_add(&f, &pv, &tmp1);
        apf_mul(&tmp2, &fv, &compound);
        apf_add(&f, &f, &tmp2);
        
        /* Check convergence */
        apf_abs(&tmp1, &f);
        if (apf_cmp(&tmp1, &tolerance) < 0) {
            break;
        }
        
        /* Compute derivative numerically */
        {
            apf i_plus, f_plus, h;
            apf_from_str(&h, "1e-8");
            apf_add(&i_plus, &i, &h);
            
            apf_add(&one_plus_i, &one, &i_plus);
            apfx_pow(&compound, &one_plus_i, &neg_n);
            if (!apf_is_zero(&i_plus)) {
                apf_sub(&tmp1, &one, &compound);
                apf_div(&annuity, &tmp1, &i_plus);
            }
            
            if (tvm_regs.beg) {
                apf_mul(&tmp1, &pmt, &one_plus_i);
            } else {
                apf_copy(&tmp1, &pmt);
            }
            apf_mul(&tmp1, &tmp1, &annuity);
            apf_add(&f_plus, &pv, &tmp1);
            apf_mul(&tmp2, &fv, &compound);
            apf_add(&f_plus, &f_plus, &tmp2);
            
            apf_sub(&df, &f_plus, &f);
            apf_div(&df, &df, &h);
        }
        
        /* i_new = i - f/df */
        apf_div(&delta, &f, &df);
        apf_sub(&i_new, &i, &delta);
        apf_copy(&i, &i_new);
    }
    
    /* Convert to percent */
    apf_mul(result, &i, &hundred);
    return 1;
}

/* Calculate PV given N, I%, PMT, FV */
int tvm_calc_pv(apf *result) {
    apf annuity, compound, neg_n, one, one_plus_i;
    apf tmp1, tmp2;
    
    apf_from_int(&one, 1);
    apf_neg(&neg_n, &tvm_regs.n);
    
    /* (1+i)^-n */
    calc_compound(&compound, &tvm_regs.i, &neg_n);
    
    /* Annuity factor */
    calc_annuity(&annuity, &tvm_regs.i, &tvm_regs.n);
    
    /* PV = -PMT*annuity*(1+beg*i) - FV*compound */
    if (tvm_regs.beg) {
        apf_add(&one_plus_i, &one, &tvm_regs.i);
        apf_mul(&tmp1, &tvm_regs.pmt, &one_plus_i);
    } else {
        apf_copy(&tmp1, &tvm_regs.pmt);
    }
    apf_mul(&tmp1, &tmp1, &annuity);
    apf_neg(&tmp1, &tmp1);
    
    apf_mul(&tmp2, &tvm_regs.fv, &compound);
    apf_neg(&tmp2, &tmp2);
    
    apf_add(result, &tmp1, &tmp2);
    return 1;
}

/* Calculate PMT given N, I%, PV, FV */
int tvm_calc_pmt(apf *result) {
    apf annuity, compound, neg_n, one, one_plus_i;
    apf numer, denom, tmp;
    
    apf_from_int(&one, 1);
    apf_neg(&neg_n, &tvm_regs.n);
    
    /* (1+i)^-n */
    calc_compound(&compound, &tvm_regs.i, &neg_n);
    
    /* Annuity factor */
    calc_annuity(&annuity, &tvm_regs.i, &tvm_regs.n);
    
    /* PMT = -(PV + FV*compound) / (annuity * (1+beg*i)) */
    apf_mul(&tmp, &tvm_regs.fv, &compound);
    apf_add(&numer, &tvm_regs.pv, &tmp);
    apf_neg(&numer, &numer);
    
    if (tvm_regs.beg) {
        apf_add(&one_plus_i, &one, &tvm_regs.i);
        apf_mul(&denom, &annuity, &one_plus_i);
    } else {
        apf_copy(&denom, &annuity);
    }
    
    apf_div(result, &numer, &denom);
    return 1;
}

/* Calculate FV given N, I%, PV, PMT */
int tvm_calc_fv(apf *result) {
    apf annuity, compound, neg_n, one, one_plus_i;
    apf tmp1, tmp2;
    
    apf_from_int(&one, 1);
    apf_neg(&neg_n, &tvm_regs.n);
    
    /* (1+i)^-n = 1 / (1+i)^n */
    calc_compound(&compound, &tvm_regs.i, &neg_n);
    
    /* Annuity factor */
    calc_annuity(&annuity, &tvm_regs.i, &tvm_regs.n);
    
    /* FV = -(PV + PMT*annuity*(1+beg*i)) / compound */
    if (tvm_regs.beg) {
        apf_add(&one_plus_i, &one, &tvm_regs.i);
        apf_mul(&tmp1, &tvm_regs.pmt, &one_plus_i);
    } else {
        apf_copy(&tmp1, &tvm_regs.pmt);
    }
    apf_mul(&tmp1, &tmp1, &annuity);
    apf_add(&tmp2, &tvm_regs.pv, &tmp1);
    apf_neg(&tmp2, &tmp2);
    
    apf_div(result, &tmp2, &compound);
    return 1;
}

void tvm_show(void) {
    char buf[64];
    apf tmp, hundred;
    int digits = display_digits ? display_digits : 6;
    
    apf_from_int(&hundred, 100);
    
    printf("TVM Registers (%s mode):\n", tvm_regs.beg ? "BEGIN" : "END");
    
    apf_to_str(buf, sizeof(buf), &tvm_regs.n, digits);
    printf("  N   = %s\n", buf);
    
    apf_mul(&tmp, &tvm_regs.i, &hundred);
    apf_to_str(buf, sizeof(buf), &tmp, digits);
    printf("  I%%  = %s%%\n", buf);
    
    apf_to_str(buf, sizeof(buf), &tvm_regs.pv, digits);
    printf("  PV  = %s\n", buf);
    
    apf_to_str(buf, sizeof(buf), &tvm_regs.pmt, digits);
    printf("  PMT = %s\n", buf);
    
    apf_to_str(buf, sizeof(buf), &tvm_regs.fv, digits);
    printf("  FV  = %s\n", buf);
}

void cmd_tvm(const char *args) {
    char buf[64];
    apf val, result;
    const char *p = args;
    int digits = display_digits ? display_digits : 6;
    
    /* Skip whitespace */
    while (*p && isspace((unsigned char)*p)) p++;
    
    if (*p == '\0') {
        tvm_show();
        return;
    }
    
    if (strcmp(p, "clr") == 0 || strcmp(p, "clear") == 0) {
        tvm_clear();
        printf("TVM cleared.\n");
        return;
    }
    
    if (strcmp(p, "beg") == 0 || strcmp(p, "begin") == 0) {
        tvm_set_beg(1);
        printf("Begin mode.\n");
        return;
    }
    
    if (strcmp(p, "end") == 0) {
        tvm_set_beg(0);
        printf("End mode.\n");
        return;
    }
    
    /* Parse "var value" or "var" (to calculate) */
    if (p[0] == 'n' && (p[1] == ' ' || p[1] == '\0')) {
        if (p[1] == '\0') {
            tvm_calc_n(&result);
            tvm_set_n(&result);
            apf_to_str(buf, sizeof(buf), &result, digits);
            printf("N = %s\n", buf);
        } else {
            apf_from_str(&val, p + 2);
            tvm_set_n(&val);
        }
        return;
    }
    
    if ((p[0] == 'i' || p[0] == 'I') && (p[1] == ' ' || p[1] == '\0' || p[1] == '%')) {
        const char *vp = p + 1;
        if (*vp == '%') vp++;
        while (*vp == ' ') vp++;
        if (*vp == '\0') {
            apf hundred;
            tvm_calc_i(&result);
            apf_from_int(&hundred, 100);
            apf_div(&val, &result, &hundred);
            apf_copy(&tvm_regs.i, &val);
            apf_to_str(buf, sizeof(buf), &result, digits);
            printf("I%% = %s%%\n", buf);
        } else {
            apf_from_str(&val, vp);
            tvm_set_i(&val);
        }
        return;
    }
    
    if ((p[0] == 'p' || p[0] == 'P') && (p[1] == 'v' || p[1] == 'V')) {
        const char *vp = p + 2;
        while (*vp == ' ') vp++;
        if (*vp == '\0') {
            tvm_calc_pv(&result);
            tvm_set_pv(&result);
            apf_to_str(buf, sizeof(buf), &result, digits);
            printf("PV = %s\n", buf);
        } else {
            apf_from_str(&val, vp);
            tvm_set_pv(&val);
        }
        return;
    }
    
    if ((p[0] == 'p' || p[0] == 'P') && (p[1] == 'm' || p[1] == 'M') && (p[2] == 't' || p[2] == 'T')) {
        const char *vp = p + 3;
        while (*vp == ' ') vp++;
        if (*vp == '\0') {
            tvm_calc_pmt(&result);
            tvm_set_pmt(&result);
            apf_to_str(buf, sizeof(buf), &result, digits);
            printf("PMT = %s\n", buf);
        } else {
            apf_from_str(&val, vp);
            tvm_set_pmt(&val);
        }
        return;
    }
    
    if ((p[0] == 'f' || p[0] == 'F') && (p[1] == 'v' || p[1] == 'V')) {
        const char *vp = p + 2;
        while (*vp == ' ') vp++;
        if (*vp == '\0') {
            tvm_calc_fv(&result);
            tvm_set_fv(&result);
            apf_to_str(buf, sizeof(buf), &result, digits);
            printf("FV = %s\n", buf);
        } else {
            apf_from_str(&val, vp);
            tvm_set_fv(&val);
        }
        return;
    }
    
    printf("Usage: tvm [n|i|pv|pmt|fv] [value]\n");
    printf("  tvm n 360      - set N to 360\n");
    printf("  tvm i 5        - set I%% to 5%%\n");
    printf("  tvm pv 100000  - set PV to 100000\n");
    printf("  tvm pmt        - calculate PMT\n");
    printf("  tvm fv 0       - set FV to 0\n");
    printf("  tvm beg/end    - set payment timing\n");
    printf("  tvm clr        - clear all\n");
}

#endif /* HAVE_TVM */
