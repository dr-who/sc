/*
 * scalar_funcs.c - Table-driven scalar function dispatch
 * 
 * C89 compliant for Watcom C / DOS
 */
#include <stdio.h>
#include <string.h>
#include <limits.h>
#include "scalar_funcs.h"
#include "apfx.h"
#include "apf_native.h"
#include "mathx.h"
#include "stats.h"

/* ==========================================================================
 * WRAPPER FUNCTIONS for functions not in apfx.h
 * ========================================================================== */

/* log10(x) = ln(x) / ln(10) */
static void sf_log10(apf *r, const apf *x) {
    apf ln10, ten;
    apfx_log(r, x);
    apf_from_int(&ten, 10);
    apfx_log(&ln10, &ten);
    apf_div(r, r, &ln10);
}

/* log2(x) = ln(x) / ln(2) */
static void sf_log2(apf *r, const apf *x) {
    apf ln2, two;
    apfx_log(r, x);
    apf_from_int(&two, 2);
    apfx_log(&ln2, &two);
    apf_div(r, r, &ln2);
}

/* exp2(x) = 2^x */
static void sf_exp2(apf *r, const apf *x) {
    apf two;
    apf_from_int(&two, 2);
    apfx_pow(r, &two, x);
}

/* exp10(x) = 10^x */
static void sf_exp10(apf *r, const apf *x) {
    apf ten;
    apf_from_int(&ten, 10);
    apfx_pow(r, &ten, x);
}

/* cbrt(x) = x^(1/3), sign preserving */
static void sf_cbrt(apf *r, const apf *x) {
    apf third, abs_x;
    int neg;
    apf_from_str(&third, "0.333333333333333333333333333333333");
    neg = apf_cmp_int(x, 0) < 0;
    if (neg) {
        apf_neg(&abs_x, x);
        apfx_pow(r, &abs_x, &third);
        apf_neg(r, r);
    } else {
        apfx_pow(r, x, &third);
    }
}

/* log1p(x) = ln(1+x), accurate for small x */
static void sf_log1p(apf *r, const apf *x) {
    apf one_plus_x, one;
    apf_from_int(&one, 1);
    apf_add(&one_plus_x, &one, x);
    apfx_log(r, &one_plus_x);
}

/* expm1(x) = e^x - 1, accurate for small x */
static void sf_expm1(apf *r, const apf *x) {
    apf one;
    apfx_exp(r, x);
    apf_from_int(&one, 1);
    apf_sub(r, r, &one);
}

/* deg2rad(x) = x * pi / 180 */
static void sf_deg2rad(apf *r, const apf *x) {
    apf pi, d180;
    apfx_pi(&pi);
    apf_from_int(&d180, 180);
    apf_mul(r, x, &pi);
    apf_div(r, r, &d180);
}

/* rad2deg(x) = x * 180 / pi */
static void sf_rad2deg(apf *r, const apf *x) {
    apf pi, d180;
    apfx_pi(&pi);
    apf_from_int(&d180, 180);
    apf_mul(r, x, &d180);
    apf_div(r, r, &pi);
}

/* sign(x) = -1, 0, or 1 */
static void sf_sign(apf *r, const apf *x) {
    int cmp = apf_cmp_int(x, 0);
    if (cmp < 0) apf_from_int(r, -1);
    else if (cmp > 0) apf_from_int(r, 1);
    else apf_from_int(r, 0);
}

/* factorial(n) = n! */
static void sf_fact(apf *r, const apf *x) {
    long n = apf_to_long(x);
    if (n < 0) n = 0;  /* Handle negative by returning 0 */
    if (n > 100000) n = 100000;  /* Cap at max */
    apfx_fact(r, n);
}

/* round(x) = nearest integer */
static void sf_round(apf *r, const apf *x) {
    apf half;
    apf_from_str(&half, "0.5");
    if (apf_cmp_int(x, 0) >= 0) {
        apf_add(r, x, &half);
    } else {
        apf_sub(r, x, &half);
    }
    apf_trunc(r, r);
}

/* sech(x) = 1/cosh(x) */
static void sf_sech(apf *r, const apf *x) {
    apf cosh_x, one;
    apfx_cosh(&cosh_x, x);
    apf_from_int(&one, 1);
    apf_div(r, &one, &cosh_x);
}

/* csch(x) = 1/sinh(x) */
static void sf_csch(apf *r, const apf *x) {
    apf sinh_x, one;
    apfx_sinh(&sinh_x, x);
    apf_from_int(&one, 1);
    apf_div(r, &one, &sinh_x);
}

/* coth(x) = cosh(x)/sinh(x) */
static void sf_coth(apf *r, const apf *x) {
    apf sinh_x, cosh_x;
    apfx_sinh(&sinh_x, x);
    apfx_cosh(&cosh_x, x);
    apf_div(r, &cosh_x, &sinh_x);
}

/* asech(x) = acosh(1/x) */
static void sf_asech(apf *r, const apf *x) {
    apf inv, one;
    apf_from_int(&one, 1);
    apf_div(&inv, &one, x);
    apfx_acosh(r, &inv);
}

/* acsch(x) = asinh(1/x) */
static void sf_acsch(apf *r, const apf *x) {
    apf inv, one;
    apf_from_int(&one, 1);
    apf_div(&inv, &one, x);
    apfx_asinh(r, &inv);
}

/* acoth(x) = atanh(1/x) */
static void sf_acoth(apf *r, const apf *x) {
    apf inv, one;
    apf_from_int(&one, 1);
    apf_div(&inv, &one, x);
    apfx_atanh(r, &inv);
}

/* wrapToPi: wrap angle to [-pi, pi] */
static void sf_wrapToPi(apf *r, const apf *x) {
    apf two_pi, pi_val;
    apfx_pi(&pi_val);
    apf_mul(&two_pi, &pi_val, &pi_val);
    apf_from_int(&two_pi, 2);
    apf_mul(&two_pi, &two_pi, &pi_val);
    apf_mod(r, x, &two_pi);
    /* Adjust to [-pi, pi] */
    if (apf_cmp(r, &pi_val) > 0) {
        apf_sub(r, r, &two_pi);
    } else {
        apf neg_pi;
        apf_neg(&neg_pi, &pi_val);
        if (apf_cmp(r, &neg_pi) < 0) {
            apf_add(r, r, &two_pi);
        }
    }
}

/* wrapTo180: wrap angle in degrees to [-180, 180] */
static void sf_wrapTo180(apf *r, const apf *x) {
    apf d360, d180, neg180;
    apf_from_int(&d360, 360);
    apf_from_int(&d180, 180);
    apf_neg(&neg180, &d180);
    apf_mod(r, x, &d360);
    if (apf_cmp(r, &d180) > 0) {
        apf_sub(r, r, &d360);
    } else if (apf_cmp(r, &neg180) < 0) {
        apf_add(r, r, &d360);
    }
}

/* wrapTo360: wrap angle in degrees to [0, 360) */
static void sf_wrapTo360(apf *r, const apf *x) {
    apf d360, zero;
    apf_from_int(&d360, 360);
    apf_zero(&zero);
    apf_mod(r, x, &d360);
    if (apf_cmp(r, &zero) < 0) {
        apf_add(r, r, &d360);
    }
}

/* nextpow2: next power of 2 >= |x| */
static void sf_nextpow2(apf *r, const apf *x) {
    apf abs_x, two, log2x;
    long n;
    apf_abs(&abs_x, x);
    if (apf_is_zero(&abs_x)) {
        apf_from_int(r, 0);
        return;
    }
    apf_from_int(&two, 2);
    /* log2(abs_x) */
    apfx_log(&log2x, &abs_x);
    {
        apf ln2;
        apfx_log(&ln2, &two);
        apf_div(&log2x, &log2x, &ln2);
    }
    apf_ceil(&log2x, &log2x);
    n = apf_to_long(&log2x);
    apf_from_int(r, n);
}

/* realsqrt - returns NaN for negative, never complex */
static void sf_realsqrt(apf *r, const apf *x) {
    if (apf_is_neg(x)) {
        apf_set_nan(r);
    } else {
        apf_sqrt(r, x);
    }
}

/* reallog - returns NaN for non-positive, never complex */
static void sf_reallog(apf *r, const apf *x) {
    if (apf_is_neg(x) || apf_is_zero(x)) {
        apf_set_nan(r);
    } else {
        apfx_log(r, x);
    }
}

/* sinc(x) = sin(pi*x)/(pi*x), with sinc(0) = 1 */
static void sf_sinc(apf *r, const apf *x) {
    if (apf_is_zero(x)) {
        apf_from_int(r, 1);
    } else {
        apf pi_x, sin_val, pi_val;
        apfx_pi(&pi_val);
        apf_mul(&pi_x, &pi_val, x);
        apfx_sin(&sin_val, &pi_x);
        apf_div(r, &sin_val, &pi_x);
    }
}

/* sigmoid(x) = 1 / (1 + exp(-x)) */
static void sf_sigmoid(apf *r, const apf *x) {
    apf neg_x, exp_val, one, denom;
    apf_neg(&neg_x, x);
    apfx_exp(&exp_val, &neg_x);
    apf_from_int(&one, 1);
    apf_add(&denom, &one, &exp_val);
    apf_div(r, &one, &denom);
}

/* softplus(x) = ln(1 + exp(x)) */
static void sf_softplus(apf *r, const apf *x) {
    apf exp_val, one, sum_val;
    apfx_exp(&exp_val, x);
    apf_from_int(&one, 1);
    apf_add(&sum_val, &one, &exp_val);
    apfx_log(r, &sum_val);
}

/* Heaviside step: 0 if x<0, 0.5 if x=0, 1 if x>0 */
static void sf_heaviside(apf *r, const apf *x) {
    if (apf_is_zero(x)) {
        apf_from_str(r, "0.5");
    } else if (apf_is_neg(x)) {
        apf_zero(r);
    } else {
        apf_from_int(r, 1);
    }
}

/* rect(x) = 1 if |x| < 0.5, 0.5 if |x| = 0.5, 0 otherwise */
static void sf_rect(apf *r, const apf *x) {
    apf abs_x, half;
    int cmp;
    apf_abs(&abs_x, x);
    apf_from_str(&half, "0.5");
    cmp = apf_cmp(&abs_x, &half);
    if (cmp < 0) {
        apf_from_int(r, 1);
    } else if (cmp == 0) {
        apf_from_str(r, "0.5");
    } else {
        apf_zero(r);
    }
}

/* tri(x) = max(0, 1 - |x|) - triangle function */
static void sf_tri(apf *r, const apf *x) {
    apf abs_x, one, diff;
    apf_abs(&abs_x, x);
    apf_from_int(&one, 1);
    apf_sub(&diff, &one, &abs_x);
    if (apf_is_neg(&diff)) {
        apf_zero(r);
    } else {
        apf_copy(r, &diff);
    }
}

/* hypot(x, y) = sqrt(x^2 + y^2) */
static void sf_hypot(apf *r, const apf *x, const apf *y) {
    apf x2, y2, sum;
    apf_mul(&x2, x, x);
    apf_mul(&y2, y, y);
    apf_add(&sum, &x2, &y2);
    apf_sqrt(r, &sum);
}

/* beta(a, b) = gamma(a)*gamma(b)/gamma(a+b) */
static void sf_beta(apf *r, const apf *a, const apf *b) {
    apf ga, gb, gab, ab, num;
    apfx_tgamma(&ga, a);
    apfx_tgamma(&gb, b);
    apf_add(&ab, a, b);
    apfx_tgamma(&gab, &ab);
    apf_mul(&num, &ga, &gb);
    apf_div(r, &num, &gab);
}

/* normpdf(x) = (1/sqrt(2*pi)) * exp(-x^2/2) */
static void sf_normpdf(apf *r, const apf *x) {
    apf x2, neg_half_x2, exp_val, two_pi, sqrt_2pi, coeff;
    apf two;
    
    /* x^2 */
    apf_mul(&x2, x, x);
    
    /* -x^2/2 */
    apf_from_int(&two, 2);
    apf_div(&neg_half_x2, &x2, &two);
    apf_neg(&neg_half_x2, &neg_half_x2);
    
    /* exp(-x^2/2) */
    apfx_exp(&exp_val, &neg_half_x2);
    
    /* 1/sqrt(2*pi) */
    apfx_pi(&two_pi);
    apf_mul(&two_pi, &two_pi, &two);
    apf_sqrt(&sqrt_2pi, &two_pi);
    apf_from_int(&coeff, 1);
    apf_div(&coeff, &coeff, &sqrt_2pi);
    
    apf_mul(r, &coeff, &exp_val);
}

/* ============== PROBABILITY DISTRIBUTION WRAPPERS ============== */
/* These wrap the double-based stat_ functions from stats.c */

#ifdef HAVE_STATS
/* t-distribution PDF */
static void sf_tpdf(apf *r, const apf *x, const apf *df) {
    double xd = apf_to_double(x);
    double dfd = apf_to_double(df);
    apf_from_double(r, stat_tpdf(xd, dfd));
}

/* t-distribution CDF */
static void sf_tcdf(apf *r, const apf *x, const apf *df) {
    double xd = apf_to_double(x);
    double dfd = apf_to_double(df);
    apf_from_double(r, stat_tcdf(xd, dfd));
}

/* t-distribution inverse CDF */
static void sf_tinv(apf *r, const apf *p, const apf *df) {
    double pd = apf_to_double(p);
    double dfd = apf_to_double(df);
    apf_from_double(r, stat_tinv(pd, dfd));
}

/* Chi-square PDF */
static void sf_chi2pdf(apf *r, const apf *x, const apf *df) {
    double xd = apf_to_double(x);
    double dfd = apf_to_double(df);
    apf_from_double(r, stat_chi2pdf(xd, dfd));
}

/* Chi-square CDF */
static void sf_chi2cdf(apf *r, const apf *x, const apf *df) {
    double xd = apf_to_double(x);
    double dfd = apf_to_double(df);
    apf_from_double(r, stat_chi2cdf(xd, dfd));
}

/* Chi-square inverse CDF */
static void sf_chi2inv(apf *r, const apf *p, const apf *df) {
    double pd = apf_to_double(p);
    double dfd = apf_to_double(df);
    apf_from_double(r, stat_chi2inv(pd, dfd));
}

/* F-distribution and Binomial are handled in parser for 3+ args */

/* Poisson PDF */
static void sf_poisspdf(apf *r, const apf *k, const apf *lambda) {
    int ki = (int)apf_to_long(k);
    double ld = apf_to_double(lambda);
    apf_from_double(r, stat_poisspdf(ki, ld));
}

/* Poisson CDF */
static void sf_poisscdf(apf *r, const apf *k, const apf *lambda) {
    int ki = (int)apf_to_long(k);
    double ld = apf_to_double(lambda);
    apf_from_double(r, stat_poisscdf(ki, ld));
}

/* Exponential PDF */
static void sf_exppdf(apf *r, const apf *x, const apf *lambda) {
    double xd = apf_to_double(x);
    double ld = apf_to_double(lambda);
    apf_from_double(r, stat_exppdf(xd, ld));
}

/* Exponential CDF */
static void sf_expcdf(apf *r, const apf *x, const apf *lambda) {
    double xd = apf_to_double(x);
    double ld = apf_to_double(lambda);
    apf_from_double(r, stat_expcdf(xd, ld));
}

/* Exponential inverse CDF */
static void sf_expinv(apf *r, const apf *p, const apf *lambda) {
    double pd = apf_to_double(p);
    double ld = apf_to_double(lambda);
    apf_from_double(r, stat_expinv(pd, ld));
}

/* Standard error of mean (need to wrap for vectors - placeholder) */
static void sf_sem(apf *r, const apf *std, const apf *n) {
    apf sqn;
    apf_sqrt(&sqn, n);
    apf_div(r, std, &sqn);
}
#endif /* HAVE_STATS */

/* nthroot(x, n) - real nth root preserving sign for odd n */
static void sf_nthroot(apf *r, const apf *x, const apf *n) {
    apf one, inv_n, abs_x;
    long n_long = apf_to_long(n);
    
    apf_from_int(&one, 1);
    apf_div(&inv_n, &one, n);
    
    if (apf_is_neg(x) && (n_long % 2) != 0) {
        /* Odd root of negative: -|x|^(1/n) */
        apf_abs(&abs_x, x);
        apfx_pow(r, &abs_x, &inv_n);
        apf_neg(r, r);
    } else {
        apfx_pow(r, x, &inv_n);
    }
}

/* rem(a, b) = a - trunc(a/b) * b */
static void sf_rem(apf *r, const apf *a, const apf *b) {
    apf quot, truncated;
    apf_div(&quot, a, b);
    apf_trunc(&truncated, &quot);
    apf_mul(&truncated, &truncated, b);
    apf_sub(r, a, &truncated);
}

/* copysign(x, y) - magnitude of x with sign of y */
static void sf_copysign(apf *r, const apf *x, const apf *y) {
    apf_abs(r, x);
    if (apf_is_neg(y)) {
        apf_neg(r, r);
    }
}

/* fdim(x, y) = max(x-y, 0) */
static void sf_fdim(apf *r, const apf *x, const apf *y) {
    apf diff;
    apf_sub(&diff, x, y);
    if (apf_is_neg(&diff)) {
        apf_zero(r);
    } else {
        apf_copy(r, &diff);
    }
}

/* fmod - same as rem */
static void sf_fmod(apf *r, const apf *a, const apf *b) {
    sf_rem(r, a, b);
}

/* isequal(a, b) - test equality, returns 1 or 0 */
static void sf_isequal(apf *r, const apf *a, const apf *b) {
    apf_from_int(r, apf_cmp(a, b) == 0 ? 1 : 0);
}

/* approxeq(a, b) - approximate equality within tolerance */
static void sf_approxeq(apf *r, const apf *a, const apf *b) {
    apf diff, tol;
    apf_sub(&diff, a, b);
    apf_abs(&diff, &diff);
    apf_from_str(&tol, "1e-10");
    apf_from_int(r, apf_cmp(&diff, &tol) <= 0 ? 1 : 0);
}

/* min2(a, b) - minimum of two values */
static void sf_min2(apf *r, const apf *a, const apf *b) {
    apf_copy(r, apf_cmp(a, b) <= 0 ? a : b);
}

/* ============================================================
 * Integer function wrappers (long -> long converted to apf)
 * ============================================================ */

/* digit sum wrapper */
static void sf_digsum(apf *r, long n) {
    apf_from_int(r, digsum_long(n < 0 ? -n : n));
}

/* number of digits wrapper */
static void sf_numdigits(apf *r, long n) {
    apf_from_int(r, numdigits_long(n < 0 ? -n : n));
}

/* digital root wrapper */
static void sf_digitalroot(apf *r, long n) {
    apf_from_int(r, digitalroot_long(n < 0 ? -n : n));
}

/* reverse digits wrapper */
static void sf_reversedigits(apf *r, long n) {
    apf_from_int(r, reverse_long(n < 0 ? -n : n));
}

/* ispalindrome wrapper */
static void sf_ispalindrome(apf *r, long n) {
    apf_from_int(r, ispalindrome_long(n < 0 ? -n : n));
}

/* isperfect wrapper */
static void sf_isperfect(apf *r, long n) {
    apf_from_int(r, isperfect_long(n));
}

/* issquarefree wrapper */
static void sf_issquarefree(apf *r, long n) {
    apf_from_int(r, issquarefree_long(n < 0 ? -n : n));
}

/* totient/eulerphi wrapper */
static void sf_totient(apf *r, long n) {
    long res, temp_n, p;
    if (n < 1) n = 1;
    res = n;
    temp_n = n;
    for (p = 2; p * p <= temp_n; p++) {
        if (temp_n % p == 0) {
            while (temp_n % p == 0) temp_n /= p;
            res -= res / p;
        }
    }
    if (temp_n > 1) res -= res / temp_n;
    apf_from_int(r, res);
}

/* divisorsum wrapper */
static void sf_divisorsum(apf *r, long n) {
    long sum = 0, i;
    if (n < 1) n = 1;
    for (i = 1; i * i <= n; i++) {
        if (n % i == 0) {
            sum += i;
            if (i != n / i) sum += n / i;
        }
    }
    apf_from_int(r, sum);
}

/* derangements/subfactorial wrapper */
static void sf_derangements(apf *r, long n) {
    /* D(n) = (n-1) * (D(n-1) + D(n-2)), D(0)=1, D(1)=0 */
    long d0 = 1, d1 = 0, d2, i;
    if (n <= 0) { apf_from_int(r, 1); return; }
    if (n == 1) { apf_from_int(r, 0); return; }
    for (i = 2; i <= n; i++) {
        d2 = (i - 1) * (d0 + d1);
        d0 = d1;
        d1 = d2;
    }
    apf_from_int(r, d1);
}

/* even - returns 1 if even, 0 otherwise */
static void sf_even(apf *r, long n) {
    apf_from_int(r, (n % 2 == 0) ? 1 : 0);
}

/* odd - returns 1 if odd, 0 otherwise */
static void sf_odd(apf *r, long n) {
    apf_from_int(r, (n % 2 != 0) ? 1 : 0);
}

/* hexagonal number: H(n) = n(2n-1) */
static void sf_hexagonal(apf *r, long n) {
    apf_from_int(r, n * (2 * n - 1));
}

/* triangular number: T(n) = n(n+1)/2 */
static void sf_triangular(apf *r, long n) {
    apf_from_int(r, n * (n + 1) / 2);
}

/* pentagonal number: P(n) = n(3n-1)/2 */
static void sf_pentagonal(apf *r, long n) {
    apf_from_int(r, n * (3 * n - 1) / 2);
}

/* square number: S(n) = n^2 */
static void sf_squarenum(apf *r, long n) {
    apf_from_int(r, n * n);
}

/* cube number: C(n) = n^3 */
static void sf_cubenum(apf *r, long n) {
    apf_from_int(r, n * n * n);
}

/* generalized harmonic: H_{n,m} = sum 1/k^m */
static void sf_genharmonic(apf *r, long n, long m) {
    apfx_harmonic_gen(r, n, m);
}

/* divisorsum with power: sigma_k(n) = sum d^k for d|n */
static void sf_sigmak(apf *r, long n, long k) {
    long sum = 0, i;
    if (n < 1) n = 1;
    for (i = 1; i * i <= n; i++) {
        if (n % i == 0) {
            long pk = 1, j;
            for (j = 0; j < k; j++) pk *= i;
            sum += pk;
            if (i != n / i) {
                long qi = n / i;
                pk = 1;
                for (j = 0; j < k; j++) pk *= qi;
                sum += pk;
            }
        }
    }
    apf_from_int(r, sum);
}

/* frac - fractional part */
static void sf_frac(apf *r, const apf *x) {
    apf truncated;
    apf_trunc(&truncated, x);
    apf_sub(r, x, &truncated);
}

/* isnan - test if NaN */
static void sf_isnan(apf *r, const apf *x) {
    apf_from_int(r, (x->cls == APF_CLASS_NAN) ? 1 : 0);
}

/* isinf - test if infinite */
static void sf_isinf(apf *r, const apf *x) {
    apf_from_int(r, (x->cls == APF_CLASS_INF) ? 1 : 0);
}

/* isfinite - test if normal or zero */
static void sf_isfinite(apf *r, const apf *x) {
    apf_from_int(r, (x->cls == APF_CLASS_NORMAL || x->cls == APF_CLASS_ZERO) ? 1 : 0);
}

/* isnumeric - always true */
static void sf_isnumeric(apf *r, const apf *x) {
    (void)x;
    apf_from_int(r, 1);
}

/* isinteger - test if integer */
static void sf_isinteger(apf *r, const apf *x) {
    int is_int = 0;
    if (x->cls == APF_CLASS_NORMAL) {
        apf truncated, diff;
        apf_trunc(&truncated, x);
        apf_sub(&diff, x, &truncated);
        is_int = apf_is_zero(&diff);
    } else if (x->cls == APF_CLASS_ZERO) {
        is_int = 1;
    }
    apf_from_int(r, is_int);
}

/* ispositive - test if positive */
static void sf_ispositive(apf *r, const apf *x) {
    apf_from_int(r, (!apf_is_neg(x) && !apf_is_zero(x)) ? 1 : 0);
}

/* isnegative - test if negative */
static void sf_isnegative(apf *r, const apf *x) {
    apf_from_int(r, apf_is_neg(x) ? 1 : 0);
}

/* lsl/shl/bitshift - left shift wrapper */
static void sf_lsl(apf *r, const apf *x, long n) {
    apf_lsl(r, x, (int)n);
}

/* lsr/shr - right shift wrapper */
static void sf_lsr(apf *r, const apf *x, long n) {
    apf_lsr(r, x, (int)n);
}

/* sinpi(x) = sin(pi*x) with exact values at integers */
static void sf_sinpi(apf *r, const apf *x) {
    apf pi_x, trunc_x, frac_x;
    apf_trunc(&trunc_x, x);
    apf_sub(&frac_x, x, &trunc_x);
    if (apf_is_zero(&frac_x)) {
        apf_zero(r);
    } else {
        apfx_pi(&pi_x);
        apf_mul(&pi_x, &pi_x, x);
        apfx_sin(r, &pi_x);
    }
}

/* cospi(x) = cos(pi*x) */
static void sf_cospi(apf *r, const apf *x) {
    apf pi_x;
    apfx_pi(&pi_x);
    apf_mul(&pi_x, &pi_x, x);
    apfx_cos(r, &pi_x);
}

/* max2(a, b) - maximum of two values */
static void sf_max2(apf *r, const apf *a, const apf *b) {
    apf_copy(r, apf_cmp(a, b) >= 0 ? a : b);
}

/* Duration conversions to seconds */
static void sf_seconds(apf *r, const apf *x) { apf_copy(r, x); }
static void sf_minutes(apf *r, const apf *x) { 
    apf c; apf_from_int(&c, 60); apf_mul(r, x, &c); 
}
static void sf_hours(apf *r, const apf *x) { 
    apf c; apf_from_int(&c, 3600); apf_mul(r, x, &c); 
}
static void sf_days(apf *r, const apf *x) { 
    apf c; apf_from_int(&c, 86400); apf_mul(r, x, &c); 
}
static void sf_milliseconds(apf *r, const apf *x) { 
    apf c; apf_from_int(&c, 1000); apf_div(r, x, &c); 
}
static void sf_weeks(apf *r, const apf *x) { 
    apf c; apf_from_int(&c, 604800); apf_mul(r, x, &c); 
}

/* fibonacci(n) - Fibonacci number F(n) */
static void sf_fibonacci(apf *r, const apf *x) {
    long n, i;
    apf a, b, tmp;
    n = apf_to_long(x);
    if (n < 0) n = 0;
    if (n > 10000) n = 10000;
    if (n == 0) {
        apf_from_int(r, 0);
        return;
    }
    apf_from_int(&a, 0);
    apf_from_int(&b, 1);
    for (i = 2; i <= n; i++) {
        apf_add(&tmp, &a, &b);
        a = b;
        b = tmp;
    }
    apf_copy(r, &b);
}

/* lucas(n) - Lucas number L(n): L(0)=2, L(1)=1, L(n)=L(n-1)+L(n-2) */
static void sf_lucas(apf *r, const apf *x) {
    long n, i;
    apf a, b, tmp;
    n = apf_to_long(x);
    if (n < 0) n = 0;
    if (n > 10000) n = 10000;
    if (n == 0) {
        apf_from_int(r, 2);
        return;
    }
    if (n == 1) {
        apf_from_int(r, 1);
        return;
    }
    apf_from_int(&a, 2);
    apf_from_int(&b, 1);
    for (i = 2; i <= n; i++) {
        apf_add(&tmp, &a, &b);
        a = b;
        b = tmp;
    }
    apf_copy(r, &b);
}

/* catalan(n) - Catalan number C(n) = (2n)! / ((n+1)! * n!) */
static void sf_catalan(apf *r, const apf *x) {
    long n, i;
    apf den, two_n_fact, n_fact, n1_fact, tmp;
    n = apf_to_long(x);
    if (n < 0) n = 0;
    if (n > 500) n = 500;
    
    /* (2n)! */
    apf_from_int(&two_n_fact, 1);
    for (i = 2; i <= 2*n; i++) {
        apf_from_int(&tmp, i);
        apf_mul(&two_n_fact, &two_n_fact, &tmp);
    }
    /* n! */
    apf_from_int(&n_fact, 1);
    for (i = 2; i <= n; i++) {
        apf_from_int(&tmp, i);
        apf_mul(&n_fact, &n_fact, &tmp);
    }
    /* (n+1)! */
    apf_from_int(&tmp, n + 1);
    apf_mul(&n1_fact, &n_fact, &tmp);
    
    /* C(n) = (2n)! / ((n+1)! * n!) */
    apf_mul(&den, &n1_fact, &n_fact);
    apf_div(r, &two_n_fact, &den);
}

/* factorial2(n) - double factorial n!! */
static void sf_factorial2(apf *r, const apf *x) {
    long n, i;
    apf tmp;
    n = apf_to_long(x);
    if (n < 0) n = 0;
    if (n > 500) n = 500;
    apf_from_int(r, 1);
    for (i = n; i > 1; i -= 2) {
        apf_from_int(&tmp, i);
        apf_mul(r, r, &tmp);
    }
}

/* digamma(x) = psi(x) = d/dx ln(Gamma(x)) */
static void sf_digamma(apf *r, const apf *x) {
    apfx_digamma(r, x);
}

/* zeta(s) - Riemann zeta function */
static void sf_zeta(apf *r, const apf *x) {
    apfx_zeta(r, x);
}

/* harmonic(n) - nth harmonic number H_n = 1 + 1/2 + ... + 1/n */
static void sf_harmonic(apf *r, const apf *x) {
    long n = apf_to_long(x);
    if (n < 0) n = 0;
    apfx_harmonic(r, n);
}

/* ==========================================================================
 * HASH TABLE
 * ========================================================================== */

/* Hash table for O(1) lookup - sized for ~2000 entries */
#define SF_HASH_SIZE 509
static ScalarFuncEntry *sf_hash[SF_HASH_SIZE];

/* Static storage for entries - room for 8x growth */
#define SF_MAX_ENTRIES 2048
static ScalarFuncEntry sf_entries[SF_MAX_ENTRIES];
static int sf_count = 0;

/* External: angle mode from runtime */
extern int angle_mode;
#define ANGLE_RAD 0
#define ANGLE_DEG 1
#define ANGLE_GRAD 2

/* Helper: convert angle to radians */
static void angle_to_rad_local(apf *rad, const apf *angle, int mode) {
    if (mode == ANGLE_DEG) {
        apf pi_val, d180;
        apfx_pi(&pi_val);
        apf_from_int(&d180, 180);
        apf_mul(rad, angle, &pi_val);
        apf_div(rad, rad, &d180);
    } else if (mode == ANGLE_GRAD) {
        apf pi_val, d200;
        apfx_pi(&pi_val);
        apf_from_int(&d200, 200);
        apf_mul(rad, angle, &pi_val);
        apf_div(rad, rad, &d200);
    } else {
        apf_copy(rad, angle);
    }
}

/* Helper: convert radians to angle */
static void rad_to_angle_local(apf *angle, const apf *rad, int mode) {
    if (mode == ANGLE_DEG) {
        apf pi_val, d180;
        apfx_pi(&pi_val);
        apf_from_int(&d180, 180);
        apf_mul(angle, rad, &d180);
        apf_div(angle, angle, &pi_val);
    } else if (mode == ANGLE_GRAD) {
        apf pi_val, d200;
        apfx_pi(&pi_val);
        apf_from_int(&d200, 200);
        apf_mul(angle, rad, &d200);
        apf_div(angle, angle, &pi_val);
    } else {
        apf_copy(angle, rad);
    }
}

/* Simple string hash */
static unsigned sf_hash_name(const char *name) {
    unsigned h = 0;
    while (*name) {
        char c = *name++;
        if (c >= 'A' && c <= 'Z') c += 32;
        h = h * 31 + (unsigned char)c;
    }
    return h % SF_HASH_SIZE;
}

/* Case-insensitive compare */
static int sf_streq(const char *a, const char *b) {
    /* Case-sensitive string comparison */
    while (*a && *b) {
        if (*a++ != *b++) return 0;
    }
    return *a == *b;
}

/* Register a real->real function */
static void reg_real1(const char *name, apf_fn_1 fn) {
    unsigned h;
    ScalarFuncEntry *e;
    if (sf_count >= SF_MAX_ENTRIES) return;
    e = &sf_entries[sf_count++];
    e->name = name;
    e->type = SF_REAL_1;
    e->fn.f1 = fn;
    h = sf_hash_name(name);
    e->next = sf_hash[h];
    sf_hash[h] = e;
}

/* Register integer sequence function: f(r, n) where n is long */
static void reg_int1(const char *name, apf_int_fn fn, long minv, long maxv) {
    unsigned h;
    ScalarFuncEntry *e;
    if (sf_count >= SF_MAX_ENTRIES) return;
    e = &sf_entries[sf_count++];
    e->name = name;
    e->type = SF_INT_1;
    e->fn.fi = fn;
    e->min_val = minv;
    e->max_val = maxv;
    h = sf_hash_name(name);
    e->next = sf_hash[h];
    sf_hash[h] = e;
}

/* Register 2-arg integer function: f(r, n, k) where n,k are long */
static void reg_int2(const char *name, apf_int2_fn fn) {
    unsigned h;
    ScalarFuncEntry *e;
    if (sf_count >= SF_MAX_ENTRIES) return;
    e = &sf_entries[sf_count++];
    e->name = name;
    e->type = SF_INT_2;
    e->fn.fi2 = fn;
    h = sf_hash_name(name);
    e->next = sf_hash[h];
    sf_hash[h] = e;
}

/* Register mixed apf + int function: f(r, x, n) where x is apf, n is long */
static void reg_apf_int(const char *name, apf_apf_int_fn fn) {
    unsigned h;
    ScalarFuncEntry *e;
    if (sf_count >= SF_MAX_ENTRIES) return;
    e = &sf_entries[sf_count++];
    e->name = name;
    e->type = SF_APF_INT;
    e->fn.fai = fn;
    h = sf_hash_name(name);
    e->next = sf_hash[h];
    sf_hash[h] = e;
}

static void reg_int_apf(const char *name, apf_int_apf_fn fn) {
    unsigned h;
    ScalarFuncEntry *e;
    if (sf_count >= SF_MAX_ENTRIES) return;
    e = &sf_entries[sf_count++];
    e->name = name;
    e->type = SF_INT_APF;
    e->fn.fia = fn;
    h = sf_hash_name(name);
    e->next = sf_hash[h];
    sf_hash[h] = e;
}

/* Register the complex(re, im) function */
static void reg_make_complex(const char *name) {
    unsigned h;
    ScalarFuncEntry *e;
    if (sf_count >= SF_MAX_ENTRIES) return;
    e = &sf_entries[sf_count++];
    e->name = name;
    e->type = SF_MAKE_COMPLEX;
    h = sf_hash_name(name);
    e->next = sf_hash[h];
    sf_hash[h] = e;
}

/* Register a 2-arg real->real function */
static void reg_real2(const char *name, apf_fn_2 fn) {
    unsigned h;
    ScalarFuncEntry *e;
    if (sf_count >= SF_MAX_ENTRIES) return;
    e = &sf_entries[sf_count++];
    e->name = name;
    e->type = SF_REAL_2;
    e->fn.f2 = fn;
    h = sf_hash_name(name);
    e->next = sf_hash[h];
    sf_hash[h] = e;
}

/* Register a complex->complex function */
static void reg_complex1(const char *name, apfc_fn_1 fn) {
    unsigned h;
    ScalarFuncEntry *e;
    if (sf_count >= SF_MAX_ENTRIES) return;
    e = &sf_entries[sf_count++];
    e->name = name;
    e->type = SF_COMPLEX_1;
    e->fn.c1 = fn;
    h = sf_hash_name(name);
    e->next = sf_hash[h];
    sf_hash[h] = e;
}

/* Register a 2-arg complex->complex function */
static void reg_complex2(const char *name, apfc_fn_2 fn) {
    unsigned h;
    ScalarFuncEntry *e;
    if (sf_count >= SF_MAX_ENTRIES) return;
    e = &sf_entries[sf_count++];
    e->name = name;
    e->type = SF_COMPLEX_2;
    e->fn.c2 = fn;
    h = sf_hash_name(name);
    e->next = sf_hash[h];
    sf_hash[h] = e;
}

/* Register a trig function (uses angle_mode) */
static void reg_trig(const char *name, apfc_fn_1 fn) {
    unsigned h;
    ScalarFuncEntry *e;
    if (sf_count >= SF_MAX_ENTRIES) return;
    e = &sf_entries[sf_count++];
    e->name = name;
    e->type = SF_TRIG;
    e->fn.c1 = fn;
    h = sf_hash_name(name);
    e->next = sf_hash[h];
    sf_hash[h] = e;
}

/* Register an inverse trig function (outputs in angle_mode) */
static void reg_trig_inv(const char *name, apf_fn_1 fn) {
    unsigned h;
    ScalarFuncEntry *e;
    if (sf_count >= SF_MAX_ENTRIES) return;
    e = &sf_entries[sf_count++];
    e->name = name;
    e->type = SF_TRIG_INV;
    e->fn.f1 = fn;
    h = sf_hash_name(name);
    e->next = sf_hash[h];
    sf_hash[h] = e;
}

/* Register a degree trig function (always degrees in) */
static void reg_trig_deg(const char *name, apfc_fn_1 fn) {
    unsigned h;
    ScalarFuncEntry *e;
    if (sf_count >= SF_MAX_ENTRIES) return;
    e = &sf_entries[sf_count++];
    e->name = name;
    e->type = SF_TRIG_DEG;
    e->fn.c1 = fn;
    h = sf_hash_name(name);
    e->next = sf_hash[h];
    sf_hash[h] = e;
}

/* Register a real trig function (uses angle_mode, real only) */
static void reg_trig_real(const char *name, apf_fn_1 fn) {
    unsigned h;
    ScalarFuncEntry *e;
    if (sf_count >= SF_MAX_ENTRIES) return;
    e = &sf_entries[sf_count++];
    e->name = name;
    e->type = SF_TRIG_REAL;
    e->fn.f1 = fn;
    h = sf_hash_name(name);
    e->next = sf_hash[h];
    sf_hash[h] = e;
}

/* Register a degree trig function (real only) */
static void reg_trig_deg_real(const char *name, apf_fn_1 fn) {
    unsigned h;
    ScalarFuncEntry *e;
    if (sf_count >= SF_MAX_ENTRIES) return;
    e = &sf_entries[sf_count++];
    e->name = name;
    e->type = SF_TRIG_DEG_REAL;
    e->fn.f1 = fn;
    h = sf_hash_name(name);
    e->next = sf_hash[h];
    sf_hash[h] = e;
}

/* Register an inverse trig returning degrees */
static void reg_trig_inv_deg(const char *name, apf_fn_1 fn) {
    unsigned h;
    ScalarFuncEntry *e;
    if (sf_count >= SF_MAX_ENTRIES) return;
    e = &sf_entries[sf_count++];
    e->name = name;
    e->type = SF_TRIG_INV_DEG;
    e->fn.f1 = fn;
    h = sf_hash_name(name);
    e->next = sf_hash[h];
    sf_hash[h] = e;
}

/* Register a 2-arg inverse trig (atan2) */
static void reg_trig_inv_2(const char *name, apf_fn_2 fn) {
    unsigned h;
    ScalarFuncEntry *e;
    if (sf_count >= SF_MAX_ENTRIES) return;
    e = &sf_entries[sf_count++];
    e->name = name;
    e->type = SF_TRIG_INV_2;
    e->fn.f2 = fn;
    h = sf_hash_name(name);
    e->next = sf_hash[h];
    sf_hash[h] = e;
}

static void reg_trig_inv_2_deg(const char *name, apf_fn_2 fn) {
    unsigned h;
    ScalarFuncEntry *e;
    if (sf_count >= SF_MAX_ENTRIES) return;
    e = &sf_entries[sf_count++];
    e->name = name;
    e->type = SF_TRIG_INV_2_DEG;
    e->fn.f2 = fn;
    h = sf_hash_name(name);
    e->next = sf_hash[h];
    sf_hash[h] = e;
}

/* Register a special type */
static void reg_special(const char *name, ScalarFuncType type) {
    unsigned h;
    ScalarFuncEntry *e;
    if (sf_count >= SF_MAX_ENTRIES) return;
    e = &sf_entries[sf_count++];
    e->name = name;
    e->type = type;
    e->fn.f1 = NULL;
    h = sf_hash_name(name);
    e->next = sf_hash[h];
    sf_hash[h] = e;
}

/* === NUMBER THEORY WRAPPERS (for table dispatch) === */

/* omega(n) - number of distinct prime factors */
static void sf_omega(apf *r, long n) {
    apf_from_int(r, omega_long(n));
}

/* bigomega(n) - number of prime factors with multiplicity */
static void sf_bigomega(apf *r, long n) {
    apf_from_int(r, bigomega_long(n));
}

/* moebius(n) - Möbius mu function */
static void sf_moebius(apf *r, long n) {
    apf_from_int(r, moebius_long(n));
}

/* sigma(n) - sum of divisors */
static void sf_sigma(apf *r, long n) {
    apf_from_int(r, sigma_long(n, 1));
}

/* sigma0(n) - number of divisors */
static void sf_sigma0(apf *r, long n) {
    apf_from_int(r, sigma_long(n, 0));
}

/* isreal(z) - check if imaginary part is zero */
static void sf_isreal(apfc *r, const apfc *z) {
    apf_from_int(&r->re, apf_is_zero(&z->im) ? 1 : 0);
    apf_zero(&r->im);
}

/* === COMPARISON FUNCTIONS === */
static void sf_eq(apf *r, const apf *a, const apf *b) {
    apf_from_int(r, apf_cmp(a, b) == 0 ? 1 : 0);
}

static void sf_ne(apf *r, const apf *a, const apf *b) {
    apf_from_int(r, apf_cmp(a, b) != 0 ? 1 : 0);
}

static void sf_lt(apf *r, const apf *a, const apf *b) {
    apf_from_int(r, apf_cmp(a, b) < 0 ? 1 : 0);
}

static void sf_le(apf *r, const apf *a, const apf *b) {
    apf_from_int(r, apf_cmp(a, b) <= 0 ? 1 : 0);
}

static void sf_gt(apf *r, const apf *a, const apf *b) {
    apf_from_int(r, apf_cmp(a, b) > 0 ? 1 : 0);
}

static void sf_ge(apf *r, const apf *a, const apf *b) {
    apf_from_int(r, apf_cmp(a, b) >= 0 ? 1 : 0);
}

/* === NUMBER THEORY FUNCTIONS === */
static void sf_isprime(apf *r, const apf *x) {
    long n, d, rr, a, xx, i, j;
    int witnesses[] = {2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37};
    int is_prime = 1;
    n = apf_to_long(x);
    if (n < 2) { is_prime = 0; }
    else if (n == 2 || n == 3) { is_prime = 1; }
    else if (n % 2 == 0) { is_prime = 0; }
    else {
        d = n - 1; rr = 0;
        while (d % 2 == 0) { d /= 2; rr++; }
        for (i = 0; i < 12 && is_prime; i++) {
            a = witnesses[i];
            if (a >= n) continue;
            xx = 1;
            { long base = a % n, exp = d;
              while (exp > 0) {
                  if (exp % 2 == 1) xx = (xx * base) % n;
                  exp /= 2; base = (base * base) % n;
              }
            }
            if (xx == 1 || xx == n - 1) continue;
            for (j = 0; j < rr - 1; j++) {
                xx = (xx * xx) % n;
                if (xx == n - 1) break;
            }
            if (xx != n - 1) is_prime = 0;
        }
    }
    apf_from_int(r, is_prime);
}

static void sf_iseven(apf *r, const apf *x) {
    apf_from_int(r, (apf_to_long(x) % 2 == 0) ? 1 : 0);
}

static void sf_isodd(apf *r, const apf *x) {
    apf_from_int(r, (apf_to_long(x) % 2 != 0) ? 1 : 0);
}

static void sf_nextprime(apf *r, const apf *x) {
    long n = apf_to_long(x);
    apf test;
    if (n <= 2) { apf_from_int(r, 2); return; }
    if (n % 2 == 0) n++;
    while (1) {
        apf_from_int(&test, n);
        sf_isprime(r, &test);
        if (apf_to_long(r) == 1) { apf_from_int(r, n); return; }
        n += 2;
    }
}

static void sf_prevprime(apf *r, const apf *x) {
    long n = apf_to_long(x);
    apf test;
    if (n <= 2) { apf_set_nan(r); return; }
    n--;
    if (n == 2) { apf_from_int(r, 2); return; }
    if (n % 2 == 0) n--;
    while (n > 1) {
        apf_from_int(&test, n);
        sf_isprime(r, &test);
        if (apf_to_long(r) == 1) { apf_from_int(r, n); return; }
        n -= 2;
    }
    apf_set_nan(r);
}

static void sf_modinv(apf *r, const apf *a, const apf *m) {
    long aa = apf_to_long(a), mm = apf_to_long(m);
    long m0 = mm, t, q, x0 = 0, x1 = 1;
    if (mm == 1) { apf_zero(r); return; }
    while (aa > 1) {
        q = aa / mm; t = mm; mm = aa % mm; aa = t;
        t = x0; x0 = x1 - q * x0; x1 = t;
    }
    if (x1 < 0) x1 += m0;
    apf_from_int(r, x1);
}

/* === LOGICAL FUNCTIONS === */
static void sf_nand(apf *r, const apf *a, const apf *b) {
    int a_val = !apf_is_zero(a);
    int b_val = !apf_is_zero(b);
    apf_from_int(r, !(a_val && b_val));
}

static void sf_nor(apf *r, const apf *a, const apf *b) {
    int a_val = !apf_is_zero(a);
    int b_val = !apf_is_zero(b);
    apf_from_int(r, !(a_val || b_val));
}

static void sf_xnor(apf *r, const apf *a, const apf *b) {
    int a_val = !apf_is_zero(a);
    int b_val = !apf_is_zero(b);
    apf_from_int(r, !(a_val ^ b_val));
}

static void sf_implies(apf *r, const apf *a, const apf *b) {
    int a_val = !apf_is_zero(a);
    int b_val = !apf_is_zero(b);
    apf_from_int(r, !a_val || b_val);
}

void scalar_funcs_init(void) {
    int i;
    for (i = 0; i < SF_HASH_SIZE; i++) sf_hash[i] = NULL;
    sf_count = 0;
    
    /* === TRIG (angle_mode aware, complex capable) === */
    reg_trig("sin", apfc_sin);
    reg_trig("cos", apfc_cos);
    reg_trig("tan", apfc_tan);
    
    /* === RECIPROCAL TRIG (angle_mode aware, real only) === */
    reg_trig_real("sec", apfx_sec);
    reg_trig_real("csc", apfx_csc);
    reg_trig_real("cosec", apfx_csc);
    reg_trig_real("cot", apfx_cot);
    reg_trig_real("cotan", apfx_cot);
    
    /* === INVERSE TRIG (output in angle_mode) === */
    reg_trig_inv("asin", apfx_asin);
    reg_trig_inv("acos", apfx_acos);
    reg_trig_inv("atan", apfx_atan);
    reg_trig_inv("arcsin", apfx_asin);
    reg_trig_inv("arccos", apfx_acos);
    reg_trig_inv("arctan", apfx_atan);
    reg_trig_inv("asec", apfx_asec);
    reg_trig_inv("acsc", apfx_acsc);
    reg_trig_inv("acot", apfx_acot);
    reg_trig_inv("arcsec", apfx_asec);
    reg_trig_inv("arccosec", apfx_acsc);
    reg_trig_inv("arccsc", apfx_acsc);
    reg_trig_inv("arccotan", apfx_acot);
    reg_trig_inv("arccot", apfx_acot);
    
    /* === DEGREE TRIG (always degrees in, complex capable) === */
    reg_trig_deg("sind", apfc_sin);
    reg_trig_deg("cosd", apfc_cos);
    reg_trig_deg("tand", apfc_tan);
    
    /* === DEGREE RECIPROCAL TRIG (always degrees in, real only) === */
    reg_trig_deg_real("secd", apfx_sec);
    reg_trig_deg_real("cscd", apfx_csc);
    reg_trig_deg_real("cotd", apfx_cot);
    
    /* === INVERSE TRIG DEGREES (output always degrees) === */
    reg_trig_inv_deg("asind", apfx_asin);
    reg_trig_inv_deg("acosd", apfx_acos);
    reg_trig_inv_deg("atand", apfx_atan);
    reg_trig_inv_deg("asecd", apfx_asec);
    reg_trig_inv_deg("acscd", apfx_acsc);
    reg_trig_inv_deg("acotd", apfx_acot);
    
    /* === 2-ARG INVERSE TRIG === */
    reg_trig_inv_2("atan2", apfx_atan2);
    reg_trig_inv_2_deg("atan2d", apfx_atan2);
    
    /* === 2-ARG REAL FUNCTIONS === */
    reg_real2("gcd", apf_gcd);
    reg_real2("lcm", apf_lcm);
    reg_real2("mod", apf_mod);
    reg_real2("pow", apfx_pow);
    reg_real2("logb", apfx_logb);
    reg_real2("log_", apfx_logb);
    
    /* === 2-ARG COMPLEX FUNCTIONS === */
    reg_complex2("cpow", apfc_pow);
    
    /* Note: bitwise functions (band, bor, bxor, bnot) are declared 
       in apf.h but not implemented yet */
    
    /* === HYPERBOLIC === */
    reg_complex1("sinh", apfc_sinh);
    reg_complex1("cosh", apfc_cosh);
    reg_real1("tanh", apfx_tanh);
    reg_real1("asinh", apfx_asinh);
    reg_real1("acosh", apfx_acosh);
    reg_real1("atanh", apfx_atanh);
    reg_real1("arcsinh", apfx_asinh);
    reg_real1("arccosh", apfx_acosh);
    reg_real1("arctanh", apfx_atanh);
    
    /* === HYPERBOLIC RECIPROCALS === */
    reg_real1("sech", sf_sech);
    reg_real1("csch", sf_csch);
    reg_real1("coth", sf_coth);
    reg_real1("asech", sf_asech);
    reg_real1("acsch", sf_acsch);
    reg_real1("acoth", sf_acoth);
    
    /* === EXP/LOG (complex) === */
    reg_complex1("exp", apfc_exp);
    reg_complex1("log", apfc_log);
    reg_complex1("ln", apfc_log);
    reg_complex1("sqrt", apfc_sqrt);
    
    /* === EXP/LOG VARIANTS (real) === */
    reg_real1("log10", sf_log10);
    reg_real1("log2", sf_log2);
    reg_real1("exp2", sf_exp2);
    reg_real1("pow2", sf_exp2);
    reg_real1("exp10", sf_exp10);
    reg_real1("pow10", sf_exp10);
    reg_real1("log1p", sf_log1p);
    reg_real1("expm1", sf_expm1);
    reg_real1("cbrt", sf_cbrt);
    
    /* === ROUNDING/TRUNCATION === */
    reg_real1("floor", apf_floor);
    reg_real1("ceil", apf_ceil);
    reg_real1("trunc", apf_trunc);
    reg_real1("fix", apf_trunc);
    reg_real1("fabs", apf_abs);
    reg_real1("int", apf_trunc);
    reg_real1("sign", sf_sign);
    reg_real1("sgn", sf_sign);
    reg_real1("round", sf_round);
    
    /* === COMBINATORICS === */
    reg_real1("fact", sf_fact);
    reg_real1("factorial", sf_fact);
    
    /* === ANGLE CONVERSIONS === */
    reg_real1("deg2rad", sf_deg2rad);
    reg_real1("rad2deg", sf_rad2deg);
    reg_real1("wrapToPi", sf_wrapToPi);
    reg_real1("wrapTo180", sf_wrapTo180);
    reg_real1("wrap180", sf_wrapTo180);
    reg_real1("wrapTo360", sf_wrapTo360);
    reg_real1("wrap360", sf_wrapTo360);
    reg_real1("wrapTo360", sf_wrapTo360);
    
    /* === UTILITY === */
    reg_real1("nextpow2", sf_nextpow2);
    reg_real1("realsqrt", sf_realsqrt);
    reg_real1("reallog", sf_reallog);
    reg_real1("sinc", sf_sinc);
    reg_real1("sigmoid", sf_sigmoid);
    reg_real1("logistic", sf_sigmoid);
    reg_real1("softplus", sf_softplus);
    reg_real1("step", sf_heaviside);
    reg_real1("heaviside", sf_heaviside);
    reg_real1("rect", sf_rect);
    reg_real1("rectangle", sf_rect);
    reg_real1("tri", sf_tri);
    reg_real1("triangle", sf_tri);
    reg_real1("signum", sf_sign);  /* alias for sign */
    
    /* === INTEGER SEQUENCE FUNCTIONS === */
    reg_int1("fibonacci", apf_fibonacci, 0, 10000);
    reg_int1("fib", apf_fibonacci, 0, 10000);
    reg_int1("lucas", apf_lucas, 0, 10000);
    reg_int1("catalan", apf_catalan, 0, 500);
    reg_int1("bell", apf_bell, 0, 100);
    reg_int1("partition", apf_partition, 0, 1000);
    reg_int1("partitions", apf_partition, 0, 1000);
    reg_int1("subfactorial", apf_subfactorial, 0, 1000);
    reg_int1("derangements", sf_derangements, 0, 1000);
    reg_int1("digsum", sf_digsum, 0, LONG_MAX);
    reg_int1("digitsum", sf_digsum, 0, LONG_MAX);
    reg_int1("numdigits", sf_numdigits, 0, LONG_MAX);
    reg_int1("ndigits", sf_numdigits, 0, LONG_MAX);
    reg_int1("digitalroot", sf_digitalroot, 0, LONG_MAX);
    
    /* Number theory functions */
    /* NOTE: omega/Omega must be in reverse order because hash table prepends */
    reg_int1("Omega", sf_bigomega, 1, LONG_MAX);
    reg_int1("bigomega", sf_bigomega, 1, LONG_MAX);
    reg_int1("omega", sf_omega, 1, LONG_MAX);
    reg_int1("moebius", sf_moebius, 1, LONG_MAX);
    reg_int1("mobius", sf_moebius, 1, LONG_MAX);
    reg_int1("mu", sf_moebius, 1, LONG_MAX);
    reg_int1("sigma", sf_sigma, 1, LONG_MAX);
    reg_int1("sigma0", sf_sigma0, 1, LONG_MAX);
    reg_int1("numdivisors", sf_sigma0, 1, LONG_MAX);
    reg_int1("tau", sf_sigma0, 1, LONG_MAX);
    reg_int1("digroot", sf_digitalroot, 0, LONG_MAX);
    reg_int1("reversedigits", sf_reversedigits, 0, LONG_MAX);
    reg_int1("revdigits", sf_reversedigits, 0, LONG_MAX);
    reg_int1("ispalindrome", sf_ispalindrome, 0, LONG_MAX);
    reg_int1("isperfect", sf_isperfect, 1, LONG_MAX);
    reg_int1("issquarefree", sf_issquarefree, 1, LONG_MAX);
    reg_int1("totient", sf_totient, 1, LONG_MAX);
    reg_int1("eulerphi", sf_totient, 1, LONG_MAX);
    reg_int1("divisorsum", sf_divisorsum, 1, LONG_MAX);
    reg_int1("even", sf_even, LONG_MIN, LONG_MAX);
    reg_int1("odd", sf_odd, LONG_MIN, LONG_MAX);
    reg_int1("hexagonal", sf_hexagonal, 0, LONG_MAX);
    reg_int1("hex_num", sf_hexagonal, 0, LONG_MAX);
    reg_int1("triangular", sf_triangular, 0, LONG_MAX);
    reg_int1("tri_num", sf_triangular, 0, LONG_MAX);
    reg_int1("pentagonal", sf_pentagonal, 0, LONG_MAX);
    reg_int1("pent_num", sf_pentagonal, 0, LONG_MAX);
    reg_int1("square_num", sf_squarenum, 0, LONG_MAX);
    reg_int1("cube_num", sf_cubenum, 0, LONG_MAX);
    reg_int2("genharmonic", sf_genharmonic);
    reg_int2("harmonic2", sf_genharmonic);
    reg_int2("sigmak", sf_sigmak);
    reg_int2("divisorsum_k", sf_sigmak);
    
    /* === FRACTIONAL PART === */
    reg_real1("frac", sf_frac);
    
    /* === BOOLEAN TEST FUNCTIONS === */
    reg_real1("isnan", sf_isnan);
    reg_real1("isinf", sf_isinf);
    reg_real1("isfinite", sf_isfinite);
    reg_real1("isnumeric", sf_isnumeric);
    reg_real1("isinteger", sf_isinteger);
    reg_real1("ispositive", sf_ispositive);
    reg_real1("isnegative", sf_isnegative);
    
    /* === BITSHIFT FUNCTIONS === */
    reg_apf_int("lsl", sf_lsl);
    reg_apf_int("shl", sf_lsl);
    reg_apf_int("bitshift", sf_lsl);
    reg_apf_int("lsr", sf_lsr);
    reg_apf_int("shr", sf_lsr);
    
    /* === 2-ARG INTEGER FUNCTIONS === */
    reg_int2("npr", apf_npr);
    reg_int2("perm", apf_npr);
    reg_int2("ncr", apf_ncr);
    reg_int2("comb", apf_ncr);
    reg_int2("choose", apf_ncr);
    reg_int2("stirling2", apf_stirling2);
    
    /* === 2-ARG REAL FUNCTIONS === */
    reg_real2("hypot", sf_hypot);
    
    /* === DURATION CONVERSIONS (to seconds) === */
    reg_real1("seconds", sf_seconds);
    reg_real1("minutes", sf_minutes);
    reg_real1("hours", sf_hours);
    reg_real1("hr", sf_hours);
    reg_real1("days", sf_days);
    reg_real1("milliseconds", sf_milliseconds);
    reg_real1("ms", sf_milliseconds);
    reg_real1("weeks", sf_weeks);
    
    /* === STATISTICAL DISTRIBUTIONS === */
    reg_apf_int("tcdf", apfx_tcdf);
    reg_apf_int("student_cdf", apfx_tcdf);
    reg_apf_int("chi2cdf", apfx_chi2cdf);
    reg_apf_int("chicdf", apfx_chi2cdf);
    
    /* === BESSEL FUNCTIONS (int n, apf x) === */
    reg_int_apf("besselj", apfx_besselj);
    reg_int_apf("besselJ", apfx_besselj);
    reg_int_apf("bessely", apfx_bessely);
    reg_int_apf("besselY", apfx_bessely);
    
    /* === MORE 2-ARG REAL FUNCTIONS === */
    reg_real2("beta", sf_beta);
    reg_real2("nthroot", sf_nthroot);
    reg_real2("rem", sf_rem);
    reg_real2("fmod", sf_fmod);
    reg_real2("copysign", sf_copysign);
    reg_real2("fdim", sf_fdim);
    reg_real2("isequal", sf_isequal);
    reg_real2("approxeq", sf_approxeq);
    reg_real2("min2", sf_min2);
    reg_real2("max2", sf_max2);
    
    /* === 1-ARG PROBABILITY FUNCTIONS === */
    reg_real1("normpdf", sf_normpdf);
    
    /* === BITWISE FUNCTIONS === */
    reg_real2("and", apf_and);
    reg_real2("band", apf_and);
    reg_real2("bitand", apf_and);
    reg_real2("or", apf_or);
    reg_real2("bor", apf_or);
    reg_real2("bitor", apf_or);
    reg_real2("xor", apf_xor);
    reg_real2("bxor", apf_xor);
    reg_real2("bitxor", apf_xor);
    reg_real1("not", apf_not);
    reg_real1("bnot", apf_not);
    
    /* === TRIGPI FUNCTIONS === */
    reg_real1("sinpi", sf_sinpi);
    reg_real1("cospi", sf_cospi);
    
    /* === SPECIAL VALUES/CONSTANTS === */
    /* Note: INF, NAN, omega are handled specially in the parser */
    
    /* === COMPARISON FUNCTIONS === */
    reg_real2("eq", sf_eq);
    reg_real2("ne", sf_ne);
    reg_real2("lt", sf_lt);
    reg_real2("le", sf_le);
    reg_real2("gt", sf_gt);
    reg_real2("ge", sf_ge);
    
    /* === NUMBER THEORY FUNCTIONS === */
    reg_real1("isprime", sf_isprime);
    reg_real1("iseven", sf_iseven);
    reg_real1("isodd", sf_isodd);
    reg_real1("nextprime", sf_nextprime);
    reg_real1("prevprime", sf_prevprime);
    reg_real2("modinv", sf_modinv);
    
    /* === LOGICAL FUNCTIONS === */
    reg_real2("nand", sf_nand);
    reg_real2("nor", sf_nor);
    reg_real2("xnor", sf_xnor);
    reg_real2("implies", sf_implies);
    reg_real2("imply", sf_implies);
    
    /* === COMPLEX CONSTRUCTION === */
    reg_make_complex("complex");
    
    /* === COMPLEX PREDICATES === */
    reg_complex1("isreal", sf_isreal);
    reg_real1("fibonacci", sf_fibonacci);
    reg_real1("fib", sf_fibonacci);
    reg_real1("lucas", sf_lucas);
    reg_real1("catalan", sf_catalan);
    reg_real1("factorial2", sf_factorial2);
    reg_real1("fact2", sf_factorial2);
    reg_real1("digamma", sf_digamma);
    reg_real1("psi", sf_digamma);
    reg_real1("zeta", sf_zeta);
    reg_real1("riemann_zeta", sf_zeta);
    reg_real1("harmonic", sf_harmonic);
    /* Note: isprime handled in parser (needs result_is_boolean) */
    
    /* === SPECIAL FUNCTIONS === */
    reg_real1("erf", apfx_erf);
    reg_real1("erfc", apfx_erfc);
    reg_real1("normcdf", apfx_normcdf);
    reg_real1("norminv", apfx_norminv);
    reg_real1("lgamma", apfx_lgamma);
    reg_real1("tgamma", apfx_tgamma);
    reg_real1("gamma", apfx_tgamma);
    
#ifdef HAVE_STATS
    /* === PROBABILITY DISTRIBUTIONS === */
    reg_real2("tpdf", sf_tpdf);
    reg_real2("tcdf", sf_tcdf);
    reg_real2("tinv", sf_tinv);
    reg_real2("chi2pdf", sf_chi2pdf);
    reg_real2("chi2cdf", sf_chi2cdf);
    reg_real2("chi2inv", sf_chi2inv);
    reg_real2("poisspdf", sf_poisspdf);
    reg_real2("poisscdf", sf_poisscdf);
    reg_real2("exppdf", sf_exppdf);
    reg_real2("expcdf", sf_expcdf);
    reg_real2("expinv", sf_expinv);
    reg_real2("sem", sf_sem);
#endif
    
    /* === BESSEL FUNCTIONS === */
    /* Note: j0, j1, y0, y1 etc. don't work because the lexer 
       treats them as variable 'j' followed by number '0'. 
       Use besselJ(n, x) and besselY(n, x) instead. */
    
    /* abs needs special handling for complex */
    reg_special("abs", SF_ABS_COMPLEX);
    
    /* === COMPLEX PARTS === */
    reg_special("real", SF_REAL_PART);
    reg_special("re", SF_REAL_PART);
    reg_special("imag", SF_IMAG_PART);
    reg_special("im", SF_IMAG_PART);
    reg_special("conj", SF_CONJ);
    reg_special("arg", SF_ARG);
    reg_special("angle", SF_ARG);
    reg_special("phase", SF_ARG);
}

ScalarFuncEntry *scalar_func_lookup(const char *name) {
    unsigned h = sf_hash_name(name);
    ScalarFuncEntry *e = sf_hash[h];
    while (e) {
        if (sf_streq(e->name, name)) return e;
        e = e->next;
    }
    return NULL;
}

int scalar_func_eval(ScalarFuncEntry *entry, apfc *result, const apfc *arg) {
    switch (entry->type) {
        case SF_REAL_1:
            entry->fn.f1(&result->re, &arg->re);
            apf_zero(&result->im);
            return 1;
            
        case SF_COMPLEX_1:
            entry->fn.c1(result, arg);
            return 1;
            
        case SF_TRIG: {
            /* Convert input based on angle_mode */
            if (angle_mode != ANGLE_RAD && apfc_is_real(arg)) {
                apfc rad_arg;
                angle_to_rad_local(&rad_arg.re, &arg->re, angle_mode);
                apf_zero(&rad_arg.im);
                entry->fn.c1(result, &rad_arg);
            } else {
                entry->fn.c1(result, arg);
            }
            return 1;
        }
        
        case SF_TRIG_INV: {
            /* Output converted based on angle_mode */
            apf rad_result;
            /* Accept NaN or real numbers */
            if (arg->re.cls == APF_CLASS_NAN) {
                apf_set_nan(&result->re);
                apf_zero(&result->im);
                return 1;
            }
            if (!apfc_is_real(arg)) {
                printf("Error: inverse trig requires real argument\n");
                return 0;
            }
            entry->fn.f1(&rad_result, &arg->re);
            rad_to_angle_local(&result->re, &rad_result, angle_mode);
            apf_zero(&result->im);
            return 1;
        }
        
        case SF_TRIG_DEG: {
            /* Input always in degrees */
            apfc rad_arg;
            apf pi_val, d180;
            apfx_pi(&pi_val);
            apf_from_int(&d180, 180);
            apf_mul(&rad_arg.re, &arg->re, &pi_val);
            apf_div(&rad_arg.re, &rad_arg.re, &d180);
            apf_zero(&rad_arg.im);
            entry->fn.c1(result, &rad_arg);
            return 1;
        }
        
        case SF_TRIG_REAL: {
            /* Real trig with angle mode conversion */
            if (!apfc_is_real(arg)) {
                printf("Error: function requires real argument\n");
                return 0;
            }
            if (angle_mode != ANGLE_RAD) {
                apf rad_arg;
                angle_to_rad_local(&rad_arg, &arg->re, angle_mode);
                entry->fn.f1(&result->re, &rad_arg);
            } else {
                entry->fn.f1(&result->re, &arg->re);
            }
            apf_zero(&result->im);
            return 1;
        }
        
        case SF_TRIG_DEG_REAL: {
            /* Real trig, input always in degrees */
            apf rad_arg, pi_val, d180;
            if (!apfc_is_real(arg)) {
                printf("Error: function requires real argument\n");
                return 0;
            }
            apfx_pi(&pi_val);
            apf_from_int(&d180, 180);
            apf_mul(&rad_arg, &arg->re, &pi_val);
            apf_div(&rad_arg, &rad_arg, &d180);
            entry->fn.f1(&result->re, &rad_arg);
            apf_zero(&result->im);
            return 1;
        }
        
        case SF_TRIG_INV_DEG: {
            /* Output always in degrees */
            apf rad_result, pi_val, d180;
            /* Accept NaN */
            if (arg->re.cls == APF_CLASS_NAN) {
                apf_set_nan(&result->re);
                apf_zero(&result->im);
                return 1;
            }
            if (!apfc_is_real(arg)) {
                printf("Error: inverse trig requires real argument\n");
                return 0;
            }
            entry->fn.f1(&rad_result, &arg->re);
            apfx_pi(&pi_val);
            apf_from_int(&d180, 180);
            apf_mul(&result->re, &rad_result, &d180);
            apf_div(&result->re, &result->re, &pi_val);
            apf_zero(&result->im);
            return 1;
        }
        
        case SF_REAL_PART:
            apf_copy(&result->re, &arg->re);
            apf_zero(&result->im);
            return 1;
            
        case SF_IMAG_PART:
            apf_copy(&result->re, &arg->im);
            apf_zero(&result->im);
            return 1;
            
        case SF_CONJ:
            apf_copy(&result->re, &arg->re);
            apf_neg(&result->im, &arg->im);
            return 1;
            
        case SF_ARG: {
            apfx_atan2(&result->re, &arg->im, &arg->re);
            /* Convert to angle_mode */
            if (angle_mode != ANGLE_RAD) {
                apf temp;
                apf_copy(&temp, &result->re);
                rad_to_angle_local(&result->re, &temp, angle_mode);
            }
            apf_zero(&result->im);
            return 1;
        }
        
        case SF_ABS_COMPLEX: {
            /* Complex magnitude: sqrt(re^2 + im^2) */
            apf re_sq, im_sq, sum;
            apf_mul(&re_sq, &arg->re, &arg->re);
            apf_mul(&im_sq, &arg->im, &arg->im);
            apf_add(&sum, &re_sq, &im_sq);
            apf_sqrt(&result->re, &sum);
            apf_zero(&result->im);
            return 1;
        }
        
        case SF_INT_1: {
            /* Integer sequence function: convert input to long */
            long n;
            if (!apfc_is_real(arg)) {
                printf("Error: function requires real integer\n");
                return 0;
            }
            n = apf_to_long(&arg->re);
            if (n < entry->min_val || n > entry->max_val) {
                printf("Error: argument out of range (%ld to %ld)\n", 
                       entry->min_val, entry->max_val);
                return 0;
            }
            entry->fn.fi(&result->re, n);
            apf_zero(&result->im);
            return 1;
        }
        
        default:
            return 0;
    }
}

int scalar_func_eval2(ScalarFuncEntry *entry, apfc *result, const apfc *arg1, const apfc *arg2) {
    if (entry->type == SF_REAL_2) {
        entry->fn.f2(&result->re, &arg1->re, &arg2->re);
        apf_zero(&result->im);
        return 1;
    } else if (entry->type == SF_COMPLEX_2) {
        entry->fn.c2(result, arg1, arg2);
        return 1;
    } else if (entry->type == SF_TRIG_INV_2) {
        /* 2-arg inverse trig (atan2) with angle mode output */
        entry->fn.f2(&result->re, &arg1->re, &arg2->re);
        if (angle_mode != ANGLE_RAD) {
            apf tmp;
            rad_to_angle_local(&tmp, &result->re, angle_mode);
            apf_copy(&result->re, &tmp);
        }
        apf_zero(&result->im);
        return 1;
    } else if (entry->type == SF_TRIG_INV_2_DEG) {
        /* 2-arg inverse trig always outputting degrees */
        apf tmp;
        entry->fn.f2(&result->re, &arg1->re, &arg2->re);
        rad_to_angle_local(&tmp, &result->re, ANGLE_DEG);
        apf_copy(&result->re, &tmp);
        apf_zero(&result->im);
        return 1;
    } else if (entry->type == SF_INT_2) {
        /* 2-arg integer function: both args converted to long */
        long n = apf_to_long(&arg1->re);
        long k = apf_to_long(&arg2->re);
        entry->fn.fi2(&result->re, n, k);
        apf_zero(&result->im);
        return 1;
    } else if (entry->type == SF_APF_INT) {
        /* Mixed apf + int function */
        long n = apf_to_long(&arg2->re);
        entry->fn.fai(&result->re, &arg1->re, n);
        apf_zero(&result->im);
        return 1;
    } else if (entry->type == SF_INT_APF) {
        /* Int first, apf second (bessel functions) */
        int n = (int)apf_to_long(&arg1->re);
        entry->fn.fia(&result->re, n, &arg2->re);
        apf_zero(&result->im);
        return 1;
    } else if (entry->type == SF_MAKE_COMPLEX) {
        /* complex(re, im) -> result */
        apf_copy(&result->re, &arg1->re);
        apf_copy(&result->im, &arg2->re);
        return 1;
    }
    return 0;
}
