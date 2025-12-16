/* apfx.c - Extended math functions for APF
 * Pure C89, standalone library - no external dependencies
 */
#include "apfx.h"

/* Include config.h when building as part of sc for feature flags */
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

/* If no feature flags defined, enable all (standalone library mode) */
#ifndef HAVE_GAMMA
#define HAVE_GAMMA 1
#endif
#ifndef HAVE_ERF
#define HAVE_ERF 1
#endif
#ifndef HAVE_BESSEL
#define HAVE_BESSEL 1
#endif
#ifndef HAVE_ELLIPTIC
#define HAVE_ELLIPTIC 1
#endif
#ifndef HAVE_LAMBERTW
#define HAVE_LAMBERTW 1
#endif
#ifndef HAVE_DISTRIBUTIONS
#define HAVE_DISTRIBUTIONS 1
#endif

/* Null pointer constant */
#ifndef NULL
#define NULL ((void *)0)
#endif

/* String constants - shorter versions for memory-constrained platforms */
#if defined(SCALC_MEDIUM) || defined(SCALC_TINY) || defined(SCALC_MINIMAL) || defined(SCALC_VIC20)
  /* 40 digits - enough for 128-bit (38 decimal digits) */
  #define STR_LN2      "0.6931471805599453094172321214581765680755"
  #define STR_EULER    "0.5772156649015328606065120900824024310421"
  #define STR_SQRT2    "1.4142135623730950488016887242096980785696"
  #define STR_2SQRTPI  "1.1283791670955125738961589031215451716881"
  #define STR_HALFLN2PI "0.9189385332046727417803297364056176398614"
#else
  /* 100+ digits - for 256-bit and beyond */
  #define STR_LN2      "0.69314718055994530941723212145817656807550013436025525412068000949339362196969471560586332699641868754"
  #define STR_EULER    "0.57721566490153286060651209008240243104215933593992359880576723488486772677766467093694706329174674951"
  #define STR_SQRT2    "1.41421356237309504880168872420969807856967187537694807317667973799073247846210703885038753432764157273"
  #define STR_2SQRTPI  "1.12837916709551257389615890312154517168810125865799771368817144342128493688298682897348732040421472688"
  #define STR_HALFLN2PI "0.91893853320467274178032973640561763986139747363778341281715154048276569592726039769474329863595419762"
#endif

/* exp(x) using Taylor series with aggressive argument reduction
 * Reduce to |x| < 0.5 for fast convergence (only ~25 terms needed)
 * exp(x) = exp(x/2^k)^(2^k)
 */
void apfx_exp(apf *r, const apf *x)
{
    apf sum, term, n_apf, temp, x_reduced;
    INT16 n, max_iter = 30;  /* With |x| < 0.5, 40 iterations is plenty */
    INT16 reductions = 0;
    apf half;
    
    /* Handle special values */
    if (x->cls == APF_CLASS_NAN) {
        apf_set_nan(r);
        return;
    }
    if (x->cls == APF_CLASS_INF) {
        if (x->sign) {
            /* exp(-Inf) = 0 */
            apf_zero(r);
        } else {
            /* exp(+Inf) = +Inf */
            apf_set_inf(r, 0);
        }
        return;
    }
    
    /* Handle zero */
    if (apf_is_zero(x)) {
        apf_from_int(r, 1);
        return;
    }
    
    apf_from_str(&half, "0.5");
    apf_copy(&x_reduced, x);
    
    /* Argument reduction: while |x| > 0.5, halve it */
    while (1) {
        apf abs_x;
        apf_abs(&abs_x, &x_reduced);
        if (apf_le(&abs_x, &half)) break;
        if (reductions >= 40) break;  /* safety limit */
        
        /* x_reduced /= 2 using exponent adjustment (fast!) */
        x_reduced.exp--;
        reductions++;
    }
    
    /* Taylor series: exp(x) = 1 + x + x^2/2! + x^3/3! + ... */
    apf_from_int(&sum, 1);
    apf_copy(&term, &x_reduced);
    apf_add(&temp, &sum, &term);  /* sum = 1 + x */
    sum = temp;
    
    for (n = 2; n < max_iter; n++) {
        apf_mul(&temp, &term, &x_reduced);
        apf_from_int(&n_apf, (INT32)n);  /* n fits in INT16, cast for clarity */
        apf_div(&term, &temp, &n_apf);
        apf_add(&temp, &sum, &term);
        sum = temp;
        if (apf_is_zero(&term)) break;
    }
    
    /* Undo reductions: square the result 'reductions' times */
    for (n = 0; n < reductions; n++) {
        apf_mul(&temp, &sum, &sum);
        sum = temp;
    }
    
    apf_copy(r, &sum);
}

/* Cached ln(2) constant - length depends on platform */
#if defined(SCALC_MEDIUM) || defined(SCALC_TINY) || defined(SCALC_MINIMAL) || defined(SCALC_VIC20)
static const char *ln2_str = "0.6931471805599453094172321214581765680755";
#else
static const char *ln2_str = "0.69314718055994530941723212145817656807550013436025525412068000949339362196969471560586332699641868754200148102057068573368552023575813";
#endif
static apf cached_ln2;
static int ln2_initialized = 0;

/* Forward declarations for cached constants */
static apf cached_pi;
static int pi_initialized;
static apf cached_e;
static int e_initialized;

/* Initialize all static variables - MUST be called at program start on DOS */
void apfx_init(void)
{
    ln2_initialized = 0;
    pi_initialized = 0;
    e_initialized = 0;
}

static void get_ln2(apf *r) {
    if (!ln2_initialized) {
        apf_from_str(&cached_ln2, ln2_str);
        ln2_initialized = 1;
    }
    apf_copy(r, &cached_ln2);
}

/* log(x) using argument reduction and fast series
 * For x = m * 2^k where m in [1, 2), log(x) = log(m) + k*ln(2)
 * For m in [1, 2): use log(m) = 2*atanh((m-1)/(m+1))
 * atanh(y) = y + y^3/3 + y^5/5 + ... converges fast when |y| < 1/3
 */
void apfx_log(apf *r, const apf *x)
{
    apf m, y, y2, term, sum, temp, k_apf, ln2;
    INT16 i, max_iter = 30;  /* With |y| < 1/3, series converges fast */
    INT32 k, n;  /* k can be large (binary exponent), n is small but use INT32 for safety */
    
    /* Handle special values */
    if (x->cls == APF_CLASS_NAN) {
        apf_set_nan(r);
        return;
    }
    if (x->cls == APF_CLASS_INF) {
        if (x->sign) {
            /* ln(-Inf) = NaN (or complex, but return NaN for real) */
            apf_set_nan(r);
        } else {
            /* ln(+Inf) = +Inf */
            apf_set_inf(r, 0);
        }
        return;
    }
    if (apf_is_zero(x)) {
        /* ln(0) = -Inf */
        apf_set_inf(r, 1);
        return;
    }
    if (x->sign) {
        /* ln(negative) = NaN for real (complex path handles this) */
        apf_set_nan(r);
        return;
    }
    
    /* Argument reduction: extract k such that m = x / 2^k is in [1, 2) */
    apf_copy(&m, x);
    k = (INT32)(AP_BITS - 1) + m.exp;
    m.exp = -(AP_BITS - 1);  /* Now m is in [1, 2) */
    
    /* Compute y = (m-1)/(m+1), which is in [0, 1/3) for m in [1, 2) */
    {
        apf one, num, den;
        apf_from_int(&one, 1);
        apf_sub(&num, &m, &one);
        apf_add(&den, &m, &one);
        apf_div(&y, &num, &den);
    }
    
    /* log(m) = 2 * atanh(y) = 2y * (1 + y^2/3 + y^4/5 + y^6/7 + ...) */
    apf_mul(&y2, &y, &y);
    
    apf_from_int(&sum, 1);
    apf_copy(&term, &y2);
    
    for (n = 3, i = 0; i < max_iter; n += 2, i++) {
        apf n_apf, contrib;
        apf_from_int(&n_apf, n);
        apf_div(&contrib, &term, &n_apf);
        apf_add(&temp, &sum, &contrib);
        sum = temp;
        
        apf_mul(&temp, &term, &y2);
        term = temp;
        
        if (apf_is_zero(&contrib)) break;
    }
    
    /* log(m) = 2 * y * sum */
    apf_mul(&temp, &y, &sum);
    temp.exp++;  /* multiply by 2 via exponent */
    
    /* Add k * ln(2) */
    if (k != 0) {
        get_ln2(&ln2);
        apf_from_int(&k_apf, k);
        apf_mul(&sum, &k_apf, &ln2);
        apf_add(r, &temp, &sum);
    } else {
        apf_copy(r, &temp);
    }
}

/* sin(x) using Taylor series */
/* sin(x) using Taylor series with argument reduction
 * Reduce x to [-pi/4, pi/4] for faster convergence
 */
void apfx_sin(apf *r, const apf *x)
{
    apf sum, term, x2, temp, n_apf, arg;
    INT16 n, max_iter = (AP_BITS + 3) / 4;  /* Scale iterations with precision */
    INT16 neg = 0;
    INT16 quad = 0;
    apf pi, pi2, pi_half, pi_quarter;
    
    /* Get pi values */
    apfx_pi(&pi);
    apf_from_int(&temp, 2);
    apf_mul(&pi2, &pi, &temp);
    apf_div(&pi_half, &pi, &temp);
    apf_from_int(&temp, 4);
    apf_div(&pi_quarter, &pi, &temp);
    
    /* Reduce to [0, 2*pi) */
    apf_copy(&arg, x);
    if (!apf_is_zero(&arg)) {
        apf_div(&temp, &arg, &pi2);
        {
            INT32 int_part = apf_to_long(&temp);
            apf int_apf;
            apf_from_int(&int_apf, int_part);
            apf_sub(&temp, &temp, &int_apf);
            apf_mul(&arg, &temp, &pi2);
        }
        /* Handle negative */
        if (arg.sign) {
            apf_add(&arg, &arg, &pi2);
        }
    }
    
    /* Determine quadrant and reduce to [0, pi/2) */
    /* quad: 0=[0,pi/2), 1=[pi/2,pi), 2=[pi,3pi/2), 3=[3pi/2,2pi) */
    if (apf_cmp(&arg, &pi_half) < 0) {
        quad = 0;
    } else if (apf_cmp(&arg, &pi) < 0) {
        quad = 1;
        apf_sub(&arg, &pi, &arg);
    } else {
        apf three_pi_half;
        apf_from_int(&temp, 3);
        apf_mul(&three_pi_half, &pi_half, &temp);
        if (apf_cmp(&arg, &three_pi_half) < 0) {
            quad = 2;
            apf_sub(&arg, &arg, &pi);
        } else {
            quad = 3;
            apf_sub(&arg, &pi2, &arg);
        }
    }
    
    /* Taylor series: sin(x) = x - x^3/3! + x^5/5! - ... */
    apf_copy(&sum, &arg);
    apf_copy(&term, &arg);
    apf_mul(&x2, &arg, &arg);
    
    for (n = 1; n < max_iter; n++) {
        /* NOTE: (2*n)*(2*n+1) overflows INT16 at n=91. Use INT32 arithmetic. */
        INT32 divisor = (2L * n) * (2L * n + 1);
        apf_mul(&temp, &term, &x2);
        apf_from_int(&n_apf, divisor);
        apf_div(&term, &temp, &n_apf);
        neg = !neg;
        if (neg) {
            apf_sub(&temp, &sum, &term);
        } else {
            apf_add(&temp, &sum, &term);
        }
        sum = temp;
        if (apf_is_zero(&term)) break;
    }
    
    /* Apply quadrant correction */
    if (quad == 2 || quad == 3) {
        sum.sign = !sum.sign;
    }
    
    apf_copy(r, &sum);
}

/* cos(x) using Taylor series with argument reduction */
void apfx_cos(apf *r, const apf *x)
{
    apf sum, term, x2, temp, n_apf, arg;
    INT16 n, max_iter = 30;
    INT16 neg = 0;
    INT16 quad = 0;
    apf pi, pi2, pi_half;
    
    /* Get pi values */
    apfx_pi(&pi);
    apf_from_int(&temp, 2);
    apf_mul(&pi2, &pi, &temp);
    apf_div(&pi_half, &pi, &temp);
    
    /* Reduce to [0, 2*pi) */
    apf_copy(&arg, x);
    if (!apf_is_zero(&arg)) {
        apf_div(&temp, &arg, &pi2);
        {
            INT32 int_part = apf_to_long(&temp);
            apf int_apf;
            apf_from_int(&int_apf, int_part);
            apf_sub(&temp, &temp, &int_apf);
            apf_mul(&arg, &temp, &pi2);
        }
        if (arg.sign) {
            apf_add(&arg, &arg, &pi2);
        }
    }
    
    /* Determine quadrant */
    if (apf_cmp(&arg, &pi_half) < 0) {
        quad = 0;
    } else if (apf_cmp(&arg, &pi) < 0) {
        quad = 1;
        apf_sub(&arg, &pi, &arg);
    } else {
        apf three_pi_half;
        apf_from_int(&temp, 3);
        apf_mul(&three_pi_half, &pi_half, &temp);
        if (apf_cmp(&arg, &three_pi_half) < 0) {
            quad = 2;
            apf_sub(&arg, &arg, &pi);
        } else {
            quad = 3;
            apf_sub(&arg, &pi2, &arg);
        }
    }
    
    /* Taylor series: cos(x) = 1 - x^2/2! + x^4/4! - ... */
    apf_from_int(&sum, 1);
    apf_from_int(&term, 1);
    apf_mul(&x2, &arg, &arg);
    
    for (n = 1; n < max_iter; n++) {
        /* NOTE: (2n-1)*(2n) overflows INT16 at n=91. Use INT32 arithmetic. */
        INT32 divisor = (2L * n - 1) * (2L * n);
        apf_mul(&temp, &term, &x2);
        apf_from_int(&n_apf, divisor);
        apf_div(&term, &temp, &n_apf);
        neg = !neg;
        if (neg) {
            apf_sub(&temp, &sum, &term);
        } else {
            apf_add(&temp, &sum, &term);
        }
        sum = temp;
        if (apf_is_zero(&term)) break;
    }
    
    /* Apply quadrant correction */
    if (quad == 1 || quad == 2) {
        sum.sign = !sum.sign;
    }
    
    apf_copy(r, &sum);
}

/* Combined sin and cos - shares argument reduction */
void apfx_sincos(apf *sin_r, apf *cos_r, const apf *x)
{
    apf sin_sum, cos_sum, sin_term, cos_term, x2, temp, n_apf, arg;
    INT16 n, max_iter = 30;  /* With arg in [0,pi/2], converges in ~35 terms */
    INT16 sin_neg = 0, cos_neg = 0;
    INT16 quad = 0;
    apf pi, pi2, pi_half;
    
    /* Get pi values */
    apfx_pi(&pi);
    apf_from_int(&temp, 2);
    apf_mul(&pi2, &pi, &temp);
    apf_div(&pi_half, &pi, &temp);
    
    /* Reduce to [0, 2*pi) */
    apf_copy(&arg, x);
    if (!apf_is_zero(&arg)) {
        apf_div(&temp, &arg, &pi2);
        {
            INT32 int_part = apf_to_long(&temp);
            apf int_apf;
            apf_from_int(&int_apf, int_part);
            apf_sub(&temp, &temp, &int_apf);
            apf_mul(&arg, &temp, &pi2);
        }
        if (arg.sign) {
            apf_add(&arg, &arg, &pi2);
        }
    }
    
    /* Determine quadrant */
    if (apf_cmp(&arg, &pi_half) < 0) {
        quad = 0;
    } else if (apf_cmp(&arg, &pi) < 0) {
        quad = 1;
        apf_sub(&arg, &pi, &arg);
    } else {
        apf three_pi_half;
        apf_from_int(&temp, 3);
        apf_mul(&three_pi_half, &pi_half, &temp);
        if (apf_cmp(&arg, &three_pi_half) < 0) {
            quad = 2;
            apf_sub(&arg, &arg, &pi);
        } else {
            quad = 3;
            apf_sub(&arg, &pi2, &arg);
        }
    }
    
    apf_mul(&x2, &arg, &arg);
    
    /* sin Taylor series: sin(x) = x - x^3/3! + x^5/5! - ... */
    apf_copy(&sin_sum, &arg);
    apf_copy(&sin_term, &arg);
    
    /* cos Taylor series: cos(x) = 1 - x^2/2! + x^4/4! - ... */
    apf_from_int(&cos_sum, 1);
    apf_from_int(&cos_term, 1);
    
    for (n = 1; n < max_iter; n++) {
        /* NOTE: These products overflow INT16 at n~91. Use INT32 arithmetic. */
        INT32 sin_div = (2L * n) * (2L * n + 1);
        INT32 cos_div = (2L * n - 1) * (2L * n);
        
        /* sin term */
        apf_mul(&temp, &sin_term, &x2);
        apf_from_int(&n_apf, sin_div);
        apf_div(&sin_term, &temp, &n_apf);
        sin_neg = !sin_neg;
        if (sin_neg) {
            apf_sub(&temp, &sin_sum, &sin_term);
        } else {
            apf_add(&temp, &sin_sum, &sin_term);
        }
        sin_sum = temp;
        
        /* cos term */
        apf_mul(&temp, &cos_term, &x2);
        apf_from_int(&n_apf, cos_div);
        apf_div(&cos_term, &temp, &n_apf);
        cos_neg = !cos_neg;
        if (cos_neg) {
            apf_sub(&temp, &cos_sum, &cos_term);
        } else {
            apf_add(&temp, &cos_sum, &cos_term);
        }
        cos_sum = temp;
        
        if (apf_is_zero(&sin_term) && apf_is_zero(&cos_term)) break;
    }
    
    /* Apply quadrant corrections */
    /* sin: negate in quads 2,3 */
    if (quad == 2 || quad == 3) {
        sin_sum.sign = !sin_sum.sign;
    }
    /* cos: negate in quads 1,2 */
    if (quad == 1 || quad == 2) {
        cos_sum.sign = !cos_sum.sign;
    }
    
    if (sin_r) apf_copy(sin_r, &sin_sum);
    if (cos_r) apf_copy(cos_r, &cos_sum);
}

/* tan(x) = sin(x) / cos(x) */
void apfx_tan(apf *r, const apf *x)
{
    apf s, c;
    apfx_sincos(&s, &c, x);  /* Combined function - one argument reduction */
    apf_div(r, &s, &c);
}

/* Combined sinh and cosh - only 1 exp call using exp(-x) = 1/exp(x) */
void apfx_sinhcosh(apf *sinh_r, apf *cosh_r, const apf *x)
{
    apf ex, emx, two, temp;
    
    /* Compute exp(x) */
    apfx_exp(&ex, x);
    
    /* exp(-x) = 1/exp(x) - much faster than second exp call */
    {
        apf one;
        apf_from_int(&one, 1);
        apf_div(&emx, &one, &ex);
    }
    
    apf_from_int(&two, 2);
    
    if (sinh_r) {
        apf_sub(&temp, &ex, &emx);
        apf_div(sinh_r, &temp, &two);
    }
    if (cosh_r) {
        apf_add(&temp, &ex, &emx);
        apf_div(cosh_r, &temp, &two);
    }
}

/* sinh(x) = (exp(x) - exp(-x)) / 2 */
void apfx_sinh(apf *r, const apf *x)
{
    apfx_sinhcosh(r, NULL, x);
}

/* cosh(x) = (exp(x) + exp(-x)) / 2 */
void apfx_cosh(apf *r, const apf *x)
{
    apfx_sinhcosh(NULL, r, x);
}

/* tanh(x) = sinh(x) / cosh(x) */
void apfx_tanh(apf *r, const apf *x)
{
    apf sinh_val, cosh_val;
    
    /* Handle special values */
    if (x->cls == APF_CLASS_NAN) {
        apf_set_nan(r);
        return;
    }
    if (x->cls == APF_CLASS_INF) {
        /* tanh(±Inf) = ±1 */
        apf_from_int(r, 1);
        r->sign = x->sign;
        return;
    }
    if (apf_is_zero(x)) {
        apf_zero(r);
        return;
    }
    
    apfx_sinhcosh(&sinh_val, &cosh_val, x);
    apf_div(r, &sinh_val, &cosh_val);
}

/* asin(x) = atan(x / sqrt(1 - x^2)) for |x| <= 1
 * Returns NaN for |x| > 1
 */
void apfx_asin(apf *r, const apf *x)
{
    apf abs_x, one, temp, x2, sqrt_term;
    
    /* Handle special values */
    if (x->cls == APF_CLASS_NAN) {
        apf_set_nan(r);
        return;
    }
    if (x->cls == APF_CLASS_INF) {
        apf_set_nan(r);  /* asin(Inf) = NaN */
        return;
    }
    if (apf_is_zero(x)) {
        apf_zero(r);
        return;
    }
    
    apf_from_int(&one, 1);
    apf_abs(&abs_x, x);
    
    /* Check |x| > 1 */
    if (apf_cmp(&abs_x, &one) > 0) {
        apf_set_nan(r);  /* Domain error */
        return;
    }
    
    /* asin(1) = pi/2, asin(-1) = -pi/2 */
    if (apf_cmp(&abs_x, &one) == 0) {
        apf pi;
        apfx_pi(&pi);
        apf_from_int(&temp, 2);
        apf_div(r, &pi, &temp);
        r->sign = x->sign;
        return;
    }
    
    /* asin(x) = atan(x / sqrt(1 - x^2)) */
    apf_mul(&x2, x, x);
    apf_sub(&temp, &one, &x2);
    apf_sqrt(&sqrt_term, &temp);
    apf_div(&temp, x, &sqrt_term);
    apfx_atan(r, &temp);
}

/* acos(x) = pi/2 - asin(x) for |x| <= 1
 * Returns NaN for |x| > 1
 */
void apfx_acos(apf *r, const apf *x)
{
    apf asin_val, pi_half, pi, two;
    
    /* Handle special values */
    if (x->cls == APF_CLASS_NAN) {
        apf_set_nan(r);
        return;
    }
    if (x->cls == APF_CLASS_INF) {
        apf_set_nan(r);  /* acos(Inf) = NaN */
        return;
    }
    
    /* acos(1) = 0 */
    {
        apf one, abs_x;
        apf_from_int(&one, 1);
        apf_abs(&abs_x, x);
        
        if (apf_cmp(&abs_x, &one) > 0) {
            apf_set_nan(r);  /* Domain error */
            return;
        }
        
        if (!x->sign && apf_cmp(x, &one) == 0) {
            apf_zero(r);
            return;
        }
    }
    
    /* acos(x) = pi/2 - asin(x) */
    apfx_asin(&asin_val, x);
    apfx_pi(&pi);
    apf_from_int(&two, 2);
    apf_div(&pi_half, &pi, &two);
    apf_sub(r, &pi_half, &asin_val);
}

/* atan(x) using Taylor series with argument reduction
 * Uses half-angle formula: atan(x) = 2*atan(x/(1+sqrt(1+x²)))
 * This reduces the argument, making Taylor series converge faster
 */
void apfx_atan(apf *r, const apf *x)
{
    apf sum, term, x2, temp, n_apf, one, half;
    apf x_work;
    INT16 n, max_iter = 30;  /* With |x| < 0.5, converges quickly */
    INT16 neg = 0;
    INT16 reductions = 0;
    
    /* Handle special values */
    if (x->cls == APF_CLASS_NAN) {
        apf_set_nan(r);
        return;
    }
    
    if (apf_is_zero(x)) {
        apf_zero(r);
        return;
    }
    
    /* atan(±Inf) = ±pi/2 */
    if (x->cls == APF_CLASS_INF) {
        apf pi;
        apfx_pi(&pi);
        pi.exp -= 1;  /* pi/2 */
        if (x->sign) {
            pi.sign = 1;  /* -pi/2 */
        }
        apf_copy(r, &pi);
        return;
    }
    
    apf_from_int(&one, 1);
    apf_from_str(&half, "0.5");
    apf_copy(&x_work, x);
    
    /* Handle negative values */
    if (x_work.sign) {
        x_work.sign = 0;
        neg = 1;
    }
    
    /* Reduce argument using half-angle formula until |x| < 0.5
     * atan(x) = 2*atan(x/(1+sqrt(1+x²)))
     */
    while (1) {
        apf abs_x;
        apf_abs(&abs_x, &x_work);
        if (apf_le(&abs_x, &half)) break;
        if (reductions >= 20) break;
        
        {
            apf x2_local, sqrt_term, denom;
            apf_mul(&x2_local, &x_work, &x_work);
            apf_add(&temp, &one, &x2_local);
            apf_sqrt(&sqrt_term, &temp);
            apf_add(&denom, &one, &sqrt_term);
            apf_div(&x_work, &x_work, &denom);
        }
        reductions++;
    }
    
    /* Taylor series: atan(x) = x - x³/3 + x⁵/5 - x⁷/7 + ... */
    apf_copy(&sum, &x_work);
    apf_copy(&term, &x_work);
    apf_mul(&x2, &x_work, &x_work);
    
    for (n = 1; n < max_iter; n++) {
        /* 2*n+1 is safe: max 61 for n=30, well under INT16 limit */
        INT16 sign_flip = (n & 1);
        apf_mul(&temp, &term, &x2);
        term = temp;
        apf_from_int(&n_apf, (INT32)(2 * n + 1));
        apf_div(&temp, &term, &n_apf);
        if (sign_flip) {
            apf_sub(&sum, &sum, &temp);
        } else {
            apf_add(&sum, &sum, &temp);
        }
        if (apf_is_zero(&temp)) break;
    }
    
    /* Undo reductions: multiply by 2^reductions using exponent */
    sum.exp += reductions;
    
    /* Apply sign */
    if (neg) {
        sum.sign = 1;
    }
    apf_copy(r, &sum);
}

/* atan2(y, x) */
void apfx_atan2(apf *r, const apf *y, const apf *x)
{
    apf pi, pi2, ratio;
    
    apfx_pi(&pi);
    apf_from_int(&pi2, 2);
    apf_div(&pi2, &pi, &pi2);
    
    if (apf_is_zero(x)) {
        if (apf_is_zero(y)) {
            apf_zero(r);
        } else if (y->sign) {
            apf_neg(r, &pi2);
        } else {
            apf_copy(r, &pi2);
        }
        return;
    }
    
    apf_div(&ratio, y, x);
    apfx_atan(r, &ratio);
    
    if (x->sign) {
        if (y->sign) {
            apf_sub(r, r, &pi);
        } else {
            apf_add(r, r, &pi);
        }
    }
}

/* pow(base, exp) = exp(exp * log(base)) */
void apfx_pow(apf *r, const apf *base, const apf *exp)
{
    apf log_base, temp;
    apf one;
    apf abs_base, abs_exp;
    int cmp;
    
    apf_from_int(&one, 1);
    
    /* IEEE 754-2008 Section 9.2.1: x^0 = 1 for ANY x, including NaN
     * This must be checked BEFORE NaN propagation */
    if (apf_is_zero(exp)) {
        apf_from_int(r, 1);
        return;
    }
    
    /* IEEE 754-2008 Section 9.2.1: 1^y = 1 for ANY y, including NaN and Inf
     * This must be checked BEFORE NaN propagation */
    if (base->cls == APF_CLASS_NORMAL || base->cls == APF_CLASS_ZERO) {
        apf_abs(&abs_base, base);
        if (apf_cmp(&abs_base, &one) == 0 && !base->sign) {
            /* base is exactly +1 */
            apf_from_int(r, 1);
            return;
        }
    }
    
    /* Handle NaN (after IEEE special cases) */
    if (base->cls == APF_CLASS_NAN || exp->cls == APF_CLASS_NAN) {
        apf_set_nan(r);
        return;
    }
    
    /* 0^x: positive exp -> 0, negative exp -> Inf */
    if (apf_is_zero(base)) {
        if (exp->sign) {
            apf_set_inf(r, 0);
        } else {
            apf_zero(r);
        }
        return;
    }
    
    /* Handle Inf exponent */
    if (exp->cls == APF_CLASS_INF) {
        apf_abs(&abs_base, base);
        cmp = apf_cmp(&abs_base, &one);
        
        if (cmp > 0) {
            /* |base| > 1: Inf exp -> Inf, -Inf exp -> 0 */
            if (exp->sign) {
                apf_zero(r);
            } else {
                apf_set_inf(r, 0);
            }
        } else if (cmp < 0) {
            /* |base| < 1: Inf exp -> 0, -Inf exp -> Inf */
            if (exp->sign) {
                apf_set_inf(r, 0);
            } else {
                apf_zero(r);
            }
        } else {
            /* |base| == 1: IEEE 754 says 1^Inf = 1 */
            apf_from_int(r, 1);
        }
        return;
    }
    
    /* Handle Inf base */
    if (base->cls == APF_CLASS_INF) {
        if (exp->sign) {
            /* Inf^(-x) = 0 for positive x */
            apf_zero(r);
        } else {
            /* Inf^x = Inf for positive x
             * But (-Inf)^n = -Inf if n is an odd integer */
            int neg_result = 0;
            if (base->sign) {
                /* Check if exponent is an odd integer */
                long n = apf_to_long(exp);
                apf n_apf, diff;
                apf_from_int(&n_apf, n);
                apf_sub(&diff, exp, &n_apf);
                if (apf_is_zero(&diff) && (n & 1)) {
                    neg_result = 1;
                }
            }
            apf_set_inf(r, neg_result);
        }
        return;
    }
    
    /* Early check for extreme exponents beyond INT32 range.
     * If |exp| > ~10^11, the result will definitely overflow to Inf or underflow to 0.
     * For normalized APF: value ≈ 2^(stored_exp + AP_BITS-1)
     * Threshold: stored_exp >= -(AP_BITS-1) + 37 means |exp| >= 2^37 ≈ 10^11
     * This catches exponents like 10^30 but allows INT32-range values to compute.
     */
    apf_abs(&abs_exp, exp);
    if (abs_exp.cls == APF_CLASS_NORMAL && abs_exp.exp >= -(AP_BITS - 1) + 37) {
        /* Exponent is astronomically large (> ~10^9) */
        apf_abs(&abs_base, base);
        
        if (apf_cmp(&abs_base, &one) > 0) {
            /* |base| > 1: large positive exp -> Inf, large negative exp -> 0 */
            if (exp->sign) {
                apf_zero(r);
            } else {
                apf_set_inf(r, 0);
            }
        } else if (apf_cmp(&abs_base, &one) < 0) {
            /* |base| < 1: large positive exp -> 0, large negative exp -> Inf */
            if (exp->sign) {
                apf_set_inf(r, 0);
            } else {
                apf_zero(r);
            }
        } else {
            /* |base| == 1: result is 1 (or -1 for negative base, but that's complex) */
            apf_from_int(r, 1);
        }
        return;
    }
    
    /* Check for integer exponent */
    {
        INT32 n = apf_to_long(exp);
        apf n_apf, diff;
        apf_from_int(&n_apf, n);
        apf_sub(&diff, exp, &n_apf);
        
        /* Check if it's an integer - use binary exponentiation */
        if (apf_is_zero(&diff)) {
            apf result, b;
            INT16 neg_exp = 0;
            UINT32 abs_n;  /* Use unsigned to handle INT32_MIN correctly */
            
            if (n < 0) {
                neg_exp = 1;
                /* Handle INT32_MIN carefully: -(-2147483648) overflows INT32 */
                /* Cast to unsigned first, then negate works correctly */
                abs_n = (UINT32)(-(n + 1)) + 1U;
            } else {
                abs_n = (UINT32)n;
            }
            
            /* Integer power - use binary exponentiation O(log n) */
            apf_from_int(&result, 1);
            apf_copy(&b, base);
            
            while (abs_n > 0) {
                if (abs_n & 1) {
                    apf_mul(&temp, &result, &b);
                    result = temp;
                    /* Check for overflow/underflow */
                    if (result.cls == APF_CLASS_INF) {
                        if (neg_exp) {
                            apf_zero(r);
                        } else {
                            *r = result;
                        }
                        return;
                    }
                    if (result.cls == APF_CLASS_ZERO) {
                        if (neg_exp) {
                            r->cls = APF_CLASS_INF;
                            r->sign = 0;
                        } else {
                            apf_zero(r);
                        }
                        return;
                    }
                }
                abs_n >>= 1;
                if (abs_n > 0) {
                    apf_mul(&temp, &b, &b);
                    b = temp;
                    /* Check for overflow in base */
                    if (b.cls == APF_CLASS_INF) {
                        /* Result depends on remaining bits of exponent */
                        if (neg_exp) {
                            apf_zero(r);
                        } else {
                            r->cls = APF_CLASS_INF;
                            r->sign = 0;
                        }
                        return;
                    }
                    if (b.cls == APF_CLASS_ZERO) {
                        if (neg_exp) {
                            r->cls = APF_CLASS_INF;
                            r->sign = 0;
                        } else {
                            apf_zero(r);
                        }
                        return;
                    }
                }
            }
            
            if (neg_exp) {
                apf_div(r, &one, &result);
            } else {
                apf_copy(r, &result);
            }
            return;
        }
    }
    
    /* General case: exp(exp * log(base)) */
    apfx_log(&log_base, base);
    apf_mul(&temp, exp, &log_base);
    apfx_exp(r, &temp);
}

/* Cached pi constant - 40 decimal digits */
#if defined(SCALC_MEDIUM) || defined(SCALC_TINY) || defined(SCALC_MINIMAL) || defined(SCALC_VIC20)
static const char *pi_str = "3.1415926535897932384626433832795028841971";
#else
static const char *pi_str = "3.14159265358979323846264338327950288419716939937510582097494459230781640628620899862803482534211706798214808651328230664709384460955058";
#endif

/* pi - returns cached value */
void apfx_pi(apf *r)
{
    if (!pi_initialized) {
        apf_from_str(&cached_pi, pi_str);
        pi_initialized = 1;
    }
    apf_copy(r, &cached_pi);
}

/* Cached e constant - 40 decimal digits */
#if defined(SCALC_MEDIUM) || defined(SCALC_TINY) || defined(SCALC_MINIMAL) || defined(SCALC_VIC20)
static const char *e_str = "2.7182818284590452353602874713526624977572";
#else
static const char *e_str = "2.71828182845904523536028747135266249775724709369995957496696762772407663035354759457138217852516642742746639193200305992181741359662904";
#endif

/* e - returns cached value */
void apfx_e(apf *r)
{
    if (!e_initialized) {
        apf_from_str(&cached_e, e_str);
        e_initialized = 1;
    }
    apf_copy(r, &cached_e);
}

/* factorial */
void apfx_fact(apf *r, long n)
{
    apf result, i_apf, temp;
    long i;
    
    apf_from_int(&result, 1);
    
    for (i = 2; i <= n; i++) {
        apf_from_int(&i_apf, i);
        apf_mul(&temp, &result, &i_apf);
        result = temp;
    }
    
    apf_copy(r, &result);
}

/* ========== Inverse Hyperbolic Functions ========== */

/* asinh(x) = ln(x + sqrt(x^2 + 1)) */
void apfx_asinh(apf *r, const apf *x)
{
    apf x2, tmp, one;
    
    apf_from_int(&one, 1);
    apf_mul(&x2, x, x);        /* x^2 */
    apf_add(&tmp, &x2, &one);  /* x^2 + 1 */
    apf_sqrt(&tmp, &tmp);     /* sqrt(x^2 + 1) */
    apf_add(&tmp, x, &tmp);    /* x + sqrt(x^2 + 1) */
    apfx_log(r, &tmp);         /* ln(...) */
}

/* acosh(x) = ln(x + sqrt(x^2 - 1)) for x >= 1 */
void apfx_acosh(apf *r, const apf *x)
{
    apf x2, tmp, one;
    
    /* Check domain: x >= 1 */
    apf_from_int(&one, 1);
    if (apf_cmp(x, &one) < 0) {
        apf_set_nan(r);
        return;
    }
    
    apf_mul(&x2, x, x);        /* x^2 */
    apf_sub(&tmp, &x2, &one);  /* x^2 - 1 */
    apf_sqrt(&tmp, &tmp);     /* sqrt(x^2 - 1) */
    apf_add(&tmp, x, &tmp);    /* x + sqrt(x^2 - 1) */
    apfx_log(r, &tmp);         /* ln(...) */
}

/* atanh(x) = 0.5 * ln((1 + x) / (1 - x)) for |x| < 1 */
void apfx_atanh(apf *r, const apf *x)
{
    apf one, numer, denom, tmp, half;
    
    apf_from_int(&one, 1);
    
    /* Check domain: |x| < 1 */
    apf_abs(&tmp, x);
    if (apf_cmp(&tmp, &one) >= 0) {
        if (apf_cmp(&tmp, &one) == 0) {
            /* atanh(±1) = ±Inf */
            apf_set_inf(r, x->sign);
        } else {
            apf_set_nan(r);
        }
        return;
    }
    
    apf_add(&numer, &one, x);    /* 1 + x */
    apf_sub(&denom, &one, x);    /* 1 - x */
    apf_div(&tmp, &numer, &denom);  /* (1 + x) / (1 - x) */
    apfx_log(&tmp, &tmp);        /* ln(...) */
    
    /* Multiply by 0.5 */
    apf_from_int(&half, 1);
    apf_from_int(&one, 2);
    apf_div(&half, &half, &one);  /* 0.5 */
    apf_mul(r, &tmp, &half);
}

/* ================================================================
 * PHASE 2: Scientific Essentials
 * ================================================================ */

#ifdef HAVE_GAMMA
/* lgamma(x) - Log-gamma function using Stirling's approximation
 * with Bernoulli number corrections for high precision
 * Returns ln(|Gamma(x)|)
 */
void apfx_lgamma(apf *r, const apf *x)
{
    apf z, result, tmp, tmp2;
    int shift = 0;
    
    /* Handle special cases */
    if (x->cls == APF_CLASS_NAN) { apf_set_nan(r); return; }
    if (x->cls == APF_CLASS_INF) {
        apf_set_inf(r, 0);  /* lgamma(+Inf) = +Inf */
        return;
    }
    
    /* For x <= 0, use reflection formula */
    if (x->sign || apf_is_zero(x)) {
        apf one, pi_val, pi_x, sin_pix, abs_sin, ln_pi, lg_1mx;
        
        if (apf_is_zero(x)) {
            apf_set_inf(r, 0);
            return;
        }
        
        apf_from_int(&one, 1);
        apfx_pi(&pi_val);
        apf_mul(&pi_x, &pi_val, x);
        apfx_sin(&sin_pix, &pi_x);
        apf_abs(&abs_sin, &sin_pix);
        
        if (apf_is_zero(&abs_sin)) {
            apf_set_inf(r, 0);  /* Pole at negative integers */
            return;
        }
        
        apfx_log(&ln_pi, &pi_val);
        apf_sub(&tmp, &one, x);  /* 1 - x */
        apfx_lgamma(&lg_1mx, &tmp);
        apfx_log(&tmp, &abs_sin);
        
        /* lgamma(x) = ln(pi) - ln|sin(pi*x)| - lgamma(1-x) */
        apf_sub(r, &ln_pi, &tmp);
        apf_sub(r, r, &lg_1mx);
        return;
    }
    
    apf_copy(&z, x);
    
    /* Use recurrence Gamma(x+1) = x*Gamma(x) to shift x to larger value
     * where Stirling approximation is more accurate.
     * lgamma(x) = lgamma(x+n) - ln(x*(x+1)*...*(x+n-1))
     */
    {
        apf eight;
        apf_from_int(&eight, 8);
        while (apf_cmp(&z, &eight) < 0) {
            shift++;
            {
                apf one;
                apf_from_int(&one, 1);
                apf_add(&z, &z, &one);
            }
        }
    }
    
    /* Stirling's approximation:
     * lgamma(z) = (z - 0.5)*ln(z) - z + 0.5*ln(2*pi)
     *           + 1/(12*z) - 1/(360*z^3) + 1/(1260*z^5) - ...
     */
    {
        apf half, ln_z, z2, z3, z5, z7;
        
        apf_from_str(&half, "0.5");
        apfx_log(&ln_z, &z);
        
        /* (z - 0.5) * ln(z) */
        apf_sub(&tmp, &z, &half);
        apf_mul(&result, &tmp, &ln_z);
        
        /* - z */
        apf_sub(&result, &result, &z);
        
        /* + 0.5 * ln(2*pi) = 0.9189385332046727... */
        apf_from_str(&tmp, STR_HALFLN2PI);
        apf_add(&result, &result, &tmp);
        
        /* Bernoulli corrections */
        apf_mul(&z2, &z, &z);     /* z^2 */
        apf_mul(&z3, &z2, &z);    /* z^3 */
        apf_mul(&z5, &z3, &z2);   /* z^5 */
        apf_mul(&z7, &z5, &z2);   /* z^7 */
        
        /* + 1/(12*z) */
        apf_from_int(&tmp, 12);
        apf_mul(&tmp, &tmp, &z);
        apf_from_int(&tmp2, 1);
        apf_div(&tmp, &tmp2, &tmp);
        apf_add(&result, &result, &tmp);
        
        /* - 1/(360*z^3) */
        apf_from_int(&tmp, 360);
        apf_mul(&tmp, &tmp, &z3);
        apf_div(&tmp, &tmp2, &tmp);
        apf_sub(&result, &result, &tmp);
        
        /* + 1/(1260*z^5) */
        apf_from_int(&tmp, 1260);
        apf_mul(&tmp, &tmp, &z5);
        apf_div(&tmp, &tmp2, &tmp);
        apf_add(&result, &result, &tmp);
        
        /* - 1/(1680*z^7) */
        apf_from_int(&tmp, 1680);
        apf_mul(&tmp, &tmp, &z7);
        apf_div(&tmp, &tmp2, &tmp);
        apf_sub(&result, &result, &tmp);
    }
    
    /* Undo the shift: subtract ln(x*(x+1)*...*(x+n-1)) */
    if (shift > 0) {
        apf prod, i_apf;
        int i;
        
        apf_copy(&prod, x);
        for (i = 1; i < shift; i++) {
            apf_from_int(&i_apf, i);
            apf_add(&tmp, x, &i_apf);
            apf_mul(&prod, &prod, &tmp);
        }
        apfx_log(&tmp, &prod);
        apf_sub(&result, &result, &tmp);
    }
    
    apf_copy(r, &result);
}

/* tgamma(x) - Gamma function */
void apfx_tgamma(apf *r, const apf *x)
{
    apf lg;
    
    /* Handle special cases */
    if (x->cls == APF_CLASS_NAN) { apf_set_nan(r); return; }
    if (x->cls == APF_CLASS_INF) {
        if (x->sign) apf_set_nan(r);  /* gamma(-Inf) = NaN */
        else apf_set_inf(r, 0);       /* gamma(+Inf) = +Inf */
        return;
    }
    
    /* Check for negative integers (poles) */
    if (x->sign) {
        long n = apf_to_long(x);
        apf n_apf, diff;
        apf_from_int(&n_apf, n);
        apf_sub(&diff, x, &n_apf);
        if (apf_is_zero(&diff) && n <= 0) {
            apf_set_inf(r, 0);  /* Pole at non-positive integers */
            return;
        }
    }
    
    /* gamma(x) = exp(lgamma(x)) */
    apfx_lgamma(&lg, x);
    apfx_exp(r, &lg);
    
    /* Handle sign for negative x using reflection */
    if (x->sign) {
        apf pi_val, pi_x, sin_pix;
        apfx_pi(&pi_val);
        apf_mul(&pi_x, &pi_val, x);
        apfx_sin(&sin_pix, &pi_x);
        if (sin_pix.sign) {
            r->sign = 1;
            sin_pix.sign = 0;
        }
    }
}
#endif /* HAVE_GAMMA */

#ifdef HAVE_ERF
/* erf(x) - Error function using Taylor series and asymptotic approximation
 * erf(x) = (2/sqrt(pi)) * integral from 0 to x of exp(-t^2) dt
 */
void apfx_erf(apf *r, const apf *x)
{
    apf ax, result, one;
    int neg;
    
    if (x->cls == APF_CLASS_NAN) { apf_set_nan(r); return; }
    if (x->cls == APF_CLASS_INF) {
        apf_from_int(r, x->sign ? -1 : 1);
        return;
    }
    
    neg = x->sign;
    apf_abs(&ax, x);
    apf_from_int(&one, 1);
    
    /* For all x, use Taylor series with sufficient terms
     * erf(x) = (2/sqrt(pi)) * sum_{n=0}^inf (-1)^n * x^(2n+1) / (n! * (2n+1))
     */
    {
        apf sum, term, x2, x_pow, two_sqrt_pi, fact_n, tmp;
        int n;
        
        apf_mul(&x2, &ax, &ax);  /* x^2 */
        apf_copy(&sum, &ax);     /* First term: x */
        apf_copy(&x_pow, &ax);   /* x^1 */
        apf_from_int(&fact_n, 1);
        
        for (n = 1; n < 100; n++) {
            apf n2p1, n_apf;
            apf_mul(&x_pow, &x_pow, &x2);  /* x^(2n+1) */
            apf_neg(&x_pow, &x_pow);       /* alternating sign */
            apf_from_int(&n_apf, n);
            apf_mul(&fact_n, &fact_n, &n_apf);  /* n! */
            apf_from_int(&n2p1, 2*n + 1);
            apf_mul(&term, &fact_n, &n2p1);
            apf_div(&term, &x_pow, &term);
            apf_add(&sum, &sum, &term);
            
            /* Check convergence */
            apf_abs(&tmp, &term);
            if (tmp.exp < sum.exp - AP_BITS + 5) break;
        }
        
        /* Multiply by 2/sqrt(pi) = 1.1283791670955125739... */
        apf_from_str(&two_sqrt_pi, STR_2SQRTPI);
        apf_mul(&result, &sum, &two_sqrt_pi);
        
        /* Clamp to [-1, 1] */
        if (apf_cmp(&result, &one) > 0) apf_copy(&result, &one);
        {
            apf neg_one;
            apf_from_int(&neg_one, -1);
            if (apf_cmp(&result, &neg_one) < 0) apf_copy(&result, &neg_one);
        }
    }
    
    if (neg) apf_neg(r, &result);
    else apf_copy(r, &result);
}

/* erfc(x) - Complementary error function: 1 - erf(x) */
void apfx_erfc(apf *r, const apf *x)
{
    apf erf_val, one;
    
    if (x->cls == APF_CLASS_NAN) { apf_set_nan(r); return; }
    if (x->cls == APF_CLASS_INF) {
        apf_from_int(r, x->sign ? 2 : 0);
        return;
    }
    
    apfx_erf(&erf_val, x);
    apf_from_int(&one, 1);
    apf_sub(r, &one, &erf_val);
}

/* normcdf(x) - Standard normal CDF: Phi(x) = 0.5 * (1 + erf(x/sqrt(2))) */
void apfx_normcdf(apf *r, const apf *x)
{
    apf sqrt2, x_scaled, erf_val, one, half, two;
    
    if (x->cls == APF_CLASS_NAN) { apf_set_nan(r); return; }
    if (x->cls == APF_CLASS_INF) {
        apf_from_int(r, x->sign ? 0 : 1);
        return;
    }
    
    /* sqrt(2) = 1.4142135623730950488... */
    apf_from_str(&sqrt2, STR_SQRT2);
    apf_div(&x_scaled, x, &sqrt2);  /* x / sqrt(2) */
    
    apfx_erf(&erf_val, &x_scaled);
    
    apf_from_int(&one, 1);
    apf_from_int(&two, 2);
    apf_div(&half, &one, &two);
    
    apf_add(r, &one, &erf_val);
    apf_mul(r, r, &half);
}

/* norminv(p) - Inverse standard normal CDF (probit function)
 * Uses Acklam's rational approximation for p in (0,1)
 */
void apfx_norminv(apf *r, const apf *p)
{
    apf zero, one, half, low, high, tmp;
    
    apf_from_int(&zero, 0);
    apf_from_int(&one, 1);
    
    /* Check bounds */
    if (p->cls == APF_CLASS_NAN || apf_cmp(p, &zero) <= 0 || apf_cmp(p, &one) >= 0) {
        if (apf_cmp(p, &zero) == 0) apf_set_inf(r, 1);  /* -Inf */
        else if (apf_cmp(p, &one) == 0) apf_set_inf(r, 0);  /* +Inf */
        else apf_set_nan(r);
        return;
    }
    
    /* Constants for Acklam's approximation */
    apf_from_int(&tmp, 2);
    apf_div(&half, &one, &tmp);
    apf_from_str(&low, "0.02425");
    apf_from_str(&high, "0.97575");
    
    /* Central region */
    if (apf_cmp(p, &low) > 0 && apf_cmp(p, &high) < 0) {
        /* Rational approximation for central region */
        static const char *a[] = {
            "-3.969683028665376e+01", "2.209460984245205e+02",
            "-2.759285104469687e+02", "1.383577518672690e+02",
            "-3.066479806614716e+01", "2.506628277459239e+00"
        };
        static const char *b[] = {
            "-5.447609879822406e+01", "1.615858368580409e+02",
            "-1.556989798598866e+02", "6.680131188771972e+01",
            "-1.328068155288572e+01"
        };
        
        apf q, q2, numer, denom;
        int i;
        
        apf_sub(&q, p, &half);  /* q = p - 0.5 */
        apf_mul(&q2, &q, &q);   /* q^2 */
        
        /* Horner's method for numerator */
        apf_from_str(&numer, a[0]);
        for (i = 1; i < 6; i++) {
            apf_mul(&numer, &numer, &q2);
            apf_from_str(&tmp, a[i]);
            apf_add(&numer, &numer, &tmp);
        }
        apf_mul(&numer, &numer, &q);
        
        /* Horner's method for denominator */
        apf_from_str(&denom, b[0]);
        for (i = 1; i < 5; i++) {
            apf_mul(&denom, &denom, &q2);
            apf_from_str(&tmp, b[i]);
            apf_add(&denom, &denom, &tmp);
        }
        apf_mul(&denom, &denom, &q2);
        apf_add(&denom, &denom, &one);
        
        apf_div(r, &numer, &denom);
    } else {
        /* Tail regions */
        static const char *c[] = {
            "-7.784894002430293e-03", "-3.223964580411365e-01",
            "-2.400758277161838e+00", "-2.549732539343734e+00",
            "4.374664141464968e+00", "2.938163982698783e+00"
        };
        static const char *d[] = {
            "7.784695709041462e-03", "3.224671290700398e-01",
            "2.445134137142996e+00", "3.754408661907416e+00"
        };
        
        apf q, numer, denom;
        int i, is_low;
        
        is_low = (apf_cmp(p, &low) <= 0);
        
        if (is_low) {
            apfx_log(&q, p);
            apf_neg(&q, &q);
            apf_from_int(&tmp, 2);
            apf_mul(&q, &q, &tmp);
            apf_sqrt(&q, &q);  /* q = sqrt(-2*ln(p)) */
        } else {
            apf one_minus_p;
            apf_sub(&one_minus_p, &one, p);
            apfx_log(&q, &one_minus_p);
            apf_neg(&q, &q);
            apf_from_int(&tmp, 2);
            apf_mul(&q, &q, &tmp);
            apf_sqrt(&q, &q);  /* q = sqrt(-2*ln(1-p)) */
        }
        
        /* Horner's method for numerator */
        apf_from_str(&numer, c[0]);
        for (i = 1; i < 6; i++) {
            apf_mul(&numer, &numer, &q);
            apf_from_str(&tmp, c[i]);
            apf_add(&numer, &numer, &tmp);
        }
        
        /* Horner's method for denominator */
        apf_from_str(&denom, d[0]);
        for (i = 1; i < 4; i++) {
            apf_mul(&denom, &denom, &q);
            apf_from_str(&tmp, d[i]);
            apf_add(&denom, &denom, &tmp);
        }
        apf_mul(&denom, &denom, &q);
        apf_add(&denom, &denom, &one);
        
        apf_div(r, &numer, &denom);
        apf_sub(r, &q, r);
        if (!is_low) apf_neg(r, r);
    }
}
#endif /* HAVE_ERF */

#ifdef HAVE_BESSEL
/* Bessel J0(x) - Bessel function of first kind, order 0 */
void apfx_j0(apf *r, const apf *x)
{
    apf ax, x2, sum, term, tmp;
    int n;
    
    if (x->cls == APF_CLASS_NAN) { apf_set_nan(r); return; }
    if (x->cls == APF_CLASS_INF) { apf_zero(r); return; }
    
    apf_abs(&ax, x);
    apf_mul(&x2, &ax, &ax);  /* x^2 */
    
    /* For small x, use Taylor series:
     * J0(x) = sum_{n=0}^inf (-1)^n * (x/2)^(2n) / (n!)^2 */
    {
        apf threshold;
        apf_from_int(&threshold, 8);
        
        if (apf_cmp(&ax, &threshold) <= 0) {
            apf x_half_sq, fact_n_sq;
            
            /* (x/2)^2 = x^2/4 */
            apf_from_int(&tmp, 4);
            apf_div(&x_half_sq, &x2, &tmp);
            
            apf_from_int(&sum, 1);  /* n=0 term */
            apf_from_int(&term, 1);
            apf_from_int(&fact_n_sq, 1);
            
            for (n = 1; n < 50; n++) {
                apf n_apf;
                apf_from_int(&n_apf, n);
                apf_mul(&fact_n_sq, &fact_n_sq, &n_apf);
                apf_mul(&fact_n_sq, &fact_n_sq, &n_apf);  /* (n!)^2 */
                
                apf_mul(&term, &term, &x_half_sq);
                apf_neg(&term, &term);  /* alternating sign */
                
                apf_div(&tmp, &term, &fact_n_sq);
                apf_add(&sum, &sum, &tmp);
                
                /* Check convergence */
                apf_abs(&tmp, &tmp);
                if (tmp.exp < sum.exp - AP_BITS + 10) break;
            }
            
            apf_copy(r, &sum);
            return;
        }
    }
    
    /* For large x, use asymptotic expansion:
     * J0(x) ~ sqrt(2/(pi*x)) * cos(x - pi/4) */
    {
        apf pi_val, two, sqrt_2_pi_x, phase, cos_phase;
        
        apfx_pi(&pi_val);
        apf_from_int(&two, 2);
        
        /* sqrt(2/(pi*x)) */
        apf_mul(&tmp, &pi_val, &ax);
        apf_div(&tmp, &two, &tmp);
        apf_sqrt(&sqrt_2_pi_x, &tmp);
        
        /* phase = x - pi/4 */
        apf_from_int(&tmp, 4);
        apf_div(&tmp, &pi_val, &tmp);
        apf_sub(&phase, &ax, &tmp);
        
        apfx_cos(&cos_phase, &phase);
        apf_mul(r, &sqrt_2_pi_x, &cos_phase);
    }
}

/* Bessel J1(x) - Bessel function of first kind, order 1 */
void apfx_j1(apf *r, const apf *x)
{
    apf ax, x2, sum, term, tmp;
    int n, neg;
    
    if (x->cls == APF_CLASS_NAN) { apf_set_nan(r); return; }
    if (x->cls == APF_CLASS_INF) { apf_zero(r); return; }
    
    neg = x->sign;
    apf_abs(&ax, x);
    apf_mul(&x2, &ax, &ax);  /* x^2 */
    
    /* For small x, use Taylor series:
     * J1(x) = (x/2) * sum_{n=0}^inf (-1)^n * (x/2)^(2n) / (n! * (n+1)!) */
    {
        apf threshold;
        apf_from_int(&threshold, 8);
        
        if (apf_cmp(&ax, &threshold) <= 0) {
            apf x_half, x_half_sq, fact_n, fact_np1;
            
            apf_from_int(&tmp, 2);
            apf_div(&x_half, &ax, &tmp);  /* x/2 */
            apf_mul(&x_half_sq, &x_half, &x_half);  /* (x/2)^2 */
            
            apf_from_int(&sum, 1);  /* n=0 term coefficient */
            apf_from_int(&term, 1);
            apf_from_int(&fact_n, 1);
            apf_from_int(&fact_np1, 1);
            
            for (n = 1; n < 50; n++) {
                apf n_apf, np1_apf;
                apf_from_int(&n_apf, n);
                apf_from_int(&np1_apf, n + 1);
                apf_mul(&fact_n, &fact_n, &n_apf);      /* n! */
                apf_mul(&fact_np1, &fact_np1, &np1_apf); /* (n+1)! */
                
                apf_mul(&term, &term, &x_half_sq);
                apf_neg(&term, &term);  /* alternating sign */
                
                {
                    apf denom;
                    apf_mul(&denom, &fact_n, &fact_np1);
                    apf_div(&tmp, &term, &denom);
                }
                apf_add(&sum, &sum, &tmp);
                
                /* Check convergence */
                apf_abs(&tmp, &tmp);
                if (tmp.exp < sum.exp - AP_BITS + 10) break;
            }
            
            apf_mul(r, &x_half, &sum);
            if (neg) apf_neg(r, r);
            return;
        }
    }
    
    /* For large x, use asymptotic expansion:
     * J1(x) ~ sqrt(2/(pi*x)) * cos(x - 3*pi/4) */
    {
        apf pi_val, two, three, sqrt_2_pi_x, phase, cos_phase;
        
        apfx_pi(&pi_val);
        apf_from_int(&two, 2);
        apf_from_int(&three, 3);
        
        /* sqrt(2/(pi*x)) */
        apf_mul(&tmp, &pi_val, &ax);
        apf_div(&tmp, &two, &tmp);
        apf_sqrt(&sqrt_2_pi_x, &tmp);
        
        /* phase = x - 3*pi/4 */
        apf_mul(&tmp, &three, &pi_val);
        {
            apf four;
            apf_from_int(&four, 4);
            apf_div(&tmp, &tmp, &four);
        }
        apf_sub(&phase, &ax, &tmp);
        
        apfx_cos(&cos_phase, &phase);
        apf_mul(r, &sqrt_2_pi_x, &cos_phase);
        if (neg) apf_neg(r, r);
    }
}

/* Bessel function J_n(x) for integer order n using recurrence */
void apfx_besselj(apf *r, int n, const apf *x)
{
    apf jn, jn1, jn2, coef;
    int k, neg;
    
    if (x->cls == APF_CLASS_NAN) { apf_set_nan(r); return; }
    
    neg = 0;
    if (n < 0) {
        n = -n;
        neg = (n % 2);  /* J_{-n}(x) = (-1)^n * J_n(x) */
    }
    
    if (n == 0) {
        apfx_j0(r, x);
        if (neg) apf_neg(r, r);
        return;
    }
    if (n == 1) {
        apfx_j1(r, x);
        if (neg) apf_neg(r, r);
        return;
    }
    
    /* Use upward recurrence: J_{n+1} = (2n/x)*J_n - J_{n-1} */
    apfx_j0(&jn2, x);
    apfx_j1(&jn1, x);
    
    for (k = 1; k < n; k++) {
        /* jn = (2*k/x) * jn1 - jn2 */
        apf_from_int(&coef, 2 * k);
        apf_div(&coef, &coef, x);
        apf_mul(&jn, &coef, &jn1);
        apf_sub(&jn, &jn, &jn2);
        
        apf_copy(&jn2, &jn1);
        apf_copy(&jn1, &jn);
    }
    
    apf_copy(r, &jn1);
    if (neg) apf_neg(r, r);
}

/* Bessel function Y_n(x) for integer order n using recurrence */
void apfx_bessely(apf *r, int n, const apf *x)
{
    apf yn, yn1, yn2, coef;
    int k, neg;
    
    if (x->cls == APF_CLASS_NAN) { apf_set_nan(r); return; }
    
    /* Y_n is only defined for x > 0 */
    if (x->cls == APF_CLASS_ZERO || x->sign < 0) {
        apf_set_nan(r);
        return;
    }
    
    neg = 0;
    if (n < 0) {
        n = -n;
        neg = (n % 2);  /* Y_{-n}(x) = (-1)^n * Y_n(x) */
    }
    
    if (n == 0) {
        apfx_y0(r, x);
        if (neg) apf_neg(r, r);
        return;
    }
    if (n == 1) {
        apfx_y1(r, x);
        if (neg) apf_neg(r, r);
        return;
    }
    
    /* Use upward recurrence: Y_{n+1} = (2n/x)*Y_n - Y_{n-1} */
    apfx_y0(&yn2, x);
    apfx_y1(&yn1, x);
    
    for (k = 1; k < n; k++) {
        /* yn = (2*k/x) * yn1 - yn2 */
        apf_from_int(&coef, 2 * k);
        apf_div(&coef, &coef, x);
        apf_mul(&yn, &coef, &yn1);
        apf_sub(&yn, &yn, &yn2);
        
        apf_copy(&yn2, &yn1);
        apf_copy(&yn1, &yn);
    }
    
    apf_copy(r, &yn1);
    if (neg) apf_neg(r, r);
}
#endif /* HAVE_BESSEL */

#ifdef HAVE_ELLIPTIC
/* Complete elliptic integral of the first kind K(m)
 * Using AGM (Arithmetic-Geometric Mean) method
 * K(m) = pi / (2 * AGM(1, sqrt(1-m)))
 */
void apfx_ellipk(apf *r, const apf *m)
{
    apf one, a, b, tmp, pi_val;
    int iter;
    
    if (m->cls == APF_CLASS_NAN) { apf_set_nan(r); return; }
    
    apf_from_int(&one, 1);
    
    /* Check bounds: m must be < 1 */
    if (apf_cmp(m, &one) >= 0) {
        apf_set_inf(r, 0);  /* K(1) = +Inf */
        return;
    }
    if (m->sign) {
        /* For negative m, use transformation */
        apf_set_nan(r);  /* TODO: extend to negative m */
        return;
    }
    
    /* AGM iteration: a_0 = 1, b_0 = sqrt(1-m) */
    apf_copy(&a, &one);
    apf_sub(&tmp, &one, m);
    apf_sqrt(&b, &tmp);
    
    for (iter = 0; iter < 50; iter++) {
        apf a_new, b_new, diff;
        
        /* a_{n+1} = (a_n + b_n) / 2 */
        apf_add(&a_new, &a, &b);
        {
            apf two;
            apf_from_int(&two, 2);
            apf_div(&a_new, &a_new, &two);
        }
        
        /* b_{n+1} = sqrt(a_n * b_n) */
        apf_mul(&tmp, &a, &b);
        apf_sqrt(&b_new, &tmp);
        
        /* Check convergence */
        apf_sub(&diff, &a_new, &b_new);
        apf_abs(&diff, &diff);
        if (diff.exp < a_new.exp - AP_BITS + 10) {
            apf_copy(&a, &a_new);
            break;
        }
        
        apf_copy(&a, &a_new);
        apf_copy(&b, &b_new);
    }
    
    /* K(m) = pi / (2 * AGM) */
    apfx_pi(&pi_val);
    {
        apf two;
        apf_from_int(&two, 2);
        apf_mul(&tmp, &two, &a);
    }
    apf_div(r, &pi_val, &tmp);
}

/* Complete elliptic integral of the second kind E(m)
 * E(m) = K(m) * (1 - sum of c_n^2 * 2^(n-1))
 * where c_n are differences from AGM
 */
void apfx_ellipe(apf *r, const apf *m)
{
    apf one, a, b, c, tmp, pi_val, sum_c2;
    apf two_pow;
    int iter;
    
    if (m->cls == APF_CLASS_NAN) { apf_set_nan(r); return; }
    
    apf_from_int(&one, 1);
    
    /* E(0) = pi/2 */
    if (apf_is_zero(m)) {
        apfx_pi(&pi_val);
        apf_from_int(&tmp, 2);
        apf_div(r, &pi_val, &tmp);
        return;
    }
    
    /* E(1) = 1 */
    if (apf_cmp(m, &one) == 0) {
        apf_copy(r, &one);
        return;
    }
    
    if (apf_cmp(m, &one) > 0 || m->sign) {
        apf_set_nan(r);  /* TODO: extend domain */
        return;
    }
    
    /* AGM with c_n tracking */
    apf_copy(&a, &one);
    apf_sub(&tmp, &one, m);
    apf_sqrt(&b, &tmp);
    apf_sub(&c, &a, &b);  /* c_0 = 1 - sqrt(1-m) */
    
    apf_zero(&sum_c2);
    apf_from_int(&two_pow, 1);  /* 2^0 = 1 for first iteration */
    
    for (iter = 0; iter < 50; iter++) {
        apf a_new, b_new, c_new, c2_term;
        
        /* Add c_n^2 * 2^(n-1) to sum */
        apf_mul(&c2_term, &c, &c);
        if (iter > 0) {
            apf_mul(&c2_term, &c2_term, &two_pow);
        }
        apf_add(&sum_c2, &sum_c2, &c2_term);
        
        /* Update 2^n for next iteration */
        {
            apf two;
            apf_from_int(&two, 2);
            apf_mul(&two_pow, &two_pow, &two);
        }
        
        /* AGM step */
        apf_add(&a_new, &a, &b);
        {
            apf two;
            apf_from_int(&two, 2);
            apf_div(&a_new, &a_new, &two);
        }
        apf_mul(&tmp, &a, &b);
        apf_sqrt(&b_new, &tmp);
        apf_sub(&c_new, &a_new, &b_new);
        
        /* Check convergence */
        apf_abs(&tmp, &c_new);
        if (tmp.exp < a_new.exp - AP_BITS + 10) break;
        
        apf_copy(&a, &a_new);
        apf_copy(&b, &b_new);
        apf_copy(&c, &c_new);
    }
    
    /* E(m) = K(m) * (1 - sum_c2/2) = (pi/(2*AGM)) * (1 - sum_c2/2) */
    apfx_pi(&pi_val);
    {
        apf two, denom, factor;
        apf_from_int(&two, 2);
        apf_mul(&denom, &two, &a);
        apf_div(&tmp, &pi_val, &denom);  /* K(m) = pi/(2*AGM) */
        
        apf_div(&sum_c2, &sum_c2, &two);
        apf_sub(&factor, &one, &sum_c2);
        apf_mul(r, &tmp, &factor);
    }
}
#endif /* HAVE_ELLIPTIC */

/* ================================================================
 * PHASE 4: Advanced Functions
 * ================================================================ */

#ifdef HAVE_LAMBERTW
/* Lambert W function W_0(x) - principal branch
 * Solves W * exp(W) = x using Halley's iteration
 */
void apfx_lambertw(apf *r, const apf *x)
{
    apf w, e_w, f, f_prime, f_dprime, tmp, delta;
    apf neg_inv_e;
    int iter;
    
    if (x->cls == APF_CLASS_NAN) { apf_set_nan(r); return; }
    if (x->cls == APF_CLASS_INF) {
        if (x->sign) apf_set_nan(r);  /* W(-Inf) undefined */
        else apf_set_inf(r, 0);       /* W(+Inf) = +Inf */
        return;
    }
    
    /* W(0) = 0 */
    if (apf_is_zero(x)) {
        apf_zero(r);
        return;
    }
    
    /* Domain check: W_0(x) defined for x >= -1/e */
    apf_from_str(&neg_inv_e, "-0.36787944117144232159552377016146086744581113103177");
    if (apf_cmp(x, &neg_inv_e) < 0) {
        apf_set_nan(r);  /* x < -1/e: no real solution */
        return;
    }
    
    /* Initial guess */
    if (x->sign) {
        /* For -1/e <= x < 0, use series: W(x) ~ x - x^2 + 1.5*x^3 */
        apf x2;
        apf_mul(&x2, x, x);
        apf_sub(&w, x, &x2);
    } else {
        /* For x > 0, use log approximation: W(x) ~ ln(1+x) */
        apf one;
        apf_from_int(&one, 1);
        apf_add(&tmp, &one, x);
        apfx_log(&w, &tmp);
    }
    
    /* Halley's iteration:
     * w_{n+1} = w_n - f/(f' - f*f''/(2*f'))
     * where f = w*e^w - x, f' = e^w*(1+w), f'' = e^w*(2+w)
     */
    for (iter = 0; iter < 50; iter++) {
        apfx_exp(&e_w, &w);  /* e^w */
        
        apf_mul(&f, &w, &e_w);
        apf_sub(&f, &f, x);  /* f = w*e^w - x */
        
        {
            apf one, one_plus_w, two_plus_w;
            apf_from_int(&one, 1);
            apf_add(&one_plus_w, &one, &w);
            apf_mul(&f_prime, &e_w, &one_plus_w);  /* f' = e^w*(1+w) */
            
            {
                apf two;
                apf_from_int(&two, 2);
                apf_add(&two_plus_w, &two, &w);
            }
            apf_mul(&f_dprime, &e_w, &two_plus_w);  /* f'' = e^w*(2+w) */
        }
        
        /* Halley correction */
        {
            apf numer, denom, two;
            apf_mul(&numer, &f, &f_dprime);
            apf_from_int(&two, 2);
            apf_mul(&denom, &two, &f_prime);
            apf_div(&tmp, &numer, &denom);
            apf_sub(&denom, &f_prime, &tmp);
            apf_div(&delta, &f, &denom);
        }
        
        apf_sub(&w, &w, &delta);
        
        /* Check convergence */
        apf_abs(&tmp, &delta);
        if (tmp.exp < w.exp - AP_BITS + 10) break;
    }
    
    apf_copy(r, &w);
}
#endif /* HAVE_LAMBERTW */

#ifdef HAVE_DISTRIBUTIONS
/* Gauss-Legendre numerical integration
 * Integrates f(x) from a to b using n-point quadrature
 * This is a helper - actual integration needs function pointer
 * Returns sum of weights * f(nodes), caller transforms interval
 */

/* 5-point Gauss-Legendre nodes and weights (on [-1,1]) */
static const char *gl5_nodes[] = {
    "0",
    "-0.5384693101056830910363144207002",
    "0.5384693101056830910363144207002",
    "-0.9061798459386639927976268782993",
    "0.9061798459386639927976268782993"
};
static const char *gl5_weights[] = {
    "0.5688888888888888888888888888889",
    "0.4786286704993664680412915148356",
    "0.4786286704993664680412915148356",
    "0.2369268850561890875142640407200",
    "0.2369268850561890875142640407200"
};

/* Get Gauss-Legendre node (index 0..4 for 5-point) */
void apfx_gl_node(apf *r, int n, int i)
{
    (void)n;  /* For now, only 5-point implemented */
    if (i < 0 || i >= 5) {
        apf_zero(r);
        return;
    }
    apf_from_str(r, gl5_nodes[i]);
}

/* Get Gauss-Legendre weight */
void apfx_gl_weight(apf *r, int n, int i)
{
    (void)n;
    if (i < 0 || i >= 5) {
        apf_zero(r);
        return;
    }
    apf_from_str(r, gl5_weights[i]);
}

/* ================================================================
 * Beta Function and Incomplete Beta
 * ================================================================ */

/* Beta(a,b) = Gamma(a)*Gamma(b)/Gamma(a+b) */
void apfx_beta(apf *r, const apf *a, const apf *b)
{
    apf ga, gb, gab, sum;
    
    if (a->cls == APF_CLASS_NAN || b->cls == APF_CLASS_NAN) {
        apf_set_nan(r);
        return;
    }
    
    /* Use log-gamma for numerical stability: Beta = exp(lgamma(a) + lgamma(b) - lgamma(a+b)) */
    apfx_lgamma(&ga, a);
    apfx_lgamma(&gb, b);
    apf_add(&sum, a, b);
    apfx_lgamma(&gab, &sum);
    
    apf_add(&sum, &ga, &gb);
    apf_sub(&sum, &sum, &gab);
    apfx_exp(r, &sum);
}

/* Regularized incomplete gamma P(a,x) = gamma(a,x)/Gamma(a)
 * Uses series for x < a+1, continued fraction otherwise
 */
void apfx_gammainc(apf *r, const apf *a, const apf *x)
{
    apf sum, term, aa, tmp, ga;
    int n;
    
    if (x->cls == APF_CLASS_NAN || a->cls == APF_CLASS_NAN) {
        apf_set_nan(r);
        return;
    }
    if (x->sign) {
        apf_set_nan(r);  /* x must be >= 0 */
        return;
    }
    if (apf_is_zero(x)) {
        apf_zero(r);
        return;
    }
    
    /* Series representation: P(a,x) = (x^a * e^(-x) / Gamma(a)) * sum_{n=0}^inf x^n / (a)_n
     * where (a)_n = a*(a+1)*...*(a+n-1) is rising factorial */
    apfx_pow(&term, x, a);
    apf_neg(&tmp, x);
    apfx_exp(&tmp, &tmp);
    apf_mul(&term, &term, &tmp);
    
    apfx_tgamma(&ga, a);
    apf_div(&term, &term, &ga);
    
    /* Sum: 1 + x/(a+1) + x^2/((a+1)(a+2)) + ... */
    apf_from_int(&sum, 1);
    apf_from_int(&tmp, 1);  /* Current term coefficient */
    apf_copy(&aa, a);
    
    for (n = 1; n < 200; n++) {
        apf one;
        apf_from_int(&one, 1);
        apf_add(&aa, &aa, &one);  /* aa = a + n */
        apf_mul(&tmp, &tmp, x);
        apf_div(&tmp, &tmp, &aa);
        apf_add(&sum, &sum, &tmp);
        
        /* Check convergence */
        if (tmp.cls == APF_CLASS_ZERO || tmp.exp < sum.exp - AP_BITS + 5) break;
    }
    
    apf_mul(r, &term, &sum);
    
    /* Clamp to [0,1] */
    {
        apf one;
        apf_from_int(&one, 1);
        if (apf_cmp(r, &one) > 0) apf_copy(r, &one);
    }
    if (r->sign) apf_zero(r);
}

/* Regularized incomplete beta I_x(a,b)
 * Uses continued fraction representation
 */
void apfx_betainc(apf *r, const apf *x, const apf *a, const apf *b)
{
    apf front, f, c, d, tmp, m_apf, one, two;
    int m;
    
    if (x->cls == APF_CLASS_NAN || a->cls == APF_CLASS_NAN || b->cls == APF_CLASS_NAN) {
        apf_set_nan(r);
        return;
    }
    
    /* Bounds check */
    apf_from_int(&one, 1);
    if (x->sign || apf_cmp(x, &one) > 0) {
        apf_set_nan(r);
        return;
    }
    if (apf_is_zero(x)) {
        apf_zero(r);
        return;
    }
    if (apf_cmp(x, &one) == 0) {
        apf_copy(r, &one);
        return;
    }
    
    /* front = x^a * (1-x)^b / (a * Beta(a,b)) */
    apfx_pow(&front, x, a);
    apf_sub(&tmp, &one, x);
    apfx_pow(&tmp, &tmp, b);
    apf_mul(&front, &front, &tmp);
    apfx_beta(&tmp, a, b);
    apf_mul(&tmp, &tmp, a);
    apf_div(&front, &front, &tmp);
    
    /* Lentz's continued fraction algorithm */
    apf_from_str(&f, "1e-30");
    apf_copy(&c, &f);
    apf_from_int(&d, 0);
    apf_from_int(&two, 2);
    
    for (m = 1; m <= 100; m++) {
        apf d_term, delta;
        
        /* Compute d_m term in continued fraction */
        apf_from_int(&m_apf, m);
        
        {
            apf num, den, amb;
            /* d_{2m} = m(b-m)x / ((a+2m-1)(a+2m)) */
            apf m2;
            apf_mul(&m2, &two, &m_apf);  /* 2m */
            
            apf_sub(&tmp, b, &m_apf);   /* b-m */
            apf_mul(&num, &m_apf, &tmp);
            apf_mul(&num, &num, x);
            
            apf_add(&amb, a, &m2);
            apf_sub(&tmp, &amb, &one);
            apf_mul(&den, &tmp, &amb);
            
            apf_div(&d_term, &num, &den);
        }
        
        /* Update c and d */
        apf_add(&tmp, &one, &d_term);
        apf_div(&tmp, &tmp, &c);
        apf_copy(&c, &tmp);
        
        apf_add(&tmp, &one, &d_term);
        apf_mul(&tmp, &d, &tmp);
        apf_add(&tmp, &tmp, &one);
        apf_div(&d, &one, &tmp);
        
        apf_mul(&delta, &c, &d);
        apf_mul(&f, &f, &delta);
        
        /* Check convergence */
        apf_sub(&tmp, &delta, &one);
        apf_abs(&tmp, &tmp);
        if (tmp.exp < -AP_BITS + 10) break;
    }
    
    apf_mul(r, &front, &f);
    
    /* Clamp to [0,1] */
    if (apf_cmp(r, &one) > 0) apf_copy(r, &one);
    if (r->sign) apf_zero(r);
}

/* ================================================================
 * Statistical Distributions
 * ================================================================ */

/* Student's t CDF using incomplete beta */
void apfx_tcdf(apf *r, const apf *t, long df)
{
    apf x, half, df_apf, t2, tmp, one;
    
    if (t->cls == APF_CLASS_NAN || df < 1) {
        apf_set_nan(r);
        return;
    }
    
    apf_from_int(&one, 1);
    apf_from_str(&half, "0.5");
    apf_from_int(&df_apf, df);
    
    /* x = df / (df + t^2) */
    apf_mul(&t2, t, t);
    apf_add(&tmp, &df_apf, &t2);
    apf_div(&x, &df_apf, &tmp);
    
    /* CDF = 1 - 0.5 * I_x(df/2, 0.5) for t >= 0 */
    {
        apf a, b, ibeta;
        apf_mul(&a, &df_apf, &half);  /* df/2 */
        apf_copy(&b, &half);           /* 0.5 */
        
        apfx_betainc(&ibeta, &x, &a, &b);
        apf_mul(&ibeta, &ibeta, &half);
        
        if (t->sign) {
            apf_copy(r, &ibeta);  /* t < 0: CDF = 0.5 * I_x */
        } else {
            apf_sub(r, &one, &ibeta);  /* t >= 0: CDF = 1 - 0.5 * I_x */
        }
    }
}

/* Inverse Student's t using Newton-Raphson */
void apfx_tinv(apf *r, const apf *p, long df)
{
    apf t, ft, delta, one, half, pdf_val;
    int iter;
    
    if (p->cls == APF_CLASS_NAN || df < 1) {
        apf_set_nan(r);
        return;
    }
    
    apf_from_int(&one, 1);
    apf_from_str(&half, "0.5");
    
    /* Check bounds */
    if (p->sign || apf_cmp(p, &one) > 0) {
        apf_set_nan(r);
        return;
    }
    if (apf_is_zero(p)) {
        apf_set_inf(r, 1);  /* -Inf */
        return;
    }
    if (apf_cmp(p, &one) == 0) {
        apf_set_inf(r, 0);  /* +Inf */
        return;
    }
    
    /* Initial guess from normal approximation */
    apfx_norminv(&t, p);
    
    /* Newton-Raphson iteration */
    for (iter = 0; iter < 50; iter++) {
        apfx_tcdf(&ft, &t, df);
        apf_sub(&ft, &ft, p);  /* f(t) = tcdf(t) - p */
        
        /* PDF of t-distribution for derivative */
        {
            apf t2, df_apf, tmp;
            apf_mul(&t2, &t, &t);
            apf_from_int(&df_apf, df);
            apf_div(&tmp, &t2, &df_apf);
            apf_add(&tmp, &tmp, &one);
            
            /* pdf = (1 + t^2/df)^(-(df+1)/2) / (sqrt(df) * Beta(df/2, 0.5)) */
            {
                apf exp_val, two;
                apf_from_int(&exp_val, df + 1);
                apf_neg(&exp_val, &exp_val);
                apf_from_int(&two, 2);
                apf_div(&exp_val, &exp_val, &two);
                apfx_pow(&pdf_val, &tmp, &exp_val);
            }
            apf_sqrt(&tmp, &df_apf);
            {
                apf beta_val, a, b;
                apf_mul(&a, &df_apf, &half);
                apf_copy(&b, &half);
                apfx_beta(&beta_val, &a, &b);
                apf_mul(&tmp, &tmp, &beta_val);
            }
            apf_div(&pdf_val, &pdf_val, &tmp);
        }
        
        apf_div(&delta, &ft, &pdf_val);
        apf_sub(&t, &t, &delta);
        
        apf_abs(&delta, &delta);
        if (delta.exp < t.exp - AP_BITS + 10) break;
    }
    
    apf_copy(r, &t);
}

/* Chi-square CDF: P(X <= x) where X ~ chi2(df) */
void apfx_chi2cdf(apf *r, const apf *x, long df)
{
    apf half_x, half_df, half;
    
    if (x->cls == APF_CLASS_NAN || df < 1) {
        apf_set_nan(r);
        return;
    }
    if (x->sign) {
        apf_zero(r);  /* chi2 always positive */
        return;
    }
    if (apf_is_zero(x)) {
        apf_zero(r);
        return;
    }
    
    /* chi2 CDF = gammainc(df/2, x/2) */
    apf_from_str(&half, "0.5");
    apf_mul(&half_x, x, &half);
    {
        apf df_apf;
        apf_from_int(&df_apf, df);
        apf_mul(&half_df, &df_apf, &half);
    }
    
    apfx_gammainc(r, &half_df, &half_x);
}

/* Inverse chi-square using Newton-Raphson */
void apfx_chi2inv(apf *r, const apf *p, long df)
{
    apf x, fx, pdf_val, delta, one;
    int iter;
    
    if (p->cls == APF_CLASS_NAN || df < 1) {
        apf_set_nan(r);
        return;
    }
    
    apf_from_int(&one, 1);
    if (p->sign || apf_cmp(p, &one) > 0) {
        apf_set_nan(r);
        return;
    }
    if (apf_is_zero(p)) {
        apf_zero(r);
        return;
    }
    if (apf_cmp(p, &one) == 0) {
        apf_set_inf(r, 0);
        return;
    }
    
    /* Initial guess: x ~ df for moderate df */
    apf_from_int(&x, df);
    
    for (iter = 0; iter < 50; iter++) {
        apfx_chi2cdf(&fx, &x, df);
        apf_sub(&fx, &fx, p);
        
        /* PDF for derivative: x^(df/2-1) * e^(-x/2) / (2^(df/2) * Gamma(df/2)) */
        {
            apf half, df_apf, exp_val, tmp, coef;
            apf_from_str(&half, "0.5");
            apf_from_int(&df_apf, df);
            
            apf_mul(&exp_val, &df_apf, &half);
            apf_sub(&exp_val, &exp_val, &one);
            apfx_pow(&pdf_val, &x, &exp_val);
            
            apf_mul(&tmp, &x, &half);
            apf_neg(&tmp, &tmp);
            apfx_exp(&tmp, &tmp);
            apf_mul(&pdf_val, &pdf_val, &tmp);
            
            {
                apf two, gam_val;
                apf_from_int(&two, 2);
                apf_mul(&tmp, &df_apf, &half);
                apfx_pow(&coef, &two, &tmp);
                apfx_tgamma(&gam_val, &tmp);
                apf_mul(&coef, &coef, &gam_val);
            }
            apf_div(&pdf_val, &pdf_val, &coef);
        }
        
        apf_div(&delta, &fx, &pdf_val);
        apf_sub(&x, &x, &delta);
        
        /* Keep x positive */
        if (x.sign) {
            apf_from_str(&x, "0.001");
        }
        
        apf_abs(&delta, &delta);
        if (delta.exp < x.exp - AP_BITS + 10) break;
    }
    
    apf_copy(r, &x);
}

/* F-distribution CDF */
void apfx_fcdf(apf *r, const apf *x, long df1, long df2)
{
    apf tmp, d1_apf, d2_apf, half, a, b;
    
    if (x->cls == APF_CLASS_NAN || df1 < 1 || df2 < 1) {
        apf_set_nan(r);
        return;
    }
    if (x->sign || apf_is_zero(x)) {
        apf_zero(r);
        return;
    }
    
    /* F CDF = I_x(d1/2, d2/2) where x = d1*F/(d1*F + d2) */
    apf_from_int(&d1_apf, df1);
    apf_from_int(&d2_apf, df2);
    apf_from_str(&half, "0.5");
    
    apf_mul(&tmp, &d1_apf, x);
    {
        apf denom;
        apf_add(&denom, &tmp, &d2_apf);
        apf_div(&tmp, &tmp, &denom);  /* x = d1*F/(d1*F+d2) */
    }
    
    apf_mul(&a, &d1_apf, &half);
    apf_mul(&b, &d2_apf, &half);
    
    apfx_betainc(r, &tmp, &a, &b);
}

/* Binomial PMF: P(X = k) where X ~ Binom(n, p) */
void apfx_binompdf(apf *r, long k, long n, const apf *p)
{
    apf coef, pk, qnk, q, one;
    long i;
    
    if (p->cls == APF_CLASS_NAN || k < 0 || k > n || n < 0) {
        apf_set_nan(r);
        return;
    }
    
    apf_from_int(&one, 1);
    if (p->sign || apf_cmp(p, &one) > 0) {
        apf_set_nan(r);
        return;
    }
    
    /* Special cases */
    if (apf_is_zero(p)) {
        if (k == 0) apf_from_int(r, 1);
        else apf_zero(r);
        return;
    }
    if (apf_cmp(p, &one) == 0) {
        if (k == n) apf_from_int(r, 1);
        else apf_zero(r);
        return;
    }
    
    /* C(n,k) using multiplicative formula */
    apf_from_int(&coef, 1);
    for (i = 0; i < k; i++) {
        apf num, den;
        apf_from_int(&num, n - i);
        apf_from_int(&den, i + 1);
        apf_mul(&coef, &coef, &num);
        apf_div(&coef, &coef, &den);
    }
    
    /* p^k */
    {
        apf k_apf;
        apf_from_int(&k_apf, k);
        apfx_pow(&pk, p, &k_apf);
    }
    
    /* (1-p)^(n-k) */
    apf_sub(&q, &one, p);
    {
        apf nk_apf;
        apf_from_int(&nk_apf, n - k);
        apfx_pow(&qnk, &q, &nk_apf);
    }
    
    apf_mul(r, &coef, &pk);
    apf_mul(r, r, &qnk);
}

/* Binomial CDF: P(X <= k) */
void apfx_binomcdf(apf *r, long k, long n, const apf *p)
{
    apf sum, term;
    long i;
    
    if (p->cls == APF_CLASS_NAN || k < 0 || n < 0) {
        apf_set_nan(r);
        return;
    }
    if (k >= n) {
        apf_from_int(r, 1);
        return;
    }
    
    apf_zero(&sum);
    for (i = 0; i <= k; i++) {
        apfx_binompdf(&term, i, n, p);
        apf_add(&sum, &sum, &term);
    }
    
    apf_copy(r, &sum);
}

/* Poisson PMF: P(X = k) where X ~ Poisson(lambda) */
void apfx_poisspdf(apf *r, long k, const apf *lambda)
{
    apf term, exp_neg_lam;
    long i;
    
    if (lambda->cls == APF_CLASS_NAN || k < 0 || lambda->sign) {
        apf_set_nan(r);
        return;
    }
    if (apf_is_zero(lambda)) {
        if (k == 0) apf_from_int(r, 1);
        else apf_zero(r);
        return;
    }
    
    /* P(X=k) = lambda^k * e^(-lambda) / k! */
    {
        apf k_apf;
        apf_from_int(&k_apf, k);
        apfx_pow(&term, lambda, &k_apf);
    }
    
    {
        apf neg_lam;
        apf_neg(&neg_lam, lambda);
        apfx_exp(&exp_neg_lam, &neg_lam);
    }
    apf_mul(&term, &term, &exp_neg_lam);
    
    /* Divide by k! */
    for (i = 2; i <= k; i++) {
        apf div;
        apf_from_int(&div, i);
        apf_div(&term, &term, &div);
    }
    
    apf_copy(r, &term);
}

/* Poisson CDF: P(X <= k) */
void apfx_poisscdf(apf *r, long k, const apf *lambda)
{
    apf sum, term;
    long i;
    
    if (lambda->cls == APF_CLASS_NAN || k < 0) {
        apf_set_nan(r);
        return;
    }
    
    apf_zero(&sum);
    for (i = 0; i <= k; i++) {
        apfx_poisspdf(&term, i, lambda);
        apf_add(&sum, &sum, &term);
    }
    
    apf_copy(r, &sum);
}
#endif /* HAVE_DISTRIBUTIONS */

#ifdef HAVE_BESSEL
/* ================================================================
 * Additional Bessel Functions
 * ================================================================ */

/* Bessel Y_0(x) - Second kind, order 0 */
void apfx_y0(apf *r, const apf *x)
{
    apf j0_val, ln_x, euler, pi_val, term, sum;
    int k;
    
    if (x->cls == APF_CLASS_NAN) { apf_set_nan(r); return; }
    if (x->sign || apf_is_zero(x)) {
        apf_set_inf(r, 1);  /* Y0(x<=0) = -Inf */
        return;
    }
    
    /* For small x: Y0(x) = (2/pi)*(ln(x/2) + gamma)*J0(x) + series correction */
    /* For large x: asymptotic expansion */
    
    apfx_j0(&j0_val, x);
    apfx_log(&ln_x, x);
    {
        apf ln2;
        apf_from_str(&ln2, STR_LN2);
        apf_sub(&ln_x, &ln_x, &ln2);  /* ln(x/2) */
    }
    
    apf_from_str(&euler, STR_EULER);
    apf_add(&term, &ln_x, &euler);
    
    apfx_pi(&pi_val);
    {
        apf two_pi, coef;
        apf_from_int(&two_pi, 2);
        apf_div(&coef, &two_pi, &pi_val);
        apf_mul(&term, &term, &coef);
        apf_mul(r, &term, &j0_val);
    }
    
    /* Add series correction for better accuracy */
    /* Y0(x) = (2/pi)*[J0(x)*(ln(x/2)+gamma) + sum_{k=1}^inf (-1)^(k+1)*H_k*(x/2)^(2k)/(k!)^2] */
    {
        apf x2, x2_pow, hk, fact_sq, one, neg_one;
        apf_from_int(&one, 1);
        apf_from_int(&neg_one, -1);
        apf_mul(&x2, x, x);
        {
            apf four;
            apf_from_int(&four, 4);
            apf_div(&x2, &x2, &four);  /* (x/2)^2 */
        }
        apf_copy(&x2_pow, &x2);
        apf_from_int(&hk, 1);  /* H_1 = 1 */
        apf_from_int(&fact_sq, 1);
        
        apf_zero(&sum);
        for (k = 1; k <= 30; k++) {
            apf term_k, sign;
            
            /* Update factorial squared */
            {
                apf k_apf;
                apf_from_int(&k_apf, k);
                apf_mul(&fact_sq, &fact_sq, &k_apf);
                apf_mul(&fact_sq, &fact_sq, &k_apf);
            }
            
            /* sign = (-1)^(k+1) */
            if ((k + 1) % 2 == 0) apf_copy(&sign, &one);
            else apf_copy(&sign, &neg_one);
            
            apf_mul(&term_k, &sign, &hk);
            apf_mul(&term_k, &term_k, &x2_pow);
            apf_div(&term_k, &term_k, &fact_sq);
            apf_add(&sum, &sum, &term_k);
            
            /* Update H_k */
            {
                apf k1_apf;
                apf_from_int(&k1_apf, k + 1);
                apf_div(&term_k, &one, &k1_apf);
                apf_add(&hk, &hk, &term_k);
            }
            
            /* Update power */
            apf_mul(&x2_pow, &x2_pow, &x2);
            
            if (term_k.exp < sum.exp - AP_BITS) break;
        }
        
        {
            apf two_pi, coef;
            apf_from_int(&two_pi, 2);
            apf_div(&coef, &two_pi, &pi_val);
            apf_mul(&sum, &sum, &coef);
        }
        apf_add(r, r, &sum);
    }
}

/* Bessel Y_1(x) - Second kind, order 1 */
void apfx_y1(apf *r, const apf *x)
{
    apf j1_val, ln_x, term, pi_val;
    
    if (x->cls == APF_CLASS_NAN) { apf_set_nan(r); return; }
    if (x->sign || apf_is_zero(x)) {
        apf_set_inf(r, 1);
        return;
    }
    
    /* Y1(x) = (2/pi)*[J1(x)*ln(x/2) - 1/x + ...] */
    apfx_j1(&j1_val, x);
    apfx_log(&ln_x, x);
    {
        apf ln2;
        apf_from_str(&ln2, STR_LN2);
        apf_sub(&ln_x, &ln_x, &ln2);
    }
    
    apf_mul(&term, &j1_val, &ln_x);
    
    {
        apf inv_x;
        apf_from_int(&inv_x, 1);
        apf_div(&inv_x, &inv_x, x);
        apf_sub(&term, &term, &inv_x);
    }
    
    apfx_pi(&pi_val);
    {
        apf two_pi;
        apf_from_int(&two_pi, 2);
        apf_div(&two_pi, &two_pi, &pi_val);
        apf_mul(r, &term, &two_pi);
    }
}

/* Modified Bessel I_0(x) - First kind, order 0 */
void apfx_i0(apf *r, const apf *x)
{
    apf sum, term, x2, one;
    int k;
    
    if (x->cls == APF_CLASS_NAN) { apf_set_nan(r); return; }
    if (x->cls == APF_CLASS_INF) { apf_set_inf(r, 0); return; }
    if (apf_is_zero(x)) { apf_from_int(r, 1); return; }
    
    /* I_0(x) = sum_{k=0}^inf (x^2/4)^k / (k!)^2 */
    apf_from_int(&one, 1);
    apf_mul(&x2, x, x);
    {
        apf four;
        apf_from_int(&four, 4);
        apf_div(&x2, &x2, &four);  /* x^2/4 */
    }
    
    apf_from_int(&sum, 1);
    apf_copy(&term, &x2);
    
    for (k = 1; k <= 100; k++) {
        apf k_apf, k_sq;
        apf_from_int(&k_apf, k);
        apf_mul(&k_sq, &k_apf, &k_apf);
        apf_div(&term, &term, &k_sq);
        apf_add(&sum, &sum, &term);
        
        if (term.exp < sum.exp - AP_BITS) break;
        
        apf_mul(&term, &term, &x2);
    }
    
    apf_copy(r, &sum);
}

/* Modified Bessel I_1(x) - First kind, order 1 */
void apfx_i1(apf *r, const apf *x)
{
    apf sum, term, x2, half_x, one;
    int k;
    
    if (x->cls == APF_CLASS_NAN) { apf_set_nan(r); return; }
    if (x->cls == APF_CLASS_INF) { apf_set_inf(r, x->sign); return; }
    if (apf_is_zero(x)) { apf_zero(r); return; }
    
    /* I_1(x) = (x/2) * sum_{k=0}^inf (x^2/4)^k / (k!(k+1)!) */
    apf_from_int(&one, 1);
    {
        apf two;
        apf_from_int(&two, 2);
        apf_div(&half_x, x, &two);
    }
    
    apf_mul(&x2, x, x);
    {
        apf four;
        apf_from_int(&four, 4);
        apf_div(&x2, &x2, &four);
    }
    
    apf_from_int(&sum, 1);
    apf_copy(&term, &x2);
    
    for (k = 1; k <= 100; k++) {
        apf k_apf, k1_apf, denom;
        apf_from_int(&k_apf, k);
        apf_from_int(&k1_apf, k + 1);
        apf_mul(&denom, &k_apf, &k1_apf);
        apf_div(&term, &term, &denom);
        apf_add(&sum, &sum, &term);
        
        if (term.exp < sum.exp - AP_BITS) break;
        
        apf_mul(&term, &term, &x2);
    }
    
    apf_mul(r, &half_x, &sum);
}

/* Modified Bessel K_0(x) - Second kind, order 0 */
void apfx_k0(apf *r, const apf *x)
{
    apf i0_val, ln_x, euler, term;
    
    if (x->cls == APF_CLASS_NAN) { apf_set_nan(r); return; }
    if (x->sign || apf_is_zero(x)) { apf_set_inf(r, 0); return; }
    if (x->cls == APF_CLASS_INF) { apf_zero(r); return; }
    
    /* K_0(x) = -[ln(x/2) + gamma] * I_0(x) + series */
    apfx_i0(&i0_val, x);
    apfx_log(&ln_x, x);
    {
        apf ln2;
        apf_from_str(&ln2, STR_LN2);
        apf_sub(&ln_x, &ln_x, &ln2);
    }
    apf_from_str(&euler, STR_EULER);
    apf_add(&term, &ln_x, &euler);
    apf_neg(&term, &term);
    apf_mul(r, &term, &i0_val);
}

/* Modified Bessel K_1(x) - Second kind, order 1 */
void apfx_k1(apf *r, const apf *x)
{
    apf i1_val, ln_x, term, inv_x;
    
    if (x->cls == APF_CLASS_NAN) { apf_set_nan(r); return; }
    if (x->sign || apf_is_zero(x)) { apf_set_inf(r, 0); return; }
    if (x->cls == APF_CLASS_INF) { apf_zero(r); return; }
    
    /* K_1(x) ~ 1/x + x/2 * [...] for small x */
    apfx_i1(&i1_val, x);
    apfx_log(&ln_x, x);
    {
        apf ln2;
        apf_from_str(&ln2, STR_LN2);
        apf_sub(&ln_x, &ln_x, &ln2);
    }
    
    apf_from_int(&inv_x, 1);
    apf_div(&inv_x, &inv_x, x);
    
    /* Simplified: K_1(x) ~ 1/x for small x */
    apf_neg(&term, &ln_x);
    apf_mul(&term, &term, &i1_val);
    apf_add(r, &term, &inv_x);
}
#endif /* HAVE_BESSEL */

/* ========== Additional Special Functions ========== */

/* Digamma function psi(x) = d/dx ln(Gamma(x)) = Gamma'(x)/Gamma(x)
 * Uses asymptotic expansion for large x, reflection for negative x
 */
void apfx_digamma(apf *r, const apf *x)
{
    apf sum, term, x_work, one, eight, inv_x, inv_x2, x2;
    
    if (x->cls == APF_CLASS_NAN) { apf_set_nan(r); return; }
    if (x->cls == APF_CLASS_INF) { 
        if (x->sign) apf_set_nan(r);
        else apf_set_inf(r, 0);
        return;
    }
    
    apf_from_int(&one, 1);
    
    /* For negative x, use reflection: psi(1-x) - psi(x) = pi*cot(pi*x) */
    if (x->sign) {
        apf one_minus_x, psi_pos, pi_x, cot_val, pi;
        apf_sub(&one_minus_x, &one, x);
        apfx_digamma(&psi_pos, &one_minus_x);
        apfx_pi(&pi);
        apf_mul(&pi_x, &pi, x);
        apfx_cos(&term, &pi_x);
        apfx_sin(&cot_val, &pi_x);
        apf_div(&cot_val, &term, &cot_val);
        apf_mul(&cot_val, &cot_val, &pi);
        apf_sub(r, &psi_pos, &cot_val);
        return;
    }
    
    /* Use recurrence to shift x to large value */
    apf_copy(&x_work, x);
    apf_zero(&sum);
    apf_from_int(&eight, 8);
    while (apf_cmp(&x_work, &eight) < 0) {
        apf_from_int(&one, 1);
        apf_div(&term, &one, &x_work);
        apf_sub(&sum, &sum, &term);
        apf_add(&x_work, &x_work, &one);
    }
    
    /* Asymptotic expansion: psi(x) ~ ln(x) - 1/(2x) - sum B_{2k}/(2k*x^{2k}) */
    apfx_log(&term, &x_work);
    apf_add(&sum, &sum, &term);
    
    apf_from_int(&one, 1);
    apf_div(&inv_x, &one, &x_work);
    apf_mul(&inv_x2, &inv_x, &inv_x);
    
    /* -1/(2x) */
    {
        apf two;
        apf_from_int(&two, 2);
        apf_mul(&term, &two, &x_work);
        apf_from_int(&one, 1);
        apf_div(&term, &one, &term);
        apf_sub(&sum, &sum, &term);
    }
    
    /* Bernoulli number terms: B_2/2 = 1/12, B_4/4 = -1/120, B_6/6 = 1/252, ... */
    apf_copy(&x2, &inv_x2);
    
    /* 1/12 / x^2 */
    {
        apf b2;
        apf_from_str(&b2, "0.0833333333333333333333333333333333333");
        apf_mul(&term, &b2, &x2);
        apf_sub(&sum, &sum, &term);
    }
    apf_mul(&x2, &x2, &inv_x2);
    
    /* 1/120 / x^4 */
    {
        apf b4;
        apf_from_str(&b4, "0.0083333333333333333333333333333333333");
        apf_mul(&term, &b4, &x2);
        apf_add(&sum, &sum, &term);
    }
    apf_mul(&x2, &x2, &inv_x2);
    
    /* 1/252 / x^6 */
    {
        apf b6;
        apf_from_str(&b6, "0.0039682539682539682539682539682539683");
        apf_mul(&term, &b6, &x2);
        apf_sub(&sum, &sum, &term);
    }
    
    apf_copy(r, &sum);
}

/* Riemann zeta function for s > 1
 * Uses direct summation with acceleration for small s
 */
void apfx_zeta(apf *r, const apf *s)
{
    apf sum, term, one, n_apf, prev_sum;
    long n;
    int max_iter = 1000;
    
    if (s->cls == APF_CLASS_NAN) { apf_set_nan(r); return; }
    if (s->cls == APF_CLASS_INF) {
        apf_from_int(r, s->sign ? 0 : 1);
        return;
    }
    
    apf_from_int(&one, 1);
    
    /* Check s > 1 */
    if (apf_cmp(s, &one) <= 0) {
        /* For s = 1, pole (return inf) */
        if (apf_cmp(s, &one) == 0) {
            apf_set_inf(r, 0);
        } else {
            /* For s <= 0, would need analytic continuation - not implemented */
            apf_set_nan(r);
        }
        return;
    }
    
    /* Direct sum: zeta(s) = sum_{n=1}^inf 1/n^s */
    apf_zero(&sum);
    for (n = 1; n <= max_iter; n++) {
        apf_copy(&prev_sum, &sum);
        apf_from_int(&n_apf, n);
        apfx_pow(&term, &n_apf, s);
        apf_from_int(&one, 1);
        apf_div(&term, &one, &term);
        apf_add(&sum, &sum, &term);
        
        /* Convergence check */
        if (n > 10) {
            apf diff;
            apf_sub(&diff, &sum, &prev_sum);
            if (diff.exp < sum.exp - AP_BITS + 10) break;
        }
    }
    
    apf_copy(r, &sum);
}

/* Harmonic number H_n = 1 + 1/2 + 1/3 + ... + 1/n */
void apfx_harmonic(apf *r, long n)
{
    apf sum, term, one;
    long k;
    
    if (n <= 0) { apf_zero(r); return; }
    
    apf_zero(&sum);
    apf_from_int(&one, 1);
    
    for (k = 1; k <= n; k++) {
        apf k_apf;
        apf_from_int(&k_apf, k);
        apf_div(&term, &one, &k_apf);
        apf_add(&sum, &sum, &term);
    }
    
    apf_copy(r, &sum);
}

/* Generalized harmonic number H_{n,m} = sum_{k=1}^n 1/k^m */
void apfx_harmonic_gen(apf *r, long n, long m)
{
    apf sum, term, one;
    long k;
    
    if (n <= 0) { apf_zero(r); return; }
    if (m == 0) { apf_from_int(r, n); return; }
    
    apf_zero(&sum);
    apf_from_int(&one, 1);
    
    for (k = 1; k <= n; k++) {
        apf k_apf, k_pow, m_apf;
        apf_from_int(&k_apf, k);
        apf_from_int(&m_apf, m);
        apfx_pow(&k_pow, &k_apf, &m_apf);
        apf_div(&term, &one, &k_pow);
        apf_add(&sum, &sum, &term);
    }
    
    apf_copy(r, &sum);
}

/* Pochhammer symbol (rising factorial) (x)_n = x(x+1)(x+2)...(x+n-1) */
void apfx_pochhammer(apf *r, const apf *x, long n)
{
    apf prod, term, one;
    long k;
    
    if (n <= 0) { apf_from_int(r, 1); return; }
    
    apf_copy(&prod, x);
    apf_from_int(&one, 1);
    
    for (k = 1; k < n; k++) {
        apf x_plus_k;
        apf_from_int(&term, k);
        apf_add(&x_plus_k, x, &term);
        apf_mul(&prod, &prod, &x_plus_k);
    }
    
    apf_copy(r, &prod);
}

/* Falling factorial x^(n) = x(x-1)(x-2)...(x-n+1) */
void apfx_falling(apf *r, const apf *x, long n)
{
    apf prod, term, one;
    long k;
    
    if (n <= 0) { apf_from_int(r, 1); return; }
    
    apf_copy(&prod, x);
    apf_from_int(&one, 1);
    
    for (k = 1; k < n; k++) {
        apf x_minus_k;
        apf_from_int(&term, k);
        apf_sub(&x_minus_k, x, &term);
        apf_mul(&prod, &prod, &x_minus_k);
    }
    
    apf_copy(r, &prod);
}


