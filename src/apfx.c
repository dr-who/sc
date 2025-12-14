/* apfx.c - Extended math functions for APF */
#include "apfx.h"
#include <string.h>

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

/* Cached ln(2) constant - 40 decimal digits */
static const char *ln2_str = "0.6931471805599453094172321214581765680755";
static apf cached_ln2;
static int ln2_initialized = 0;

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
    INT16 n, max_iter = 30;  /* Reduced from 100 - 60 is plenty for 128 bits */
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
     * For normalized APF with 128-bit mantissa: value ≈ 2^(stored_exp + 127)
     * Threshold: stored_exp >= -127 + 37 = -90 means |exp| >= 2^37 ≈ 10^11
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
static const char *pi_str = "3.1415926535897932384626433832795028841971";
static apf cached_pi;
static int pi_initialized = 0;

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
static const char *e_str = "2.7182818284590452353602874713526624977572";
static apf cached_e;
static int e_initialized = 0;

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
