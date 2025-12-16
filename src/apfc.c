/* apfc.c - Arbitrary Precision Complex numbers
 * Pure C89, standalone library - no external dependencies
 */
#include "apfc.h"

/* Simple string helpers to avoid string.h dependency */
static void str_copy(char *dst, const char *src)
{
    while (*src) *dst++ = *src++;
    *dst = '\0';
}

static int str_eq(const char *a, const char *b)
{
    while (*a && *b && *a == *b) { a++; b++; }
    return *a == *b;
}

static char *str_append(char *dst, const char *src)
{
    while (*dst) dst++;
    while (*src) *dst++ = *src++;
    *dst = '\0';
    return dst;
}

void apfc_zero(apfc *z)
{
    apf_zero(&z->re);
    apf_zero(&z->im);
}

void apfc_from_real(apfc *z, const apf *r)
{
    apf_copy(&z->re, r);
    apf_zero(&z->im);
}

void apfc_copy(apfc *dest, const apfc *src)
{
    apf_copy(&dest->re, &src->re);
    apf_copy(&dest->im, &src->im);
}

int apfc_is_zero(const apfc *z)
{
    return apf_is_zero(&z->re) && apf_is_zero(&z->im);
}

int apfc_is_real(const apfc *z)
{
    return apf_is_zero(&z->im);
}

/* (a + bi) + (c + di) = (a+c) + (b+d)i */
void apfc_add(apfc *r, const apfc *a, const apfc *b)
{
    apf_add(&r->re, &a->re, &b->re);
    apf_add(&r->im, &a->im, &b->im);
}

/* (a + bi) - (c + di) = (a-c) + (b-d)i */
void apfc_sub(apfc *r, const apfc *a, const apfc *b)
{
    apf_sub(&r->re, &a->re, &b->re);
    apf_sub(&r->im, &a->im, &b->im);
}

/* (a + bi) * (c + di) = (ac - bd) + (ad + bc)i */
void apfc_mul(apfc *r, const apfc *a, const apfc *b)
{
    apf ac, bd, ad, bc;
    
    /* Special case: both purely real - avoid Inf*0 issues */
    if (apf_is_zero(&a->im) && apf_is_zero(&b->im)) {
        apf_mul(&r->re, &a->re, &b->re);
        apf_zero(&r->im);
        return;
    }
    
    /* Special case: both purely imaginary */
    if (apf_is_zero(&a->re) && apf_is_zero(&b->re)) {
        apf_mul(&r->re, &a->im, &b->im);
        apf_neg(&r->re, &r->re);  /* i*i = -1 */
        apf_zero(&r->im);
        return;
    }
    
    apf_mul(&ac, &a->re, &b->re);
    apf_mul(&bd, &a->im, &b->im);
    apf_mul(&ad, &a->re, &b->im);
    apf_mul(&bc, &a->im, &b->re);
    apf_sub(&r->re, &ac, &bd);
    apf_add(&r->im, &ad, &bc);
}

/* (a + bi) / (c + di) = ((ac + bd) + (bc - ad)i) / (c^2 + d^2) */
void apfc_div(apfc *r, const apfc *a, const apfc *b)
{
    apf c2, d2, denom, ac, bd, bc, ad;
    int a_re_nan = (a->re.cls == APF_CLASS_NAN);
    int a_im_nan = (a->im.cls == APF_CLASS_NAN);
    int b_re_nan = (b->re.cls == APF_CLASS_NAN);
    int b_im_nan = (b->im.cls == APF_CLASS_NAN);
    int b_zero = apf_is_zero(&b->re) && apf_is_zero(&b->im);
    int a_zero = apf_is_zero(&a->re) && apf_is_zero(&a->im);
    int b_re_inf = (b->re.cls == APF_CLASS_INF);
    int b_im_inf = (b->im.cls == APF_CLASS_INF);
    int a_re_inf = (a->re.cls == APF_CLASS_INF);
    int a_im_inf = (a->im.cls == APF_CLASS_INF);
    
    /* NaN propagation: if any input is NaN, result is NaN */
    if (a_re_nan || a_im_nan || b_re_nan || b_im_nan) {
        apf_set_nan(&r->re);
        apf_set_nan(&r->im);
        return;
    }
    
    /* Handle division by zero */
    if (b_zero) {
        if (a_zero) {
            /* 0/0 = NaN */
            apf_set_nan(&r->re);
            apf_set_nan(&r->im);
        } else {
            /* nonzero/0 = Inf (sign based on numerator) */
            apf_set_inf(&r->re, a->re.sign);
            /* Imaginary part: if purely real num, result is real Inf */
            if (apf_is_zero(&a->im)) {
                apf_zero(&r->im);
            } else {
                apf_set_inf(&r->im, a->im.sign);
            }
        }
        return;
    }
    
    /* Handle division by Inf: finite/Inf = 0, Inf/Inf = NaN */
    if (b_re_inf || b_im_inf) {
        if (a_re_inf || a_im_inf) {
            /* Inf/Inf = NaN */
            apf_set_nan(&r->re);
            apf_set_nan(&r->im);
        } else {
            /* finite/Inf = 0 */
            apf_zero(&r->re);
            apf_zero(&r->im);
        }
        return;
    }
    
    /* Handle Inf numerator with finite denominator: Inf/finite = Inf */
    if (a_re_inf || a_im_inf) {
        /* For purely real: Inf/x = Inf with sign(Inf)*sign(x) */
        if (apf_is_zero(&a->im) && apf_is_zero(&b->im)) {
            /* Real Inf / Real finite = Real Inf */
            int sign = a->re.sign ^ b->re.sign;
            apf_set_inf(&r->re, sign);
            apf_zero(&r->im);
        } else {
            /* Complex Inf case - use sign rules */
            if (a_re_inf) {
                apf_set_inf(&r->re, a->re.sign ^ b->re.sign);
            } else {
                apf_zero(&r->re);
            }
            if (a_im_inf) {
                apf_set_inf(&r->im, a->im.sign ^ b->re.sign);
            } else {
                apf_zero(&r->im);
            }
        }
        return;
    }
    
    apf_mul(&c2, &b->re, &b->re);
    apf_mul(&d2, &b->im, &b->im);
    apf_add(&denom, &c2, &d2);
    
    apf_mul(&ac, &a->re, &b->re);
    apf_mul(&bd, &a->im, &b->im);
    apf_mul(&bc, &a->im, &b->re);
    apf_mul(&ad, &a->re, &b->im);
    
    apf_add(&r->re, &ac, &bd);
    apf_div(&r->re, &r->re, &denom);
    apf_sub(&r->im, &bc, &ad);
    apf_div(&r->im, &r->im, &denom);
}

void apfc_neg(apfc *r, const apfc *a)
{
    apf_neg(&r->re, &a->re);
    apf_neg(&r->im, &a->im);
}

void apfc_conj(apfc *r, const apfc *a)
{
    apf_copy(&r->re, &a->re);
    apf_neg(&r->im, &a->im);
}

/* |z| = sqrt(re^2 + im^2) */
void apfc_abs(apf *r, const apfc *z)
{
    apf re2, im2, sum;
    apf_mul(&re2, &z->re, &z->re);
    apf_mul(&im2, &z->im, &z->im);
    apf_add(&sum, &re2, &im2);
    apf_sqrt(r, &sum);
}

/* arg(z) = atan2(im, re) */
void apfc_arg(apf *r, const apfc *z)
{
    apfx_atan2(r, &z->im, &z->re);
}

/* sqrt(z) for complex numbers */
void apfc_sqrt(apfc *r, const apfc *z)
{
    apf mod, two, re_part, im_part;
    int re_zero, im_zero;
    
    re_zero = apf_is_zero(&z->re);
    im_zero = apf_is_zero(&z->im);
    
    if (re_zero && im_zero) {
        apfc_zero(r);
        return;
    }
    
    /* For real numbers */
    if (im_zero) {
        if (!z->re.sign) {
            /* Positive real: sqrt is real */
            apf_sqrt(&r->re, &z->re);
            apf_zero(&r->im);
        } else {
            /* Negative real: sqrt is imaginary */
            apf neg_re;
            apf_neg(&neg_re, &z->re);
            apf_zero(&r->re);
            apf_sqrt(&r->im, &neg_re);
        }
        return;
    }
    
    /* General complex case: sqrt(a+bi) = sqrt((|z|+a)/2) + sign(b)*i*sqrt((|z|-a)/2) */
    apfc_abs(&mod, z);
    apf_from_int(&two, 2);
    
    {
        apf sum, diff;
        apf_add(&sum, &mod, &z->re);
        apf_div(&sum, &sum, &two);
        apf_sqrt(&re_part, &sum);
        
        apf_sub(&diff, &mod, &z->re);
        apf_div(&diff, &diff, &two);
        apf_sqrt(&im_part, &diff);
        
        if (z->im.sign) {
            apf_neg(&im_part, &im_part);
        }
    }
    
    apf_copy(&r->re, &re_part);
    apf_copy(&r->im, &im_part);
}

/* exp(a + bi) = exp(a) * (cos(b) + i*sin(b)) */
void apfc_exp(apfc *r, const apfc *z)
{
    apf ea, cos_b, sin_b;
    
    /* For purely real input, use real exp to avoid Inf*0 issue */
    if (apf_is_zero(&z->im)) {
        apfx_exp(&r->re, &z->re);
        apf_zero(&r->im);
        return;
    }
    
    apfx_exp(&ea, &z->re);
    apfx_sincos(&sin_b, &cos_b, &z->im);  /* Combined: one argument reduction */
    apf_mul(&r->re, &ea, &cos_b);
    apf_mul(&r->im, &ea, &sin_b);
}

/* log(z) = log|z| + i*arg(z) */
void apfc_log(apfc *r, const apfc *z)
{
    apf mod;
    apfc_abs(&mod, z);
    apfx_log(&r->re, &mod);
    apfc_arg(&r->im, z);
}

/* z^w = exp(w * log(z)) */
void apfc_pow(apfc *r, const apfc *base, const apfc *exp)
{
    apfc log_base, prod;
    apf one;
    
    apf_from_int(&one, 1);
    
    /* IEEE 754-2008: x^0 = 1 for ANY x, including NaN */
    if (apf_is_zero(&exp->re) && apf_is_zero(&exp->im)) {
        apf_from_int(&r->re, 1);
        apf_zero(&r->im);
        return;
    }
    
    /* IEEE 754-2008: 1^y = 1 for ANY y, including NaN and Inf
     * This applies to +1 only */
    if (apf_is_zero(&base->im)) {
        apf abs_base;
        apf_abs(&abs_base, &base->re);
        if (apf_cmp(&abs_base, &one) == 0 && !base->re.sign) {
            /* base is exactly +1 (real) */
            apf_from_int(&r->re, 1);
            apf_zero(&r->im);
            return;
        }
    }
    
    /* IEEE 754-2008: Handle Inf exponent BEFORE complex path
     * pow(x, ±Inf) depends on |x| vs 1 */
    if (apf_is_zero(&exp->im) && exp->re.cls == APF_CLASS_INF) {
        if (apf_is_zero(&base->im)) {
            /* Real base with Inf exponent */
            apf abs_base;
            int cmp;
            apf_abs(&abs_base, &base->re);
            cmp = apf_cmp(&abs_base, &one);
            
            if (cmp > 0) {
                /* |base| > 1: Inf -> Inf, -Inf -> 0 */
                if (exp->re.sign) {
                    apf_zero(&r->re);
                } else {
                    apf_set_inf(&r->re, 0);
                }
                apf_zero(&r->im);
            } else if (cmp < 0) {
                /* |base| < 1: Inf -> 0, -Inf -> Inf */
                if (exp->re.sign) {
                    apf_set_inf(&r->re, 0);
                } else {
                    apf_zero(&r->re);
                }
                apf_zero(&r->im);
            } else {
                /* |base| == 1 (including -1): result is 1 */
                apf_from_int(&r->re, 1);
                apf_zero(&r->im);
            }
            return;
        }
    }
    
    /* Special case: real base with real exponent */
    if (apf_is_zero(&base->im) && apf_is_zero(&exp->im)) {
        /* Check if base is negative */
        if (base->re.sign) {
            /* Negative real base - check if exponent is integer */
            long n = apf_to_long(&exp->re);
            apf n_apf, diff;
            apf_from_int(&n_apf, n);
            apf_sub(&diff, &exp->re, &n_apf);
            
            if (!apf_is_zero(&diff)) {
                /* Non-integer exponent with negative base: use complex path
                 * e.g., (-1)^0.5 = i, (-8)^(1/3) = 1 + 1.732i */
                apfc_log(&log_base, base);
                apfc_mul(&prod, exp, &log_base);
                apfc_exp(r, &prod);
                return;
            }
        }
        /* Positive base or integer exponent: use real pow */
        apfx_pow(&r->re, &base->re, &exp->re);
        apf_zero(&r->im);
        return;
    }
    
    /* Special case: integer real exponent - use binary exponentiation */
    if (apf_is_zero(&exp->im)) {
        long n = apf_to_long(&exp->re);
        apf n_apf, diff;
        apf_from_int(&n_apf, n);
        apf_sub(&diff, &exp->re, &n_apf);
        
        if (apf_is_zero(&diff) && n >= -10000000L && n <= 10000000L) {
            apfc result, b, temp;
            int neg_exp = 0;
            
            if (n < 0) {
                neg_exp = 1;
                n = -n;
            }
            
            apf_from_int(&result.re, 1);
            apf_zero(&result.im);
            apfc_copy(&b, base);
            
            /* Binary exponentiation: O(log n) multiplications */
            while (n > 0) {
                if (n & 1) {
                    apfc_mul(&temp, &result, &b);
                    result = temp;
                }
                apfc_mul(&temp, &b, &b);
                b = temp;
                n >>= 1;
            }
            
            if (neg_exp) {
                apfc one;
                apf_from_int(&one.re, 1);
                apf_zero(&one.im);
                apfc_div(r, &one, &result);
            } else {
                apfc_copy(r, &result);
            }
            return;
        }
    }
    
    apfc_log(&log_base, base);
    apfc_mul(&prod, exp, &log_base);
    apfc_exp(r, &prod);
}

/* sin(a + bi) = sin(a)cosh(b) + i*cos(a)sinh(b) */
void apfc_sin(apfc *r, const apfc *z)
{
    apf sin_a, cos_a, sinh_b, cosh_b;
    apfx_sincos(&sin_a, &cos_a, &z->re);  /* Combined sincos */
    apfx_sinhcosh(&sinh_b, &cosh_b, &z->im);  /* Combined sinhcosh */
    apf_mul(&r->re, &sin_a, &cosh_b);
    apf_mul(&r->im, &cos_a, &sinh_b);
}

/* cos(a + bi) = cos(a)cosh(b) - i*sin(a)sinh(b) */
void apfc_cos(apfc *r, const apfc *z)
{
    apf sin_a, cos_a, sinh_b, cosh_b, temp;
    apfx_sincos(&sin_a, &cos_a, &z->re);  /* Combined sincos */
    apfx_sinhcosh(&sinh_b, &cosh_b, &z->im);  /* Combined sinhcosh */
    apf_mul(&r->re, &cos_a, &cosh_b);
    apf_mul(&temp, &sin_a, &sinh_b);
    apf_neg(&r->im, &temp);
}

/* tan(z) = sin(z) / cos(z) */
void apfc_tan(apfc *r, const apfc *z)
{
    apfc s, c;
    apfc_sin(&s, z);
    apfc_cos(&c, z);
    apfc_div(r, &s, &c);
}

/* sinh(z) = (exp(z) - exp(-z)) / 2 */
void apfc_sinh(apfc *r, const apfc *z)
{
    apfc ez, emz, neg_z, temp, two_c;
    apfc_neg(&neg_z, z);
    apfc_exp(&ez, z);
    apfc_exp(&emz, &neg_z);
    apfc_sub(&temp, &ez, &emz);
    apf_from_int(&two_c.re, 2);
    apf_zero(&two_c.im);
    apfc_div(r, &temp, &two_c);
}

/* cosh(z) = (exp(z) + exp(-z)) / 2 */
void apfc_cosh(apfc *r, const apfc *z)
{
    apfc ez, emz, neg_z, temp, two_c;
    apfc_neg(&neg_z, z);
    apfc_exp(&ez, z);
    apfc_exp(&emz, &neg_z);
    apfc_add(&temp, &ez, &emz);
    apf_from_int(&two_c.re, 2);
    apf_zero(&two_c.im);
    apfc_div(r, &temp, &two_c);
}

/* String output */
void apfc_to_str(char *buf, int bufsize, const apfc *z, int max_frac)
{
    char re_buf[256], im_buf[256];
    int re_zero, im_zero;
    char *p;
    
    (void)bufsize;  /* silence warning; TODO: add proper bounds checking */
    
    /* Special case: if both are NaN, just show "NaN" */
    if (apf_isnan(&z->re) && apf_isnan(&z->im)) {
        str_copy(buf, "NaN");
        return;
    }
    
    re_zero = apf_is_zero(&z->re);
    im_zero = apf_is_zero(&z->im);
    
    if (re_zero && im_zero) {
        str_copy(buf, "0");
        return;
    }
    
    apf_to_str(re_buf, sizeof(re_buf), &z->re, max_frac);
    apf_to_str(im_buf, sizeof(im_buf), &z->im, max_frac);
    
    /* If imaginary is zero, just show real */
    if (im_zero) {
        str_copy(buf, re_buf);
        return;
    }
    
    /* If real is zero and imaginary is not */
    if (re_zero) {
        /* Pure imaginary */
        if (str_eq(im_buf, "1")) {
            str_copy(buf, "i");
        } else if (str_eq(im_buf, "-1")) {
            str_copy(buf, "-i");
        } else {
            /* im_buf + "i" */
            str_copy(buf, im_buf);
            str_append(buf, "i");
        }
        return;
    }
    
    /* Both parts: "re + im*i" or "re - im*i" */
    str_copy(buf, re_buf);
    if (z->im.sign) {
        p = str_append(buf, " - ");
        str_copy(p, im_buf + 1);  /* skip '-' */
        str_append(buf, "i");
    } else {
        p = str_append(buf, " + ");
        str_copy(p, im_buf);
        str_append(buf, "i");
    }
}
