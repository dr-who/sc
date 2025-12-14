/* format.c - Output formatting (decimal, hex, binary, fraction)
 * C89 compliant for Watcom C / DOS
 */
#include "sc.h"

/* ========== Fraction Output ========== */

/* Simple GCD for reducing fractions (only if mathx.c not providing it) */
#ifndef HAVE_GCD
static long gcd_long(long a, long b)
{
    if (a < 0) a = -a;
    if (b < 0) b = -b;
    while (b != 0) {
        long t = b;
        b = a % b;
        a = t;
    }
    return a;
}
#endif

/* Convert apf to best rational approximation using continued fractions */
static int apf_to_fraction(const apf *val, long *numer, long *denom, long max_denom)
{
    apf x, int_part, frac;
    long h0, h1, h2, k0, k1, k2;
    long a;
    int neg = 0;
    int iter;
    
    if (val->cls != APF_CLASS_NORMAL && val->cls != APF_CLASS_ZERO) {
        return 0;  /* Can't convert inf/nan */
    }
    
    if (val->cls == APF_CLASS_ZERO) {
        *numer = 0;
        *denom = 1;
        return 1;
    }
    
    /* Work with absolute value */
    apf_copy(&x, val);
    if (x.sign) {
        neg = 1;
        x.sign = 0;
    }
    
    /* Continued fraction convergents */
    h0 = 0; h1 = 1;
    k0 = 1; k1 = 0;
    
    for (iter = 0; iter < 50; iter++) {
        /* a = floor(x) - get integer part */
        a = apf_to_long(&x);
        
        /* Update convergents: h2/k2 = (a*h1 + h0) / (a*k1 + k0) */
        h2 = a * h1 + h0;
        k2 = a * k1 + k0;
        
        /* Check if denominator exceeded max */
        if (k2 > max_denom || k2 < 0) {  /* k2 < 0 means overflow */
            break;
        }
        
        h0 = h1; h1 = h2;
        k0 = k1; k1 = k2;
        
        /* frac = x - a */
        apf_from_int(&int_part, a);
        apf_sub(&frac, &x, &int_part);
        
        if (apf_is_zero(&frac)) {
            break;  /* Exact representation found */
        }
        
        /* x = 1 / frac */
        {
            apf one;
            apf_from_int(&one, 1);
            apf_div(&x, &one, &frac);
        }
        
        /* Check for overflow or very large intermediate */
        if (x.exp > 60) {  /* x > 2^60, getting too big */
            break;
        }
    }
    
    *numer = neg ? -h1 : h1;
    *denom = k1;
    
    /* Reduce by GCD just in case */
    {
        long g = gcd_long(h1, k1);
        if (g > 1) {
            *numer = (neg ? -h1 : h1) / g;
            *denom = k1 / g;
        }
    }
    
    return 1;
}

/* Print a value as a fraction if possible */
static void print_fraction(const apf *val)
{
    long numer, denom;
    
    if (apf_to_fraction(val, &numer, &denom, 1000000L)) {
        if (denom == 1) {
            printf("%ld", numer);
        } else {
            printf("%ld/%ld", numer, denom);
        }
    } else {
        /* Fall back to decimal */
        static char buf[256];
        apf_to_str(buf, sizeof(buf), val, display_digits);
        printf("%s", buf);
    }
}

/* ========== Result Output ========== */

void print_result(const apfc *val)
{
    static char buf[512];
    
    if (current_mode == MODE_FRACTION && apfc_is_real(val)) {
        printf("= ");
        print_fraction(&val->re);
        printf("\n");
        return;
    }
    
    apfc_to_str(buf, sizeof(buf), val, display_digits);
    printf("= %s\n", buf);
}

void print_value(const value_t *val)
{
    if (val->type == VAL_SCALAR) {
        print_result(&val->v.scalar);
    } else {
        printf("=\n");
        mat_print(&val->v.matrix);
    }
}

void value_to_scalar(apfc *r, const value_t *val)
{
    if (val->type == VAL_SCALAR) {
        *r = val->v.scalar;
    } else if (val->v.matrix.rows == 1 && val->v.matrix.cols == 1) {
        *r = MAT_AT(&val->v.matrix, 0, 0);
    } else {
        printf("Error: cannot convert matrix to scalar\n");
        apf_zero(&r->re);
        apf_zero(&r->im);
    }
}

void value_from_scalar(value_t *r, const apfc *val)
{
    r->type = VAL_SCALAR;
    r->v.scalar = *val;
}

void value_from_matrix(value_t *r, const matrix_t *m)
{
    r->type = VAL_MATRIX;
    mat_copy(&r->v.matrix, m);
}

/* ========== Hex Output ========== */

void apf_to_hex_str(const apf *val)
{
    static char hbuf[256];
    char *p;
    long lval;
    
    if (val->cls == APF_CLASS_ZERO) {
        printf("0x0\n");
        return;
    }
    if (val->cls == APF_CLASS_INF) {
        printf("%sInf\n", val->sign ? "-" : "");
        return;
    }
    if (val->cls == APF_CLASS_NAN) {
        printf("NaN\n");
        return;
    }
    
    /* For numbers that fit in long, use simple conversion */
    lval = apf_to_long(val);
    
    /* Check if number is too large for exact representation */
    if (val->exp > 30) {
        printf("~2^%ld (too large for exact hex)\n", val->exp);
        return;
    }
    
    p = hbuf;
    if (lval < 0) {
        *p++ = '-';
        lval = -lval;
    }
    sprintf(p, "0x%lX", lval);
    printf("%s\n", hbuf);
}

/* ========== Binary Output ========== */

void apf_to_bin_str(const apf *val)
{
    static char bbuf[256];
    char *p;
    long lval;
    int i, started = 0;
    
    if (val->cls == APF_CLASS_ZERO) {
        printf("0b0\n");
        return;
    }
    if (val->cls == APF_CLASS_INF) {
        printf("%sInf\n", val->sign ? "-" : "");
        return;
    }
    if (val->cls == APF_CLASS_NAN) {
        printf("NaN\n");
        return;
    }
    
    /* For numbers that fit in long, use simple conversion */
    lval = apf_to_long(val);
    
    /* Check if number is too large for exact representation */
    if (val->exp > 30) {
        printf("~2^%ld (too large for exact binary)\n", val->exp);
        return;
    }
    
    p = bbuf;
    if (lval < 0) {
        *p++ = '-';
        lval = -lval;
    }
    *p++ = '0';
    *p++ = 'b';
    
    /* Print bits from MSB to LSB */
    for (i = 30; i >= 0; i--) {
        int bit = (lval >> i) & 1;
        if (bit || started) {
            *p++ = '0' + bit;
            started = 1;
        }
    }
    if (!started) {
        *p++ = '0';
    }
    *p = '\0';
    printf("%s\n", bbuf);
}
