/* format.c - Output formatting (decimal, hex, binary, fraction)
 * C89 compliant for Watcom C / DOS
 */
#include <stdio.h>
#include <string.h>
#include "sc.h"

/* Forward declarations */
static void format_complex_mode(char *out, int outsize, const apfc *val,
                                calc_mode_t mode, int fixed_places);

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
    
    switch (current_mode) {
        case MODE_FRACTION:
            if (apfc_is_real(val)) {
                printf("= ");
                print_fraction(&val->re);
                printf("\n");
                return;
            }
            break;
        case MODE_HEX:
            if (apfc_is_real(val)) {
                printf("= ");
                apf_to_hex_str(&val->re);
                return;
            }
            break;
        case MODE_BIN:
            if (apfc_is_real(val)) {
                printf("= ");
                apf_to_bin_str(&val->re);
                return;
            }
            break;
        case MODE_IEEE:
            apfc_ieee_display(val);
            return;
        case MODE_FIX:
        case MODE_SCI:
        case MODE_ENG:
            format_complex_mode(buf, sizeof(buf), val, current_mode, display_fixed_places);
            printf("= %s\n", buf);
            return;
        default:
            break;
    }
    
    apfc_to_str(buf, sizeof(buf), val, display_digits);
    printf("= %s\n", buf);
}

void print_value_ex(const value_t *val, int quiet)
{
    static char buf[256];
    if (val->type == VAL_SCALAR) {
        switch (current_mode) {
            case MODE_HEX:
                if (!quiet) printf("= ");
                apf_to_hex_str(&val->v.scalar.re);
                return;
            case MODE_BIN:
                if (!quiet) printf("= ");
                apf_to_bin_str(&val->v.scalar.re);
                return;
            case MODE_IEEE:
                apfc_ieee_display(&val->v.scalar);
                return;
            case MODE_FRACTION:
                if (apfc_is_real(&val->v.scalar)) {
                    if (!quiet) printf("= ");
                    print_fraction(&val->v.scalar.re);
                    printf("\n");
                    return;
                }
                /* Fall through to decimal for complex */
                break;
            case MODE_FIX:
            case MODE_SCI:
            case MODE_ENG:
                format_complex_mode(buf, sizeof(buf), &val->v.scalar, 
                                   current_mode, display_fixed_places);
                if (quiet) {
                    printf("%s\n", buf);
                } else {
                    printf("= %s\n", buf);
                }
                return;
            default:
                break;
        }
        apfc_to_str(buf, sizeof(buf), &val->v.scalar, display_digits);
        if (quiet) {
            printf("%s\n", buf);
        } else {
            printf("= %s\n", buf);
        }
    } else {
        if (!quiet) printf("=\n");
        mat_print(&val->v.matrix);
    }
}

void print_value(const value_t *val)
{
    print_value_ex(val, 0);
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

/* ========== IEEE Educational Display ========== */

/* Helper: print binary with grouping */
static void print_binary_grouped(unsigned long val, int bits, int group_size)
{
    int i;
    for (i = bits - 1; i >= 0; i--) {
        printf("%d", (int)((val >> i) & 1));
        if (i > 0 && i % group_size == 0) {
            printf(" ");
        }
    }
}

/* Helper: print mantissa bits from limbs (skip implicit leading 1) */
static void print_mantissa_binary(const apf *val, int max_bits)
{
    int bit_count = 0;
    int limb_idx, bit_idx;
    int first_limb = 1;
    
    /* Mantissa is stored with MSB in mant[AP_LIMBS-1], bit 15 */
    /* The MSB is the implicit leading 1, so skip it */
    for (limb_idx = AP_LIMBS - 1; limb_idx >= 0 && bit_count < max_bits; limb_idx--) {
        int start_bit = (first_limb) ? 14 : 15;  /* Skip bit 15 on first limb (implicit 1) */
        first_limb = 0;
        for (bit_idx = start_bit; bit_idx >= 0 && bit_count < max_bits; bit_idx--) {
            int bit = (val->mant[limb_idx] >> bit_idx) & 1;
            printf("%d", bit);
            bit_count++;
            if (bit_count % 4 == 0 && bit_count < max_bits) {
                printf(" ");
            }
        }
    }
}

/* Get mantissa as apf value (1.xxxxx in binary) */
static void get_mantissa_apf(apf *result, const apf *val)
{
    apf frac, bit_val;
    int limb_idx, bit_idx;
    int bit_count = 0;
    
    /* Start with implicit leading 1 */
    apf_from_int(result, 1);
    apf_from_str(&frac, "0.5");
    
    /* Add fractional bits after the leading 1 */
    for (limb_idx = AP_LIMBS - 1; limb_idx >= 0 && bit_count < AP_BITS; limb_idx--) {
        int start_bit = (limb_idx == AP_LIMBS - 1) ? 14 : 15;  /* Skip MSB on first limb */
        for (bit_idx = start_bit; bit_idx >= 0 && bit_count < AP_BITS; bit_idx--) {
            int bit = (val->mant[limb_idx] >> bit_idx) & 1;
            if (bit) {
                apf_add(result, result, &frac);
            }
            /* frac = frac / 2 */
            apf_from_int(&bit_val, 2);
            apf_div(&frac, &frac, &bit_val);
            bit_count++;
        }
    }
}

/* Compute base^exp where exp is a long integer using binary exponentiation */
static void apf_pow_long(apf *r, const apf *base, long exp)
{
    apf result, b, temp;
    int neg_exp = 0;
    
    if (exp == 0) {
        apf_from_int(r, 1);
        return;
    }
    
    if (exp < 0) {
        neg_exp = 1;
        exp = -exp;
    }
    
    apf_from_int(&result, 1);
    apf_copy(&b, base);
    
    while (exp > 0) {
        if (exp & 1) {
            apf_mul(&temp, &result, &b);
            apf_copy(&result, &temp);
        }
        apf_mul(&temp, &b, &b);
        apf_copy(&b, &temp);
        exp >>= 1;
    }
    
    if (neg_exp) {
        apf one;
        apf_from_int(&one, 1);
        apf_div(r, &one, &result);
    } else {
        apf_copy(r, &result);
    }
}

/* Educational IEEE-style display */
void apf_ieee_display(const apf *val)
{
    apf mantissa, two, power, evaluated, neg_one;
    char mantissa_str[128];
    char power_str[128];
    char eval_str[128];
    char orig_str[128];
    char expr_buf[512];
    long exp_val;
    long actual_exp;
    int is_correct;
    
    printf("\n");
    printf("══════════════════════════════════════════════════════════════════\n");
    printf("           FLOATING POINT BREAKDOWN: How the number is stored\n");
    printf("══════════════════════════════════════════════════════════════════\n\n");
    
    /* Handle special cases */
    if (val->cls == APF_CLASS_ZERO) {
        printf("CLASS: ZERO\n\n");
        printf("SIGN\n");
        printf("  %d (%s)\n\n", val->sign, val->sign ? "-" : "+");
        printf("VALUE: 0\n");
        printf("\n══════════════════════════════════════════════════════════════════\n");
        return;
    }
    
    if (val->cls == APF_CLASS_INF) {
        printf("CLASS: INFINITY\n\n");
        printf("SIGN\n");
        printf("  %d (%s)\n\n", val->sign, val->sign ? "-" : "+");
        printf("VALUE: %sInfinity\n", val->sign ? "-" : "+");
        printf("\n══════════════════════════════════════════════════════════════════\n");
        return;
    }
    
    if (val->cls == APF_CLASS_NAN) {
        printf("CLASS: NaN (Not a Number)\n\n");
        printf("VALUE: NaN\n");
        printf("\n══════════════════════════════════════════════════════════════════\n");
        return;
    }
    
    /* Normal number */
    exp_val = val->exp;
    actual_exp = exp_val + (AP_BITS - 1);  /* Convert from internal to actual exponent */
    
    printf("CLASS: NORMAL\n\n");
    
    /* Sign bit */
    printf("──────────────────────────────────────────────────────────────────\n");
    printf("SIGN BIT\n");
    printf("──────────────────────────────────────────────────────────────────\n");
    printf("  %d (%s)\n\n", val->sign, val->sign ? "negative" : "positive");
    
    /* Exponent */
    printf("──────────────────────────────────────────────────────────────────\n");
    printf("EXPONENT\n");
    printf("──────────────────────────────────────────────────────────────────\n");
    printf("  ");
    if (actual_exp >= 0 && actual_exp < 256) {
        print_binary_grouped((unsigned long)actual_exp, 8, 4);
    } else if (actual_exp < 0 && actual_exp > -256) {
        printf("-");
        print_binary_grouped((unsigned long)(-actual_exp), 8, 4);
    } else {
        /* Large exponent - show more bits */
        if (actual_exp >= 0) {
            print_binary_grouped((unsigned long)actual_exp, 16, 4);
        } else {
            printf("-");
            print_binary_grouped((unsigned long)(-actual_exp), 16, 4);
        }
    }
    printf(" (%ld)\n\n", actual_exp);
    
    /* Mantissa - compute using apf */
    get_mantissa_apf(&mantissa, val);
    apf_to_str(mantissa_str, sizeof(mantissa_str), &mantissa, 15);
    
    printf("──────────────────────────────────────────────────────────────────\n");
    printf("MANTISSA (with implicit leading 1)\n");
    printf("──────────────────────────────────────────────────────────────────\n");
    printf("  1.");
    print_mantissa_binary(val, 24);
    printf("... (%s)\n\n", mantissa_str);
    
    /* Compute 2^exponent using apf */
    apf_from_int(&two, 2);
    apf_pow_long(&power, &two, actual_exp);
    apf_to_str(power_str, sizeof(power_str), &power, 15);
    
    /* Show the formula and calculation */
    printf("──────────────────────────────────────────────────────────────────\n");
    printf("CALCULATION\n");
    printf("──────────────────────────────────────────────────────────────────\n");
    printf("  Formula: (-1)^sign * mantissa * 2^exponent\n\n");
    
    printf("  Step 1: (-1)^%d = %s\n", val->sign, val->sign ? "-1" : "+1");
    printf("  Step 2: mantissa = %s\n", mantissa_str);
    printf("  Step 3: 2^%ld = %s\n", actual_exp, power_str);
    
    /* Build combining expression */
    sprintf(expr_buf, "%s * %s * %s", 
            val->sign ? "-1" : "1", mantissa_str, power_str);
    printf("\n  Combining: %s\n", expr_buf);
    
    /* Evaluate: sign * mantissa * 2^exp using apf */
    apf_mul(&evaluated, &mantissa, &power);
    if (val->sign) {
        apf_from_int(&neg_one, -1);
        apf_mul(&evaluated, &evaluated, &neg_one);
    }
    
    /* Get strings for comparison */
    apf_to_str(eval_str, sizeof(eval_str), &evaluated, 15);
    apf_to_str(orig_str, sizeof(orig_str), val, 15);
    
    /* Check if evaluation matches original */
    is_correct = apf_eq(&evaluated, val);
    
    printf("\n");
    printf("══════════════════════════════════════════════════════════════════\n");
    printf("  Evaluated: %s (%s)\n", eval_str, is_correct ? "correct" : "MISMATCH");
    printf("══════════════════════════════════════════════════════════════════\n");
}

/* Complex number IEEE display */
void apfc_ieee_display(const apfc *val)
{
    if (!apf_is_zero(&val->im)) {
        printf("\n*** REAL PART ***\n");
        apf_ieee_display(&val->re);
        printf("\n*** IMAGINARY PART ***\n");
        apf_ieee_display(&val->im);
    } else {
        apf_ieee_display(&val->re);
    }
}

/* ========== FIX/SCI/ENG Display Modes ========== */

/* Format a number in fixed-point notation with N decimal places */
static void format_fix(char *out, int outsize, const apf *val, int places)
{
    char buf[256];
    char *p, *dot;
    int i;
    
    if (val->cls == APF_CLASS_NAN) { strcpy(out, "NaN"); return; }
    if (val->cls == APF_CLASS_INF) {
        strcpy(out, val->sign ? "-Inf" : "Inf");
        return;
    }
    if (val->cls == APF_CLASS_ZERO) {
        strcpy(out, "0");
        if (places > 0) {
            strcat(out, ".");
            for (i = 0; i < places && i < outsize - 10; i++) strcat(out, "0");
        }
        return;
    }
    
    /* Get full precision string */
    apf_to_str(buf, sizeof(buf), val, places + 20);
    
    /* Find decimal point */
    p = out;
    dot = strchr(buf, '.');
    
    if (buf[0] == '-') {
        *p++ = '-';
        memmove(buf, buf + 1, strlen(buf));
        if (dot) dot--;
    }
    
    /* Handle scientific notation in input */
    {
        char *e = strchr(buf, 'e');
        if (e) {
            /* Convert from scientific notation */
            int exp_val = atoi(e + 1);
            *e = '\0';
            dot = strchr(buf, '.');
            
            if (exp_val >= 0) {
                /* Move decimal right */
                char mantissa[256];
                int len;
                strcpy(mantissa, buf);
                if (dot) {
                    memmove(dot, dot + 1, strlen(dot));
                }
                len = (int)strlen(mantissa);
                /* Pad with zeros if needed */
                while (len < exp_val + 1) {
                    mantissa[len++] = '0';
                }
                mantissa[len] = '\0';
                /* Insert new decimal point */
                if (places > 0) {
                    int insert_pos = exp_val + 1;
                    memmove(mantissa + insert_pos + 1, mantissa + insert_pos, len - insert_pos + 1);
                    mantissa[insert_pos] = '.';
                }
                strcpy(buf, mantissa);
            } else {
                /* Move decimal left - prepend zeros */
                char mantissa[256];
                int zeros = -exp_val - 1;
                int j;
                strcpy(mantissa, "0.");
                for (j = 0; j < zeros; j++) strcat(mantissa, "0");
                if (dot) memmove(dot, dot + 1, strlen(dot));
                strcat(mantissa, buf);
                strcpy(buf, mantissa);
            }
        }
    }
    
    /* Now buf has the number, format to places decimal digits */
    dot = strchr(buf, '.');
    if (dot) {
        int frac_len = (int)strlen(dot + 1);
        if (frac_len > places) {
            /* Round */
            dot[places + 1] = '\0';
        } else {
            /* Pad with zeros */
            while (frac_len < places) {
                strcat(buf, "0");
                frac_len++;
            }
        }
    } else if (places > 0) {
        strcat(buf, ".");
        for (i = 0; i < places; i++) strcat(buf, "0");
    }
    
    strcpy(p, buf);
    (void)outsize;
}

/* Format a number in scientific notation with N decimal places */
static void format_sci(char *out, int outsize, const apf *val, int places)
{
    char buf[256], mantissa[64];
    int exp_val = 0;
    int i, neg = 0;
    char *e, *dot;
    
    if (val->cls == APF_CLASS_NAN) { strcpy(out, "NaN"); return; }
    if (val->cls == APF_CLASS_INF) {
        strcpy(out, val->sign ? "-Inf" : "Inf");
        return;
    }
    if (val->cls == APF_CLASS_ZERO) {
        strcpy(out, "0");
        if (places > 0) {
            strcat(out, ".");
            for (i = 0; i < places && i < outsize - 10; i++) strcat(out, "0");
        }
        strcat(out, "e+0");
        return;
    }
    
    /* Get string representation */
    apf_to_str(buf, sizeof(buf), val, places + 5);
    
    if (buf[0] == '-') {
        neg = 1;
        memmove(buf, buf + 1, strlen(buf));
    }
    
    /* Already in scientific notation? */
    e = strchr(buf, 'e');
    if (e) {
        exp_val = atoi(e + 1);
        *e = '\0';
    }
    
    /* Get mantissa digits */
    dot = strchr(buf, '.');
    if (dot) {
        /* Normalize: first digit before decimal */
        int first_nonzero = -1;
        char *p = buf;
        
        /* Find first non-zero digit */
        while (*p) {
            if (*p >= '1' && *p <= '9') {
                first_nonzero = (int)(p - buf);
                break;
            }
            p++;
        }
        
        if (first_nonzero >= 0) {
            int shift;
            dot = strchr(buf, '.');
            if (first_nonzero < (dot - buf)) {
                shift = first_nonzero;
            } else {
                shift = first_nonzero - 1;  /* Account for decimal */
            }
            
            if (shift != 0 || e == NULL) {
                /* Rebuild mantissa */
                char digits[128];
                int j = 0;
                for (i = 0; buf[i]; i++) {
                    if (buf[i] != '.') digits[j++] = buf[i];
                }
                digits[j] = '\0';
                
                /* Skip leading zeros */
                p = digits;
                while (*p == '0') { p++; exp_val--; }
                
                if (*p == '\0') {
                    /* All zeros */
                    strcpy(buf, "0");
                    exp_val = 0;
                } else {
                    /* Format as D.DDDD */
                    mantissa[0] = *p++;
                    mantissa[1] = '.';
                    strncpy(mantissa + 2, p, places);
                    mantissa[2 + places] = '\0';
                    /* Pad if needed */
                    while ((int)strlen(mantissa) < places + 2) strcat(mantissa, "0");
                    strcpy(buf, mantissa);
                }
            }
        }
    } else {
        /* Integer - convert to scientific */
        int len = (int)strlen(buf);
        if (len > 1) {
            mantissa[0] = buf[0];
            mantissa[1] = '.';
            strncpy(mantissa + 2, buf + 1, places);
            mantissa[2 + places] = '\0';
            while ((int)strlen(mantissa) < places + 2) strcat(mantissa, "0");
            exp_val += len - 1;
            strcpy(buf, mantissa);
        }
    }
    
    /* Build output */
    out[0] = '\0';
    if (neg) strcat(out, "-");
    strcat(out, buf);
    sprintf(out + strlen(out), "e%+d", exp_val);
    (void)outsize;
}

/* Format a number in engineering notation (exponent multiple of 3) */
static void format_eng(char *out, int outsize, const apf *val, int places)
{
    char buf[256];
    int exp_val = 0;
    int i, neg = 0, shift;
    char *e;
    
    if (val->cls == APF_CLASS_NAN) { strcpy(out, "NaN"); return; }
    if (val->cls == APF_CLASS_INF) {
        strcpy(out, val->sign ? "-Inf" : "Inf");
        return;
    }
    if (val->cls == APF_CLASS_ZERO) {
        strcpy(out, "0");
        if (places > 0) {
            strcat(out, ".");
            for (i = 0; i < places && i < outsize - 10; i++) strcat(out, "0");
        }
        strcat(out, "e+0");
        return;
    }
    
    /* First get scientific notation */
    format_sci(buf, sizeof(buf), val, places + 3);
    
    if (buf[0] == '-') {
        neg = 1;
        memmove(buf, buf + 1, strlen(buf));
    }
    
    e = strchr(buf, 'e');
    if (e) {
        exp_val = atoi(e + 1);
        *e = '\0';
    }
    
    /* Adjust to make exponent multiple of 3 */
    shift = exp_val % 3;
    if (shift < 0) shift += 3;
    
    if (shift != 0) {
        /* Move decimal point right by 'shift' positions */
        char *dot = strchr(buf, '.');
        if (dot) {
            char digits[128];
            int j = 0;
            for (i = 0; buf[i]; i++) {
                if (buf[i] != '.') digits[j++] = buf[i];
            }
            digits[j] = '\0';
            
            /* Insert decimal at new position */
            memmove(digits + shift + 2, digits + shift + 1, strlen(digits + shift + 1) + 1);
            memmove(digits + 1, digits, shift + 1);
            digits[shift + 1] = '.';
            /* Ensure proper length */
            digits[shift + 2 + places] = '\0';
            strcpy(buf, digits);
            
            exp_val -= shift;
        }
    }
    
    /* Ensure decimal places */
    {
        char *dot = strchr(buf, '.');
        if (dot) {
            int frac_len = (int)strlen(dot + 1);
            while (frac_len < places) {
                strcat(buf, "0");
                frac_len++;
            }
            dot[places + 1] = '\0';
        }
    }
    
    out[0] = '\0';
    if (neg) strcat(out, "-");
    strcat(out, buf);
    sprintf(out + strlen(out), "e%+d", exp_val);
    (void)outsize;
}

/* Format complex number in FIX/SCI/ENG mode */
static void format_complex_mode(char *out, int outsize, const apfc *val, 
                                calc_mode_t mode, int places)
{
    char re_buf[128], im_buf[128];
    int re_zero = apf_is_zero(&val->re);
    int im_zero = apf_is_zero(&val->im);
    
    if (re_zero && im_zero) {
        switch (mode) {
            case MODE_FIX: format_fix(out, outsize, &val->re, places); break;
            case MODE_SCI: format_sci(out, outsize, &val->re, places); break;
            case MODE_ENG: format_eng(out, outsize, &val->re, places); break;
            default: strcpy(out, "0"); break;
        }
        return;
    }
    
    switch (mode) {
        case MODE_FIX:
            format_fix(re_buf, sizeof(re_buf), &val->re, places);
            format_fix(im_buf, sizeof(im_buf), &val->im, places);
            break;
        case MODE_SCI:
            format_sci(re_buf, sizeof(re_buf), &val->re, places);
            format_sci(im_buf, sizeof(im_buf), &val->im, places);
            break;
        case MODE_ENG:
            format_eng(re_buf, sizeof(re_buf), &val->re, places);
            format_eng(im_buf, sizeof(im_buf), &val->im, places);
            break;
        default:
            apfc_to_str(out, outsize, val, places);
            return;
    }
    
    if (im_zero) {
        strcpy(out, re_buf);
    } else if (re_zero) {
        sprintf(out, "%si", im_buf);
    } else {
        if (val->im.sign) {
            sprintf(out, "%s - %si", re_buf, im_buf + 1);  /* Skip minus sign */
        } else {
            sprintf(out, "%s + %si", re_buf, im_buf);
        }
    }
}
