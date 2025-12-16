/* mathx.c - Extended math functions
 * C89 portable for DOS, Linux, VIC-20
 */

#include "config.h"
#include "mathx.h"
#include "apf.h"
#include "apfx.h"
#include "sc.h"  /* For angle_mode */

/* ========== Combinatorics ========== */
#ifdef HAVE_COMB

/* nPr = n! / (n-r)! = n * (n-1) * ... * (n-r+1) */
void apf_npr(apf *r, long n, long k) {
    long i;
    apf tmp, mul;
    
    if (k < 0 || n < 0 || k > n) {
        apf_zero(r);
        return;
    }
    if (k == 0) {
        apf_from_int(r, 1);
        return;
    }
    
    apf_from_int(r, n);
    for (i = 1; i < k; i++) {
        apf_from_int(&mul, n - i);
        apf_mul(&tmp, r, &mul);
        apf_copy(r, &tmp);
    }
}

/* nCr = n! / (r! * (n-r)!) = nPr / r! */
void apf_ncr(apf *r, long n, long k) {
    long i;
    apf tmp, div;
    
    if (k < 0 || n < 0 || k > n) {
        apf_zero(r);
        return;
    }
    
    /* Use symmetry: C(n,k) = C(n, n-k) */
    if (k > n - k) {
        k = n - k;
    }
    
    if (k == 0) {
        apf_from_int(r, 1);
        return;
    }
    
    /* Calculate n * (n-1) * ... * (n-k+1) / k! iteratively */
    apf_from_int(r, n);
    for (i = 1; i < k; i++) {
        apf_from_int(&tmp, n - i);
        apf_mul(&div, r, &tmp);
        apf_from_int(&tmp, i + 1);
        apf_div(r, &div, &tmp);
    }
}

#endif /* HAVE_COMB */

/* ========== Number Theory ========== */
#ifdef HAVE_GCD

long gcd_long(long a, long b) {
    long t;
    if (a < 0) a = -a;
    if (b < 0) b = -b;
    while (b != 0) {
        t = b;
        b = a % b;
        a = t;
    }
    return a;
}

void apf_gcd(apf *r, const apf *a, const apf *b) {
    long ai = apf_to_long(a);
    long bi = apf_to_long(b);
    apf_from_int(r, gcd_long(ai, bi));
}

long lcm_long(long a, long b) {
    long g;
    if (a < 0) a = -a;
    if (b < 0) b = -b;
    if (a == 0 || b == 0) return 0;
    g = gcd_long(a, b);
    return (a / g) * b;  /* Avoid overflow by dividing first */
}

void apf_lcm(apf *r, const apf *a, const apf *b) {
    long ai = apf_to_long(a);
    long bi = apf_to_long(b);
    apf_from_int(r, lcm_long(ai, bi));
}

int is_prime_long(long n) {
    long i;
    if (n < 2) return 0;
    if (n == 2) return 1;
    if (n % 2 == 0) return 0;
    for (i = 3; i * i <= n; i += 2) {
        if (n % i == 0) return 0;
    }
    return 1;
}

int prime_factors(long n, long *factors, int max_factors) {
    int count = 0;
    long d = 2;
    
    if (n < 0) n = -n;
    if (n <= 1) return 0;
    
    while (d * d <= n && count < max_factors) {
        while (n % d == 0) {
            factors[count++] = d;
            n /= d;
            if (count >= max_factors) return count;
        }
        d++;
    }
    if (n > 1 && count < max_factors) {
        factors[count++] = n;
    }
    return count;
}

#endif /* HAVE_GCD */

/* ========== Bitwise Operations ========== */
#ifdef HAVE_BITWISE

void apf_and(apf *r, const apf *a, const apf *b) {
    long ai = apf_to_long(a);
    long bi = apf_to_long(b);
    apf_from_int(r, ai & bi);
}

void apf_or(apf *r, const apf *a, const apf *b) {
    long ai = apf_to_long(a);
    long bi = apf_to_long(b);
    apf_from_int(r, ai | bi);
}

void apf_xor(apf *r, const apf *a, const apf *b) {
    long ai = apf_to_long(a);
    long bi = apf_to_long(b);
    apf_from_int(r, ai ^ bi);
}

void apf_not(apf *r, const apf *a) {
    long ai = apf_to_long(a);
    apf_from_int(r, ~ai);
}

void apf_lsl(apf *r, const apf *a, int bits) {
    long ai = apf_to_long(a);
    if (bits < 0) bits = 0;
    if (bits > 63) bits = 63;
    apf_from_int(r, ai << bits);
}

void apf_lsr(apf *r, const apf *a, int bits) {
    unsigned long ai = (unsigned long)apf_to_long(a);
    if (bits < 0) bits = 0;
    if (bits > 63) bits = 63;
    apf_from_int(r, (long)(ai >> bits));
}

#endif /* HAVE_BITWISE */

/* ========== Angle Mode Conversions ========== */

/* pi/180 for deg->rad conversion */
static void get_deg_to_rad_factor(apf *r) {
    apf pi, d180;
    apfx_pi(&pi);
    apf_from_int(&d180, 180);
    apf_div(r, &pi, &d180);
}

/* 180/pi for rad->deg conversion */
static void get_rad_to_deg_factor(apf *r) {
    apf pi, d180;
    apfx_pi(&pi);
    apf_from_int(&d180, 180);
    apf_div(r, &d180, &pi);
}

/* pi/200 for grad->rad conversion */
static void get_grad_to_rad_factor(apf *r) {
    apf pi, d200;
    apfx_pi(&pi);
    apf_from_int(&d200, 200);
    apf_div(r, &pi, &d200);
}

/* 200/pi for rad->grad conversion */
static void get_rad_to_grad_factor(apf *r) {
    apf pi, d200;
    apfx_pi(&pi);
    apf_from_int(&d200, 200);
    apf_div(r, &d200, &pi);
}

void deg_to_rad(apf *r, const apf *deg) {
    apf factor;
    get_deg_to_rad_factor(&factor);
    apf_mul(r, deg, &factor);
}

void rad_to_deg(apf *r, const apf *rad) {
    apf factor;
    get_rad_to_deg_factor(&factor);
    apf_mul(r, rad, &factor);
}

void grad_to_rad(apf *r, const apf *grad) {
    apf factor;
    get_grad_to_rad_factor(&factor);
    apf_mul(r, grad, &factor);
}

void rad_to_grad(apf *r, const apf *rad) {
    apf factor;
    get_rad_to_grad_factor(&factor);
    apf_mul(r, rad, &factor);
}

/* Convert angle to radians based on mode (0=rad, 1=deg, 2=grad) */
void angle_to_rad(apf *r, const apf *angle, int mode) {
    switch (mode) {
        case 1: deg_to_rad(r, angle); break;
        case 2: grad_to_rad(r, angle); break;
        default: apf_copy(r, angle); break;
    }
}

/* Convert radians to angle based on mode */
void rad_to_angle(apf *r, const apf *rad, int mode) {
    switch (mode) {
        case 1: rad_to_deg(r, rad); break;
        case 2: rad_to_grad(r, rad); break;
        default: apf_copy(r, rad); break;
    }
}

/* ========== Reciprocal Trig Functions ========== */

void apfx_sec(apf *r, const apf *x) {
    apf cos_x, one;
    apfx_cos(&cos_x, x);
    apf_from_int(&one, 1);
    apf_div(r, &one, &cos_x);
}

void apfx_csc(apf *r, const apf *x) {
    apf sin_x, one;
    apfx_sin(&sin_x, x);
    apf_from_int(&one, 1);
    apf_div(r, &one, &sin_x);
}

void apfx_cot(apf *r, const apf *x) {
    apf tan_x, one;
    apfx_tan(&tan_x, x);
    apf_from_int(&one, 1);
    apf_div(r, &one, &tan_x);
}

void apfx_asec(apf *r, const apf *x) {
    /* arcsec(x) = acos(1/x) */
    apf inv, one;
    apf_from_int(&one, 1);
    apf_div(&inv, &one, x);
    apfx_acos(r, &inv);
}

void apfx_acsc(apf *r, const apf *x) {
    /* arccsc(x) = asin(1/x) */
    apf inv, one;
    apf_from_int(&one, 1);
    apf_div(&inv, &one, x);
    apfx_asin(r, &inv);
}

void apfx_acot(apf *r, const apf *x) {
    /* arccot(x) = atan(1/x) for x > 0, pi + atan(1/x) for x < 0 */
    apf inv, one;
    apf_from_int(&one, 1);
    apf_div(&inv, &one, x);
    apfx_atan(r, &inv);
    /* Adjust for negative x */
    if (x->sign && !apf_is_zero(x)) {
        apf pi;
        apfx_pi(&pi);
        apf_add(r, r, &pi);
    }
}

/* ========== Log with Arbitrary Base ========== */

void apfx_logb(apf *r, const apf *x, const apf *base) {
    /* log_b(x) = ln(x) / ln(b) */
    apf log_x, log_base;
    apfx_log(&log_x, x);
    apfx_log(&log_base, base);
    apf_div(r, &log_x, &log_base);
}
