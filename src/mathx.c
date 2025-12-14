/* mathx.c - Extended math functions
 * C89 portable for DOS, Linux, VIC-20
 */

#include "config.h"
#include "mathx.h"
#include "apf.h"
#include "apfx.h"

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
