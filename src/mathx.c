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

/* nth prime number (1-indexed: nthprime(1)=2, nthprime(2)=3, ...) */
long nthprime_long(long n) {
    long count = 0, candidate = 1;
    if (n <= 0) return 0;
    if (n == 1) return 2;
    while (count < n) {
        candidate++;
        if (is_prime_long(candidate)) count++;
    }
    return candidate;
}

/* Smallest prime >= n */
long nextprime_long(long n) {
    if (n <= 2) return 2;
    if (n % 2 == 0) n++;
    while (!is_prime_long(n)) n += 2;
    return n;
}

/* Largest prime <= n, or 0 if n < 2 */
long prevprime_long(long n) {
    if (n < 2) return 0;
    if (n == 2) return 2;
    if (n % 2 == 0) n--;
    while (n >= 2 && !is_prime_long(n)) n -= 2;
    return (n >= 2) ? n : 0;
}

/* Count of primes <= n */
long primepi_long(long n) {
    long count = 0, i;
    if (n < 2) return 0;
    count = 1;  /* 2 is prime */
    for (i = 3; i <= n; i += 2) {
        if (is_prime_long(i)) count++;
    }
    return count;
}

/* Product of distinct prime factors (radical) */
long radical_long(long n) {
    long result = 1, d = 2;
    if (n < 0) n = -n;
    if (n <= 1) return n;
    while (d * d <= n) {
        if (n % d == 0) {
            result *= d;
            while (n % d == 0) n /= d;
        }
        d++;
    }
    if (n > 1) result *= n;
    return result;
}

/* Count of distinct prime factors (omega) */
long omega_long(long n) {
    long count = 0, d = 2;
    if (n < 0) n = -n;
    if (n <= 1) return 0;
    while (d * d <= n) {
        if (n % d == 0) {
            count++;
            while (n % d == 0) n /= d;
        }
        d++;
    }
    if (n > 1) count++;
    return count;
}

/* Count of prime factors with multiplicity (bigomega/Omega) */
long bigomega_long(long n) {
    long count = 0, d = 2;
    if (n < 0) n = -n;
    if (n <= 1) return 0;
    while (d * d <= n) {
        while (n % d == 0) {
            count++;
            n /= d;
        }
        d++;
    }
    if (n > 1) count++;
    return count;
}

/* Is squarefree (no repeated prime factors)? */
int issquarefree_long(long n) {
    long d = 2;
    if (n < 0) n = -n;
    if (n <= 1) return (n == 1);
    while (d * d <= n) {
        if (n % d == 0) {
            n /= d;
            if (n % d == 0) return 0;  /* d^2 divides original n */
        }
        d++;
    }
    return 1;
}

/* Möbius function: 0 if n has squared prime factor, (-1)^k if k distinct primes */
int moebius_long(long n) {
    long count = 0, d = 2;
    if (n < 0) n = -n;
    if (n == 1) return 1;
    while (d * d <= n) {
        if (n % d == 0) {
            n /= d;
            if (n % d == 0) return 0;  /* Squared factor */
            count++;
        }
        d++;
    }
    if (n > 1) count++;
    return (count % 2 == 0) ? 1 : -1;
}

/* Euler's totient function */
long totient_long(long n) {
    long result = n, p = 2;
    if (n <= 0) return 0;
    while (p * p <= n) {
        if (n % p == 0) {
            while (n % p == 0) n /= p;
            result -= result / p;
        }
        p++;
    }
    if (n > 1) result -= result / n;
    return result;
}

/* Sum of divisors (sigma function) */
long sigma_long(long n, int k) {
    long result = 0, d, pk;
    if (n <= 0) return 0;
    for (d = 1; d * d <= n; d++) {
        if (n % d == 0) {
            /* Compute d^k */
            pk = 1;
            if (k > 0) { long i; for (i = 0; i < k; i++) pk *= d; }
            result += pk;
            if (d != n / d) {
                pk = 1;
                if (k > 0) { long i, dd = n/d; for (i = 0; i < k; i++) pk *= dd; }
                result += pk;
            }
        }
    }
    return result;
}

/* Number of divisors (sigma_0) */
long numdivisors_long(long n) {
    return sigma_long(n, 0);
}

/* Sum of digits in base 10 */
long digsum_long(long n) {
    long sum = 0;
    if (n < 0) n = -n;
    while (n > 0) {
        sum += n % 10;
        n /= 10;
    }
    return sum;
}

/* Number of digits in base 10 */
long numdigits_long(long n) {
    long count = 0;
    if (n == 0) return 1;
    if (n < 0) n = -n;
    while (n > 0) {
        count++;
        n /= 10;
    }
    return count;
}

/* Digital root (repeated digit sum until single digit) */
long digitalroot_long(long n) {
    if (n < 0) n = -n;
    if (n == 0) return 0;
    return 1 + (n - 1) % 9;
}

/* Reverse digits */
long reverse_long(long n) {
    long rev = 0, sign = 1;
    if (n < 0) { sign = -1; n = -n; }
    while (n > 0) {
        rev = rev * 10 + n % 10;
        n /= 10;
    }
    return rev * sign;
}

/* Is palindrome? */
int ispalindrome_long(long n) {
    if (n < 0) return 0;
    return n == reverse_long(n);
}

/* Is perfect number? (sum of proper divisors == n) */
int isperfect_long(long n) {
    return (n > 0 && sigma_long(n, 1) - n == n);
}

/* Is abundant? (sum of proper divisors > n) */
int isabundant_long(long n) {
    return (n > 0 && sigma_long(n, 1) - n > n);
}

/* Is deficient? (sum of proper divisors < n) */
int isdeficient_long(long n) {
    return (n > 0 && sigma_long(n, 1) - n < n);
}

/* nth triangular number = n*(n+1)/2 */
long triangular_long(long n) {
    return (n >= 0) ? n * (n + 1) / 2 : 0;
}

/* Is triangular? */
int istriangular_long(long n) {
    /* n is triangular if 8n+1 is a perfect square */
    long t, s;
    if (n < 0) return 0;
    t = 8 * n + 1;
    s = 1;
    while (s * s < t) s++;
    return (s * s == t);
}

/* nth square number = n^2 */
long square_long(long n) {
    return n * n;
}

/* Is perfect square? */
int issquare_long(long n) {
    long s;
    if (n < 0) return 0;
    s = 1;
    while (s * s < n) s++;
    return (s * s == n);
}

/* nth pentagonal number = n*(3n-1)/2 */
long pentagonal_long(long n) {
    return (n >= 0) ? n * (3 * n - 1) / 2 : 0;
}

/* nth hexagonal number = n*(2n-1) */
long hexagonal_long(long n) {
    return (n >= 0) ? n * (2 * n - 1) : 0;
}

/* Is perfect power? (n = a^b for some a, b >= 2) */
int ispower_long(long n) {
    long b, a, t;
    if (n < 0) n = -n;
    if (n <= 1) return 0;
    for (b = 2; b <= 40; b++) {  /* 2^40 > 10^12 */
        a = 2;
        t = 1;
        while (t < n) {
            long i;
            t = 1;
            for (i = 0; i < b; i++) t *= a;
            if (t == n) return 1;
            a++;
        }
    }
    return 0;
}

/* Subfactorial / derangements !n */
void apf_subfactorial(apf *r, long n) {
    /* !n = n! * sum(k=0 to n) of (-1)^k / k!
     * Recurrence: !n = (n-1) * (!(n-1) + !(n-2)), !0 = 1, !1 = 0 */
    long i;
    apf a, b, tmp, tmp2, n_minus_1;
    if (n < 0) { apf_zero(r); return; }
    if (n == 0) { apf_from_int(r, 1); return; }
    if (n == 1) { apf_from_int(r, 0); return; }
    
    apf_from_int(&a, 1);  /* !0 = 1 */
    apf_from_int(&b, 0);  /* !1 = 0 */
    for (i = 2; i <= n; i++) {
        apf_from_int(&n_minus_1, i - 1);
        apf_add(&tmp, &a, &b);
        apf_mul(&tmp2, &n_minus_1, &tmp);
        a = b;
        b = tmp2;
    }
    *r = b;
}

/* Bell numbers B(n) = number of partitions of set of n elements */
void apf_bell(apf *r, long n) {
    /* Bell triangle method */
    apf row[101];  /* Support up to n=100 */
    long i, j;
    if (n < 0 || n > 100) { apf_from_int(r, 0); return; }
    if (n == 0) { apf_from_int(r, 1); return; }
    
    apf_from_int(&row[0], 1);
    for (i = 1; i <= n; i++) {
        /* row[0] = previous row's last element */
        apf prev_last = row[i-1];
        row[0] = prev_last;
        for (j = 1; j <= i; j++) {
            apf_add(&row[j], &row[j-1], &row[j-1]);
            /* Actually: row[j] = row[j-1] + old_row[j-1] 
             * Need proper Bell triangle calculation */
        }
    }
    /* Simplified: use recurrence B(n) = sum(k=0..n-1) C(n-1,k) * B(k) */
    {
        apf bells[101], sum, comb, tmp;
        apf_from_int(&bells[0], 1);
        for (i = 1; i <= n; i++) {
            apf_from_int(&sum, 0);
            for (j = 0; j < i; j++) {
                apf_ncr(&comb, i - 1, j);
                apf_mul(&tmp, &comb, &bells[j]);
                apf_add(&sum, &sum, &tmp);
            }
            bells[i] = sum;
        }
        *r = bells[n];
    }
}

/* Stirling numbers of the second kind S(n,k) */
void apf_stirling2(apf *r, long n, long k) {
    /* S(n,k) = k*S(n-1,k) + S(n-1,k-1) */
    /* S(n,0) = 0 for n>0, S(0,0) = 1, S(n,n) = 1 */
    apf S[101][101];
    long i, j;
    if (n < 0 || k < 0 || k > n || n > 100) { apf_from_int(r, 0); return; }
    
    for (i = 0; i <= n; i++) {
        for (j = 0; j <= k; j++) {
            if (i == 0 && j == 0) {
                apf_from_int(&S[i][j], 1);
            } else if (j == 0) {
                apf_from_int(&S[i][j], 0);
            } else if (i == j) {
                apf_from_int(&S[i][j], 1);
            } else if (j > i) {
                apf_from_int(&S[i][j], 0);
            } else {
                apf tmp1, tmp2, kval;
                apf_from_int(&kval, j);
                apf_mul(&tmp1, &kval, &S[i-1][j]);
                apf_add(&tmp2, &tmp1, &S[i-1][j-1]);
                S[i][j] = tmp2;
            }
        }
    }
    *r = S[n][k];
}

/* Partition function p(n) - number of ways to write n as sum of positive integers */
void apf_partition(apf *r, long n) {
    /* Using pentagonal number theorem recurrence */
    apf p[1001];
    long i, k, pent;
    int sign;
    if (n < 0 || n > 1000) { apf_from_int(r, 0); return; }
    
    apf_from_int(&p[0], 1);
    for (i = 1; i <= n; i++) {
        apf sum;
        apf_from_int(&sum, 0);
        for (k = 1; ; k++) {
            /* Generalized pentagonal: k(3k-1)/2 and k(3k+1)/2 */
            pent = k * (3*k - 1) / 2;
            if (pent > i) break;
            sign = (k % 2 == 1) ? 1 : -1;
            if (sign > 0) {
                apf_add(&sum, &sum, &p[i - pent]);
            } else {
                apf_sub(&sum, &sum, &p[i - pent]);
            }
            
            pent = k * (3*k + 1) / 2;
            if (pent > i) break;
            if (sign > 0) {
                apf_add(&sum, &sum, &p[i - pent]);
            } else {
                apf_sub(&sum, &sum, &p[i - pent]);
            }
        }
        p[i] = sum;
    }
    *r = p[n];
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
