/* mathx.h - Extended math functions
 * C89 portable for DOS, Linux, VIC-20
 * 
 * Combinatorics: nPr, nCr
 * Number theory: GCD, LCM
 * Bitwise: AND, OR, XOR, NOT, shifts
 */

#ifndef MATHX_H
#define MATHX_H

#include "config.h"
#include "apf.h"

/* ========== Combinatorics ========== */
#ifdef HAVE_COMB

/* nPr = n! / (n-r)! - permutations */
void apf_npr(apf *r, long n, long k);

/* nCr = n! / (r! * (n-r)!) - combinations */
void apf_ncr(apf *r, long n, long k);

#endif /* HAVE_COMB */

/* ========== Number Theory ========== */
#ifdef HAVE_GCD

/* Greatest Common Divisor */
long gcd_long(long a, long b);
void apf_gcd(apf *r, const apf *a, const apf *b);

/* Least Common Multiple */
long lcm_long(long a, long b);
void apf_lcm(apf *r, const apf *a, const apf *b);

/* Prime test (simple trial division) */
int is_prime_long(long n);

/* Prime factorization (stores factors in array, returns count) */
int prime_factors(long n, long *factors, int max_factors);

/* nth prime number (nthprime(1)=2, nthprime(2)=3, ...) */
long nthprime_long(long n);

/* Next prime >= n */
long nextprime_long(long n);

/* Previous prime <= n */
long prevprime_long(long n);

/* Count of primes <= n (prime counting function) */
long primepi_long(long n);

/* Product of distinct prime factors */
long radical_long(long n);

/* Count of distinct prime factors (omega) */
long omega_long(long n);

/* Count of prime factors with multiplicity (bigomega/Omega) */
long bigomega_long(long n);

/* Is squarefree? */
int issquarefree_long(long n);

/* Möbius function */
int moebius_long(long n);

/* Euler's totient */
long totient_long(long n);

/* Sum of divisors to power k (sigma_k) */
long sigma_long(long n, int k);

/* Number of divisors */
long numdivisors_long(long n);

/* Sum of digits */
long digsum_long(long n);

/* Number of digits */
long numdigits_long(long n);

/* Digital root */
long digitalroot_long(long n);

/* Reverse digits */
long reverse_long(long n);

/* Is palindrome? */
int ispalindrome_long(long n);

/* Perfect, abundant, deficient tests */
int isperfect_long(long n);
int isabundant_long(long n);
int isdeficient_long(long n);

/* Figurate numbers */
long triangular_long(long n);
int istriangular_long(long n);
long square_long(long n);
int issquare_long(long n);
long pentagonal_long(long n);
long hexagonal_long(long n);

/* Is perfect power? */
int ispower_long(long n);

/* Subfactorial / derangements */
void apf_subfactorial(apf *r, long n);

/* Bell numbers */
void apf_bell(apf *r, long n);
void apf_fibonacci(apf *r, long n);
void apf_lucas(apf *r, long n);
void apf_catalan(apf *r, long n);

/* Stirling numbers of second kind */
void apf_stirling2(apf *r, long n, long k);

/* Partition function */
void apf_partition(apf *r, long n);

#endif /* HAVE_GCD */

/* ========== Bitwise Operations ========== */
#ifdef HAVE_BITWISE

/* Bitwise AND */
void apf_and(apf *r, const apf *a, const apf *b);

/* Bitwise OR */
void apf_or(apf *r, const apf *a, const apf *b);

/* Bitwise XOR */
void apf_xor(apf *r, const apf *a, const apf *b);

/* Bitwise NOT (ones complement) */
void apf_not(apf *r, const apf *a);

/* Left shift */
void apf_lsl(apf *r, const apf *a, int bits);

/* Right shift */
void apf_lsr(apf *r, const apf *a, int bits);

#endif /* HAVE_BITWISE */

/* ========== Angle Mode Conversions ========== */
/* Convert angle to radians based on current mode */
void angle_to_rad(apf *r, const apf *angle, int mode);

/* Convert radians to angle based on current mode */
void rad_to_angle(apf *r, const apf *rad, int mode);

/* Degree/radian/gradian conversions */
void deg_to_rad(apf *r, const apf *deg);
void rad_to_deg(apf *r, const apf *rad);
void grad_to_rad(apf *r, const apf *grad);
void rad_to_grad(apf *r, const apf *rad);

/* Reciprocal trig functions */
void apfx_sec(apf *r, const apf *x);   /* sec(x) = 1/cos(x) */
void apfx_csc(apf *r, const apf *x);   /* csc(x) = 1/sin(x) */
void apfx_cot(apf *r, const apf *x);   /* cot(x) = 1/tan(x) */
void apfx_asec(apf *r, const apf *x);  /* arcsec(x) = acos(1/x) */
void apfx_acsc(apf *r, const apf *x);  /* arccsc(x) = asin(1/x) */
void apfx_acot(apf *r, const apf *x);  /* arccot(x) = atan(1/x) */

/* Log with arbitrary base */
void apfx_logb(apf *r, const apf *x, const apf *base);

#endif /* MATHX_H */
