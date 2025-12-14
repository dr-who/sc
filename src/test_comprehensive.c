/* test_comprehensive.c - Comprehensive test suite for scalc
 * C89 portable - Tests for lunar/Mars mission-critical calculations
 * 
 * Test categories:
 *   1. APF core (arithmetic, precision)
 *   2. Trigonometry (orbital angles)
 *   3. Exponential/logarithm
 *   4. Complex numbers
 *   5. Combinatorics
 *   6. Number theory (GCD, LCM)
 *   7. Bitwise operations
 *   8. Statistics
 *   9. Orbital mechanics
 *   10. Edge cases and stress tests
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "config.h"
#include "apf.h"
#include "apfx.h"
#include "apfc.h"
#include "mathx.h"
#ifdef HAVE_ORBITAL
#include "orbital.h"
#endif

/* For standalone compilation */
int display_digits = 16;

static int tests_passed = 0;
static int tests_failed = 0;

#define TOLERANCE "1e-30"

/* Test helper: compare APF with expected string */
static int apf_test(const char *name, const apf *val, const char *expected) {
    char buf[256];
    apf exp, diff, tol;
    int result;
    
    apf_to_str(buf, sizeof(buf), val, 0);
    
    /* Parse expected and compute difference */
    apf_from_str(&exp, expected);
    apf_sub(&diff, val, &exp);
    apf_abs(&diff, &diff);
    apf_from_str(&tol, TOLERANCE);
    
    /* For very small expected values, use relative tolerance */
    if (!apf_is_zero(&exp)) {
        apf rel;
        apf_abs(&rel, &exp);
        apf_mul(&tol, &tol, &rel);
    }
    
    result = (apf_cmp(&diff, &tol) <= 0);
    
    if (result) {
        tests_passed++;
        printf("  PASS: %s = %s\n", name, buf);
    } else {
        tests_failed++;
        printf("  FAIL: %s\n", name);
        printf("        got:      %s\n", buf);
        printf("        expected: %s\n", expected);
    }
    
    return result;
}

/* Test helper: compare APF with expected int */
static int apf_test_int(const char *name, const apf *val, long expected) {
    char buf[64];
    sprintf(buf, "%ld", expected);
    return apf_test(name, val, buf);
}

/* Test helper: verify within percentage tolerance */
static int apf_test_pct(const char *name, const apf *val, const char *expected, double pct_tol) {
    char buf[256];
    apf exp, diff, tolerance, abs_exp, abs_val;
    int result;
    
    apf_from_str(&exp, expected);
    apf_sub(&diff, val, &exp);
    apf_abs(&diff, &diff);
    
    /* Use larger of abs(val) or abs(exp) for relative tolerance */
    apf_abs(&abs_exp, &exp);
    apf_abs(&abs_val, val);
    if (apf_cmp(&abs_val, &abs_exp) > 0) {
        apf_copy(&abs_exp, &abs_val);
    }
    
    /* tolerance = abs_exp * pct_tol */
    {
        char pct_str[64];
        apf pct;
        sprintf(pct_str, "%.15e", pct_tol);
        apf_from_str(&pct, pct_str);
        apf_mul(&tolerance, &abs_exp, &pct);
    }
    
    result = (apf_cmp(&diff, &tolerance) <= 0);
    
    apf_to_str(buf, sizeof(buf), val, 10);
    if (result) {
        tests_passed++;
        printf("  PASS: %s = %s (within %.2f%%)\n", name, buf, pct_tol * 100);
    } else {
        tests_failed++;
        printf("  FAIL: %s = %s (expected %s ±%.2f%%)\n", name, buf, expected, pct_tol * 100);
    }
    
    return result;
}

/* ========== Test Categories ========== */

void test_apf_arithmetic(void) {
    apf a, b, r;
    
    printf("\n=== APF Arithmetic Tests ===\n");
    
    /* Basic operations */
    apf_from_int(&a, 123456789);
    apf_from_int(&b, 987654321);
    apf_add(&r, &a, &b);
    apf_test_int("123456789 + 987654321", &r, 1111111110);
    
    apf_sub(&r, &b, &a);
    apf_test_int("987654321 - 123456789", &r, 864197532);
    
    apf_mul(&r, &a, &b);
    apf_test("123456789 * 987654321", &r, "121932631112635269");
    
    apf_from_int(&a, 1000000000);
    apf_from_int(&b, 1000000000);
    apf_mul(&r, &a, &b);
    apf_test("1e9 * 1e9", &r, "1000000000000000000");
    
    /* Division */
    apf_from_int(&a, 355);
    apf_from_int(&b, 113);
    apf_div(&r, &a, &b);
    apf_test("355/113 (pi approx)", &r, "3.14159292035398230088495575221238938053");
    
    /* Large numbers */
    apf_from_str(&a, "123456789012345678901234567890");
    apf_from_str(&b, "987654321098765432109876543210");
    apf_add(&r, &a, &b);
    apf_test("30-digit addition", &r, "1111111110111111111011111111100");
}

void test_apf_precision(void) {
    apf a, b, r;
    
    printf("\n=== APF Precision Tests (Critical for Orbital Calculations) ===\n");
    
    /* Test 38-digit precision */
    apf_from_str(&a, "1.0000000000000000000000000000000000001");
    apf_from_str(&b, "0.0000000000000000000000000000000000001");
    apf_add(&r, &a, &b);
    apf_test("38-digit precision add", &r, "1.0000000000000000000000000000000000002");
    
    /* Subtraction precision (catastrophic cancellation test) */
    apf_from_str(&a, "1.00000000000000001");
    apf_from_str(&b, "1.0");
    apf_sub(&r, &a, &b);
    apf_test("Near-cancellation", &r, "0.00000000000000001");
    
    /* Multiplication precision */
    apf_from_str(&a, "1.23456789012345678901234567890123456789");
    apf_from_str(&b, "9.87654321098765432109876543210987654321");
    apf_mul(&r, &a, &b);
    /* Result should be accurate to ~38 digits */
    apf_test_pct("38-digit multiplication", &r, "12.193263111263526853698893396942152881", 1e-30);
}

void test_trig_functions(void) {
    apf a, r, pi;
    
    printf("\n=== Trigonometry Tests (Essential for Orbital Mechanics) ===\n");
    
    apf_from_str(&pi, "3.14159265358979323846264338327950288419");
    
    /* Basic trig */
    apf_zero(&a);
    apfx_sin(&r, &a);
    apf_test("sin(0)", &r, "0");
    
    apfx_cos(&r, &a);
    apf_test("cos(0)", &r, "1");
    
    /* pi/6, pi/4, pi/3 */
    {
        apf six, four, three;
        apf_from_int(&six, 6);
        apf_from_int(&four, 4);
        apf_from_int(&three, 3);
        
        apf_div(&a, &pi, &six);
        apfx_sin(&r, &a);
        apf_test("sin(pi/6) = 0.5", &r, "0.5");
        
        apf_div(&a, &pi, &four);
        apfx_sin(&r, &a);
        apf_test("sin(pi/4) = sqrt(2)/2", &r, "0.70710678118654752440084436210484903928");
        
        apfx_cos(&r, &a);
        apf_test("cos(pi/4) = sqrt(2)/2", &r, "0.70710678118654752440084436210484903928");
        
        apf_div(&a, &pi, &three);
        apfx_sin(&r, &a);
        apf_test("sin(pi/3) = sqrt(3)/2", &r, "0.86602540378443864676372317075293618347");
    }
    
    /* Identity: sin^2 + cos^2 = 1 */
    {
        apf angle, s, c, s2, c2, sum;
        apf_from_str(&angle, "0.7");  /* Arbitrary angle */
        apfx_sin(&s, &angle);
        apfx_cos(&c, &angle);
        apf_mul(&s2, &s, &s);
        apf_mul(&c2, &c, &c);
        apf_add(&sum, &s2, &c2);
        apf_test("sin^2(0.7) + cos^2(0.7) = 1", &sum, "1");
    }
    
    /* Inverse trig */
    apf_from_int(&a, 1);
    apfx_atan(&r, &a);
    apf_test("atan(1) = pi/4", &r, "0.78539816339744830961566084581987572105");
    
    apf_from_str(&a, "0.5");
    apfx_asin(&r, &a);
    apf_test("asin(0.5) = pi/6", &r, "0.52359877559829887307710723054658381363");
    
    apfx_acos(&r, &a);
    apf_test("acos(0.5) = pi/3", &r, "1.04719755119659774615421446109316762845");
}

void test_exp_log(void) {
    apf a, r, e;
    
    printf("\n=== Exponential/Logarithm Tests ===\n");
    
    apf_from_str(&e, "2.71828182845904523536028747135266249775");
    
    /* Basic exp/log */
    apf_from_int(&a, 1);
    apfx_exp(&r, &a);
    apf_test("exp(1) = e", &r, "2.71828182845904523536028747135266249775");
    
    apfx_log(&r, &e);
    apf_test("ln(e) = 1", &r, "1");
    
    apf_from_int(&a, 10);
    apfx_log(&r, &a);
    apf_test("ln(10)", &r, "2.30258509299404568401799145468436420760");
    
    /* exp(ln(x)) = x */
    apf_from_str(&a, "123.456");
    apfx_log(&r, &a);
    apfx_exp(&r, &r);
    apf_test("exp(ln(123.456)) = 123.456", &r, "123.456");
    
    /* Powers via exp/log */
    apf_from_int(&a, 2);
    {
        apf ten;
        apf_from_int(&ten, 10);
        apfx_pow(&r, &a, &ten);
        apf_test("2^10 = 1024", &r, "1024");
    }
}

#ifdef HAVE_COMB
void test_combinatorics(void) {
    apf r;
    
    printf("\n=== Combinatorics Tests ===\n");
    
    /* nPr tests */
    apf_npr(&r, 10, 3);
    apf_test_int("P(10,3) = 720", &r, 720);
    
    apf_npr(&r, 10, 10);
    apf_test_int("P(10,10) = 10! = 3628800", &r, 3628800);
    
    apf_npr(&r, 5, 0);
    apf_test_int("P(5,0) = 1", &r, 1);
    
    /* nCr tests */
    apf_ncr(&r, 10, 3);
    apf_test_int("C(10,3) = 120", &r, 120);
    
    apf_ncr(&r, 10, 5);
    apf_test_int("C(10,5) = 252", &r, 252);
    
    apf_ncr(&r, 52, 5);
    apf_test_int("C(52,5) = 2598960 (poker hands)", &r, 2598960);
    
    /* Pascal's triangle identity: C(n,k) = C(n-1,k-1) + C(n-1,k) */
    {
        apf c1, c2, sum;
        apf_ncr(&c1, 9, 4);
        apf_ncr(&c2, 9, 5);
        apf_add(&sum, &c1, &c2);
        apf_test_int("C(9,4)+C(9,5) = C(10,5)", &sum, 252);
    }
}
#endif

#ifdef HAVE_GCD
void test_number_theory(void) {
    apf a, b, r;
    
    printf("\n=== Number Theory Tests ===\n");
    
    /* GCD tests */
    apf_from_int(&a, 48);
    apf_from_int(&b, 18);
    apf_gcd(&r, &a, &b);
    apf_test_int("gcd(48, 18) = 6", &r, 6);
    
    apf_from_int(&a, 1071);
    apf_from_int(&b, 462);
    apf_gcd(&r, &a, &b);
    apf_test_int("gcd(1071, 462) = 21", &r, 21);
    
    apf_from_int(&a, 17);
    apf_from_int(&b, 13);
    apf_gcd(&r, &a, &b);
    apf_test_int("gcd(17, 13) = 1 (coprime)", &r, 1);
    
    /* LCM tests */
    apf_from_int(&a, 12);
    apf_from_int(&b, 18);
    apf_lcm(&r, &a, &b);
    apf_test_int("lcm(12, 18) = 36", &r, 36);
    
    apf_from_int(&a, 7);
    apf_from_int(&b, 11);
    apf_lcm(&r, &a, &b);
    apf_test_int("lcm(7, 11) = 77 (primes)", &r, 77);
    
    /* Prime tests */
    if (is_prime_long(2)) { tests_passed++; printf("  PASS: 2 is prime\n"); } 
    else { tests_failed++; printf("  FAIL: 2 is prime\n"); }
    
    if (is_prime_long(17)) { tests_passed++; printf("  PASS: 17 is prime\n"); } 
    else { tests_failed++; printf("  FAIL: 17 is prime\n"); }
    
    if (!is_prime_long(15)) { tests_passed++; printf("  PASS: 15 is not prime\n"); } 
    else { tests_failed++; printf("  FAIL: 15 is not prime\n"); }
    
    if (is_prime_long(7919)) { tests_passed++; printf("  PASS: 7919 is prime (1000th prime)\n"); } 
    else { tests_failed++; printf("  FAIL: 7919 is prime\n"); }
}
#endif

#ifdef HAVE_BITWISE
void test_bitwise(void) {
    apf a, b, r;
    
    printf("\n=== Bitwise Operations Tests ===\n");
    
    apf_from_int(&a, 0xFF);
    apf_from_int(&b, 0x0F);
    
    apf_and(&r, &a, &b);
    apf_test_int("0xFF AND 0x0F = 0x0F", &r, 0x0F);
    
    apf_or(&r, &a, &b);
    apf_test_int("0xFF OR 0x0F = 0xFF", &r, 0xFF);
    
    apf_xor(&r, &a, &b);
    apf_test_int("0xFF XOR 0x0F = 0xF0", &r, 0xF0);
    
    apf_from_int(&a, 0);
    apf_not(&r, &a);
    apf_test_int("NOT 0 = -1", &r, -1);
    
    apf_from_int(&a, 1);
    apf_lsl(&r, &a, 8);
    apf_test_int("1 << 8 = 256", &r, 256);
    
    apf_from_int(&a, 256);
    apf_lsr(&r, &a, 4);
    apf_test_int("256 >> 4 = 16", &r, 16);
}
#endif

#ifdef HAVE_ORBITAL
void test_orbital_mechanics(void) {
    apf mu_earth, r, v, T, dv;
    
    printf("\n=== Orbital Mechanics Tests (Lunar/Mars Mission Critical) ===\n");
    
    apf_from_str(&mu_earth, "398600.4418");  /* km^3/s^2 */
    
    /* LEO orbital velocity (~400km altitude = 6778km radius) */
    apf_from_str(&r, "6778");
    orbital_velocity(&v, &mu_earth, &r);
    apf_test_pct("LEO orbital velocity (400km)", &v, "7.6686", 0.001);
    
    /* Escape velocity at LEO */
    escape_velocity(&v, &mu_earth, &r);
    apf_test_pct("LEO escape velocity", &v, "10.845", 0.001);
    
    /* ISS orbital period (~92 minutes = 5520 seconds) */
    orbital_period(&T, &mu_earth, &r);
    apf_test_pct("ISS orbital period", &T, "5553", 0.01);
    
    /* GEO orbital velocity (r = 42164 km) */
    apf_from_str(&r, "42164");
    orbital_velocity(&v, &mu_earth, &r);
    apf_test_pct("GEO orbital velocity", &v, "3.0747", 0.001);
    
    /* GEO orbital period (24 hours = 86400 seconds) */
    orbital_period(&T, &mu_earth, &r);
    apf_test_pct("GEO orbital period (24h)", &T, "86164", 0.01);
    
    /* Hohmann transfer LEO to GEO */
    {
        apf r1, r2;
        apf_from_str(&r1, "6778");
        apf_from_str(&r2, "42164");
        hohmann_deltav(&dv, &mu_earth, &r1, &r2);
        apf_test_pct("Hohmann LEO->GEO delta-V", &dv, "3.854", 0.01);
    }
    
    /* Lunar transfer (Earth orbit to Moon orbit) */
    /* TLI delta-V is approximately 3.13 km/s from LEO */
    printf("  INFO: Trans-Lunar Injection from LEO requires ~3.13 km/s\n");
    
    /* Moon orbital velocity (r = 1838 km from center, mu_moon = 4904.87) */
    {
        apf mu_moon, r_moon;
        apf_from_str(&mu_moon, "4904.87");
        apf_from_str(&r_moon, "1838");  /* 100km altitude */
        orbital_velocity(&v, &mu_moon, &r_moon);
        apf_test_pct("Low lunar orbit velocity", &v, "1.634", 0.01);
    }
}

void test_kepler_equation(void) {
    apf M, E, e, nu, r, a;
    
    printf("\n=== Kepler's Equation Tests ===\n");
    
    /* Circular orbit (e=0): E = M */
    apf_from_str(&M, "1.5");
    apf_from_str(&e, "0");
    eccentric_anomaly(&E, &M, &e);
    apf_test("E = M for circular orbit", &E, "1.5");
    
    /* Low eccentricity (e=0.1) - verify M = E - e*sin(E) */
    apf_from_str(&M, "1.0");
    apf_from_str(&e, "0.1");
    eccentric_anomaly(&E, &M, &e);
    apf_test_pct("E for e=0.1, M=1.0", &E, "1.0890", 0.01);  /* Verified correct */
    
    /* Verify by computing M back from E */
    {
        apf M_check;
        mean_anomaly(&M_check, &E, &e);
        apf_test_pct("M computed from E (e=0.1)", &M_check, "1.0", 0.01);
    }
    
    /* High eccentricity (e=0.9) - solve fresh */
    apf_from_str(&M, "1.0");
    apf_from_str(&e, "0.9");
    eccentric_anomaly(&E, &M, &e);  /* Re-solve for e=0.9 */
    {
        apf M_check;
        mean_anomaly(&M_check, &E, &e);
        apf_test_pct("M computed from E (e=0.9)", &M_check, "1.0", 0.01);
    }
    
    /* True anomaly from eccentric anomaly */
    apf_from_str(&E, "1.0");
    apf_from_str(&e, "0.5");
    true_anomaly(&nu, &E, &e);
    apf_test_pct("True anomaly for E=1.0, e=0.5", &nu, "1.5155", 0.01);  /* Corrected value */
    
    /* Radius from true anomaly */
    apf_from_str(&a, "10000");
    apf_from_str(&e, "0.5");
    apf_from_str(&nu, "0");  /* Periapsis */
    radius_from_anomaly(&r, &a, &e, &nu);
    apf_test("Periapsis radius (nu=0)", &r, "5000");
    
    {
        apf pi;
        apf_from_str(&pi, "3.14159265358979323846264338327950288419");
        radius_from_anomaly(&r, &a, &e, &pi);  /* Apoapsis */
        apf_test("Apoapsis radius (nu=pi)", &r, "15000");
    }
}
#endif

void test_special_values(void) {
    apf a, b, r;
    
    printf("\n=== Special Values and Edge Cases ===\n");
    
    /* Infinity handling */
    apf_set_inf(&a, 0);  /* +Inf */
    apf_from_int(&b, 1);
    apf_add(&r, &a, &b);
    if (r.cls == APF_CLASS_INF && r.sign == 0) {
        tests_passed++;
        printf("  PASS: Inf + 1 = Inf\n");
    } else {
        tests_failed++;
        printf("  FAIL: Inf + 1 = Inf\n");
    }
    
    /* NaN propagation */
    apf_set_nan(&a);
    apf_from_int(&b, 5);
    apf_mul(&r, &a, &b);
    if (r.cls == APF_CLASS_NAN) {
        tests_passed++;
        printf("  PASS: NaN * 5 = NaN\n");
    } else {
        tests_failed++;
        printf("  FAIL: NaN * 5 = NaN\n");
    }
    
    /* Zero handling */
    apf_zero(&a);
    apf_zero(&b);
    apf_div(&r, &a, &b);
    if (r.cls == APF_CLASS_NAN) {
        tests_passed++;
        printf("  PASS: 0/0 = NaN\n");
    } else {
        tests_failed++;
        printf("  FAIL: 0/0 = NaN\n");
    }
    
    apf_from_int(&b, 0);
    apf_from_int(&a, 1);
    apf_div(&r, &a, &b);
    if (r.cls == APF_CLASS_INF) {
        tests_passed++;
        printf("  PASS: 1/0 = Inf\n");
    } else {
        tests_failed++;
        printf("  FAIL: 1/0 = Inf\n");
    }
}

void test_stress(void) {
    apf a, r, two;
    int i;
    
    printf("\n=== Stress Tests ===\n");
    
    /* Large exponent */
    apf_from_int(&two, 2);
    apf_from_int(&a, 1000);
    apfx_pow(&r, &two, &a);
    apf_test_pct("2^1000", &r, "1.0715086071862673209484250490600018106e301", 1e-8);
    
    /* Many iterations (convergence test) */
    apf_from_int(&a, 2);
    for (i = 0; i < 100; i++) {
        apf_sqrt(&r, &a);
        apf_mul(&a, &r, &r);
    }
    apf_test("sqrt precision after 100 iterations", &a, "2");
    
    /* Accumulated error test */
    apf_from_str(&a, "0.1");
    for (i = 0; i < 10; i++) {
        apf_add(&r, &a, &a);
        apf_copy(&a, &r);
        apf_add(&r, &a, &a);
        apf_copy(&a, &r);
    }
    /* 0.1 * 2^20 = 104857.6 */
    apf_test("0.1 * 2^20 (accumulated adds)", &a, "104857.6");
}

void test_mission_critical(void) {
    printf("\n=== Mission-Critical Calculations ===\n");
    printf("  These tests verify calculations essential for lunar/Mars missions\n\n");
    
#ifdef HAVE_ORBITAL
    {
        apf mu, r1, r2, dv, v1, v2;
        
        /* LEO to GEO - classic Hohmann transfer */
        printf("  LEO to GEO Hohmann Transfer:\n");
        apf_from_str(&mu, "398600.4418");
        apf_from_str(&r1, "6678");   /* 300km LEO */
        apf_from_str(&r2, "42164");  /* GEO */
        hohmann_deltav(&dv, &mu, &r1, &r2);
        apf_test_pct("    Total delta-V", &dv, "3.935", 0.02);
        
        /* Orbital velocities for verification */
        orbital_velocity(&v1, &mu, &r1);
        orbital_velocity(&v2, &mu, &r2);
        printf("    LEO velocity: ");
        {
            char buf[64];
            apf_to_str(buf, sizeof(buf), &v1, 6);
            printf("%s km/s\n", buf);
        }
        printf("    GEO velocity: ");
        {
            char buf[64];
            apf_to_str(buf, sizeof(buf), &v2, 6);
            printf("%s km/s\n", buf);
        }
        
        /* Escape velocity from LEO */
        printf("\n  Earth Escape from LEO:\n");
        escape_velocity(&v1, &mu, &r1);
        {
            apf v_circ;
            orbital_velocity(&v_circ, &mu, &r1);
            apf_sub(&dv, &v1, &v_circ);
            apf_test_pct("    Delta-V above circular", &dv, "3.22", 0.02);
        }
        
        /* Moon escape velocity */
        printf("\n  Lunar Surface Operations:\n");
        {
            apf mu_moon, r_moon;
            apf_from_str(&mu_moon, "4902.8");   /* Moon GM */
            apf_from_str(&r_moon, "1737.4");    /* Moon radius */
            escape_velocity(&v1, &mu_moon, &r_moon);
            apf_test_pct("    Lunar escape velocity", &v1, "2.38", 0.02);
            
            orbital_velocity(&v2, &mu_moon, &r_moon);
            apf_test_pct("    Lunar surface orbital v", &v2, "1.68", 0.02);
        }
        
        /* Mars operations */
        printf("\n  Mars Surface Operations:\n");
        {
            apf mu_mars, r_mars;
            apf_from_str(&mu_mars, "42828.37");  /* Mars GM */
            apf_from_str(&r_mars, "3389.5");     /* Mars radius */
            escape_velocity(&v1, &mu_mars, &r_mars);
            apf_test_pct("    Mars escape velocity", &v1, "5.03", 0.02);
            
            orbital_velocity(&v2, &mu_mars, &r_mars);
            apf_test_pct("    Mars surface orbital v", &v2, "3.55", 0.02);
        }
    }
#else
    printf("  (Orbital module not enabled)\n");
#endif
}

int main(void) {
    printf("==============================================\n");
    printf("  scalc Comprehensive Test Suite\n");
    printf("  128-bit Precision Scientific Calculator\n");
    printf("  For Lunar/Mars Mission Calculations\n");
    printf("==============================================\n");
    
    test_apf_arithmetic();
    test_apf_precision();
    test_trig_functions();
    test_exp_log();
    
#ifdef HAVE_COMB
    test_combinatorics();
#endif

#ifdef HAVE_GCD
    test_number_theory();
#endif

#ifdef HAVE_BITWISE
    test_bitwise();
#endif

#ifdef HAVE_ORBITAL
    test_orbital_mechanics();
    test_kepler_equation();
#endif

    test_special_values();
    test_stress();
    test_mission_critical();
    
    printf("\n==============================================\n");
    printf("  RESULTS: %d passed, %d failed\n", tests_passed, tests_failed);
    printf("==============================================\n");
    
    if (tests_failed > 0) {
        printf("\n  WARNING: Some tests failed!\n");
        printf("  Calculator may not be suitable for mission-critical use.\n");
        return 1;
    } else {
        printf("\n  All tests passed!\n");
        printf("  Calculator verified for mission-critical calculations.\n");
        return 0;
    }
}
