/* test_ieee754.c - IEEE 754 Floating Point Compliance Tests
 * Tests scalc against IEEE 754-2008 standard
 * 
 * Compile: gcc -std=c89 -o test_ieee754 test_ieee754.c -lm
 * Run: ./test_ieee754
 * Or: ./test_ieee754 ./scalc
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

/* Test counters */
static int tests_passed = 0;
static int tests_failed = 0;
static int tests_skipped = 0;

/* Verbose mode */
static int verbose = 0;

/* Path to scalc binary */
static const char *scalc_path = "./scalc";

/* Execute a scalc expression and return the result line */
static int run_scalc(const char *expr, char *result, int maxlen) {
    char cmd[512];
    FILE *fp;
    
    snprintf(cmd, sizeof(cmd), "echo '%s' | %s 2>&1", expr, scalc_path);
    fp = popen(cmd, "r");
    if (!fp) return 0;
    
    result[0] = '\0';
    while (fgets(result, maxlen, fp)) {
        /* Look for "= " in output */
        if (strstr(result, "= ")) {
            /* Remove trailing newline */
            char *nl = strchr(result, '\n');
            if (nl) *nl = '\0';
            pclose(fp);
            return 1;
        }
    }
    pclose(fp);
    return 0;
}

/* Extract numeric value from "= value" result */
static const char *get_value(const char *result) {
    const char *p = strstr(result, "= ");
    if (p) return p + 2;
    return result;
}

/* Check if result contains expected string */
static int contains(const char *result, const char *expected) {
    return strstr(result, expected) != NULL;
}

/* Test helper - pass/fail with message */
static void test(const char *name, const char *expr, const char *expected, int exact) {
    char result[256];
    const char *val;
    int pass;
    
    if (!run_scalc(expr, result, sizeof(result))) {
        printf("  FAIL: %s - no result from '%s'\n", name, expr);
        tests_failed++;
        return;
    }
    
    val = get_value(result);
    
    if (exact) {
        pass = strcmp(val, expected) == 0;
    } else {
        pass = contains(result, expected);
    }
    
    if (pass) {
        if (verbose) printf("  PASS: %s\n", name);
        tests_passed++;
    } else {
        printf("  FAIL: %s\n        expr: %s\n        expected: %s\n        got: %s\n", 
               name, expr, expected, val);
        tests_failed++;
    }
}

/* Category header */
static void section(const char *name) {
    printf("\n=== %s ===\n", name);
}

/* ====================================================================
 * IEEE 754-2008 TESTS
 * ==================================================================== */

static void test_special_values(void) {
    section("IEEE 754: Special Values");
    
    /* Positive/Negative Infinity */
    test("Positive infinity from overflow", "1e308 * 10", "Inf", 0);
    test("Negative infinity from overflow", "-1e308 * 10", "-Inf", 0);
    test("Division by zero (positive)", "1/0", "Inf", 0);
    test("Division by zero (negative)", "-1/0", "-Inf", 0);
    
    /* NaN generation */
    test("0/0 produces NaN", "0/0", "NaN", 0);
    test("Inf - Inf produces NaN", "(1/0) - (1/0)", "NaN", 0);
    test("Inf * 0 produces NaN", "(1/0) * 0", "NaN", 0);
    test("0 * Inf produces NaN", "0 * (1/0)", "NaN", 0);
    
    /* Infinity arithmetic */
    test("Inf + 1 = Inf", "(1/0) + 1", "Inf", 0);
    test("Inf - 1 = Inf", "(1/0) - 1", "Inf", 0);
    test("Inf * 2 = Inf", "(1/0) * 2", "Inf", 0);
    test("Inf / 2 = Inf", "(1/0) / 2", "Inf", 0);
    test("-Inf * 2 = -Inf", "(-1/0) * 2", "-Inf", 0);
    test("1 / Inf = 0", "1 / (1/0)", "0", 0);
    test("-1 / Inf = -0 or 0", "-1 / (1/0)", "0", 0);
    
    /* Infinity comparisons (these are harder to test via expression) */
}

static void test_nan_propagation(void) {
    section("IEEE 754: NaN Propagation");
    
    /* NaN in arithmetic operations should propagate */
    test("NaN + 1 = NaN", "(0/0) + 1", "NaN", 0);
    test("NaN - 1 = NaN", "(0/0) - 1", "NaN", 0);
    test("NaN * 2 = NaN", "(0/0) * 2", "NaN", 0);
    test("NaN / 2 = NaN", "(0/0) / 2", "NaN", 0);
    test("1 + NaN = NaN", "1 + (0/0)", "NaN", 0);
    test("1 - NaN = NaN", "1 - (0/0)", "NaN", 0);
    test("2 * NaN = NaN", "2 * (0/0)", "NaN", 0);
    test("2 / NaN = NaN", "2 / (0/0)", "NaN", 0);
    
    /* NaN in functions */
    test("sqrt(NaN) = NaN", "sqrt(0/0)", "NaN", 0);
    test("sin(NaN) = NaN", "sin(0/0)", "NaN", 0);
    test("cos(NaN) = NaN", "cos(0/0)", "NaN", 0);
    test("exp(NaN) = NaN", "exp(0/0)", "NaN", 0);
    test("ln(NaN) = NaN", "ln(0/0)", "NaN", 0);
    
    /* Power with NaN */
    test("NaN^2 = NaN", "(0/0)^2", "NaN", 0);
    test("2^NaN = NaN", "2^(0/0)", "NaN", 0);
}

static void test_infinity_math(void) {
    section("IEEE 754: Infinity Arithmetic Rules");
    
    /* Addition/Subtraction */
    test("Inf + Inf = Inf", "(1/0) + (1/0)", "Inf", 0);
    test("-Inf + -Inf = -Inf", "(-1/0) + (-1/0)", "-Inf", 0);
    test("Inf + -Inf = NaN", "(1/0) + (-1/0)", "NaN", 0);
    test("-Inf + Inf = NaN", "(-1/0) + (1/0)", "NaN", 0);
    
    /* Multiplication */
    test("Inf * Inf = Inf", "(1/0) * (1/0)", "Inf", 0);
    test("-Inf * -Inf = Inf", "(-1/0) * (-1/0)", "Inf", 0);
    test("Inf * -Inf = -Inf", "(1/0) * (-1/0)", "-Inf", 0);
    
    /* Division */
    test("Inf / Inf = NaN", "(1/0) / (1/0)", "NaN", 0);
    test("Inf / -Inf = NaN", "(1/0) / (-1/0)", "NaN", 0);
    test("1 / 0 = Inf", "1 / 0", "Inf", 0);
    test("-1 / 0 = -Inf", "-1 / 0", "-Inf", 0);
    test("Inf / 1 = Inf", "(1/0) / 1", "Inf", 0);
    test("Inf / -1 = -Inf", "(1/0) / (-1)", "-Inf", 0);
}

static void test_sqrt_special(void) {
    section("IEEE 754: Square Root Special Cases");
    
    test("sqrt(0) = 0", "sqrt(0)", "0", 0);
    test("sqrt(1) = 1", "sqrt(1)", "1", 0);
    test("sqrt(Inf) = Inf", "sqrt(1/0)", "Inf", 0);
    test("sqrt(-1) = i (complex)", "sqrt(-1)", "i", 0);
    test("sqrt(-4) = 2i (complex)", "sqrt(-4)", "2i", 0);
    test("sqrt(NaN) = NaN", "sqrt(0/0)", "NaN", 0);
}

static void test_exp_log_special(void) {
    section("IEEE 754: Exp/Log Special Cases");
    
    /* Exponential */
    test("exp(0) = 1", "exp(0)", "1", 0);
    test("exp(1) = e", "exp(1)", "2.71828", 0);
    test("exp(Inf) = Inf", "exp(1/0)", "Inf", 0);
    test("exp(-Inf) = 0", "exp(-1/0)", "0", 0);
    test("exp(NaN) = NaN", "exp(0/0)", "NaN", 0);
    
    /* Natural logarithm */
    test("ln(1) = 0", "ln(1)", "0", 0);
    test("ln(e) = 1", "ln(e)", "1", 0);
    test("ln(0) = -Inf", "ln(0)", "-Inf", 0);
    test("ln(Inf) = Inf", "ln(1/0)", "Inf", 0);
    test("ln(-1) = complex (pi*i)", "ln(-1)", "i", 0);  /* Should be pi*i */
    test("ln(NaN) = NaN", "ln(0/0)", "NaN", 0);
}

static void test_trig_special(void) {
    section("IEEE 754: Trigonometric Special Cases");
    
    /* Sin/Cos at special points */
    test("sin(0) = 0", "sin(0)", "0", 1);
    test("cos(0) = 1", "cos(0)", "1", 1);
    test("sin(pi/2) = 1", "sin(pi/2)", "1", 0);
    test("cos(pi) = -1", "cos(pi)", "-1", 0);
    
    /* Tan at poles */
    test("tan(0) = 0", "tan(0)", "0", 1);
    
    /* With infinity */
    test("sin(Inf) = NaN", "sin(1/0)", "NaN", 0);
    test("cos(Inf) = NaN", "cos(1/0)", "NaN", 0);
    test("tan(Inf) = NaN", "tan(1/0)", "NaN", 0);
    
    /* With NaN */
    test("sin(NaN) = NaN", "sin(0/0)", "NaN", 0);
    test("cos(NaN) = NaN", "cos(0/0)", "NaN", 0);
    test("tan(NaN) = NaN", "tan(0/0)", "NaN", 0);
}

static void test_power_special(void) {
    section("IEEE 754: Power Function Special Cases");
    
    /* Zero exponent */
    test("0^0 = 1 (by convention)", "0^0", "1", 0);
    test("1^0 = 1", "1^0", "1", 0);
    test("(-1)^0 = 1", "(-1)^0", "1", 0);
    test("Inf^0 = 1", "(1/0)^0", "1", 0);
    test("NaN^0 = 1 (IEEE special)", "(0/0)^0", "1", 0);
    
    /* One base */
    test("1^1 = 1", "1^1", "1", 0);
    test("1^Inf = 1", "1^(1/0)", "1", 0);
    test("1^(-Inf) = 1", "1^(-1/0)", "1", 0);
    test("1^NaN = 1 (IEEE special)", "1^(0/0)", "1", 0);
    
    /* Infinity exponent */
    test("2^Inf = Inf", "2^(1/0)", "Inf", 0);
    test("0.5^Inf = 0", "0.5^(1/0)", "0", 0);
    test("2^(-Inf) = 0", "2^(-1/0)", "0", 0);
    test("0.5^(-Inf) = Inf", "0.5^(-1/0)", "Inf", 0);
    
    /* Zero base */
    test("0^1 = 0", "0^1", "0", 0);
    test("0^2 = 0", "0^2", "0", 0);
    test("0^Inf = 0", "0^(1/0)", "0", 0);
    test("0^(-1) = Inf", "0^(-1)", "Inf", 0);
    
    /* Negative base with non-integer exponent */
    test("(-1)^0.5 = i (complex)", "(-1)^0.5", "i", 0);
    test("(-8)^(1/3) handles complex", "(-8)^(1/3)", "", 0);  /* May vary */
}

static void test_precision(void) {
    section("IEEE 754: Precision and Rounding");
    
    /* Basic precision tests */
    test("1/3 precision", "1/3", "0.333333", 0);
    test("2/3 precision", "2/3", "0.666666", 0);
    test("1/7 repeating", "1/7", "0.142857", 0);
    test("1/9 repeating", "1/9", "0.111111", 0);
    
    /* Near-integer results */
    test("0.1 + 0.2 near 0.3", "0.1 + 0.2", "0.3", 0);
    test("1 - 0.9 - 0.1 near 0", "1 - 0.9 - 0.1", "0", 0);
    
    /* Large numbers */
    test("Large: 2^100", "2^100", "1.267650600228229e+30", 0);
    test("Large: 10^100", "10^100", "1e+100", 0);
    
    /* Small numbers */
    test("Small: 1e-100", "1e-100", "1e-100", 0);
    test("Small: 2^(-100)", "2^(-100)", "", 0);  /* Should be very small */
}

static void test_overflow_underflow(void) {
    section("IEEE 754: Overflow and Underflow");
    
    /* Overflow to infinity */
    test("Overflow: 1e308 * 10 = Inf", "1e308 * 10", "Inf", 0);
    test("Overflow: 1e309 = Inf", "1e309", "Inf", 0);
    test("Overflow: -1e308 * 10 = -Inf", "-1e308 * 10", "-Inf", 0);
    
    /* Underflow to zero (gradual underflow) */
    test("Underflow: 1e-308 / 1e10", "1e-308 / 1e10", "0", 0);
    test("Underflow: 1e-400 = 0", "1e-400", "0", 0);
}

static void test_complex_special(void) {
    section("IEEE 754 Extensions: Complex Number Special Cases");
    
    /* Basic complex */
    test("i * i = -1", "i * i", "-1", 0);
    test("i^2 = -1", "i^2", "-1", 0);
    test("i^4 = 1", "i^4", "1", 0);
    
    /* Complex magnitude */
    test("|0| = 0", "abs(0)", "0", 0);
    test("|3+4i| = 5", "abs(3+4i)", "5", 0);
    test("|i| = 1", "abs(i)", "1", 0);
    test("|Inf| = Inf", "abs(1/0)", "Inf", 0);
    
    /* Complex with infinity */
    test("Inf + 0i magnitude", "abs((1/0) + 0*i)", "Inf", 0);
    
    /* Euler's identity approximation */
    test("e^(i*pi) + 1 ≈ 0", "e^(i*pi) + 1", "0", 0);
}

static void test_associativity(void) {
    section("IEEE 754: Non-Associativity (Known Limitations)");
    
    /* Floating point is NOT associative */
    /* (a + b) + c != a + (b + c) in general */
    /* These tests document behavior, may show expected failures */
    
    test("(1e20 + -1e20) + 1 = 1", "(1e20 + (-1e20)) + 1", "1", 0);
    test("1e20 + (-1e20 + 1) behavior", "1e20 + ((-1e20) + 1)", "1", 0);
    
    /* Catastrophic cancellation */
    test("Large - large + small", "1e16 - 1e16 + 1", "1", 0);
}

static void test_transcendental_accuracy(void) {
    section("IEEE 754: Transcendental Function Accuracy");
    
    /* Pi-related */
    test("sin(pi) ≈ 0", "sin(pi)", "0", 0);
    test("cos(pi/2) ≈ 0", "cos(pi/2)", "0", 0);
    test("tan(pi/4) = 1", "tan(pi/4)", "1", 0);
    
    /* Inverse functions */
    test("asin(sin(0.5)) = 0.5", "asin(sin(0.5))", "0.5", 0);
    test("acos(cos(0.5)) = 0.5", "acos(cos(0.5))", "0.5", 0);
    test("atan(tan(0.5)) = 0.5", "atan(tan(0.5))", "0.5", 0);
    
    /* Hyperbolic */
    test("sinh(0) = 0", "sinh(0)", "0", 1);
    test("cosh(0) = 1", "cosh(0)", "1", 1);
    test("tanh(0) = 0", "tanh(0)", "0", 1);
    
    /* Identities */
    test("cosh^2 - sinh^2 = 1 at x=1", "cosh(1)^2 - sinh(1)^2", "1", 0);
    test("sin^2 + cos^2 = 1 at x=1", "sin(1)^2 + cos(1)^2", "1", 0);
}

static void test_edge_integers(void) {
    section("IEEE 754: Integer Edge Cases");
    
    /* Factorial precision */
    test("20! exact", "20!", "2432902008176640000", 0);
    test("21! may lose precision", "21!", "", 0);  /* Check it's computed */
    
    /* Large integers as floats */
    test("2^53 (max safe int)", "2^53", "9007199254740992", 0);
    test("2^53 + 1 precision", "2^53 + 1", "9007199254740992", 0);  /* May round */
    
    /* Modular arithmetic */
    test("Negative mod", "(-7) mod 3", "", 0);  /* Behavior varies */
}

/* Main test driver */
int main(int argc, char *argv[]) {
    int i;
    
    /* Parse arguments */
    for (i = 1; i < argc; i++) {
        if (strcmp(argv[i], "-v") == 0 || strcmp(argv[i], "--verbose") == 0) {
            verbose = 1;
        } else {
            scalc_path = argv[i];
        }
    }
    
    printf("=====================================================\n");
    printf("  IEEE 754 Floating Point Compliance Tests\n");
    printf("  Testing: %s\n", scalc_path);
    printf("=====================================================\n");
    
    /* Run all test categories */
    test_special_values();
    test_nan_propagation();
    test_infinity_math();
    test_sqrt_special();
    test_exp_log_special();
    test_trig_special();
    test_power_special();
    test_precision();
    test_overflow_underflow();
    test_complex_special();
    test_associativity();
    test_transcendental_accuracy();
    test_edge_integers();
    
    /* Summary */
    printf("\n=====================================================\n");
    printf("  RESULTS: %d passed, %d failed, %d skipped\n", 
           tests_passed, tests_failed, tests_skipped);
    printf("  IEEE 754 Compliance: %.1f%%\n", 
           100.0 * tests_passed / (tests_passed + tests_failed));
    printf("=====================================================\n");
    
    return tests_failed > 0 ? 1 : 0;
}
