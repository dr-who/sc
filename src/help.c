/*
 * help.c - Function help, demo, and benchmark system
 * 
 * Uses func_registry.c for all function documentation.
 * C89 compliant.
 */
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>
#include "sc.h"
#include "help.h"
#include "func_registry.h"
#include "config.h"

/* Defined in main.c */
extern int eval_expr_line(const char *line, int quiet);

/* ========== Help System ========== */

void show_function_help(const char *name)
{
    func_help(name);
}

int function_exists(const char *name)
{
    return func_exists(name);
}

int count_help_entries(void)
{
    return func_count();
}

void list_all_functions(void)
{
    func_list();
}

/* ========== Demo System ========== */

void show_function_demo(const char *name)
{
    func_demo(name);
}

/* ========== Benchmark System ========== */

/* Benchmark a single expression */
static double bench_expr(const char *expr, int iterations)
{
    clock_t start, end;
    int i;
    extern void mat_arena_reset(void);
    
    start = clock();
    for (i = 0; i < iterations; i++) {
        mat_arena_reset();  /* Reset arena each iteration */
        eval_expr_line(expr, 2);  /* silent=2 suppresses all output */
    }
    end = clock();
    
    return (double)(end - start) / CLOCKS_PER_SEC;
}

void run_function_bench(void)
{
    int iterations = 5000;
    double elapsed;
    
    printf("\nBenchmarking core functions (%d iterations each):\n", iterations);
    printf("================================================\n\n");
    
    /* Trigonometric */
    printf("Trigonometric:\n");
    elapsed = bench_expr("sin(0.5)", iterations);
    printf("  sin(0.5):     %.4f sec (%.0f ops/sec)\n", elapsed, iterations/elapsed);
    
    elapsed = bench_expr("cos(0.5)", iterations);
    printf("  cos(0.5):     %.4f sec (%.0f ops/sec)\n", elapsed, iterations/elapsed);
    
    elapsed = bench_expr("tan(0.5)", iterations);
    printf("  tan(0.5):     %.4f sec (%.0f ops/sec)\n", elapsed, iterations/elapsed);
    
    elapsed = bench_expr("atan(0.5)", iterations);
    printf("  atan(0.5):    %.4f sec (%.0f ops/sec)\n", elapsed, iterations/elapsed);
    
    elapsed = bench_expr("sinh(1)", iterations);
    printf("  sinh(1):      %.4f sec (%.0f ops/sec)\n", elapsed, iterations/elapsed);
    
    elapsed = bench_expr("sind(45)", iterations);
    printf("  sind(45):     %.4f sec (%.0f ops/sec)\n", elapsed, iterations/elapsed);
    
    /* Exponential/logarithmic */
    printf("\nExponential/Logarithmic:\n");
    elapsed = bench_expr("exp(2)", iterations);
    printf("  exp(2):       %.4f sec (%.0f ops/sec)\n", elapsed, iterations/elapsed);
    
    elapsed = bench_expr("ln(10)", iterations);
    printf("  ln(10):       %.4f sec (%.0f ops/sec)\n", elapsed, iterations/elapsed);
    
    elapsed = bench_expr("sqrt(2)", iterations);
    printf("  sqrt(2):      %.4f sec (%.0f ops/sec)\n", elapsed, iterations/elapsed);
    
    elapsed = bench_expr("cbrt(27)", iterations);
    printf("  cbrt(27):     %.4f sec (%.0f ops/sec)\n", elapsed, iterations/elapsed);
    
    elapsed = bench_expr("2^10", iterations);
    printf("  2^10:         %.4f sec (%.0f ops/sec)\n", elapsed, iterations/elapsed);
    
    elapsed = bench_expr("log10(1000)", iterations);
    printf("  log10(1000):  %.4f sec (%.0f ops/sec)\n", elapsed, iterations/elapsed);
    
    /* Special functions */
    printf("\nSpecial Functions:\n");
    elapsed = bench_expr("gamma(5)", iterations);
    printf("  gamma(5):     %.4f sec (%.0f ops/sec)\n", elapsed, iterations/elapsed);
    
    elapsed = bench_expr("lgamma(10)", iterations);
    printf("  lgamma(10):   %.4f sec (%.0f ops/sec)\n", elapsed, iterations/elapsed);
    
    elapsed = bench_expr("erf(0.5)", iterations);
    printf("  erf(0.5):     %.4f sec (%.0f ops/sec)\n", elapsed, iterations/elapsed);
    
    elapsed = bench_expr("erfc(1)", iterations);
    printf("  erfc(1):      %.4f sec (%.0f ops/sec)\n", elapsed, iterations/elapsed);
    
    elapsed = bench_expr("beta(2,3)", iterations);
    printf("  beta(2,3):    %.4f sec (%.0f ops/sec)\n", elapsed, iterations/elapsed);
    
    elapsed = bench_expr("zeta(2)", iterations);
    printf("  zeta(2):      %.4f sec (%.0f ops/sec)\n", elapsed, iterations/elapsed);
    
    /* Number theory */
    printf("\nNumber Theory:\n");
    elapsed = bench_expr("gcd(1234, 5678)", iterations);
    printf("  gcd(1234,5678): %.4f sec (%.0f ops/sec)\n", elapsed, iterations/elapsed);
    
    elapsed = bench_expr("lcm(12, 18)", iterations);
    printf("  lcm(12,18):   %.4f sec (%.0f ops/sec)\n", elapsed, iterations/elapsed);
    
    elapsed = bench_expr("isprime(997)", iterations);
    printf("  isprime(997): %.4f sec (%.0f ops/sec)\n", elapsed, iterations/elapsed);
    
    elapsed = bench_expr("fact(20)", iterations);
    printf("  fact(20):     %.4f sec (%.0f ops/sec)\n", elapsed, iterations/elapsed);
    
    elapsed = bench_expr("fib(30)", iterations);
    printf("  fib(30):      %.4f sec (%.0f ops/sec)\n", elapsed, iterations/elapsed);
    
    elapsed = bench_expr("ncr(20,10)", iterations);
    printf("  ncr(20,10):   %.4f sec (%.0f ops/sec)\n", elapsed, iterations/elapsed);
    
    /* Rounding */
    printf("\nRounding:\n");
    elapsed = bench_expr("floor(3.7)", iterations * 2);
    printf("  floor(3.7):   %.4f sec (%.0f ops/sec)\n", elapsed, iterations*2/elapsed);
    
    elapsed = bench_expr("ceil(3.2)", iterations * 2);
    printf("  ceil(3.2):    %.4f sec (%.0f ops/sec)\n", elapsed, iterations*2/elapsed);
    
    elapsed = bench_expr("round(3.5)", iterations * 2);
    printf("  round(3.5):   %.4f sec (%.0f ops/sec)\n", elapsed, iterations*2/elapsed);
    
    elapsed = bench_expr("mod(17, 5)", iterations * 2);
    printf("  mod(17,5):    %.4f sec (%.0f ops/sec)\n", elapsed, iterations*2/elapsed);
    
    /* Arithmetic */
    printf("\nBasic Arithmetic:\n");
    elapsed = bench_expr("1+2", iterations * 10);
    printf("  1+2:          %.4f sec (%.0f ops/sec)\n", elapsed, iterations*10/elapsed);
    
    elapsed = bench_expr("3*4", iterations * 10);
    printf("  3*4:          %.4f sec (%.0f ops/sec)\n", elapsed, iterations*10/elapsed);
    
    elapsed = bench_expr("10/3", iterations * 10);
    printf("  10/3:         %.4f sec (%.0f ops/sec)\n", elapsed, iterations*10/elapsed);
    
    elapsed = bench_expr("5!", iterations * 5);
    printf("  5!:           %.4f sec (%.0f ops/sec)\n", elapsed, iterations*5/elapsed);
    
    /* Complex operations */
    printf("\nComplex Operations:\n");
    elapsed = bench_expr("(3+4i)*(1+2i)", iterations);
    printf("  (3+4i)*(1+2i): %.4f sec (%.0f ops/sec)\n", elapsed, iterations/elapsed);
    
    elapsed = bench_expr("abs(3+4i)", iterations);
    printf("  abs(3+4i):    %.4f sec (%.0f ops/sec)\n", elapsed, iterations/elapsed);
    
    elapsed = bench_expr("conj(3+4i)", iterations);
    printf("  conj(3+4i):   %.4f sec (%.0f ops/sec)\n", elapsed, iterations/elapsed);
    
    elapsed = bench_expr("arg(1+i)", iterations);
    printf("  arg(1+i):     %.4f sec (%.0f ops/sec)\n", elapsed, iterations/elapsed);
    
    /* Statistics */
    printf("\nStatistics:\n");
    elapsed = bench_expr("sum([1,2,3,4,5])", iterations);
    printf("  sum([1..5]):  %.4f sec (%.0f ops/sec)\n", elapsed, iterations/elapsed);
    
    elapsed = bench_expr("mean([1,2,3,4,5])", iterations);
    printf("  mean([1..5]): %.4f sec (%.0f ops/sec)\n", elapsed, iterations/elapsed);
    
    elapsed = bench_expr("std([1,2,3,4,5])", iterations);
    printf("  std([1..5]):  %.4f sec (%.0f ops/sec)\n", elapsed, iterations/elapsed);
    
    elapsed = bench_expr("median([1,2,3,4,5])", iterations);
    printf("  median([1..5]): %.4f sec (%.0f ops/sec)\n", elapsed, iterations/elapsed);
    
    /* Matrix operations */
    printf("\nMatrix Operations:\n");
    elapsed = bench_expr("det([1,2;3,4])", iterations);
    printf("  det(2x2):     %.4f sec (%.0f ops/sec)\n", elapsed, iterations/elapsed);
    
    elapsed = bench_expr("inv([1,2;3,4])", iterations);
    printf("  inv(2x2):     %.4f sec (%.0f ops/sec)\n", elapsed, iterations/elapsed);
    
    elapsed = bench_expr("trace([1,2;3,4])", iterations);
    printf("  trace(2x2):   %.4f sec (%.0f ops/sec)\n", elapsed, iterations/elapsed);
    
    elapsed = bench_expr("transpose([1,2;3,4])", iterations);
    printf("  transpose(2x2): %.4f sec (%.0f ops/sec)\n", elapsed, iterations/elapsed);
    
    elapsed = bench_expr("[1,2;3,4]*[5,6;7,8]", iterations);
    printf("  2x2 * 2x2:    %.4f sec (%.0f ops/sec)\n", elapsed, iterations/elapsed);
    
    elapsed = bench_expr("qr([1,2;3,4])", iterations / 2);
    printf("  qr(2x2):      %.4f sec (%.0f ops/sec)\n", elapsed, (iterations/2)/elapsed);
    
    /* Probability */
    printf("\nProbability:\n");
    elapsed = bench_expr("normcdf(1.96)", iterations);
    printf("  normcdf(1.96): %.4f sec (%.0f ops/sec)\n", elapsed, iterations/elapsed);
    
    elapsed = bench_expr("norminv(0.975)", iterations);
    printf("  norminv(0.975): %.4f sec (%.0f ops/sec)\n", elapsed, iterations/elapsed);
    
    /* Logical */
    printf("\nLogical:\n");
    elapsed = bench_expr("isnan(0/0)", iterations * 2);
    printf("  isnan(0/0):   %.4f sec (%.0f ops/sec)\n", elapsed, iterations*2/elapsed);
    
    elapsed = bench_expr("isinf(1/0)", iterations * 2);
    printf("  isinf(1/0):   %.4f sec (%.0f ops/sec)\n", elapsed, iterations*2/elapsed);
    
    elapsed = bench_expr("5>3", iterations * 2);
    printf("  5>3:          %.4f sec (%.0f ops/sec)\n", elapsed, iterations*2/elapsed);
    
    /* Angle conversion */
    printf("\nAngle Conversion:\n");
    elapsed = bench_expr("deg2rad(180)", iterations * 2);
    printf("  deg2rad(180): %.4f sec (%.0f ops/sec)\n", elapsed, iterations*2/elapsed);
    
    elapsed = bench_expr("rad2deg(pi)", iterations * 2);
    printf("  rad2deg(pi):  %.4f sec (%.0f ops/sec)\n", elapsed, iterations*2/elapsed);
    
    /* Bitwise */
    printf("\nBitwise:\n");
    elapsed = bench_expr("bitand(255, 170)", iterations * 2);
    printf("  bitand:       %.4f sec (%.0f ops/sec)\n", elapsed, iterations*2/elapsed);
    
    elapsed = bench_expr("bitxor(255, 170)", iterations * 2);
    printf("  bitxor:       %.4f sec (%.0f ops/sec)\n", elapsed, iterations*2/elapsed);
    
    printf("\n");
}

/* === SaaS Metrics === */

