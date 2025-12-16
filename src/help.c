/*
 * help.c - Function help, demo, and benchmark system
 * 
 * Provides wrappers that delegate to the func_dispatch.c documentation system.
 * C89 compliant.
 */
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>
#include "sc.h"
#include "help.h"
#include "func_defs.h"
#include "config.h"

/* Forward declarations */
extern const char *input_ptr;
extern void next_token(void);

/* Defined in main.c */
extern int eval_expr_line(const char *line, int quiet);

/* ========== Help System ========== */

void show_function_help(const char *name)
{
    func_show_help(name);
}

int function_exists(const char *name)
{
    return func_doc_exists(name);
}

int count_help_entries(void)
{
    return func_doc_count();
}

void list_all_functions(void)
{
    func_list_all();
}

/* ========== Demo System ========== */

void show_function_demo(const char *name)
{
    const FuncDoc *doc = func_find_doc(name);
    int i;
    
    if (doc == NULL) {
        printf("No demo available for '%s'\n", name);
        printf("Type 'functions' to list all available functions.\n");
        return;
    }
    
    printf("\n=== Demo: %s ===\n", doc->name);
    printf("Description: %s\n\n", doc->description);
    
    for (i = 0; i < MAX_EXAMPLES && doc->examples[i] != NULL; i++) {
        const char *ex = doc->examples[i];
        const char *eq = strchr(ex, '=');
        char expr[256];
        
        /* Extract just the expression (before =) */
        if (eq != NULL) {
            int len = (int)(eq - ex);
            if (len > 0 && len < 255) {
                strncpy(expr, ex, len);
                expr[len] = '\0';
                /* Trim trailing space */
                while (len > 0 && expr[len-1] == ' ') {
                    expr[--len] = '\0';
                }
            } else {
                strncpy(expr, ex, 255);
                expr[255] = '\0';
            }
        } else {
            strncpy(expr, ex, 255);
            expr[255] = '\0';
        }
        
        printf(">>> %s\n", expr);
        eval_expr_line(expr, 0);
        printf("\n");
    }
}

/* ========== Benchmark System ========== */

/* Benchmark a single expression */
static double bench_expr(const char *expr, int iterations)
{
    clock_t start, end;
    int i;
    
    start = clock();
    for (i = 0; i < iterations; i++) {
        eval_expr_line(expr, 1);
    }
    end = clock();
    
    return (double)(end - start) / CLOCKS_PER_SEC;
}

void run_function_bench(void)
{
    int iterations = 10000;
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
    
    /* Exponential/logarithmic */
    printf("\nExponential/Logarithmic:\n");
    elapsed = bench_expr("exp(2)", iterations);
    printf("  exp(2):       %.4f sec (%.0f ops/sec)\n", elapsed, iterations/elapsed);
    
    elapsed = bench_expr("ln(10)", iterations);
    printf("  ln(10):       %.4f sec (%.0f ops/sec)\n", elapsed, iterations/elapsed);
    
    elapsed = bench_expr("sqrt(2)", iterations);
    printf("  sqrt(2):      %.4f sec (%.0f ops/sec)\n", elapsed, iterations/elapsed);
    
    elapsed = bench_expr("2^10", iterations);
    printf("  2^10:         %.4f sec (%.0f ops/sec)\n", elapsed, iterations/elapsed);
    
    /* Special functions */
    printf("\nSpecial Functions:\n");
    elapsed = bench_expr("gamma(5)", iterations);
    printf("  gamma(5):     %.4f sec (%.0f ops/sec)\n", elapsed, iterations/elapsed);
    
    elapsed = bench_expr("erf(0.5)", iterations);
    printf("  erf(0.5):     %.4f sec (%.0f ops/sec)\n", elapsed, iterations/elapsed);
    
    /* Number theory */
    printf("\nNumber Theory:\n");
    elapsed = bench_expr("gcd(1234, 5678)", iterations);
    printf("  gcd(1234,5678): %.4f sec (%.0f ops/sec)\n", elapsed, iterations/elapsed);
    
    elapsed = bench_expr("isprime(997)", iterations);
    printf("  isprime(997): %.4f sec (%.0f ops/sec)\n", elapsed, iterations/elapsed);
    
    /* Arithmetic */
    printf("\nBasic Arithmetic:\n");
    elapsed = bench_expr("1+2", iterations * 10);
    printf("  1+2:          %.4f sec (%.0f ops/sec)\n", elapsed, iterations*10/elapsed);
    
    elapsed = bench_expr("3*4", iterations * 10);
    printf("  3*4:          %.4f sec (%.0f ops/sec)\n", elapsed, iterations*10/elapsed);
    
    elapsed = bench_expr("10/3", iterations * 10);
    printf("  10/3:         %.4f sec (%.0f ops/sec)\n", elapsed, iterations*10/elapsed);
    
    /* Complex operations */
    printf("\nComplex Operations:\n");
    elapsed = bench_expr("(3+4i)*(1+2i)", iterations);
    printf("  (3+4i)*(1+2i): %.4f sec (%.0f ops/sec)\n", elapsed, iterations/elapsed);
    
    elapsed = bench_expr("abs(3+4i)", iterations);
    printf("  abs(3+4i):    %.4f sec (%.0f ops/sec)\n", elapsed, iterations/elapsed);
    
    printf("\n");
}
