/* matrix_rand.c - Random number generation for scalc
 * C89 compliant for Watcom C / DOS
 * 16-bit clean: int is 16-bit, long is 32-bit
 * Uses Linear Congruential Generator for portability
 */
#include <stdio.h>
#include <time.h>
#include "matrix.h"
#include "apfx.h"

/* ========== LCG Random Number Generator ========== */
/* Parameters from Numerical Recipes (16-bit safe) */
static unsigned long rand_state = 1;

void mat_rand_seed(unsigned long seed)
{
    rand_state = seed;
}

/* Generate random 32-bit value */
static unsigned long rand_next(void)
{
    /* LCG: X(n+1) = (a * X(n) + c) mod m
     * Using parameters that work well with 32-bit arithmetic */
    rand_state = rand_state * 1103515245UL + 12345UL;
    return rand_state;
}

/* Generate uniform random in [0, 1) */
static void rand_uniform_internal(apf *r)
{
    /* Generate random bits for precision */
    unsigned long r1;
    apf numerator, denominator;
    
    r1 = rand_next();
    
    /* Use r1, scaled by 2^31 */
    apf_from_int(&numerator, (long)(r1 >> 1));  /* Use top 31 bits */
    apf_from_int(&denominator, 0x7FFFFFFFL);     /* 2^31 - 1 */
    apf_div(r, &numerator, &denominator);
}

/* ========== Scalar Random Functions ========== */

void rand_uniform(apfc *r)
{
    rand_uniform_internal(&r->re);
    apf_zero(&r->im);
}

void rand_uniform_range(apfc *r, const apfc *low, const apfc *high)
{
    apf u, range, scaled;
    
    rand_uniform_internal(&u);
    
    /* r = low + u * (high - low) */
    apf_sub(&range, &high->re, &low->re);
    apf_mul(&scaled, &u, &range);
    apf_add(&r->re, &low->re, &scaled);
    apf_zero(&r->im);
}

/* Box-Muller transform for normal distribution */
void rand_normal(apfc *r)
{
    apf u1, u2, two, neg_two, log_u1, tmp, sqrt_part;
    apf two_pi, angle, cos_val;
    
    /* Generate two uniform random numbers in (0, 1) */
    do {
        rand_uniform_internal(&u1);
    } while (apf_is_zero(&u1));  /* Avoid log(0) */
    
    rand_uniform_internal(&u2);
    
    /* z = sqrt(-2 * ln(u1)) * cos(2 * pi * u2) */
    apfx_log(&log_u1, &u1);
    apf_from_int(&neg_two, -2);
    apf_mul(&tmp, &neg_two, &log_u1);
    apf_sqrt(&sqrt_part, &tmp);
    
    apfx_pi(&two_pi);
    apf_from_int(&two, 2);
    apf_mul(&two_pi, &two_pi, &two);
    apf_mul(&angle, &two_pi, &u2);
    
    apfx_cos(&cos_val, &angle);
    apf_mul(&r->re, &sqrt_part, &cos_val);
    apf_zero(&r->im);
}

/* Random integer from 1 to imax */
void rand_int(apfc *r, long imax)
{
    unsigned long rval;
    long result;
    
    if (imax < 1) {
        apf_from_int(&r->re, 1);
        apf_zero(&r->im);
        return;
    }
    
    rval = rand_next();
    result = (long)(rval % (unsigned long)imax) + 1;
    
    apf_from_int(&r->re, result);
    apf_zero(&r->im);
}

/* Random integer from imin to imax (inclusive) */
void rand_int_range(apfc *r, long imin, long imax)
{
    unsigned long rval;
    long range, result;
    
    if (imax < imin) {
        long tmp = imin;
        imin = imax;
        imax = tmp;
    }
    
    range = imax - imin + 1;
    if (range <= 0) {
        apf_from_int(&r->re, imin);
        apf_zero(&r->im);
        return;
    }
    
    rval = rand_next();
    result = imin + (long)(rval % (unsigned long)range);
    
    apf_from_int(&r->re, result);
    apf_zero(&r->im);
}

/* ========== Matrix Random Functions ========== */

/* rand(rows, cols) - uniform [0, 1) */
void mat_rand(matrix_t *r, int rows, int cols)
{
    int i;
    
    if (rows < 1) rows = 1;
    if (cols < 1) cols = 1;
    if (rows > MAT_MAX_ROWS) rows = MAT_MAX_ROWS;
    if (cols > MAT_MAX_COLS) cols = MAT_MAX_COLS;
    
    r->rows = rows;
    r->cols = cols;
    
    for (i = 0; i < rows * cols; i++) {
        rand_uniform(&r->data[i]);
    }
}

/* rand(rows, cols) with range [low, high) */
void mat_rand_range(matrix_t *r, int rows, int cols,
                    const apfc *low, const apfc *high)
{
    int i;
    
    if (rows < 1) rows = 1;
    if (cols < 1) cols = 1;
    if (rows > MAT_MAX_ROWS) rows = MAT_MAX_ROWS;
    if (cols > MAT_MAX_COLS) cols = MAT_MAX_COLS;
    
    r->rows = rows;
    r->cols = cols;
    
    for (i = 0; i < rows * cols; i++) {
        rand_uniform_range(&r->data[i], low, high);
    }
}

/* randn(rows, cols) - standard normal distribution */
void mat_randn(matrix_t *r, int rows, int cols)
{
    int i;
    
    if (rows < 1) rows = 1;
    if (cols < 1) cols = 1;
    if (rows > MAT_MAX_ROWS) rows = MAT_MAX_ROWS;
    if (cols > MAT_MAX_COLS) cols = MAT_MAX_COLS;
    
    r->rows = rows;
    r->cols = cols;
    
    for (i = 0; i < rows * cols; i++) {
        rand_normal(&r->data[i]);
    }
}

/* randi(imax, rows, cols) - random integers 1 to imax */
void mat_randi(matrix_t *r, int rows, int cols, long imax)
{
    int i;
    
    if (rows < 1) rows = 1;
    if (cols < 1) cols = 1;
    if (rows > MAT_MAX_ROWS) rows = MAT_MAX_ROWS;
    if (cols > MAT_MAX_COLS) cols = MAT_MAX_COLS;
    
    r->rows = rows;
    r->cols = cols;
    
    for (i = 0; i < rows * cols; i++) {
        rand_int(&r->data[i], imax);
    }
}

/* randi([imin, imax], rows, cols) - random integers imin to imax */
void mat_randi_range(matrix_t *r, int rows, int cols, long imin, long imax)
{
    int i;
    
    if (rows < 1) rows = 1;
    if (cols < 1) cols = 1;
    if (rows > MAT_MAX_ROWS) rows = MAT_MAX_ROWS;
    if (cols > MAT_MAX_COLS) cols = MAT_MAX_COLS;
    
    r->rows = rows;
    r->cols = cols;
    
    for (i = 0; i < rows * cols; i++) {
        rand_int_range(&r->data[i], imin, imax);
    }
}

/* Auto-seed from time (call once at startup) */
void rand_auto_seed(void)
{
    mat_rand_seed((unsigned long)time(NULL));
}
