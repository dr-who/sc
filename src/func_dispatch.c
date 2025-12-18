/*
 * func_dispatch.c - Unified function dispatch implementation
 * 
 * C89 compliant for Watcom C / DOS
 */
#include <stdio.h>
#include <string.h>
#include "func_dispatch.h"
#include "apfx.h"

/* Hash table size - prime number for better distribution */
#define HASH_TABLE_SIZE 521

/* Hash table */
static FuncDispatchEntry *hash_table[HASH_TABLE_SIZE];

/* Static storage for dispatch entries */
#define MAX_DISPATCH_ENTRIES 600
static FuncDispatchEntry dispatch_entries[MAX_DISPATCH_ENTRIES];
static int num_entries = 0;

/* Case-insensitive string comparison */
static int str_eq_ci(const char *a, const char *b) {
    while (*a && *b) {
        char ca = *a, cb = *b;
        if (ca >= 'A' && ca <= 'Z') ca += 32;
        if (cb >= 'A' && cb <= 'Z') cb += 32;
        if (ca != cb) return 0;
        a++; b++;
    }
    return *a == *b;
}

/* Hash function for function names */
static unsigned int hash_name(const char *name) {
    unsigned int h = 0;
    while (*name) {
        char c = *name;
        if (c >= 'A' && c <= 'Z') c += 32;  /* Case-insensitive */
        h = h * 31 + (unsigned char)c;
        name++;
    }
    return h % HASH_TABLE_SIZE;
}

/* Register a 1-arg scalar function */
static void register_fn1(const char *name, scalar_fn_1 fn) {
    unsigned int h;
    FuncDispatchEntry *entry;
    
    if (num_entries >= MAX_DISPATCH_ENTRIES) return;
    
    entry = &dispatch_entries[num_entries++];
    entry->name = name;
    entry->category = FC_SCALAR_1;
    entry->min_args = 1;
    entry->max_args = 1;
    entry->impl.fn1 = fn;
    
    h = hash_name(name);
    entry->next = hash_table[h];
    hash_table[h] = entry;
}

/* Register a 2-arg scalar function */
static void register_fn2(const char *name, scalar_fn_2 fn) {
    unsigned int h;
    FuncDispatchEntry *entry;
    
    if (num_entries >= MAX_DISPATCH_ENTRIES) return;
    
    entry = &dispatch_entries[num_entries++];
    entry->name = name;
    entry->category = FC_SCALAR_2;
    entry->min_args = 2;
    entry->max_args = 2;
    entry->impl.fn2 = fn;
    
    h = hash_name(name);
    entry->next = hash_table[h];
    hash_table[h] = entry;
}

/* Initialize the dispatch table */
void func_dispatch_init(void) {
    int i;
    
    /* Clear hash table */
    for (i = 0; i < HASH_TABLE_SIZE; i++) {
        hash_table[i] = NULL;
    }
    num_entries = 0;
    
    /* Register 1-arg scalar functions */
    register_fn1("sin", apfx_sin);
    register_fn1("cos", apfx_cos);
    register_fn1("tan", apfx_tan);
    register_fn1("asin", apfx_asin);
    register_fn1("acos", apfx_acos);
    register_fn1("atan", apfx_atan);
    register_fn1("sinh", apfx_sinh);
    register_fn1("cosh", apfx_cosh);
    register_fn1("tanh", apfx_tanh);
    register_fn1("asinh", apfx_asinh);
    register_fn1("acosh", apfx_acosh);
    register_fn1("atanh", apfx_atanh);
    register_fn1("exp", apfx_exp);
    register_fn1("log", apfx_log);
    register_fn1("sqrt", apf_sqrt);
    register_fn1("abs", apf_abs);
    register_fn1("floor", apf_floor);
    register_fn1("ceil", apf_ceil);
    register_fn1("trunc", apf_trunc);
    
    /* Register 2-arg scalar functions */
    register_fn2("atan2", apfx_atan2);
    register_fn2("pow", apfx_pow);
    register_fn2("mod", apf_mod);
    
    /* Aliases */
    register_fn1("arcsin", apfx_asin);
    register_fn1("arccos", apfx_acos);
    register_fn1("arctan", apfx_atan);
    register_fn1("ln", apfx_log);
    register_fn1("fabs", apf_abs);
    register_fn1("int", apf_trunc);
    register_fn1("fix", apf_trunc);
}

/* Look up a function by name */
FuncDispatchEntry *func_dispatch_lookup(const char *name) {
    unsigned int h = hash_name(name);
    FuncDispatchEntry *entry = hash_table[h];
    
    while (entry) {
        if (str_eq_ci(entry->name, name)) {
            return entry;
        }
        entry = entry->next;
    }
    return NULL;
}

/* Evaluate a 1-arg scalar function */
int func_eval_scalar_1(func_value_t *result, func_value_t *args, scalar_fn_1 fn) {
    if (args[0].type == VAL_SCALAR) {
        result->type = VAL_SCALAR;
        fn(&result->v.scalar.re, &args[0].v.scalar.re);
        apf_zero(&result->v.scalar.im);
        return 1;
    } else if (args[0].type == VAL_MATRIX) {
        /* Apply element-wise */
        int i, j;
        result->type = VAL_MATRIX;
        mat_zero(&result->v.matrix, args[0].v.matrix.rows, args[0].v.matrix.cols);
        for (i = 0; i < args[0].v.matrix.rows; i++) {
            for (j = 0; j < args[0].v.matrix.cols; j++) {
                fn(&MAT_AT(&result->v.matrix, i, j).re, &MAT_AT(&args[0].v.matrix, i, j).re);
                apf_zero(&MAT_AT(&result->v.matrix, i, j).im);
            }
        }
        return 1;
    }
    return 0;
}

/* Evaluate a 2-arg scalar function */
int func_eval_scalar_2(func_value_t *result, func_value_t *args, scalar_fn_2 fn) {
    if (args[0].type == VAL_SCALAR && args[1].type == VAL_SCALAR) {
        result->type = VAL_SCALAR;
        fn(&result->v.scalar.re, &args[0].v.scalar.re, &args[1].v.scalar.re);
        apf_zero(&result->v.scalar.im);
        return 1;
    }
    /* TODO: Handle matrix cases with broadcasting */
    return 0;
}
