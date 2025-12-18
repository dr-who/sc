/*
 * func_dispatch.h - Unified function dispatch system
 * 
 * Replaces the massive str_eq chains in parser.c with a hash-table lookup
 * and generic argument parsing.
 * 
 * C89 compliant for Watcom C / DOS
 */
#ifndef FUNC_DISPATCH_H
#define FUNC_DISPATCH_H

#include "apf.h"
#include "apfc.h"
#include "matrix.h"

/* Maximum arguments for any function */
#define FUNC_MAX_ARGS 8

/* Value type (matches parser) */
typedef enum {
    VAL_SCALAR = 0,
    VAL_MATRIX = 1,
    VAL_STRING = 2
} ValueType;

typedef struct {
    ValueType type;
    union {
        apfc scalar;
        matrix_t matrix;
        char string[64];
    } v;
} func_value_t;

/* Function categories for dispatch */
typedef enum {
    FC_SCALAR_1,      /* 1-arg scalar: sin, cos, abs, etc. */
    FC_SCALAR_2,      /* 2-arg scalar: atan2, pow, mod, etc. */
    FC_SCALAR_3,      /* 3-arg scalar: clamp, lerp, etc. */
    FC_MATRIX_1,      /* 1-arg matrix: transpose, det, etc. */
    FC_MATRIX_2,      /* 2-arg matrix: dot, cross, etc. */
    FC_REDUCTION,     /* Vector reduction: sum, mean, max, etc. */
    FC_CREATION,      /* Matrix creation: zeros, ones, eye, etc. */
    FC_COMPLEX,       /* Complex number: real, imag, conj, etc. */
    FC_CUSTOM         /* Custom handling required */
} FuncCategory;

/* Scalar function pointer types */
typedef void (*scalar_fn_1)(apf *result, const apf *x);
typedef void (*scalar_fn_2)(apf *result, const apf *x, const apf *y);

/* Function dispatch entry */
typedef struct FuncDispatchEntry {
    const char *name;           /* Function name */
    FuncCategory category;      /* Dispatch category */
    int min_args;               /* Minimum arguments */
    int max_args;               /* Maximum arguments */
    union {
        scalar_fn_1 fn1;        /* For FC_SCALAR_1 */
        scalar_fn_2 fn2;        /* For FC_SCALAR_2 */
        void *custom;           /* For FC_CUSTOM */
    } impl;
    struct FuncDispatchEntry *next;  /* Hash chain */
} FuncDispatchEntry;

/* Initialize the dispatch table (call once at startup) */
void func_dispatch_init(void);

/* Look up a function by name - O(1) average case */
FuncDispatchEntry *func_dispatch_lookup(const char *name);

/* Generic argument parser
 * Returns number of arguments parsed, or -1 on error
 * Caller must call next_token() to skip past the function name first
 */
int func_parse_args(func_value_t *args, int min_args, int max_args);

/* Evaluate a scalar_1 function (handles both scalar and matrix inputs) */
int func_eval_scalar_1(func_value_t *result, func_value_t *args, scalar_fn_1 fn);

/* Evaluate a scalar_2 function */
int func_eval_scalar_2(func_value_t *result, func_value_t *args, scalar_fn_2 fn);

#endif /* FUNC_DISPATCH_H */
