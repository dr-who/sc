/*
 * matrix_funcs.h - Table-driven matrix function dispatch
 * 
 * Replaces str_eq chains in parse_value_factor with hash table lookup.
 * C89 compliant.
 */
#ifndef MATRIX_FUNCS_H
#define MATRIX_FUNCS_H

#include "matrix.h"
#include "apfc.h"

/* Forward declare value_t to avoid circular includes */
typedef struct value_s {
    int type;  /* VAL_SCALAR or VAL_MATRIX */
    union {
        apfc scalar;
        matrix_t matrix;
    } v;
} mf_value_t;

#define MF_VAL_SCALAR 0
#define MF_VAL_MATRIX 1

/* Matrix function categories */
typedef enum {
    MF_NONE = 0,
    
    /* 1-arg matrix -> matrix */
    MF_MAT_1,         /* f(A) -> B (matrix output) */
    
    /* 1-arg matrix -> scalar */
    MF_REDUCE_1,      /* f(A) -> scalar (e.g., sum, det) */
    
    /* 2-arg matrix functions */
    MF_MAT_2,         /* f(A, B) -> C */
    
    /* Matrix creation functions */
    MF_CREATE,        /* zeros(n), eye(n), etc. */
    
    /* Decomposition functions (return specific part) */
    MF_SVD_U,         /* svdu(A) -> U from SVD */
    MF_SVD_S,         /* svd(A) -> singular values as vector */
    MF_SVD_V,         /* svdv(A) -> V from SVD */
    MF_QR_Q,          /* qrq(A) -> Q from QR */
    MF_QR_R,          /* qr(A) -> R from QR */
    MF_LU_L,          /* lul(A) -> L from LU */
    MF_LU_U,          /* luu(A) -> U from LU */
    
    /* Custom handler needed */
    MF_CUSTOM
} MatFuncType;

/* Function pointer types for matrix operations */
typedef void (*mat_unary_fn)(matrix_t *r, const matrix_t *a);
typedef int  (*mat_unary_fn_err)(matrix_t *r, const matrix_t *a);
typedef void (*mat_reduce_fn)(apfc *r, const matrix_t *a);
typedef void (*mat_binary_fn)(matrix_t *r, const matrix_t *a, const matrix_t *b);

/* Matrix function table entry */
typedef struct MatFuncEntry {
    const char *name;
    MatFuncType type;
    union {
        mat_unary_fn m1;          /* For MF_MAT_1 */
        mat_unary_fn_err m1_err;  /* For MF_MAT_1 with error return */
        mat_reduce_fn red;        /* For MF_REDUCE_1 */
        mat_binary_fn m2;         /* For MF_MAT_2 */
        int custom_id;            /* For MF_CUSTOM */
    } fn;
    int can_fail;                 /* 1 if function can fail (returns int) */
    struct MatFuncEntry *next;    /* Hash chain */
} MatFuncEntry;

/* Initialize the matrix function table */
void matrix_funcs_init(void);

/* Look up a matrix function by name - returns NULL if not found */
MatFuncEntry *matrix_func_lookup(const char *name);

/* Evaluate a 1-arg matrix function
 * Returns 1 on success, 0 on error
 * arg is input, result is output
 */
int matrix_func_eval1(MatFuncEntry *entry, mf_value_t *result, mf_value_t *arg);

#endif /* MATRIX_FUNCS_H */
