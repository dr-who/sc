/*
 * scalar_funcs.h - Table-driven scalar function dispatch
 * 
 * Replaces hundreds of str_eq chains with a hash table lookup.
 * C89 compliant.
 */
#ifndef SCALAR_FUNCS_H
#define SCALAR_FUNCS_H

#include "apf.h"
#include "apfc.h"

/* Function categories for scalar dispatch */
typedef enum {
    SF_NONE = 0,
    
    /* Real -> Real, applied to real part only */
    SF_REAL_1,        /* f(x.re) -> result.re, result.im = 0 */
    SF_REAL_2,        /* f(x.re, y.re) -> result.re */
    
    /* Complex -> Complex */
    SF_COMPLEX_1,     /* f(x) -> result (complex) */
    SF_COMPLEX_2,     /* f(x, y) -> result (complex) */
    
    /* Trig with angle mode */
    SF_TRIG,          /* Uses angle_mode for input (complex capable) */
    SF_TRIG_REAL,     /* Uses angle_mode for input (real only) */
    SF_TRIG_INV,      /* Uses angle_mode for output */
    SF_TRIG_INV_2,    /* 2-arg inverse trig (atan2), output in angle_mode */
    SF_TRIG_INV_2_DEG, /* 2-arg inverse trig always outputting degrees */
    
    /* Always degrees */
    SF_TRIG_DEG,      /* Input in degrees (complex capable) */
    SF_TRIG_DEG_REAL, /* Input in degrees (real only) */
    SF_TRIG_INV_DEG,  /* Output in degrees */
    
    /* Complex number parts */
    SF_REAL_PART,     /* Returns real part */
    SF_IMAG_PART,     /* Returns imag part */
    SF_CONJ,          /* Complex conjugate */
    SF_ABS_COMPLEX,   /* Complex magnitude */
    SF_ARG,           /* Complex argument (angle) */
    
    /* Integer sequence functions (Fibonacci, Bell, etc) */
    SF_INT_1,         /* f(r, long n) - input converted to long */
    SF_INT_2,         /* f(r, long n, long k) - both inputs converted to long */
    SF_APF_INT,       /* f(r, apf x, long n) - mixed apf and int args */
    SF_INT_APF,       /* f(r, int n, apf x) - int first, apf second (bessel) */
    
    /* Complex construction */
    SF_MAKE_COMPLEX,  /* complex(re, im) -> result.re=arg1, result.im=arg2 */
    
    /* Special handling required */
    SF_CUSTOM
} ScalarFuncType;

/* Function pointer types */
typedef void (*apf_fn_1)(apf *r, const apf *x);
typedef void (*apf_fn_2)(apf *r, const apf *x, const apf *y);
typedef void (*apfc_fn_1)(apfc *r, const apfc *x);
typedef void (*apfc_fn_2)(apfc *r, const apfc *x, const apfc *y);
typedef void (*apf_int_fn)(apf *r, long n);  /* For integer sequence functions */
typedef void (*apf_int2_fn)(apf *r, long n, long k);  /* For 2-arg integer functions */
typedef void (*apf_apf_int_fn)(apf *r, const apf *x, long n);  /* Mixed apf + int */
typedef void (*apf_int_apf_fn)(apf *r, int n, const apf *x);  /* Int first, apf second */

/* Scalar function table entry */
typedef struct ScalarFuncEntry {
    const char *name;
    ScalarFuncType type;
    union {
        apf_fn_1 f1;       /* For SF_REAL_1 */
        apf_fn_2 f2;       /* For SF_REAL_2 */
        apfc_fn_1 c1;      /* For SF_COMPLEX_1, SF_TRIG, etc. */
        apfc_fn_2 c2;      /* For SF_COMPLEX_2 */
        apf_int_fn fi;     /* For SF_INT_1 */
        apf_int2_fn fi2;   /* For SF_INT_2 */
        apf_apf_int_fn fai; /* For SF_APF_INT */
        apf_int_apf_fn fia; /* For SF_INT_APF (bessel) */
    } fn;
    long min_val;          /* For SF_INT_1: minimum allowed value */
    long max_val;          /* For SF_INT_1: maximum allowed value */
    struct ScalarFuncEntry *next;  /* Hash chain */
} ScalarFuncEntry;

/* Initialize the scalar function table */
void scalar_funcs_init(void);

/* Look up a scalar function by name - returns NULL if not found */
ScalarFuncEntry *scalar_func_lookup(const char *name);

/* Evaluate a scalar function given entry, argument, and result
 * Returns 1 on success, 0 on error
 */
int scalar_func_eval(ScalarFuncEntry *entry, apfc *result, const apfc *arg);

/* Evaluate a 2-arg scalar function */
int scalar_func_eval2(ScalarFuncEntry *entry, apfc *result, const apfc *arg1, const apfc *arg2);

#endif /* SCALAR_FUNCS_H */
