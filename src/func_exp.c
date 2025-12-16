/*
 * func_exp.c - Exponential and logarithmic function documentation
 * C89 compliant
 */
#include "func_defs.h"

/* Exponential functions */
const FuncDoc func_exp_docs[] = {
    {"exp", CAT_EXP, "exp(x)",
     "Natural exponential e^x. For complex x, exp(x) = e^(re(x)) * (cos(im(x)) + i*sin(im(x))).",
     {"exp(0) = 1", "exp(1) = 2.7183", "exp(2) = 7.3891", "exp(-1) = 0.3679", "exp(ln(5)) = 5", NULL},
     "log, exp2, expm1"},
    
    {"exp2", CAT_EXP, "exp2(x)",
     "Base-2 exponential: 2^x. Useful for computing powers of 2.",
     {"exp2(0) = 1", "exp2(3) = 8", "exp2(10) = 1024", "exp2(-1) = 0.5", "exp2(0.5) = 1.4142", NULL},
     "exp, log2, pow2"},
    
    {"expm1", CAT_EXP, "expm1(x)",
     "exp(x) - 1. More accurate than exp(x)-1 for small x.",
     {"expm1(0) = 0", "expm1(1) = 1.7183", "expm1(0.001) = 0.001001", "expm1(-0.001) = -0.0009995", NULL, NULL},
     "exp, log1p"},
    
    {"pow2", CAT_EXP, "pow2(x)",
     "2^x. Alias for exp2.",
     {"pow2(8) = 256", "pow2(10) = 1024", "pow2(-3) = 0.125", NULL, NULL, NULL},
     "exp2, sqrt"},
    
    {"sqrt", CAT_EXP, "sqrt(x)",
     "Square root. For negative x, returns complex result i*sqrt(-x).",
     {"sqrt(4) = 2", "sqrt(2) = 1.4142", "sqrt(0) = 0", "sqrt(-1) = i", "sqrt(9) = 3", NULL},
     "cbrt, nthroot, exp"},
    
    {"cbrt", CAT_EXP, "cbrt(x)",
     "Cube root. Returns real result for negative x (unlike x^(1/3)).",
     {"cbrt(8) = 2", "cbrt(27) = 3", "cbrt(-8) = -2", "cbrt(1000) = 10", NULL, NULL},
     "sqrt, nthroot"},
    
    {"nthroot", CAT_EXP, "nthroot(x, n)",
     "Real nth root. For negative x with odd n, returns negative real result.",
     {"nthroot(16, 4) = 2", "nthroot(32, 5) = 2", "nthroot(-8, 3) = -2", "nthroot(81, 4) = 3", NULL, NULL},
     "sqrt, cbrt"},
    
    {"hypot", CAT_EXP, "hypot(x, y)",
     "Hypotenuse: sqrt(x^2 + y^2). Avoids overflow for large x, y.",
     {"hypot(3, 4) = 5", "hypot(5, 12) = 13", "hypot(1, 1) = 1.4142", "hypot(0, 5) = 5", NULL, NULL},
     "sqrt"},
    
    {"realsqrt", CAT_EXP, "realsqrt(x)",
     "Real-valued square root. Error for x < 0 (unlike sqrt which returns complex).",
     {"realsqrt(4) = 2", "realsqrt(2) = 1.4142", "realsqrt(0) = 0", NULL, NULL, NULL},
     "sqrt, reallog"},
    
    {"realpow", CAT_EXP, "realpow(x, y)",
     "Real-valued power. Error if result would be complex.",
     {"realpow(2, 3) = 8", "realpow(4, 0.5) = 2", "realpow(27, 1/3) = 3", NULL, NULL, NULL},
     "sqrt, nthroot"},
    
    /* Sentinel */
    {NULL, NULL, NULL, NULL, {NULL, NULL, NULL, NULL, NULL, NULL}, NULL}
};

/* Logarithmic functions */
const FuncDoc func_log_docs[] = {
    {"ln", CAT_LOG, "ln(x)",
     "Natural logarithm (base e). For complex x, returns principal value.",
     {"ln(1) = 0", "ln(e) = 1", "ln(10) = 2.3026", "ln(exp(5)) = 5", "ln(2) = 0.6931", NULL},
     "log, log10, log2, exp"},
    
    {"log", CAT_LOG, "log(x)",
     "Natural logarithm. Alias for ln.",
     {"log(1) = 0", "log(e) = 1", "log(10) = 2.3026", "log(100) = 4.6052", NULL, NULL},
     "ln, log10, log2"},
    
    {"log10", CAT_LOG, "log10(x)",
     "Base-10 logarithm. log10(x) = ln(x) / ln(10).",
     {"log10(1) = 0", "log10(10) = 1", "log10(100) = 2", "log10(1000) = 3", "log10(2) = 0.3010", NULL},
     "log, log2, exp"},
    
    {"log2", CAT_LOG, "log2(x)",
     "Base-2 logarithm. log2(x) = ln(x) / ln(2). Useful for bit counts.",
     {"log2(1) = 0", "log2(2) = 1", "log2(8) = 3", "log2(1024) = 10", "log2(256) = 8", NULL},
     "log, log10, exp2"},
    
    {"log1p", CAT_LOG, "log1p(x)",
     "ln(1 + x). More accurate than log(1+x) for small x.",
     {"log1p(0) = 0", "log1p(1) = 0.6931", "log1p(0.001) = 0.0009995", "log1p(-0.5) = -0.6931", NULL, NULL},
     "log, expm1"},
    
    {"logb", CAT_LOG, "logb(x, b)",
     "Logarithm of x with base b. logb(x, b) = ln(x) / ln(b).",
     {"logb(8, 2) = 3", "logb(1000, 10) = 3", "logb(27, 3) = 3", "logb(256, 2) = 8", NULL, NULL},
     "log, log10, log2"},
    
    {"reallog", CAT_LOG, "reallog(x)",
     "Real-valued logarithm. Error for x <= 0 (unlike log which returns complex).",
     {"reallog(1) = 0", "reallog(10) = 2.3026", "reallog(e) = 1", NULL, NULL, NULL},
     "log, realsqrt"},
    
    /* Sentinel */
    {NULL, NULL, NULL, NULL, {NULL, NULL, NULL, NULL, NULL, NULL}, NULL}
};
