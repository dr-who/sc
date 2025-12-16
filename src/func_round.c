/*
 * func_round.c - Rounding function documentation
 * C89 compliant
 */
#include "func_defs.h"

const FuncDoc func_round_docs[] = {
    {"floor", CAT_ROUND, "floor(x)",
     "Round towards negative infinity. Largest integer <= x.",
     {"floor(3.7) = 3", "floor(-2.3) = -3", "floor(5) = 5", "floor(0.9) = 0", "floor(-0.1) = -1", NULL},
     "ceil, round, trunc"},
    
    {"ceil", CAT_ROUND, "ceil(x)",
     "Round towards positive infinity. Smallest integer >= x.",
     {"ceil(3.2) = 4", "ceil(-2.7) = -2", "ceil(5) = 5", "ceil(0.1) = 1", "ceil(-0.9) = 0", NULL},
     "floor, round, trunc"},
    
    {"round", CAT_ROUND, "round(x)",
     "Round to nearest integer. Ties round away from zero.",
     {"round(3.4) = 3", "round(3.5) = 4", "round(-2.5) = -3", "round(3.6) = 4", "round(-3.4) = -3", NULL},
     "floor, ceil, trunc"},
    
    {"trunc", CAT_ROUND, "trunc(x)",
     "Round towards zero. Removes fractional part.",
     {"trunc(3.7) = 3", "trunc(-3.7) = -3", "trunc(3.2) = 3", "trunc(-0.9) = 0", NULL, NULL},
     "floor, ceil, round, fix"},
    
    {"fix", CAT_ROUND, "fix(x)",
     "Round towards zero. Alias for trunc. MATLAB compatible.",
     {"fix(3.7) = 3", "fix(-3.7) = -3", "fix(2.9) = 2", "fix(-2.9) = -2", NULL, NULL},
     "trunc, floor, ceil"},
    
    {"frac", CAT_ROUND, "frac(x)",
     "Fractional part. frac(x) = x - floor(x). Always in [0, 1).",
     {"frac(3.7) = 0.7", "frac(-2.3) = 0.7", "frac(5) = 0", "frac(0.25) = 0.25", NULL, NULL},
     "floor, trunc"},
    
    {"sign", CAT_ROUND, "sign(x)",
     "Sign of x. Returns -1 for x < 0, 0 for x = 0, 1 for x > 0.",
     {"sign(-5) = -1", "sign(0) = 0", "sign(3.14) = 1", "sign(-0.001) = -1", NULL, NULL},
     "abs, signum"},
    
    {"signum", CAT_ROUND, "signum(x)",
     "Sign function. Alias for sign.",
     {"signum(-5) = -1", "signum(0) = 0", "signum(100) = 1", NULL, NULL, NULL},
     "sign, abs"},
    
    {"mod", CAT_ROUND, "mod(x, y)",
     "Modulo. Result has same sign as y (MATLAB style). mod(x,y) = x - floor(x/y)*y.",
     {"mod(17, 5) = 2", "mod(-3, 5) = 2", "mod(10, 3) = 1", "mod(7.5, 2) = 1.5", "mod(-10, 3) = 2", NULL},
     "rem, floor"},
    
    {"rem", CAT_ROUND, "rem(x, y)",
     "Remainder. Result has same sign as x (C style). rem(x,y) = x - fix(x/y)*y.",
     {"rem(17, 5) = 2", "rem(-3, 5) = -3", "rem(10, 3) = 1", "rem(-17, 5) = -2", NULL, NULL},
     "mod, trunc"},
    
    {"clamp", CAT_ROUND, "clamp(x, lo, hi)",
     "Constrain x to range [lo, hi].",
     {"clamp(5, 0, 10) = 5", "clamp(-5, 0, 10) = 0", "clamp(15, 0, 10) = 10", "clamp(7, 5, 5) = 5", NULL, NULL},
     "min, max"},
    
    {"clip", CAT_ROUND, "clip(x, lo, hi)",
     "Alias for clamp.",
     {"clip(5, 0, 10) = 5", "clip(-5, 0, 10) = 0", "clip(15, 0, 10) = 10", NULL, NULL, NULL},
     "clamp"},
    
    /* Sentinel */
    {NULL, NULL, NULL, NULL, {NULL, NULL, NULL, NULL, NULL, NULL}, NULL}
};
