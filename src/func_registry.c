/*
 * func_registry.c - Unified function registry implementation
 * 
 * Single source of truth for all built-in functions.
 * C89 compliant
 */
#include <stdio.h>
#include <string.h>
#include "func_registry.h"

/* External evaluator for demos */
extern int eval_expr_line(const char *line, int quiet);

/* String comparison helper */
static int str_eq(const char *a, const char *b) {
    while (*a && *b) {
        char ca = *a, cb = *b;
        if (ca >= 'A' && ca <= 'Z') ca += 32;
        if (cb >= 'A' && cb <= 'Z') cb += 32;
        if (ca != cb) return 0;
        a++; b++;
    }
    return *a == *b;
}

/* ========== Master Function Registry ========== */

static const FuncEntry func_registry[] = {
    /* ===== CONSTANTS ===== */
    {"pi", FH_CONST_PI, -1, -1, FF_CONST, "Constants", "pi",
     "Mathematical constant pi = 3.14159...",
     {"pi = 3.1416", "2*pi = 6.2832", "sin(pi) = 0", NULL, NULL, NULL},
     "tau, e"},
    
    {"e", FH_CONST_E, -1, -1, FF_CONST, "Constants", "e",
     "Euler's number e = 2.71828...",
     {"e = 2.7183", "exp(1) = e", "log(e) = 1", NULL, NULL, NULL},
     "exp, log, pi"},
    
    {"i", FH_CONST_I, -1, -1, FF_CONST, "Constants", "i",
     "Imaginary unit. i^2 = -1.",
     {"i^2 = -1", "sqrt(-1) = i", "exp(i*pi) = -1", NULL, NULL, NULL},
     "j, complex, real, imag"},
    
    {"j", FH_CONST_I, -1, -1, FF_CONST, "Constants", "j",
     "Imaginary unit (alias for i). j^2 = -1. Common in electrical engineering.",
     {"j^2 = -1", "j = i", "1 + 2j", NULL, NULL, NULL},
     "i, complex, real, imag"},
    
    {"Inf", FH_CONST_INF, -1, -1, FF_CONST, "Constants", "Inf",
     "Positive infinity.",
     {"1/0 = Inf", "-1/0 = -Inf", "Inf + 1 = Inf", NULL, NULL, NULL},
     "NaN, isfinite, isinf"},
    
    {"inf", FH_CONST_INF, -1, -1, FF_CONST | FF_DEPRECATED, "Constants", "inf",
     "Positive infinity (alias for Inf).",
     {"inf = Inf", NULL, NULL, NULL, NULL, NULL},
     "Inf"},
    
    {"NaN", FH_CONST_NAN, -1, -1, FF_CONST, "Constants", "NaN",
     "Not a Number. Result of undefined operations like 0/0.",
     {"0/0 = NaN", "NaN + 1 = NaN", "NaN == NaN = false", NULL, NULL, NULL},
     "Inf, isnan, isfinite"},
    
    {"nan", FH_CONST_NAN, -1, -1, FF_CONST | FF_DEPRECATED, "Constants", "nan",
     "Not a Number (alias for NaN).",
     {"nan = NaN", NULL, NULL, NULL, NULL, NULL},
     "NaN"},
    
    {"true", FH_CONST_TRUE, -1, -1, FF_CONST, "Constants", "true",
     "Boolean true (1).",
     {"true = 1", "1 == 1 = true", NULL, NULL, NULL, NULL},
     "false"},
    
    {"false", FH_CONST_FALSE, -1, -1, FF_CONST, "Constants", "false",
     "Boolean false (0).",
     {"false = 0", "1 == 2 = false", NULL, NULL, NULL, NULL},
     "true"},
    
    {"eps", FH_CONST_EPS, -1, -1, FF_CONST, "Constants", "eps",
     "Machine epsilon (IEEE double compatible). Smallest x where 1+x != 1.",
     {"eps", "1 + eps/2 = 1", "1 + eps > 1", NULL, NULL, NULL},
     "epsp, realmin, realmax"},
    
    {"epsp", FH_CONST_EPSP, -1, -1, FF_CONST, "Constants", "epsp",
     "Native epsilon for current precision. Use 'precision' to see current bits.",
     {"epsp", "precision", NULL, NULL, NULL, NULL},
     "eps, precision"},
    
    {"epsn", FH_CONST_EPSP, -1, -1, FF_CONST, "Constants", "epsn",
     "Alias for epsp (native epsilon).",
     {"epsn", NULL, NULL, NULL, NULL, NULL},
     "epsp"},
    
    {"precision", FH_CONST_PRECISION, -1, -1, FF_CONST, "Constants", "precision",
     "Current arithmetic precision in bits (64, 128, or 256).",
     {"precision = 128", "precision = 256", NULL, NULL, NULL, NULL},
     "epsp, eps"},
    
    {"apbits", FH_CONST_PRECISION, -1, -1, FF_CONST, "Constants", "apbits",
     "Alias for precision (current bits).",
     {"apbits", NULL, NULL, NULL, NULL, NULL},
     "precision"},
    
    {"realmax", FH_CONST_REALMAX, -1, -1, FF_CONST, "Constants", "realmax",
     "Largest representable finite number.",
     {"realmax", NULL, NULL, NULL, NULL, NULL},
     "realmin, eps, Inf"},
    
    {"realmin", FH_CONST_REALMIN, -1, -1, FF_CONST, "Constants", "realmin",
     "Smallest positive normalized number.",
     {"realmin", NULL, NULL, NULL, NULL, NULL},
     "realmax, eps"},
    
    {"phi", FH_CONST_PHI, -1, -1, FF_CONST, "Constants", "phi",
     "Golden ratio phi = (1+sqrt(5))/2 = 1.618...",
     {"phi = 1.618", "phi^2 - phi - 1 = 0", NULL, NULL, NULL, NULL},
     "golden"},
    
    {"golden", FH_CONST_PHI, -1, -1, FF_CONST, "Constants", "golden",
     "Golden ratio (alias for phi).",
     {"golden = phi", NULL, NULL, NULL, NULL, NULL},
     "phi"},
    
    {"tau", FH_CONST_TAU, -1, -1, FF_CONST, "Constants", "tau",
     "Tau = 2*pi = 6.28318...",
     {"tau = 6.2832", "tau = 2*pi", "sin(tau) = 0", NULL, NULL, NULL},
     "pi"},
    
    /* ===== TRIGONOMETRIC - BASIC ===== */
    {"sin", FH_SIN, 1, 1, FF_SCALAR | FF_ANGLE | FF_COMPLEX, "Trigonometry", "sin(x)",
     "Sine of x (radians). For complex x, sin(x) = (e^(ix) - e^(-ix)) / (2i).",
     {"sin(0) = 0", "sin(pi/6) = 0.5", "sin(pi/2) = 1", "sin(pi) = 0", NULL, NULL},
     "cos, tan, asin, sind"},
    
    {"cos", FH_COS, 1, 1, FF_SCALAR | FF_ANGLE | FF_COMPLEX, "Trigonometry", "cos(x)",
     "Cosine of x (radians). For complex x, cos(x) = (e^(ix) + e^(-ix)) / 2.",
     {"cos(0) = 1", "cos(pi/3) = 0.5", "cos(pi/2) = 0", "cos(pi) = -1", NULL, NULL},
     "sin, tan, acos, cosd"},
    
    {"tan", FH_TAN, 1, 1, FF_SCALAR | FF_ANGLE | FF_COMPLEX, "Trigonometry", "tan(x)",
     "Tangent of x (radians). tan(x) = sin(x)/cos(x).",
     {"tan(0) = 0", "tan(pi/4) = 1", "tan(pi/3) = 1.7321", NULL, NULL, NULL},
     "sin, cos, atan, tand"},
    
    {"asin", FH_ASIN, 1, 1, FF_SCALAR | FF_COMPLEX, "Trigonometry", "asin(x)",
     "Arc sine (inverse sine). Returns angle in radians.",
     {"asin(0) = 0", "asin(0.5) = 0.5236", "asin(1) = 1.5708", NULL, NULL, NULL},
     "sin, acos, atan, asind"},
    
    {"acos", FH_ACOS, 1, 1, FF_SCALAR | FF_COMPLEX, "Trigonometry", "acos(x)",
     "Arc cosine (inverse cosine). Returns angle in radians.",
     {"acos(1) = 0", "acos(0.5) = 1.0472", "acos(0) = 1.5708", NULL, NULL, NULL},
     "cos, asin, atan, acosd"},
    
    {"atan", FH_ATAN, 1, 1, FF_SCALAR | FF_COMPLEX, "Trigonometry", "atan(x)",
     "Arc tangent (inverse tangent). Returns angle in radians.",
     {"atan(0) = 0", "atan(1) = 0.7854", "atan(Inf) = 1.5708", NULL, NULL, NULL},
     "tan, asin, acos, atan2, atand"},
    
    {"atan2", FH_ATAN2, 2, 2, FF_SCALAR, "Trigonometry", "atan2(y, x)",
     "Two-argument arc tangent. Returns angle in radians from -pi to pi.",
     {"atan2(0, 1) = 0", "atan2(1, 1) = 0.7854", "atan2(1, 0) = 1.5708", "atan2(-1, -1) = -2.3562", NULL, NULL},
     "atan, atan2d"},
    
    /* Hyperbolic */
    {"sinh", FH_SINH, 1, 1, FF_SCALAR | FF_COMPLEX, "Trigonometry", "sinh(x)",
     "Hyperbolic sine. sinh(x) = (e^x - e^(-x)) / 2.",
     {"sinh(0) = 0", "sinh(1) = 1.1752", "sinh(-1) = -1.1752", NULL, NULL, NULL},
     "cosh, tanh, asinh"},
    
    {"cosh", FH_COSH, 1, 1, FF_SCALAR | FF_COMPLEX, "Trigonometry", "cosh(x)",
     "Hyperbolic cosine. cosh(x) = (e^x + e^(-x)) / 2.",
     {"cosh(0) = 1", "cosh(1) = 1.5431", "cosh(-1) = 1.5431", NULL, NULL, NULL},
     "sinh, tanh, acosh"},
    
    {"tanh", FH_TANH, 1, 1, FF_SCALAR | FF_COMPLEX, "Trigonometry", "tanh(x)",
     "Hyperbolic tangent. tanh(x) = sinh(x)/cosh(x).",
     {"tanh(0) = 0", "tanh(1) = 0.7616", "tanh(Inf) = 1", NULL, NULL, NULL},
     "sinh, cosh, atanh"},
    
    {"asinh", FH_ASINH, 1, 1, FF_SCALAR | FF_COMPLEX, "Trigonometry", "asinh(x)",
     "Inverse hyperbolic sine.",
     {"asinh(0) = 0", "asinh(1) = 0.8814", "asinh(-1) = -0.8814", NULL, NULL, NULL},
     "sinh, acosh, atanh"},
    
    {"acosh", FH_ACOSH, 1, 1, FF_SCALAR | FF_COMPLEX, "Trigonometry", "acosh(x)",
     "Inverse hyperbolic cosine. x must be >= 1 for real result.",
     {"acosh(1) = 0", "acosh(2) = 1.3170", "acosh(10) = 2.9932", NULL, NULL, NULL},
     "cosh, asinh, atanh"},
    
    {"atanh", FH_ATANH, 1, 1, FF_SCALAR | FF_COMPLEX, "Trigonometry", "atanh(x)",
     "Inverse hyperbolic tangent. x must be in (-1, 1) for real result.",
     {"atanh(0) = 0", "atanh(0.5) = 0.5493", "atanh(-0.5) = -0.5493", NULL, NULL, NULL},
     "tanh, asinh, acosh"},
    
    /* Reciprocal trig */
    {"sec", FH_SEC, 1, 1, FF_SCALAR | FF_ANGLE | FF_COMPLEX, "Trigonometry", "sec(x)",
     "Secant. sec(x) = 1/cos(x).",
     {"sec(0) = 1", "sec(pi/3) = 2", "sec(pi/4) = 1.4142", NULL, NULL, NULL},
     "cos, csc, cot, asec"},
    
    {"csc", FH_CSC, 1, 1, FF_SCALAR | FF_ANGLE | FF_COMPLEX, "Trigonometry", "csc(x)",
     "Cosecant. csc(x) = 1/sin(x).",
     {"csc(pi/2) = 1", "csc(pi/6) = 2", "csc(pi/4) = 1.4142", NULL, NULL, NULL},
     "sin, sec, cot, acsc"},
    
    {"cot", FH_COT, 1, 1, FF_SCALAR | FF_ANGLE | FF_COMPLEX, "Trigonometry", "cot(x)",
     "Cotangent. cot(x) = 1/tan(x) = cos(x)/sin(x).",
     {"cot(pi/4) = 1", "cot(pi/3) = 0.5774", "cot(pi/2) = 0", NULL, NULL, NULL},
     "tan, sec, csc, acot"},
    
    {"asec", FH_ASEC, 1, 1, FF_SCALAR | FF_COMPLEX, "Trigonometry", "asec(x)",
     "Arc secant (inverse secant).",
     {"asec(1) = 0", "asec(2) = 1.0472", "asec(-1) = 3.1416", NULL, NULL, NULL},
     "sec, acos"},
    
    {"acsc", FH_ACSC, 1, 1, FF_SCALAR | FF_COMPLEX, "Trigonometry", "acsc(x)",
     "Arc cosecant (inverse cosecant).",
     {"acsc(1) = 1.5708", "acsc(2) = 0.5236", "acsc(-1) = -1.5708", NULL, NULL, NULL},
     "csc, asin"},
    
    {"acot", FH_ACOT, 1, 1, FF_SCALAR | FF_COMPLEX, "Trigonometry", "acot(x)",
     "Arc cotangent (inverse cotangent).",
     {"acot(1) = 0.7854", "acot(0) = 1.5708", "acot(-1) = -0.7854", NULL, NULL, NULL},
     "cot, atan"},
    
    /* Hyperbolic reciprocal */
    {"sech", FH_SECH, 1, 1, FF_SCALAR | FF_COMPLEX, "Trigonometry", "sech(x)",
     "Hyperbolic secant. sech(x) = 1/cosh(x).",
     {"sech(0) = 1", "sech(1) = 0.6481", "sech(-1) = 0.6481", NULL, NULL, NULL},
     "cosh, csch, coth"},
    
    {"csch", FH_CSCH, 1, 1, FF_SCALAR | FF_COMPLEX, "Trigonometry", "csch(x)",
     "Hyperbolic cosecant. csch(x) = 1/sinh(x).",
     {"csch(1) = 0.8509", "csch(-1) = -0.8509", "csch(2) = 0.2757", NULL, NULL, NULL},
     "sinh, sech, coth"},
    
    {"coth", FH_COTH, 1, 1, FF_SCALAR | FF_COMPLEX, "Trigonometry", "coth(x)",
     "Hyperbolic cotangent. coth(x) = 1/tanh(x).",
     {"coth(1) = 1.3130", "coth(-1) = -1.3130", "coth(2) = 1.0373", NULL, NULL, NULL},
     "tanh, sech, csch"},
    
    {"asech", FH_ASECH, 1, 1, FF_SCALAR | FF_COMPLEX, "Trigonometry", "asech(x)",
     "Inverse hyperbolic secant.",
     {"asech(1) = 0", "asech(0.5) = 1.3170", "asech(0.1) = 2.9932", NULL, NULL, NULL},
     "sech, acosh"},
    
    {"acsch", FH_ACSCH, 1, 1, FF_SCALAR | FF_COMPLEX, "Trigonometry", "acsch(x)",
     "Inverse hyperbolic cosecant.",
     {"acsch(1) = 0.8814", "acsch(-1) = -0.8814", "acsch(2) = 0.4812", NULL, NULL, NULL},
     "csch, asinh"},
    
    {"acoth", FH_ACOTH, 1, 1, FF_SCALAR | FF_COMPLEX, "Trigonometry", "acoth(x)",
     "Inverse hyperbolic cotangent.",
     {"acoth(2) = 0.5493", "acoth(-2) = -0.5493", "acoth(10) = 0.1003", NULL, NULL, NULL},
     "coth, atanh"},
    
    /* Degree variants */
    {"sind", FH_SIND, 1, 1, FF_SCALAR | FF_DEGREE, "Trigonometry", "sind(x)",
     "Sine of x in degrees.",
     {"sind(0) = 0", "sind(30) = 0.5", "sind(90) = 1", "sind(180) = 0", NULL, NULL},
     "sin, cosd, tand"},
    
    {"cosd", FH_COSD, 1, 1, FF_SCALAR | FF_DEGREE, "Trigonometry", "cosd(x)",
     "Cosine of x in degrees.",
     {"cosd(0) = 1", "cosd(60) = 0.5", "cosd(90) = 0", "cosd(180) = -1", NULL, NULL},
     "cos, sind, tand"},
    
    {"tand", FH_TAND, 1, 1, FF_SCALAR | FF_DEGREE, "Trigonometry", "tand(x)",
     "Tangent of x in degrees.",
     {"tand(0) = 0", "tand(45) = 1", "tand(60) = 1.7321", NULL, NULL, NULL},
     "tan, sind, cosd"},
    
    {"asind", FH_ASIND, 1, 1, FF_SCALAR | FF_DEGREE, "Trigonometry", "asind(x)",
     "Arc sine returning degrees.",
     {"asind(0) = 0", "asind(0.5) = 30", "asind(1) = 90", NULL, NULL, NULL},
     "asin, acosd, atand"},
    
    {"acosd", FH_ACOSD, 1, 1, FF_SCALAR | FF_DEGREE, "Trigonometry", "acosd(x)",
     "Arc cosine returning degrees.",
     {"acosd(1) = 0", "acosd(0.5) = 60", "acosd(0) = 90", NULL, NULL, NULL},
     "acos, asind, atand"},
    
    {"atand", FH_ATAND, 1, 1, FF_SCALAR | FF_DEGREE, "Trigonometry", "atand(x)",
     "Arc tangent returning degrees.",
     {"atand(0) = 0", "atand(1) = 45", "atand(Inf) = 90", NULL, NULL, NULL},
     "atan, asind, acosd, atan2d"},
    
    {"atan2d", FH_ATAN2D, 2, 2, FF_SCALAR | FF_DEGREE, "Trigonometry", "atan2d(y, x)",
     "Two-argument arc tangent returning degrees.",
     {"atan2d(0, 1) = 0", "atan2d(1, 1) = 45", "atan2d(1, 0) = 90", NULL, NULL, NULL},
     "atan2, atand"},
    
    /* Special trig */
    {"sinc", FH_SINC, 1, 1, FF_SCALAR, "Trigonometry", "sinc(x)",
     "Sinc function. sinc(x) = sin(pi*x)/(pi*x), sinc(0) = 1.",
     {"sinc(0) = 1", "sinc(1) = 0", "sinc(0.5) = 0.6366", NULL, NULL, NULL},
     "sin"},
    
    {"sinpi", FH_SINPI, 1, 1, FF_SCALAR, "Trigonometry", "sinpi(x)",
     "Sine of pi*x. More accurate than sin(pi*x).",
     {"sinpi(0) = 0", "sinpi(0.5) = 1", "sinpi(1) = 0", "sinpi(2) = 0", NULL, NULL},
     "sin, cospi"},
    
    {"cospi", FH_COSPI, 1, 1, FF_SCALAR, "Trigonometry", "cospi(x)",
     "Cosine of pi*x. More accurate than cos(pi*x).",
     {"cospi(0) = 1", "cospi(0.5) = 0", "cospi(1) = -1", "cospi(2) = 1", NULL, NULL},
     "cos, sinpi"},
    
    /* ===== EXPONENTIAL AND POWERS ===== */
    {"exp", FH_EXP, 1, 1, FF_SCALAR | FF_COMPLEX, "Exponential", "exp(x)",
     "Exponential function e^x.",
     {"exp(0) = 1", "exp(1) = 2.7183", "exp(-1) = 0.3679", "log(exp(x)) = x", NULL, NULL},
     "log, exp2, expm1"},
    
    {"exp2", FH_EXP2, 1, 1, FF_SCALAR | FF_COMPLEX, "Exponential", "exp2(x)",
     "Base-2 exponential. exp2(x) = 2^x.",
     {"exp2(0) = 1", "exp2(3) = 8", "exp2(10) = 1024", NULL, NULL, NULL},
     "exp, log2, pow2"},
    
    {"expm1", FH_EXPM1, 1, 1, FF_SCALAR, "Exponential", "expm1(x)",
     "exp(x) - 1. More accurate for small x.",
     {"expm1(0) = 0", "expm1(1) = 1.7183", "expm1(0.001) = 0.001001", NULL, NULL, NULL},
     "exp, log1p"},
    
    {"pow2", FH_POW2, 1, 1, FF_SCALAR, "Exponential", "pow2(x)",
     "Power of 2. pow2(x) = 2^x.",
     {"pow2(0) = 1", "pow2(8) = 256", "pow2(10) = 1024", NULL, NULL, NULL},
     "exp2, log2"},
    
    {"exp10", FH_EXP10, 1, 1, FF_SCALAR, "Exponential", "exp10(x)",
     "Power of 10. exp10(x) = 10^x.",
     {"exp10(0) = 1", "exp10(1) = 10", "exp10(2) = 100", "exp10(3) = 1000", NULL, NULL},
     "log10, pow10"},
    
    {"pow10", FH_POW10, 1, 1, FF_SCALAR, "Exponential", "pow10(x)",
     "Power of 10. Alias for exp10.",
     {"pow10(0) = 1", "pow10(2) = 100", "pow10(-1) = 0.1", NULL, NULL, NULL},
     "exp10, log10"},
    
    {"pow", FH_POW, 2, 2, FF_SCALAR | FF_COMPLEX, "Exponential", "pow(x, y)",
     "Power. pow(x, y) = x^y.",
     {"pow(2, 3) = 8", "pow(10, 2) = 100", "pow(27, 1/3) = 3", "pow(2, -1) = 0.5", NULL, NULL},
     "exp, log, sqrt"},
    
    {"cpow", FH_CPOW, 2, 2, FF_COMPLEX, "Exponential", "cpow(z, w)",
     "Complex power. cpow(z, w) = z^w for complex z and w.",
     {"cpow(2+3i, 2) = -5+12i", "cpow(e, i*pi) = -1", NULL, NULL, NULL, NULL},
     "pow, exp, log"},
    
    {"fabs", FH_FABS, 1, 1, FF_SCALAR, "Arithmetic", "fabs(x)",
     "Floating-point absolute value. Same as abs for real numbers.",
     {"fabs(-5) = 5", "fabs(3.14) = 3.14", "fabs(0) = 0", NULL, NULL, NULL},
     "abs, sign"},
    
    {"sqrt", FH_SQRT, 1, 1, FF_SCALAR | FF_COMPLEX, "Exponential", "sqrt(x)",
     "Square root. Returns complex for negative x.",
     {"sqrt(4) = 2", "sqrt(2) = 1.4142", "sqrt(-1) = i", "sqrt(0) = 0", NULL, NULL},
     "cbrt, nthroot, realsqrt"},
    
    {"cbrt", FH_CBRT, 1, 1, FF_SCALAR, "Exponential", "cbrt(x)",
     "Cube root. cbrt(x) = x^(1/3). Handles negative x.",
     {"cbrt(8) = 2", "cbrt(-8) = -2", "cbrt(27) = 3", NULL, NULL, NULL},
     "sqrt, nthroot"},
    
    {"nthroot", FH_NTHROOT, 2, 2, FF_SCALAR, "Exponential", "nthroot(x, n)",
     "Real nth root of x. For even n, x must be non-negative.",
     {"nthroot(16, 4) = 2", "nthroot(-8, 3) = -2", "nthroot(32, 5) = 2", NULL, NULL, NULL},
     "sqrt, cbrt"},
    
    {"hypot", FH_HYPOT, 2, 2, FF_SCALAR, "Exponential", "hypot(x, y)",
     "Hypotenuse. hypot(x,y) = sqrt(x^2 + y^2).",
     {"hypot(3, 4) = 5", "hypot(5, 12) = 13", "hypot(1, 1) = 1.4142", NULL, NULL, NULL},
     "sqrt, norm"},
    
    {"realsqrt", FH_REALSQRT, 1, 1, FF_SCALAR | FF_REAL_ONLY, "Exponential", "realsqrt(x)",
     "Real square root. Error for negative x.",
     {"realsqrt(4) = 2", "realsqrt(2) = 1.4142", NULL, NULL, NULL, NULL},
     "sqrt, realpow"},
    
    {"realpow", FH_REALPOW, 2, 2, FF_SCALAR | FF_REAL_ONLY, "Exponential", "realpow(x, y)",
     "Real power. Error if result would be complex.",
     {"realpow(2, 3) = 8", "realpow(27, 1/3) = 3", NULL, NULL, NULL, NULL},
     "pow, realsqrt"},
    
    /* ===== LOGARITHMIC ===== */
    {"log", FH_LOG, 1, 1, FF_SCALAR | FF_COMPLEX, "Logarithmic", "log(x)",
     "Natural logarithm (base e).",
     {"log(1) = 0", "log(e) = 1", "log(10) = 2.3026", "log(exp(x)) = x", NULL, NULL},
     "exp, log10, log2, ln"},
    
    {"ln", FH_LN, 1, 1, FF_SCALAR | FF_COMPLEX, "Logarithmic", "ln(x)",
     "Natural logarithm. Alias for log.",
     {"ln(1) = 0", "ln(e) = 1", "ln(10) = 2.3026", NULL, NULL, NULL},
     "log, exp"},
    
    {"log10", FH_LOG10, 1, 1, FF_SCALAR | FF_COMPLEX, "Logarithmic", "log10(x)",
     "Common logarithm (base 10).",
     {"log10(1) = 0", "log10(10) = 1", "log10(100) = 2", "log10(1000) = 3", NULL, NULL},
     "log, log2"},
    
    {"log2", FH_LOG2, 1, 1, FF_SCALAR | FF_COMPLEX, "Logarithmic", "log2(x)",
     "Binary logarithm (base 2).",
     {"log2(1) = 0", "log2(2) = 1", "log2(8) = 3", "log2(1024) = 10", NULL, NULL},
     "log, log10, exp2"},
    
    {"log1p", FH_LOG1P, 1, 1, FF_SCALAR, "Logarithmic", "log1p(x)",
     "log(1+x). More accurate for small x.",
     {"log1p(0) = 0", "log1p(1) = 0.6931", "log1p(0.001) = 0.0009995", NULL, NULL, NULL},
     "log, expm1"},
    
    {"logb", FH_LOGB, 2, 2, FF_SCALAR, "Logarithmic", "logb(x, b)",
     "Logarithm base b. logb(x, b) = log(x)/log(b).",
     {"logb(8, 2) = 3", "logb(1000, 10) = 3", "logb(27, 3) = 3", NULL, NULL, NULL},
     "log, log10, log2"},
    
    {"reallog", FH_REALLOG, 1, 1, FF_SCALAR | FF_REAL_ONLY, "Logarithmic", "reallog(x)",
     "Real logarithm. Error for non-positive x.",
     {"reallog(1) = 0", "reallog(e) = 1", "reallog(10) = 2.3026", NULL, NULL, NULL},
     "log, realsqrt"},
    
    /* ===== COMPLEX NUMBERS ===== */
    {"real", FH_REAL, 1, 1, FF_SCALAR | FF_COMPLEX, "Complex", "real(z)",
     "Real part of complex number.",
     {"real(3+4i) = 3", "real(5) = 5", "real(2i) = 0", NULL, NULL, NULL},
     "imag, conj, abs, re"},
    
    {"re", FH_REAL, 1, 1, FF_SCALAR | FF_COMPLEX, "Complex", "re(z)",
     "Real part. Alias for real.",
     {"re(3+4i) = 3", NULL, NULL, NULL, NULL, NULL},
     "real, im"},
    
    {"imag", FH_IMAG, 1, 1, FF_SCALAR | FF_COMPLEX, "Complex", "imag(z)",
     "Imaginary part of complex number.",
     {"imag(3+4i) = 4", "imag(5) = 0", "imag(2i) = 2", NULL, NULL, NULL},
     "real, conj, abs, im"},
    
    {"im", FH_IMAG, 1, 1, FF_SCALAR | FF_COMPLEX, "Complex", "im(z)",
     "Imaginary part. Alias for imag.",
     {"im(3+4i) = 4", NULL, NULL, NULL, NULL, NULL},
     "imag, re"},
    
    {"conj", FH_CONJ, 1, 1, FF_SCALAR | FF_COMPLEX, "Complex", "conj(z)",
     "Complex conjugate. conj(a+bi) = a-bi.",
     {"conj(3+4i) = 3-4i", "conj(5) = 5", "conj(2i) = -2i", NULL, NULL, NULL},
     "real, imag, abs"},
    
    {"abs", FH_ABS, 1, 1, FF_SCALAR | FF_COMPLEX, "Complex", "abs(z)",
     "Absolute value (magnitude). For complex z, abs(z) = sqrt(re^2 + im^2).",
     {"abs(-5) = 5", "abs(3+4i) = 5", "abs(-3-4i) = 5", NULL, NULL, NULL},
     "arg, angle, sign"},
    
    {"arg", FH_ARG, 1, 1, FF_SCALAR | FF_COMPLEX, "Complex", "arg(z)",
     "Argument (phase angle) in radians. Range is (-pi, pi].",
     {"arg(1) = 0", "arg(-1) = 3.1416", "arg(i) = 1.5708", "arg(1+i) = 0.7854", NULL, NULL},
     "abs, angle, phase"},
    
    {"angle", FH_ANGLE, 1, 1, FF_SCALAR | FF_COMPLEX, "Complex", "angle(z)",
     "Angle (phase) of complex number. Alias for arg.",
     {"angle(1) = 0", "angle(1+i) = 0.7854", NULL, NULL, NULL, NULL},
     "arg, abs, phase"},
    
    {"phase", FH_PHASE, 1, 1, FF_SCALAR | FF_COMPLEX, "Complex", "phase(z)",
     "Phase angle. Alias for arg.",
     {"phase(1) = 0", "phase(-1) = 3.1416", NULL, NULL, NULL, NULL},
     "arg, angle"},
    
    {"complex", FH_COMPLEX, 2, 2, FF_SCALAR, "Complex", "complex(re, im)",
     "Create complex number from real and imaginary parts.",
     {"complex(3, 4) = 3+4i", "complex(1, 0) = 1", "complex(0, 1) = i", NULL, NULL, NULL},
     "real, imag"},
    
    {"cart2pol", FH_CART2POL, 2, 2, FF_SCALAR, "Complex", "cart2pol(x, y)",
     "Convert Cartesian to polar coordinates. Returns [r, theta].",
     {"cart2pol(3, 4)", "cart2pol(1, 1)", NULL, NULL, NULL, NULL},
     "pol2cart"},
    
    {"pol2cart", FH_POL2CART, 2, 2, FF_SCALAR, "Complex", "pol2cart(r, theta)",
     "Convert polar to Cartesian coordinates. Returns [x, y].",
     {"pol2cart(5, 0.9273)", "pol2cart(1, pi/4)", NULL, NULL, NULL, NULL},
     "cart2pol"},
    
    /* ===== ROUNDING AND SIGN ===== */
    {"floor", FH_FLOOR, 1, 1, FF_SCALAR, "Rounding", "floor(x)",
     "Round toward negative infinity.",
     {"floor(3.7) = 3", "floor(-3.7) = -4", "floor(3) = 3", NULL, NULL, NULL},
     "ceil, round, fix, trunc"},
    
    {"ceil", FH_CEIL, 1, 1, FF_SCALAR, "Rounding", "ceil(x)",
     "Round toward positive infinity.",
     {"ceil(3.2) = 4", "ceil(-3.2) = -3", "ceil(3) = 3", NULL, NULL, NULL},
     "floor, round, fix"},
    
    {"round", FH_ROUND, 1, 1, FF_SCALAR, "Rounding", "round(x)",
     "Round to nearest integer. Ties round to even.",
     {"round(3.4) = 3", "round(3.5) = 4", "round(-3.5) = -4", "round(2.5) = 2", NULL, NULL},
     "floor, ceil, fix"},
    
    {"trunc", FH_TRUNC, 1, 1, FF_SCALAR, "Rounding", "trunc(x)",
     "Round toward zero (truncate).",
     {"trunc(3.7) = 3", "trunc(-3.7) = -3", "trunc(3) = 3", NULL, NULL, NULL},
     "floor, ceil, fix"},
    
    {"fix", FH_FIX, 1, 1, FF_SCALAR, "Rounding", "fix(x)",
     "Round toward zero. Alias for trunc.",
     {"fix(3.7) = 3", "fix(-3.7) = -3", NULL, NULL, NULL, NULL},
     "trunc, floor, ceil"},
    
    {"frac", FH_FRAC, 1, 1, FF_SCALAR, "Rounding", "frac(x)",
     "Fractional part. frac(x) = x - floor(x).",
     {"frac(3.7) = 0.7", "frac(-3.7) = 0.3", "frac(3) = 0", NULL, NULL, NULL},
     "floor, mod"},
    
    {"sign", FH_SIGN, 1, 1, FF_SCALAR, "Rounding", "sign(x)",
     "Sign function. Returns -1, 0, or 1.",
     {"sign(-5) = -1", "sign(0) = 0", "sign(5) = 1", NULL, NULL, NULL},
     "abs, signum"},
    
    {"signum", FH_SIGNUM, 1, 1, FF_SCALAR, "Rounding", "signum(x)",
     "Sign function. Alias for sign.",
     {"signum(-5) = -1", "signum(0) = 0", "signum(5) = 1", NULL, NULL, NULL},
     "sign, abs"},
    
    {"sgn", FH_SIGN, 1, 1, FF_SCALAR, "Rounding", "sgn(x)",
     "Sign function. Alias for sign.",
     {"sgn(-5) = -1", "sgn(0) = 0", "sgn(5) = 1", NULL, NULL, NULL},
     "sign, abs"},
    
    {"neg", FH_NEG, 1, 1, FF_SCALAR | FF_COMPLEX, "Arithmetic", "neg(x)",
     "Negation. neg(x) = -x.",
     {"neg(5) = -5", "neg(-3) = 3", "neg(0) = 0", NULL, NULL, NULL},
     "abs, sign"},
    
    {"copysign", FH_COPYSIGN, 2, 2, FF_SCALAR, "Rounding", "copysign(x, y)",
     "Copy sign. Returns magnitude of x with sign of y.",
     {"copysign(5, -1) = -5", "copysign(-5, 1) = 5", "copysign(3, 0) = 3", NULL, NULL, NULL},
     "sign, abs"},
    
    {"mod", FH_MOD, 2, 2, FF_SCALAR, "Rounding", "mod(x, y)",
     "Modulo operation. Result has same sign as divisor y.",
     {"mod(17, 5) = 2", "mod(-17, 5) = 3", "mod(17, -5) = -3", NULL, NULL, NULL},
     "rem, frac"},
    
    {"rem", FH_REM, 2, 2, FF_SCALAR, "Rounding", "rem(x, y)",
     "Remainder. Result has same sign as dividend x.",
     {"rem(17, 5) = 2", "rem(-17, 5) = -2", "rem(17, -5) = 2", NULL, NULL, NULL},
     "mod, frac"},
    
    {"fmod", FH_FMOD, 2, 2, FF_SCALAR, "Rounding", "fmod(x, y)",
     "Floating-point remainder. Same as rem.",
     {"fmod(5.3, 2) = 1.3", "fmod(-5.3, 2) = -1.3", NULL, NULL, NULL, NULL},
     "rem, mod"},
    
    {"fdim", FH_FDIM, 2, 2, FF_SCALAR, "Rounding", "fdim(x, y)",
     "Positive difference. Returns max(x-y, 0).",
     {"fdim(5, 3) = 2", "fdim(3, 5) = 0", "fdim(-1, -3) = 2", NULL, NULL, NULL},
     "max, abs"},
    
    {"clamp", FH_CLAMP, 3, 3, FF_SCALAR, "Rounding", "clamp(x, lo, hi)",
     "Clamp value to range [lo, hi].",
     {"clamp(5, 0, 10) = 5", "clamp(-5, 0, 10) = 0", "clamp(15, 0, 10) = 10", NULL, NULL, NULL},
     "min, max, clip"},
    
    {"clip", FH_CLIP, 3, 3, FF_SCALAR, "Rounding", "clip(x, lo, hi)",
     "Clip value to range. Alias for clamp.",
     {"clip(5, 0, 10) = 5", "clip(-5, 0, 10) = 0", NULL, NULL, NULL, NULL},
     "clamp, min, max"},
    
    /* ===== NUMBER THEORY ===== */
    {"gcd", FH_GCD, 2, 2, FF_SCALAR, "Number Theory", "gcd(a, b)",
     "Greatest common divisor.",
     {"gcd(12, 8) = 4", "gcd(17, 13) = 1", "gcd(100, 35) = 5", NULL, NULL, NULL},
     "lcm"},
    
    {"lcm", FH_LCM, 2, 2, FF_SCALAR, "Number Theory", "lcm(a, b)",
     "Least common multiple.",
     {"lcm(4, 6) = 12", "lcm(3, 5) = 15", "lcm(12, 8) = 24", NULL, NULL, NULL},
     "gcd, modinv"},
    
    {"modinv", FH_MODINV, 2, 2, FF_SCALAR, "Number Theory", "modinv(a, m)",
     "Modular multiplicative inverse. Returns x where a*x ≡ 1 (mod m).",
     {"modinv(3, 7) = 5", "modinv(2, 5) = 3", "mod(3*5, 7) = 1", NULL, NULL, NULL},
     "mod, gcd"},
    
    {"isprime", FH_ISPRIME, 1, 1, FF_SCALAR, "Number Theory", "isprime(n)",
     "Test if n is prime.",
     {"isprime(2) = true", "isprime(17) = true", "isprime(15) = false", NULL, NULL, NULL},
     "nextprime, prevprime"},
    
    {"nextprime", FH_NEXTPRIME, 1, 1, FF_SCALAR, "Number Theory", "nextprime(n)",
     "Next prime greater than n.",
     {"nextprime(10) = 11", "nextprime(13) = 17", "nextprime(100) = 101", NULL, NULL, NULL},
     "prevprime, isprime"},
    
    {"prevprime", FH_PREVPRIME, 1, 1, FF_SCALAR, "Number Theory", "prevprime(n)",
     "Previous prime less than n.",
     {"prevprime(10) = 7", "prevprime(17) = 13", "prevprime(100) = 97", NULL, NULL, NULL},
     "nextprime, isprime"},
    
    {"fact", FH_FACT, 1, 1, FF_SCALAR, "Number Theory", "fact(n)",
     "Factorial. n! = 1*2*3*...*n.",
     {"fact(0) = 1", "fact(5) = 120", "fact(10) = 3628800", NULL, NULL, NULL},
     "factorial, gamma, ncr"},
    
    {"factorial", FH_FACTORIAL, 1, 1, FF_SCALAR, "Number Theory", "factorial(n)",
     "Factorial. Alias for fact.",
     {"factorial(5) = 120", "factorial(0) = 1", NULL, NULL, NULL, NULL},
     "fact, gamma"},
    
    {"ncr", FH_NCR, 2, 2, FF_SCALAR, "Number Theory", "ncr(n, r)",
     "Binomial coefficient. n choose r.",
     {"ncr(5, 2) = 10", "ncr(10, 3) = 120", "ncr(52, 5) = 2598960", NULL, NULL, NULL},
     "npr, comb, fact"},
    
    {"npr", FH_NPR, 2, 2, FF_SCALAR, "Number Theory", "npr(n, r)",
     "Permutations. n!/(n-r)!.",
     {"npr(5, 2) = 20", "npr(10, 3) = 720", NULL, NULL, NULL, NULL},
     "ncr, perm, fact"},
    
    {"comb", FH_COMB, 2, 2, FF_SCALAR, "Number Theory", "comb(n, r)",
     "Combinations. Alias for ncr.",
     {"comb(5, 2) = 10", NULL, NULL, NULL, NULL, NULL},
     "ncr, perm"},
    
    {"perm", FH_PERM, 2, 2, FF_SCALAR, "Number Theory", "perm(n, r)",
     "Permutations. Alias for npr.",
     {"perm(5, 2) = 20", NULL, NULL, NULL, NULL, NULL},
     "npr, comb"},
    
    {"fibonacci", FH_FIBONACCI, 1, 1, FF_SCALAR, "Number Theory", "fibonacci(n)",
     "nth Fibonacci number. F(0)=0, F(1)=1, F(n)=F(n-1)+F(n-2).",
     {"fibonacci(0) = 0", "fibonacci(1) = 1", "fibonacci(10) = 55", "fibonacci(20) = 6765", NULL, NULL},
     "lucas"},
    
    {"lucas", FH_LUCAS, 1, 1, FF_SCALAR, "Number Theory", "lucas(n)",
     "nth Lucas number. L(0)=2, L(1)=1, L(n)=L(n-1)+L(n-2).",
     {"lucas(0) = 2", "lucas(1) = 1", "lucas(10) = 123", NULL, NULL, NULL},
     "fibonacci"},
    
    {"totient", FH_TOTIENT, 1, 1, FF_SCALAR, "Number Theory", "totient(n)",
     "Euler's totient function. Count of integers <= n coprime to n.",
     {"totient(1) = 1", "totient(10) = 4", "totient(12) = 4", NULL, NULL, NULL},
     "gcd"},
    
    {"even", FH_EVEN, 1, 1, FF_SCALAR, "Number Theory", "even(n)",
     "Test if n is even.",
     {"even(4) = true", "even(7) = false", "even(0) = true", NULL, NULL, NULL},
     "odd"},
    
    {"odd", FH_ODD, 1, 1, FF_SCALAR, "Number Theory", "odd(n)",
     "Test if n is odd.",
     {"odd(7) = true", "odd(4) = false", "odd(1) = true", NULL, NULL, NULL},
     "even, isodd"},
    
    {"iseven", FH_ISEVEN, 1, 1, FF_SCALAR, "Number Theory", "iseven(n)",
     "Test if n is even. Alias for even.",
     {"iseven(4) = true", "iseven(7) = false", NULL, NULL, NULL, NULL},
     "even, isodd"},
    
    {"isodd", FH_ISODD, 1, 1, FF_SCALAR, "Number Theory", "isodd(n)",
     "Test if n is odd. Alias for odd.",
     {"isodd(7) = true", "isodd(4) = false", NULL, NULL, NULL, NULL},
     "odd, iseven"},
    
    /* ===== SPECIAL FUNCTIONS ===== */
    {"gamma", FH_GAMMA, 1, 1, FF_SCALAR, "Special Functions", "gamma(x)",
     "Gamma function. gamma(n) = (n-1)! for positive integers.",
     {"gamma(5) = 24", "gamma(1) = 1", "gamma(0.5) = 1.7725", NULL, NULL, NULL},
     "lgamma, fact, beta"},
    
    {"lgamma", FH_LGAMMA, 1, 1, FF_SCALAR, "Special Functions", "lgamma(x)",
     "Natural log of gamma function. More accurate for large x.",
     {"lgamma(5) = 3.1781", "lgamma(100) = 359.13", NULL, NULL, NULL, NULL},
     "gamma"},
    
    {"beta", FH_BETA, 2, 2, FF_SCALAR, "Special Functions", "beta(a, b)",
     "Beta function. beta(a,b) = gamma(a)*gamma(b)/gamma(a+b).",
     {"beta(2, 3) = 0.0833", "beta(0.5, 0.5) = 3.1416", NULL, NULL, NULL, NULL},
     "gamma"},
    
    {"erf", FH_ERF, 1, 1, FF_SCALAR, "Special Functions", "erf(x)",
     "Error function.",
     {"erf(0) = 0", "erf(1) = 0.8427", "erf(-1) = -0.8427", NULL, NULL, NULL},
     "erfc, normcdf"},
    
    {"erfc", FH_ERFC, 1, 1, FF_SCALAR, "Special Functions", "erfc(x)",
     "Complementary error function. erfc(x) = 1 - erf(x).",
     {"erfc(0) = 1", "erfc(1) = 0.1573", "erfc(3) = 0.00002", NULL, NULL, NULL},
     "erf"},
    
    {"besselj", FH_BESSELJ, 2, 2, FF_SCALAR, "Special Functions", "besselj(n, x)",
     "Bessel function of first kind.",
     {"besselj(0, 0) = 1", "besselj(0, 1) = 0.7652", "besselj(1, 0) = 0", NULL, NULL, NULL},
     "bessely"},
    
    {"bessely", FH_BESSELY, 2, 2, FF_SCALAR, "Special Functions", "bessely(n, x)",
     "Bessel function of second kind.",
     {"bessely(0, 1) = 0.0883", "bessely(1, 1) = -0.7812", NULL, NULL, NULL, NULL},
     "besselj"},
    
    {"zeta", FH_ZETA, 1, 1, FF_SCALAR, "Special Functions", "zeta(s)",
     "Riemann zeta function.",
     {"zeta(2) = 1.6449", "zeta(4) = 1.0823", "zeta(-1) = -0.0833", NULL, NULL, NULL},
     "harmonic"},
    
    {"sigmoid", FH_SIGMOID, 1, 1, FF_SCALAR, "Special Functions", "sigmoid(x)",
     "Logistic sigmoid. sigmoid(x) = 1/(1+exp(-x)).",
     {"sigmoid(0) = 0.5", "sigmoid(1) = 0.7311", "sigmoid(-1) = 0.2689", NULL, NULL, NULL},
     "softplus"},
    
    {"heaviside", FH_HEAVISIDE, 1, 1, FF_SCALAR, "Special Functions", "heaviside(x)",
     "Heaviside step function. 0 for x<0, 0.5 for x=0, 1 for x>0.",
     {"heaviside(-1) = 0", "heaviside(0) = 0.5", "heaviside(1) = 1", NULL, NULL, NULL},
     "step, sign"},
    
    /* ===== STATISTICS ===== */
    {"sum", FH_SUM, 1, -1, FF_MATRIX | FF_VARIADIC, "Statistics", "sum(v)",
     "Sum of elements.",
     {"sum([1,2,3,4,5]) = 15", "sum([10,20,30])", "sum(1:100) = 5050", NULL, NULL, NULL},
     "prod, mean"},
    
    {"prod", FH_PROD, 1, -1, FF_MATRIX | FF_VARIADIC, "Statistics", "prod(v)",
     "Product of elements.",
     {"prod([1,2,3,4,5]) = 120", "prod(1,2,3,4) = 24", NULL, NULL, NULL, NULL},
     "sum, factorial"},
    
    {"mean", FH_MEAN, 1, -1, FF_MATRIX | FF_VARIADIC, "Statistics", "mean(v)",
     "Arithmetic mean.",
     {"mean([1,2,3,4,5]) = 3", "mean(10,20,30) = 20", NULL, NULL, NULL, NULL},
     "median, geomean, harmmean"},
    
    {"median", FH_MEDIAN, 1, -1, FF_MATRIX | FF_VARIADIC, "Statistics", "median(v)",
     "Median (middle value).",
     {"median([1,2,3,4,5]) = 3", "median([1,2,3,4]) = 2.5", NULL, NULL, NULL, NULL},
     "mean, mode"},
    
    {"mode", FH_MODE, 1, -1, FF_MATRIX | FF_VARIADIC, "Statistics", "mode(v)",
     "Mode (most frequent value).",
     {"mode([1,2,2,3]) = 2", "mode([1,1,2,2,3]) = 1", NULL, NULL, NULL, NULL},
     "mean, median"},
    
    {"std", FH_STD, 1, -1, FF_MATRIX | FF_VARIADIC, "Statistics", "std(v)",
     "Standard deviation (sample).",
     {"std([1,2,3,4,5]) = 1.5811", "std([2,4,4,4,5,5,7,9]) = 2.1381", NULL, NULL, NULL, NULL},
     "var, mean"},
    
    {"var", FH_VAR, 1, -1, FF_MATRIX | FF_VARIADIC, "Statistics", "var(v)",
     "Variance (sample).",
     {"var([1,2,3,4,5]) = 2.5", "var([2,4,4,4,5,5,7,9]) = 4.5714", NULL, NULL, NULL, NULL},
     "std, mean"},
    
    {"min", FH_MIN, 1, -1, FF_MATRIX | FF_VARIADIC, "Statistics", "min(v)",
     "Minimum value.",
     {"min([3,1,4,1,5]) = 1", "min(5,2,8,1) = 1", NULL, NULL, NULL, NULL},
     "max, range"},
    
    {"max", FH_MAX, 1, -1, FF_MATRIX | FF_VARIADIC, "Statistics", "max(v)",
     "Maximum value.",
     {"max([3,1,4,1,5]) = 5", "max(5,2,8,1) = 8", NULL, NULL, NULL, NULL},
     "min, range, max2"},
    
    {"min2", FH_MIN2, 2, 2, FF_SCALAR, "Statistics", "min2(a, b)",
     "Minimum of two values.",
     {"min2(3, 5) = 3", "min2(-1, 2) = -1", "min2(7, 7) = 7", NULL, NULL, NULL},
     "min, max2"},
    
    {"max2", FH_MAX2, 2, 2, FF_SCALAR, "Statistics", "max2(a, b)",
     "Maximum of two values.",
     {"max2(3, 5) = 5", "max2(-1, 2) = 2", "max2(7, 7) = 7", NULL, NULL, NULL},
     "max, min2"},
    
    {"sum_cols", FH_SUM_COLS, 1, 1, FF_MATRIX, "Statistics", "sum_cols(A)",
     "Sum of each column. Returns row vector.",
     {"sum_cols([1,2;3,4])", "sum_cols(ones(3,4))", NULL, NULL, NULL, NULL},
     "sum, sum_rows, mean_cols"},
    
    {"sum_rows", FH_SUM_ROWS, 1, 1, FF_MATRIX, "Statistics", "sum_rows(A)",
     "Sum of each row. Returns column vector.",
     {"sum_rows([1,2;3,4])", "sum_rows(ones(3,4))", NULL, NULL, NULL, NULL},
     "sum, sum_cols"},
    
    {"mean_cols", FH_MEAN_COLS, 1, 1, FF_MATRIX, "Statistics", "mean_cols(A)",
     "Mean of each column. Returns row vector.",
     {"mean_cols([1,2;3,4])", "mean_cols(rand(5,3))", NULL, NULL, NULL, NULL},
     "mean, std_cols, sum_cols"},
    
    {"std_cols", FH_STD_COLS, 1, 1, FF_MATRIX, "Statistics", "std_cols(A)",
     "Standard deviation of each column. Returns row vector.",
     {"std_cols([1,2;3,4])", "std_cols(randn(10,3))", NULL, NULL, NULL, NULL},
     "std, mean_cols, var"},
    
    {"range", FH_RANGE, 1, -1, FF_MATRIX | FF_VARIADIC, "Statistics", "range(v)",
     "Range (max - min).",
     {"range([1,2,3,4,5]) = 4", "range([10,20,15]) = 10", NULL, NULL, NULL, NULL},
     "min, max"},
    
    /* ===== PROBABILITY ===== */
    {"normpdf", FH_NORMPDF, 1, 3, FF_SCALAR, "Probability", "normpdf(x) or normpdf(x, mu, sigma)",
     "Normal probability density function.",
     {"normpdf(0) = 0.3989", "normpdf(0, 0, 1) = 0.3989", "normpdf(1, 0, 2) = 0.1760", NULL, NULL, NULL},
     "normcdf, norminv"},
    
    {"normcdf", FH_NORMCDF, 1, 3, FF_SCALAR, "Probability", "normcdf(x) or normcdf(x, mu, sigma)",
     "Normal cumulative distribution function.",
     {"normcdf(0) = 0.5", "normcdf(1.96) = 0.975", "normcdf(-1.96) = 0.025", NULL, NULL, NULL},
     "normpdf, norminv"},
    
    {"norminv", FH_NORMINV, 1, 3, FF_SCALAR, "Probability", "norminv(p) or norminv(p, mu, sigma)",
     "Inverse normal CDF (quantile function).",
     {"norminv(0.5) = 0", "norminv(0.975) = 1.96", "norminv(0.025) = -1.96", NULL, NULL, NULL},
     "normcdf, normpdf"},
    
    {"binopdf", FH_BINOPDF, 3, 3, FF_SCALAR, "Probability", "binopdf(k, n, p)",
     "Binomial probability mass function. P(X = k) for X ~ Binom(n, p).",
     {"binopdf(3, 10, 0.5) = 0.117", "binopdf(0, 5, 0.3) = 0.168", "binopdf(5, 5, 0.5) = 0.031", NULL, NULL, NULL},
     "binocdf, poisspdf"},
    
    {"binocdf", FH_BINOCDF, 3, 3, FF_SCALAR, "Probability", "binocdf(k, n, p)",
     "Binomial cumulative distribution. P(X <= k) for X ~ Binom(n, p).",
     {"binocdf(3, 10, 0.5) = 0.172", "binocdf(5, 10, 0.5) = 0.623", NULL, NULL, NULL, NULL},
     "binopdf, poisscdf"},
    
    {"unifpdf", FH_UNIFPDF, 1, 3, FF_SCALAR, "Probability", "unifpdf(x) or unifpdf(x, a, b)",
     "Uniform probability density on [a, b]. Default [0, 1].",
     {"unifpdf(0.5) = 1", "unifpdf(0.5, 0, 2) = 0.5", "unifpdf(3, 0, 2) = 0", NULL, NULL, NULL},
     "unifcdf, rand"},
    
    {"unifcdf", FH_UNIFCDF, 1, 3, FF_SCALAR, "Probability", "unifcdf(x) or unifcdf(x, a, b)",
     "Uniform cumulative distribution on [a, b]. Default [0, 1].",
     {"unifcdf(0.5) = 0.5", "unifcdf(1, 0, 2) = 0.5", "unifcdf(3, 0, 2) = 1", NULL, NULL, NULL},
     "unifpdf, rand"},
    
    {"fpdf", FH_FPDF, 3, 3, FF_SCALAR, "Probability", "fpdf(x, d1, d2)",
     "F-distribution probability density function. Used in ANOVA and variance tests.",
     {"fpdf(1, 5, 10) = 0.646", "fpdf(2, 5, 10) = 0.162", "fpdf(0, 5, 10) = 0", NULL, NULL, NULL},
     "fcdf, finv"},
    
    {"fcdf", FH_FCDF, 3, 3, FF_SCALAR, "Probability", "fcdf(x, d1, d2)",
     "F-distribution CDF. P(F <= x) for F ~ F(d1, d2).",
     {"fcdf(1, 5, 10) = 0.543", "fcdf(2, 5, 10) = 0.836", "fcdf(4.06, 5, 10) = 0.975", NULL, NULL, NULL},
     "fpdf, finv, anova"},
    
    {"rand", FH_RAND, 0, 2, FF_SCALAR | FF_NO_PARENS, "Probability", "rand or rand(n) or rand(m,n)",
     "Uniform random numbers in [0,1).",
     {"rand", "rand(3)", "rand(2,3)", NULL, NULL, NULL},
     "randn, randi"},
    
    {"randn", FH_RANDN, 0, 2, FF_SCALAR | FF_NO_PARENS, "Probability", "randn or randn(n) or randn(m,n)",
     "Standard normal random numbers (mean=0, std=1).",
     {"randn", "randn(3)", "randn(2,3)", NULL, NULL, NULL},
     "rand, randi"},
    
    {"randi", FH_RANDI, 1, 3, FF_SCALAR, "Probability", "randi(imax) or randi(imax, n) or randi([imin,imax], n)",
     "Random integers.",
     {"randi(10)", "randi(6, 5)", "randi([1,6], 10)", NULL, NULL, NULL},
     "rand, randn"},
    
    /* ===== LINEAR ALGEBRA ===== */
    {"det", FH_DET, 1, 1, FF_MATRIX, "Linear Algebra", "det(A)",
     "Matrix determinant.",
     {"det([1,2;3,4]) = -2", "det(eye(3)) = 1", NULL, NULL, NULL, NULL},
     "inv, trace"},
    
    {"inv", FH_INV, 1, 1, FF_MATRIX, "Linear Algebra", "inv(A)",
     "Matrix inverse.",
     {"inv([1,2;3,4])", "inv(eye(3)) = eye(3)", NULL, NULL, NULL, NULL},
     "det, pinv, linsolve"},
    
    {"trace", FH_TRACE, 1, 1, FF_MATRIX, "Linear Algebra", "trace(A)",
     "Sum of diagonal elements.",
     {"trace(eye(5)) = 5", "trace([1,2;3,4]) = 5", NULL, NULL, NULL, NULL},
     "det, diag"},
    
    {"trans", FH_TRANS, 1, 1, FF_MATRIX, "Linear Algebra", "trans(A)",
     "Matrix transpose.",
     {"trans([1,2,3])", "trans([1,2;3,4])", NULL, NULL, NULL, NULL},
     "transpose"},
    
    {"transpose", FH_TRANSPOSE, 1, 1, FF_MATRIX, "Linear Algebra", "transpose(A)",
     "Matrix transpose. Alias for trans.",
     {"transpose([1,2,3])", NULL, NULL, NULL, NULL, NULL},
     "trans, ctranspose"},
    
    {"ctranspose", FH_CTRANSPOSE, 1, 1, FF_MATRIX | FF_COMPLEX, "Linear Algebra", "ctranspose(A)",
     "Complex conjugate transpose. A' in MATLAB notation.",
     {"ctranspose([1+i,2-i])", "ctranspose([1,2;3,4])", NULL, NULL, NULL, NULL},
     "transpose, conj"},
    
    {"norm", FH_NORM, 1, 2, FF_MATRIX, "Linear Algebra", "norm(A) or norm(A, p)",
     "Matrix or vector norm.",
     {"norm([3,4]) = 5", "norm([1,2,2]) = 3", NULL, NULL, NULL, NULL},
     "norm1, norminf, cond"},
    
    {"norm1", FH_NORM1, 1, 1, FF_MATRIX, "Linear Algebra", "norm1(A)",
     "1-norm (maximum column sum for matrices, sum of abs for vectors).",
     {"norm1([1,-2,3]) = 6", "norm1([1,2;3,4]) = 6", NULL, NULL, NULL, NULL},
     "norm, norminf"},
    
    {"norminf", FH_NORMINF, 1, 1, FF_MATRIX, "Linear Algebra", "norminf(A)",
     "Infinity norm (maximum row sum for matrices, max abs for vectors).",
     {"norminf([1,-2,3]) = 3", "norminf([1,2;3,4]) = 7", NULL, NULL, NULL, NULL},
     "norm, norm1"},
    
    {"rank", FH_RANK, 1, 1, FF_MATRIX, "Linear Algebra", "rank(A)",
     "Matrix rank.",
     {"rank([1,0;0,1]) = 2", "rank([1,2;2,4]) = 1", NULL, NULL, NULL, NULL},
     "det, null"},
    
    {"eig", FH_EIG, 1, 1, FF_MATRIX, "Linear Algebra", "eig(A)",
     "Eigenvalues (and eigenvectors with [V,D]=eig(A)).",
     {"eig([1,2;2,1])", "eig(eye(3))", NULL, NULL, NULL, NULL},
     "svd"},
    
    {"svd", FH_SVD, 1, 1, FF_MATRIX, "Linear Algebra", "svd(A) or [U,S,V]=svd(A)",
     "Singular value decomposition. Returns singular values as vector, or "
     "full decomposition with [U,S,V]=svd(A) where A=U*S*V'.",
     {"svd([1,2;3,4])", "svd(eye(3))", NULL, NULL, NULL, NULL},
     "eig, qr"},
    
    {"qr", FH_QR, 1, 1, FF_MATRIX, "Linear Algebra", "qr(A) or [Q,R]=qr(A)",
     "QR decomposition. Returns R matrix, or full decomposition with "
     "[Q,R]=qr(A) where A=Q*R and Q is orthogonal.",
     {"qr([1,2;3,4])", "qr(eye(3))", NULL, NULL, NULL, NULL},
     "lu, svd"},
    
    {"lu", FH_LU, 1, 1, FF_MATRIX, "Linear Algebra", "[L,U]=lu(A) or [L,U,P]=lu(A)",
     "LU decomposition with partial pivoting. Returns L and U matrices, "
     "optionally with permutation matrix P where P*A=L*U.",
     {"lu([1,2;3,4])", "lu([4,3;6,3])", NULL, NULL, NULL, NULL},
     "qr, chol"},
    
    {"chol", FH_CHOL, 1, 1, FF_MATRIX, "Linear Algebra", "chol(A)",
     "Cholesky decomposition. Returns lower triangular L where A=L*L'.",
     {"chol([4,2;2,2])", "chol(eye(3))", NULL, NULL, NULL, NULL},
     "lu, qr"},
    
    {"pinv", FH_PINV, 1, 1, FF_MATRIX, "Linear Algebra", "pinv(A)",
     "Moore-Penrose pseudoinverse.",
     {"pinv([1,2;3,4])", "pinv([1;2;3])", NULL, NULL, NULL, NULL},
     "inv"},
    
    {"linsolve", FH_LINSOLVE, 2, 2, FF_MATRIX, "Linear Algebra", "linsolve(A, b)",
     "Solve Ax = b.",
     {"linsolve([1,0;0,1], [3;4])", NULL, NULL, NULL, NULL, NULL},
     "inv, mldivide"},
    
    {"rref", FH_RREF, 1, 1, FF_MATRIX, "Linear Algebra", "rref(A)",
     "Row Reduced Echelon Form (Gauss-Jordan elimination).",
     {"rref([1,2,3;4,5,6;7,8,9])", "rref([1,2;3,4])", "rref([2,1,-1,8;-3,-1,2,-11;-2,1,2,-3])", NULL, NULL, NULL},
     "lu, inv, rank"},
    
    {"orth", FH_ORTH, 1, 1, FF_MATRIX, "Linear Algebra", "orth(A)",
     "Orthonormal basis for column space of A (via SVD).",
     {"orth([1,0;1,0;0,1])", "orth([1,2;3,4])", NULL, NULL, NULL, NULL},
     "null, svd, rank"},
    
    {"poly", FH_POLY, 1, 1, FF_MATRIX, "Linear Algebra", "poly(A)",
     "Characteristic polynomial coefficients from square matrix eigenvalues.",
     {"poly([1,2;3,4])", "poly([1,0;0,2])", "poly(eye(3))", NULL, NULL, NULL},
     "roots, eig, det"},
    
    {"roots", FH_ROOTS, 1, 1, FF_MATRIX, "Polynomials", "roots(p)",
     "Find polynomial roots. p = [a_n, a_{n-1}, ..., a_1, a_0] for a_n*x^n + ... + a_0.",
     {"roots([1,0,-4])", "roots([1,-6,11,-6])", "roots([1,0,0,-1])", NULL, NULL, NULL},
     "poly, polyval, eig"},
    
    {"expm", FH_EXPM, 1, 1, FF_MATRIX, "Linear Algebra", "expm(A)",
     "Matrix exponential via Pade approximation. e^A for square matrix A.",
     {"expm([0,1;0,0])", "expm([1,0;0,1])", "expm(zeros(2))", NULL, NULL, NULL},
     "logm, exp, eig"},
    
    {"logm", FH_LOGM, 1, 1, FF_MATRIX, "Linear Algebra", "logm(A)",
     "Matrix logarithm. Returns B where e^B = A.",
     {"logm(eye(2))", "logm([2.7183,0;0,2.7183])", NULL, NULL, NULL, NULL},
     "expm, log, eig"},
    
    {"hess", FH_HESS, 1, 1, FF_MATRIX, "Linear Algebra", "hess(A)",
     "Hessenberg form via Householder reflections. Upper Hessenberg has zeros below first subdiagonal.",
     {"hess([1,2,3;4,5,6;7,8,9])", "hess(rand(4,4))", NULL, NULL, NULL, NULL},
     "schur, eig, qr"},
    
    {"balance", FH_BALANCE, 1, 1, FF_MATRIX, "Linear Algebra", "balance(A)",
     "Balance matrix for eigenvalue computation. Scales rows/columns for better numerical stability.",
     {"balance([1,1000;0.001,1])", "balance(rand(3,3))", NULL, NULL, NULL, NULL},
     "eig, hess, schur"},
    
    {"fft", FH_FFT, 1, 1, FF_MATRIX, "Signal Processing", "fft(x)",
     "Fast Fourier Transform (DFT). Transforms signal from time to frequency domain.",
     {"fft([1,0,0,0])", "fft([1,1,1,1])", "fft([1,2,3,4])", "ifft(fft([1,2,3,4]))", NULL, NULL},
     "ifft, abs, angle"},
    
    {"ifft", FH_IFFT, 1, 1, FF_MATRIX, "Signal Processing", "ifft(X)",
     "Inverse FFT. Transforms from frequency domain back to time domain.",
     {"ifft(fft([1,2,3,4]))", "ifft([4,0,0,0])", NULL, NULL, NULL, NULL},
     "fft"},
    
    {"xcorr", FH_XCORR, 1, 2, FF_MATRIX, "Signal Processing", "xcorr(x) or xcorr(x, y)",
     "Cross-correlation. xcorr(x) is autocorrelation. Measures similarity at different lags.",
     {"xcorr([1,2,3])", "xcorr([1,0,1],[1,1,0])", "xcorr([1,2,3,4],[1,0,0,0])", NULL, NULL, NULL},
     "xcov, conv, fft"},
    
    {"xcov", FH_XCOV, 1, 2, FF_MATRIX, "Signal Processing", "xcov(x) or xcov(x, y)",
     "Cross-covariance. Like xcorr but mean-removed (centered signals).",
     {"xcov([1,2,3])", "xcov([1,2,3],[3,2,1])", NULL, NULL, NULL, NULL},
     "xcorr, cov, corrcoef"},
    
    {"meshgrid", FH_MESHGRID, 2, 2, FF_MATRIX, "Matrix", "meshgrid(x, y)",
     "Create 2D grid matrices. Returns X where each row is a copy of x vector.",
     {"meshgrid([1,2,3],[10,20])", "meshgrid(1:3,1:2)", "meshgrid(linspace(0,1,3),linspace(0,1,2))", NULL, NULL, NULL},
     "linspace, repmat"},
    
    {"histcounts", FH_HISTCOUNTS, 1, 2, FF_MATRIX, "Statistics", "histcounts(x) or histcounts(x, nbins)",
     "Histogram bin counts. Default 10 bins from min to max.",
     {"histcounts([1,2,2,3,3,3])", "histcounts([1,2,2,3,3,3], 3)", "histcounts(randn(1,100), 20)", NULL, NULL, NULL},
     "hist, mean, std"},
    
    {"isoutlier", FH_ISOUTLIER, 1, 1, FF_MATRIX, "Statistics", "isoutlier(x)",
     "Detect outliers using MAD (median absolute deviation). Returns 1 for outliers, 0 otherwise.",
     {"isoutlier([1,2,3,100,4,5])", "isoutlier([1,1,1,1,10])", "sum(isoutlier(randn(1,100)))", NULL, NULL, NULL},
     "median, std, mad"},
    
    /* ===== MATRIX CREATION ===== */
    {"zeros", FH_ZEROS, 1, 2, FF_MATRIX, "Matrix", "zeros(n) or zeros(m, n)",
     "Matrix of zeros.",
     {"zeros(3)", "zeros(2, 3)", NULL, NULL, NULL, NULL},
     "ones, eye"},
    
    {"ones", FH_ONES, 1, 2, FF_MATRIX, "Matrix", "ones(n) or ones(m, n)",
     "Matrix of ones.",
     {"ones(3)", "ones(2, 3)", NULL, NULL, NULL, NULL},
     "zeros, eye"},
    
    {"eye", FH_EYE, 1, 2, FF_MATRIX, "Matrix", "eye(n) or eye(m, n)",
     "Identity matrix.",
     {"eye(3)", "eye(2, 3)", NULL, NULL, NULL, NULL},
     "zeros, ones, diag"},
    
    {"diag", FH_DIAG, 1, 2, FF_MATRIX, "Matrix", "diag(v) or diag(A, k)",
     "Create diagonal matrix or extract diagonal.",
     {"diag([1,2,3])", "diag([1,2;3,4])", "diag([1,2;3,4], 1)", NULL, NULL, NULL},
     "eye, trace"},
    
    {"linspace", FH_LINSPACE, 2, 3, FF_MATRIX, "Matrix", "linspace(a, b) or linspace(a, b, n)",
     "Linearly spaced vector.",
     {"linspace(0, 1, 5)", "linspace(1, 10)", NULL, NULL, NULL, NULL},
     "logspace"},
    
    {"logspace", FH_LOGSPACE, 2, 3, FF_MATRIX, "Matrix", "logspace(a, b) or logspace(a, b, n)",
     "Logarithmically spaced vector.",
     {"logspace(0, 2, 5)", "logspace(1, 3)", NULL, NULL, NULL, NULL},
     "linspace"},
    
    {"reshape", FH_RESHAPE, 3, 3, FF_MATRIX, "Matrix", "reshape(A, m, n)",
     "Reshape matrix to m x n.",
     {"reshape([1,2,3,4,5,6], 2, 3)", "reshape([1;2;3;4], 2, 2)", NULL, NULL, NULL, NULL},
     "size"},
    
    {"sort", FH_SORT, 1, 2, FF_MATRIX, "Matrix", "sort(v) or sort(A, dim)",
     "Sort elements.",
     {"sort([3,1,4,1,5])", "sort([3;1;2])", NULL, NULL, NULL, NULL},
     "unique, min, max"},
    
    {"unique", FH_UNIQUE, 1, 1, FF_MATRIX, "Matrix", "unique(v)",
     "Unique elements sorted.",
     {"unique([1,2,1,3,2,1])", "unique([3,1,4,1,5,9])", NULL, NULL, NULL, NULL},
     "sort"},
    
    {"flip", FH_FLIP, 1, 2, FF_MATRIX, "Matrix", "flip(A) or flip(A, dim)",
     "Flip array.",
     {"flip([1,2,3])", "flip([1;2;3])", NULL, NULL, NULL, NULL},
     "fliplr, flipud, rot90"},
    
    {"fliplr", FH_FLIPLR, 1, 1, FF_MATRIX, "Matrix", "fliplr(A)",
     "Flip left-right.",
     {"fliplr([1,2,3])", "fliplr([1,2;3,4])", NULL, NULL, NULL, NULL},
     "flipud, flip"},
    
    {"flipud", FH_FLIPUD, 1, 1, FF_MATRIX, "Matrix", "flipud(A)",
     "Flip up-down.",
     {"flipud([1;2;3])", "flipud([1,2;3,4])", NULL, NULL, NULL, NULL},
     "fliplr, flip"},
    
    {"rot90", FH_ROT90, 1, 2, FF_MATRIX, "Matrix", "rot90(A) or rot90(A, k)",
     "Rotate matrix 90 degrees counter-clockwise.",
     {"rot90([1,2;3,4])", "rot90([1,2;3,4], 2)", NULL, NULL, NULL, NULL},
     "flip, trans"},
    
    {"dot", FH_DOT, 2, 2, FF_MATRIX, "Matrix", "dot(a, b)",
     "Dot product of vectors.",
     {"dot([1,2,3], [4,5,6]) = 32", "dot([1,0,0], [0,1,0]) = 0", NULL, NULL, NULL, NULL},
     "cross, norm"},
    
    {"cross", FH_CROSS, 2, 2, FF_MATRIX, "Matrix", "cross(a, b)",
     "Cross product of 3-element vectors.",
     {"cross([1,0,0], [0,1,0])", "cross([1,2,3], [4,5,6])", NULL, NULL, NULL, NULL},
     "dot"},
    
    {"kron", FH_KRON, 2, 2, FF_MATRIX, "Matrix", "kron(A, B)",
     "Kronecker tensor product.",
     {"kron([1,2], [3,4])", "kron(eye(2), [1,2;3,4])", NULL, NULL, NULL, NULL},
     "repmat"},
    
    /* ===== MATRIX QUERY ===== */
    {"size", FH_SIZE, 1, 2, FF_MATRIX, "Matrix Query", "size(A) or size(A, dim)",
     "Matrix dimensions.",
     {"size([1,2,3])", "size([1,2;3,4])", "size([1,2;3,4], 1)", NULL, NULL, NULL},
     "length, numel, rows, cols"},
    
    {"length", FH_LENGTH, 1, 1, FF_MATRIX, "Matrix Query", "length(v)",
     "Length of vector (max dimension).",
     {"length([1,2,3]) = 3", "length([1,2;3,4]) = 2", NULL, NULL, NULL, NULL},
     "size, numel"},
    
    {"rows", FH_ROWS, 1, 1, FF_MATRIX, "Matrix Query", "rows(A)",
     "Number of rows.",
     {"rows([1,2,3]) = 1", "rows([1;2;3]) = 3", NULL, NULL, NULL, NULL},
     "cols, size"},
    
    {"cols", FH_COLS, 1, 1, FF_MATRIX, "Matrix Query", "cols(A)",
     "Number of columns.",
     {"cols([1,2,3]) = 3", "cols([1;2;3]) = 1", NULL, NULL, NULL, NULL},
     "rows, size"},
    
    {"numel", FH_NUMEL, 1, 1, FF_MATRIX, "Matrix Query", "numel(A)",
     "Number of elements.",
     {"numel([1,2,3]) = 3", "numel([1,2;3,4]) = 4", NULL, NULL, NULL, NULL},
     "size, length"},
    
    {"isempty", FH_ISEMPTY, 1, 1, FF_MATRIX, "Matrix Query", "isempty(A)",
     "Test if empty.",
     {"isempty([]) = true", "isempty([1]) = false", NULL, NULL, NULL, NULL},
     "isscalar, isvector"},
    
    {"isscalar", FH_ISSCALAR, 1, 1, FF_MATRIX, "Matrix Query", "isscalar(A)",
     "Test if scalar (1x1).",
     {"isscalar(5) = true", "isscalar([1,2]) = false", NULL, NULL, NULL, NULL},
     "isvector, ismatrix"},
    
    {"isvector", FH_ISVECTOR, 1, 1, FF_MATRIX, "Matrix Query", "isvector(A)",
     "Test if vector (row or column).",
     {"isvector([1,2,3]) = true", "isvector([1;2;3]) = true", "isvector([1,2;3,4]) = false", NULL, NULL, NULL},
     "isrow, iscolumn, isscalar"},
    
    {"issquare", FH_ISSQUARE, 1, 1, FF_MATRIX, "Matrix Query", "issquare(A)",
     "Test if square matrix.",
     {"issquare(eye(3)) = true", "issquare([1,2,3]) = false", NULL, NULL, NULL, NULL},
     "ismatrix, size"},
    
    {"nnz", FH_NNZ, 1, 1, FF_MATRIX, "Matrix Query", "nnz(A)",
     "Number of nonzero elements.",
     {"nnz([1,0,2,0,3]) = 3", "nnz(eye(3)) = 3", NULL, NULL, NULL, NULL},
     "find"},
    
    /* ===== SIGNAL PROCESSING ===== */
    {"diff", FH_DIFF, 1, 2, FF_MATRIX, "Signal Processing", "diff(v) or diff(v, n)",
     "Differences between adjacent elements.",
     {"diff([1,4,9,16])", "diff([1,2,4,7,11])", NULL, NULL, NULL, NULL},
     "cumsum, gradient"},
    
    {"cumsum", FH_CUMSUM, 1, 1, FF_MATRIX, "Signal Processing", "cumsum(v)",
     "Cumulative sum.",
     {"cumsum([1,2,3,4])", "cumsum([1;2;3])", NULL, NULL, NULL, NULL},
     "sum, diff, cumprod"},
    
    {"cumprod", FH_CUMPROD, 1, 1, FF_MATRIX, "Signal Processing", "cumprod(v)",
     "Cumulative product.",
     {"cumprod([1,2,3,4])", "cumprod([1;2;3])", NULL, NULL, NULL, NULL},
     "prod, cumsum"},
    
    {"conv", FH_CONV, 2, 2, FF_MATRIX, "Signal Processing", "conv(a, b)",
     "Convolution of vectors.",
     {"conv([1,2,3], [1,1])", "conv([1,0,1], [1,2,1])", NULL, NULL, NULL, NULL},
     "deconv"},
    
    {"trapz", FH_TRAPZ, 1, 2, FF_MATRIX, "Signal Processing", "trapz(y) or trapz(x, y)",
     "Trapezoidal numerical integration.",
     {"trapz([1,2,3,4]) = 7.5", "trapz([0,1,2], [0,1,4]) = 2.5", NULL, NULL, NULL, NULL},
     "sum, cumsum"},
    
    /* ===== ANGLE CONVERSION ===== */
    {"deg2rad", FH_DEG2RAD, 1, 1, FF_SCALAR, "Angle Conversion", "deg2rad(x)",
     "Convert degrees to radians.",
     {"deg2rad(180) = 3.1416", "deg2rad(90) = 1.5708", "deg2rad(360) = 6.2832", NULL, NULL, NULL},
     "rad2deg"},
    
    {"rad2deg", FH_RAD2DEG, 1, 1, FF_SCALAR, "Angle Conversion", "rad2deg(x)",
     "Convert radians to degrees.",
     {"rad2deg(pi) = 180", "rad2deg(pi/2) = 90", "rad2deg(2*pi) = 360", NULL, NULL, NULL},
     "deg2rad"},
    
    /* ===== LOGICAL ===== */
    {"any", FH_ANY, 1, 1, FF_MATRIX, "Logical", "any(v)",
     "True if any element is nonzero.",
     {"any([0,0,1,0]) = true", "any([0,0,0]) = false", NULL, NULL, NULL, NULL},
     "all"},
    
    {"all", FH_ALL, 1, 1, FF_MATRIX, "Logical", "all(v)",
     "True if all elements are nonzero.",
     {"all([1,2,3,4]) = true", "all([1,0,1]) = false", NULL, NULL, NULL, NULL},
     "any"},
    
    {"isnan", FH_ISNAN, 1, 1, FF_SCALAR, "Logical", "isnan(x)",
     "Test for NaN.",
     {"isnan(0/0) = true", "isnan(1) = false", "isnan(Inf) = false", NULL, NULL, NULL},
     "isinf, isfinite"},
    
    {"isinf", FH_ISINF, 1, 1, FF_SCALAR, "Logical", "isinf(x)",
     "Test for infinity.",
     {"isinf(1/0) = true", "isinf(-Inf) = true", "isinf(1) = false", NULL, NULL, NULL},
     "isnan, isfinite"},
    
    {"isfinite", FH_ISFINITE, 1, 1, FF_SCALAR, "Logical", "isfinite(x)",
     "Test for finite value.",
     {"isfinite(5) = true", "isfinite(Inf) = false", "isfinite(NaN) = false", NULL, NULL, NULL},
     "isnan, isinf"},
    
    {"isreal", FH_ISREAL, 1, 1, FF_SCALAR, "Logical", "isreal(x)",
     "Test if real (no imaginary part).",
     {"isreal(5) = true", "isreal(3+4i) = false", "isreal(3+0i) = true", NULL, NULL, NULL},
     "real, imag"},
    
    /* ===== BITWISE ===== */
    {"bitand", FH_BITAND, 2, 2, FF_SCALAR, "Bitwise", "bitand(a, b)",
     "Bitwise AND.",
     {"bitand(12, 10) = 8", "bitand(15, 7) = 7", NULL, NULL, NULL, NULL},
     "bitor, bitxor"},
    
    {"bitor", FH_BITOR, 2, 2, FF_SCALAR, "Bitwise", "bitor(a, b)",
     "Bitwise OR.",
     {"bitor(12, 10) = 14", "bitor(8, 4) = 12", NULL, NULL, NULL, NULL},
     "bitand, bitxor"},
    
    {"bitxor", FH_BITXOR, 2, 2, FF_SCALAR, "Bitwise", "bitxor(a, b)",
     "Bitwise XOR.",
     {"bitxor(12, 10) = 6", "bitxor(15, 15) = 0", NULL, NULL, NULL, NULL},
     "bitand, bitor"},
    
    {"bitshift", FH_BITSHIFT, 2, 2, FF_SCALAR, "Bitwise", "bitshift(a, k)",
     "Bit shift. Positive k shifts left, negative shifts right.",
     {"bitshift(1, 4) = 16", "bitshift(16, -2) = 4", NULL, NULL, NULL, NULL},
     "shl, shr"},
    
    {"shl", FH_SHL, 2, 2, FF_SCALAR, "Bitwise", "shl(a, k)",
     "Shift left.",
     {"shl(1, 8) = 256", "shl(3, 4) = 48", NULL, NULL, NULL, NULL},
     "shr, bitshift"},
    
    {"shr", FH_SHR, 2, 2, FF_SCALAR, "Bitwise", "shr(a, k)",
     "Shift right.",
     {"shr(256, 8) = 1", "shr(48, 4) = 3", NULL, NULL, NULL, NULL},
     "shl, bitshift"},
    
    /* ===== UTILITY ===== */
    {"nextpow2", FH_NEXTPOW2, 1, 1, FF_SCALAR, "Utility", "nextpow2(n)",
     "Next power of 2 exponent. 2^nextpow2(n) >= n.",
     {"nextpow2(5) = 3", "nextpow2(8) = 3", "nextpow2(9) = 4", NULL, NULL, NULL},
     "log2, pow2"},
    
    {"normalize", FH_NORMALIZE, 1, 1, FF_MATRIX, "Utility", "normalize(v)",
     "Normalize vector to unit length.",
     {"normalize([3,4])", "norm(normalize([1,2,3])) = 1", NULL, NULL, NULL, NULL},
     "norm"},
    
    {"rescale", FH_RESCALE, 3, 3, FF_MATRIX, "Utility", "rescale(x, lo, hi)",
     "Rescale values to range [lo, hi].",
     {"rescale([1,2,3,4,5], 0, 1)", "rescale([0,50,100], 0, 10)", NULL, NULL, NULL, NULL},
     "normalize, clamp"},
    
    {"disp", FH_DISP, 1, 1, FF_SCALAR | FF_MATRIX, "Utility", "disp(x)",
     "Display value.",
     {"disp(123)", "disp([1,2,3])", NULL, NULL, NULL, NULL},
     "format"},
    
    {"lerp", FH_LERP, 3, 3, FF_SCALAR, "Utility", "lerp(a, b, t)",
     "Linear interpolation. lerp(a,b,t) = a + t*(b-a).",
     {"lerp(0, 10, 0.5) = 5", "lerp(0, 100, 0.25) = 25", NULL, NULL, NULL, NULL},
     "interp1"},
    
    {"tic", FH_TIC, 0, 0, FF_CONST, "Utility", "tic",
     "Start timer. Use with toc to measure elapsed time.",
     {"tic", "tic; sum(1:10000); toc", NULL, NULL, NULL, NULL},
     "toc"},
    
    {"toc", FH_TOC, 0, 0, FF_CONST, "Utility", "toc",
     "Stop timer and display elapsed time since last tic.",
     {"toc", "tic; A=rand(100,100); B=A*A; toc", NULL, NULL, NULL, NULL},
     "tic"},
    
    {"ver", FH_VER, 0, 0, FF_CONST, "Utility", "ver",
     "Display version information.",
     {"ver", NULL, NULL, NULL, NULL, NULL},
     "precision, help"},
    
    {"clc", FH_CLC, 0, 0, FF_CONST, "Utility", "clc",
     "Clear command window.",
     {"clc", NULL, NULL, NULL, NULL, NULL},
     "clear"},
    
    {"clear", FH_CLEAR, 0, 0, FF_CONST, "Utility", "clear",
     "Clear all variables from workspace.",
     {"clear", NULL, NULL, NULL, NULL, NULL},
     "clc"},
    
    /* ===== SESSION COMMANDS ===== */
    {"help", FH_HELP, 0, 1, FF_CONST, "Commands", "help [name]",
     "Show help for a function or command. 'help' alone lists categories.",
     {"help", "help sin", "help plot", "help format", NULL, NULL},
     "demo, funcs"},
    
    {"demo", FH_DEMO, 1, 1, FF_CONST, "Commands", "demo name",
     "Run interactive demonstration of a function.",
     {"demo sin", "demo fft", "demo roots", NULL, NULL, NULL},
     "help, bench"},
    
    {"bench", FH_BENCH, 0, 0, FF_CONST, "Commands", "bench",
     "Run performance benchmarks on core functions.",
     {"bench", NULL, NULL, NULL, NULL, NULL},
     "tic, toc"},
    
    {"funcs", FH_FUNCS, 0, 0, FF_CONST, "Commands", "funcs",
     "List all available functions by category.",
     {"funcs", NULL, NULL, NULL, NULL, NULL},
     "help, list"},
    
    {"list", FH_FUNCS, 0, 0, FF_CONST, "Commands", "list",
     "Alias for funcs. List all available functions.",
     {"list", NULL, NULL, NULL, NULL, NULL},
     "funcs, vars"},
    
    {"vars", FH_VARS, 0, 0, FF_CONST, "Commands", "vars",
     "List all defined variables.",
     {"vars", "x=5; vars", NULL, NULL, NULL, NULL},
     "clear, funcs"},
    
    {"quit", FH_QUIT, 0, 0, FF_CONST, "Commands", "quit",
     "Exit the calculator.",
     {"quit", NULL, NULL, NULL, NULL, NULL},
     "exit"},
    
    {"exit", FH_EXIT, 0, 0, FF_CONST, "Commands", "exit",
     "Exit the calculator. Alias for quit.",
     {"exit", NULL, NULL, NULL, NULL, NULL},
     "quit"},
    
    {"format", FH_FORMAT, 0, 1, FF_CONST, "Commands", "format [mode]",
     "Set output format. Modes: sci, eng, fix, std, hex, bin, oct, short, long.",
     {"format sci", "format eng", "format hex", "format std", "format short", "format long"},
     "digits, precision"},
    
    {"sci", FH_FORMAT_SCI, 0, 0, FF_CONST, "Commands", "sci",
     "Scientific notation format (e.g., 1.23e+04).",
     {"sci", "12345; sci", NULL, NULL, NULL, NULL},
     "format, eng"},
    
    {"eng", FH_FORMAT_ENG, 0, 0, FF_CONST, "Commands", "eng",
     "Engineering notation (exponents are multiples of 3).",
     {"eng", "12345; eng", NULL, NULL, NULL, NULL},
     "format, sci"},
    
    {"hex", FH_FORMAT_HEX, 0, 0, FF_CONST, "Commands", "hex",
     "Hexadecimal output format.",
     {"hex", "255; hex", "hex; 0xff", NULL, NULL, NULL},
     "format, bin, oct"},
    
    {"bin", FH_FORMAT_BIN, 0, 0, FF_CONST, "Commands", "bin",
     "Binary output format.",
     {"bin", "15; bin", "bin; 0b1111", NULL, NULL, NULL},
     "format, hex, oct"},
    
    {"oct", FH_FORMAT_OCT, 0, 0, FF_CONST, "Commands", "oct",
     "Octal output format.",
     {"oct", "64; oct", "oct; 0o100", NULL, NULL, NULL},
     "format, hex, bin"},
    
    {"short", FH_FORMAT_SHORT, 0, 0, FF_CONST, "Commands", "short",
     "Short display format (fewer digits).",
     {"short", "pi; short", NULL, NULL, NULL, NULL},
     "format, long"},
    
    {"long", FH_FORMAT_LONG, 0, 0, FF_CONST, "Commands", "long",
     "Long display format (more digits).",
     {"long", "pi; long", NULL, NULL, NULL, NULL},
     "format, short"},
    
    {"digits", FH_FORMAT_DIGITS, 1, 1, FF_CONST, "Commands", "digits n",
     "Set number of display digits (1-50).",
     {"digits 10", "digits 20", "pi; digits 30", NULL, NULL, NULL},
     "format, precision"},
    
    {"rad", FH_ANGLE_RAD, 0, 0, FF_CONST, "Commands", "rad",
     "Set angle mode to radians (default).",
     {"rad", "rad; sin(pi/2)", NULL, NULL, NULL, NULL},
     "deg, grad"},
    
    {"deg", FH_ANGLE_DEG, 0, 0, FF_CONST, "Commands", "deg",
     "Set angle mode to degrees.",
     {"deg", "deg; sin(90)", NULL, NULL, NULL, NULL},
     "rad, grad"},
    
    {"grad", FH_ANGLE_GRAD, 0, 0, FF_CONST, "Commands", "grad",
     "Set angle mode to gradians (400 per circle).",
     {"grad", "grad; sin(100)", NULL, NULL, NULL, NULL},
     "rad, deg"},
    
    {"save", FH_SAVE, 1, 1, FF_CONST, "Commands", "save filename",
     "Save current session (variables and functions) to file.",
     {"save session.sc", "save mywork", NULL, NULL, NULL, NULL},
     "load"},
    
    {"load", FH_LOAD, 1, 1, FF_CONST, "Commands", "load filename",
     "Load session from file.",
     {"load session.sc", "load mywork", NULL, NULL, NULL, NULL},
     "save"},
    
    {"plot", FH_PLOT, 1, 3, FF_CONST, "Commands", "plot f(x) [xmin xmax]",
     "Plot a function in the terminal using ASCII art.",
     {"plot sin(x)", "plot sin(x) 0 2*pi", "plot x^2 -2 2", NULL, NULL, NULL},
     "hplot"},
    
    {"hplot", FH_PLOT, 1, 3, FF_CONST, "Commands", "hplot f(x) [xmin xmax]",
     "Horizontal plot (wider aspect ratio).",
     {"hplot sin(x)", "hplot cos(x) 0 2*pi", NULL, NULL, NULL, NULL},
     "plot"},
    
    {"solve", FH_SOLVE, 1, 1, FF_CONST, "Commands", "solve equation",
     "Solve quadratic or linear equations.",
     {"solve x^2-5x+6=0", "solve 2x+3=7", NULL, NULL, NULL, NULL},
     "root, fzero"},
    
    {"root", FH_SOLVE, 2, 2, FF_CONST, "Commands", "root f x0",
     "Find root of f(x)=0 using Newton-Raphson starting from x0.",
     {"f(x)=x^2-2; root f 1", "g(x)=cos(x)-x; root g 0.5", NULL, NULL, NULL, NULL},
     "solve, fzero, nsolve"},
    
    {"nsolve", FH_SOLVE, 2, 2, FF_CONST, "Commands", "nsolve f x0",
     "Numerical solver. Alias for root.",
     {"f(x)=x^3-x-1; nsolve f 1", NULL, NULL, NULL, NULL, NULL},
     "root, fzero"},
    
    {"fzero", FH_SOLVE, 2, 2, FF_CONST, "Commands", "fzero f x0",
     "Find zero of function (MATLAB compatible). Alias for root.",
     {"f(x)=x^2-2; fzero f 1", NULL, NULL, NULL, NULL, NULL},
     "root, solve"},
    
    {"const", FH_CONST_CMD, 0, 1, FF_CONST, "Commands", "const [name]",
     "List physical constants or show specific constant.",
     {"const", "const c", "const h", "const G", NULL, NULL},
     "help"},
    
    {"prec", FH_PREC, 0, 1, FF_CONST, "Commands", "prec [bits]",
     "Set or show arithmetic precision (64, 128, or 256 bits).",
     {"prec", "prec 256", "prec 64", NULL, NULL, NULL},
     "precision, digits"},
    
    {"base", FH_BASE, 1, 1, FF_CONST, "Commands", "base n",
     "Set output base (2-36).",
     {"base 16", "base 2", "base 10", NULL, NULL, NULL},
     "hex, bin, oct"},
    
    {"float", FH_FLOAT, 0, 0, FF_CONST, "Commands", "float",
     "Reset to standard floating-point display format.",
     {"float", NULL, NULL, NULL, NULL, NULL},
     "format, sci, eng"},
    
    {"eng1000", FH_ENG1000, 0, 0, FF_CONST, "Commands", "eng1000",
     "Engineering notation with SI prefixes (k, M, G, T, etc.).",
     {"eng1000", "1234567; eng1000", NULL, NULL, NULL, NULL},
     "eng, eng24"},
    
    {"eng24", FH_ENG24, 0, 0, FF_CONST, "Commands", "eng24",
     "Engineering notation optimized for 24-hour time/frequency.",
     {"eng24", NULL, NULL, NULL, NULL, NULL},
     "eng, eng1000"},
    
    {"ascii", FH_ASCII, 0, 0, FF_CONST, "Commands", "ascii",
     "Display ASCII character table.",
     {"ascii", NULL, NULL, NULL, NULL, NULL},
     "char"},
    
    {"set", FH_SET, 1, 2, FF_CONST, "Commands", "set option [value]",
     "Set calculator options (e.g., set complex on).",
     {"set complex on", "set angle deg", "set", NULL, NULL, NULL},
     "unset"},
    
    {"unset", FH_UNSET, 1, 1, FF_CONST, "Commands", "unset option",
     "Unset/reset calculator option to default.",
     {"unset complex", "unset angle", NULL, NULL, NULL, NULL},
     "set"},
    
    {"quad", FH_QUAD, 3, 3, FF_CONST, "Commands", "quad a b c",
     "Solve quadratic equation ax^2 + bx + c = 0.",
     {"quad 1 -5 6", "quad 1 0 -4", "quad 1 2 1", NULL, NULL, NULL},
     "solve, roots"},
    
    {"integrate", FH_INTEGRATE, 3, 3, FF_CONST, "Commands", "integrate f a b",
     "Numerical integration of f(x) from a to b.",
     {"f(x)=x^2; integrate f 0 1", "g(x)=sin(x); integrate g 0 pi", NULL, NULL, NULL, NULL},
     "trapz, diff"},
    
    {"orbit", FH_ORBIT, 1, 1, FF_CONST, "Commands", "orbit name",
     "Show orbital elements for a celestial body or satellite.",
     {"orbit ISS", "orbit moon", NULL, NULL, NULL, NULL},
     "tle"},
    
    {"tle", FH_TLE, 1, 1, FF_CONST, "Commands", "tle name",
     "Load Two-Line Element set for satellite tracking.",
     {"tle ISS", NULL, NULL, NULL, NULL, NULL},
     "orbit"},
    
    {"rpn", FH_RPN, 0, 0, FF_CONST, "Commands", "rpn",
     "Enter Reverse Polish Notation (RPN) mode.",
     {"rpn", NULL, NULL, NULL, NULL, NULL},
     "help"},
    
    /* ===== TRIGONOMETRIC ALIASES ===== */
    {"arcsin", FH_ASIN, 1, 1, FF_SCALAR | FF_COMPLEX, "Trigonometry", "arcsin(x)",
     "Arc sine. Alias for asin.",
     {"arcsin(0.5) = 0.5236", NULL, NULL, NULL, NULL, NULL},
     "asin"},
    
    {"arccos", FH_ACOS, 1, 1, FF_SCALAR | FF_COMPLEX, "Trigonometry", "arccos(x)",
     "Arc cosine. Alias for acos.",
     {"arccos(0.5) = 1.0472", NULL, NULL, NULL, NULL, NULL},
     "acos"},
    
    {"arctan", FH_ATAN, 1, 1, FF_SCALAR | FF_COMPLEX, "Trigonometry", "arctan(x)",
     "Arc tangent. Alias for atan.",
     {"arctan(1) = 0.7854", NULL, NULL, NULL, NULL, NULL},
     "atan"},
    
    {"arcsinh", FH_ASINH, 1, 1, FF_SCALAR | FF_COMPLEX, "Trigonometry", "arcsinh(x)",
     "Inverse hyperbolic sine. Alias for asinh.",
     {"arcsinh(1) = 0.8814", NULL, NULL, NULL, NULL, NULL},
     "asinh"},
    
    {"arccosh", FH_ACOSH, 1, 1, FF_SCALAR | FF_COMPLEX, "Trigonometry", "arccosh(x)",
     "Inverse hyperbolic cosine. Alias for acosh.",
     {"arccosh(2) = 1.3170", NULL, NULL, NULL, NULL, NULL},
     "acosh"},
    
    {"arctanh", FH_ATANH, 1, 1, FF_SCALAR | FF_COMPLEX, "Trigonometry", "arctanh(x)",
     "Inverse hyperbolic tangent. Alias for atanh.",
     {"arctanh(0.5) = 0.5493", NULL, NULL, NULL, NULL, NULL},
     "atanh"},
    
    {"arcsec", FH_ASEC, 1, 1, FF_SCALAR | FF_COMPLEX, "Trigonometry", "arcsec(x)",
     "Arc secant. Alias for asec.",
     {"arcsec(2) = 1.0472", NULL, NULL, NULL, NULL, NULL},
     "asec"},
    
    {"arccsc", FH_ACSC, 1, 1, FF_SCALAR | FF_COMPLEX, "Trigonometry", "arccsc(x)",
     "Arc cosecant. Alias for acsc.",
     {"arccsc(2) = 0.5236", NULL, NULL, NULL, NULL, NULL},
     "acsc"},
    
    {"arccot", FH_ACOT, 1, 1, FF_SCALAR | FF_COMPLEX, "Trigonometry", "arccot(x)",
     "Arc cotangent. Alias for acot.",
     {"arccot(1) = 0.7854", NULL, NULL, NULL, NULL, NULL},
     "acot"},
    
    {"arccotan", FH_ACOT, 1, 1, FF_SCALAR | FF_COMPLEX, "Trigonometry", "arccotan(x)",
     "Arc cotangent. Alias for acot.",
     {"arccotan(1) = 0.7854", NULL, NULL, NULL, NULL, NULL},
     "acot"},
    
    {"arccosec", FH_ACSC, 1, 1, FF_SCALAR | FF_COMPLEX, "Trigonometry", "arccosec(x)",
     "Arc cosecant. Alias for acsc.",
     {"arccosec(2) = 0.5236", NULL, NULL, NULL, NULL, NULL},
     "acsc"},
    
    {"cosec", FH_CSC, 1, 1, FF_SCALAR | FF_ANGLE | FF_COMPLEX, "Trigonometry", "cosec(x)",
     "Cosecant. Alias for csc.",
     {"cosec(pi/2) = 1", NULL, NULL, NULL, NULL, NULL},
     "csc"},
    
    {"cotan", FH_COT, 1, 1, FF_SCALAR | FF_ANGLE | FF_COMPLEX, "Trigonometry", "cotan(x)",
     "Cotangent. Alias for cot.",
     {"cotan(pi/4) = 1", NULL, NULL, NULL, NULL, NULL},
     "cot"},
    
    {"cotd", FH_COTD, 1, 1, FF_SCALAR | FF_DEGREE, "Trigonometry", "cotd(x)",
     "Cotangent of x in degrees.",
     {"cotd(45) = 1", "cotd(30) = 1.7321", NULL, NULL, NULL, NULL},
     "cot, tand"},
    
    {"cscd", FH_CSCD, 1, 1, FF_SCALAR | FF_DEGREE, "Trigonometry", "cscd(x)",
     "Cosecant of x in degrees.",
     {"cscd(90) = 1", "cscd(30) = 2", NULL, NULL, NULL, NULL},
     "csc, sind"},
    
    {"secd", FH_SECD, 1, 1, FF_SCALAR | FF_DEGREE, "Trigonometry", "secd(x)",
     "Secant of x in degrees.",
     {"secd(0) = 1", "secd(60) = 2", NULL, NULL, NULL, NULL},
     "sec, cosd"},
    
    {"acotd", FH_ACOTD, 1, 1, FF_SCALAR | FF_DEGREE, "Trigonometry", "acotd(x)",
     "Arc cotangent returning degrees.",
     {"acotd(1) = 45", "acotd(0) = 90", NULL, NULL, NULL, NULL},
     "acot, atand"},
    
    {"acscd", FH_ACSCD, 1, 1, FF_SCALAR | FF_DEGREE, "Trigonometry", "acscd(x)",
     "Arc cosecant returning degrees.",
     {"acscd(1) = 90", "acscd(2) = 30", NULL, NULL, NULL, NULL},
     "acsc, asind"},
    
    {"asecd", FH_ASECD, 1, 1, FF_SCALAR | FF_DEGREE, "Trigonometry", "asecd(x)",
     "Arc secant returning degrees.",
     {"asecd(1) = 0", "asecd(2) = 60", NULL, NULL, NULL, NULL},
     "asec, acosd"},
    
    /* ===== MORE NUMBER THEORY ===== */
    {"catalan", FH_CATALAN, 1, 1, FF_SCALAR, "Number Theory", "catalan(n)",
     "nth Catalan number. C(n) = (2n)! / ((n+1)! * n!).",
     {"catalan(0) = 1", "catalan(1) = 1", "catalan(5) = 42", "catalan(10) = 16796", NULL, NULL},
     "fibonacci, bell"},
    
    {"bell", FH_BELL, 1, 1, FF_SCALAR, "Number Theory", "bell(n)",
     "nth Bell number. Number of partitions of a set.",
     {"bell(0) = 1", "bell(1) = 1", "bell(5) = 52", "bell(10) = 115975", NULL, NULL},
     "catalan, stirling2"},
    
    {"stirling2", FH_STIRLING2, 2, 2, FF_SCALAR, "Number Theory", "stirling2(n, k)",
     "Stirling number of second kind. Ways to partition n elements into k subsets.",
     {"stirling2(4, 2) = 7", "stirling2(5, 3) = 25", NULL, NULL, NULL, NULL},
     "bell, ncr"},
    
    {"derangements", FH_DERANGEMENTS, 1, 1, FF_SCALAR, "Number Theory", "derangements(n)",
     "Number of derangements (permutations with no fixed points).",
     {"derangements(0) = 1", "derangements(3) = 2", "derangements(5) = 44", NULL, NULL, NULL},
     "fact, subfactorial"},
    
    {"subfactorial", FH_DERANGEMENTS, 1, 1, FF_SCALAR, "Number Theory", "subfactorial(n)",
     "Subfactorial. Alias for derangements.",
     {"subfactorial(5) = 44", NULL, NULL, NULL, NULL, NULL},
     "derangements"},
    
    {"divisorsum", FH_DIVISORSUM, 1, 1, FF_SCALAR, "Number Theory", "divisorsum(n)",
     "Sum of divisors of n (sigma function).",
     {"divisorsum(12) = 28", "divisorsum(6) = 12", "divisorsum(28) = 56", NULL, NULL, NULL},
     "divisors, totient"},
    
    {"sigma", FH_DIVISORSUM, 1, 1, FF_SCALAR, "Number Theory", "sigma(n)",
     "Sum of divisors. Alias for divisorsum.",
     {"sigma(12) = 28", NULL, NULL, NULL, NULL, NULL},
     "divisorsum"},
    
    {"divisors", FH_DIVISORS, 1, 1, FF_MATRIX, "Number Theory", "divisors(n)",
     "List of divisors of n.",
     {"divisors(12)", "divisors(28)", "length(divisors(60))", NULL, NULL, NULL},
     "divisorsum, factor"},
    
    {"factor", FH_FACTOR, 1, 1, FF_MATRIX, "Number Theory", "factor(n)",
     "Prime factorization of n.",
     {"factor(12)", "factor(100)", "factor(2310)", NULL, NULL, NULL},
     "isprime, divisors"},
    
    {"eulerphi", FH_TOTIENT, 1, 1, FF_SCALAR, "Number Theory", "eulerphi(n)",
     "Euler's totient function. Alias for totient.",
     {"eulerphi(10) = 4", NULL, NULL, NULL, NULL, NULL},
     "totient"},
    
    {"mobius", FH_MOBIUS, 1, 1, FF_SCALAR, "Number Theory", "mobius(n)",
     "Mobius function. mu(n) = (-1)^k if n is product of k distinct primes, 0 otherwise.",
     {"mobius(1) = 1", "mobius(6) = 1", "mobius(4) = 0", "mobius(30) = -1", NULL, NULL},
     "totient, factor"},
    
    {"moebius", FH_MOBIUS, 1, 1, FF_SCALAR, "Number Theory", "moebius(n)",
     "Mobius function. Alias for mobius.",
     {"moebius(6) = 1", NULL, NULL, NULL, NULL, NULL},
     "mobius"},
    
    {"mu", FH_MOBIUS, 1, 1, FF_SCALAR, "Number Theory", "mu(n)",
     "Mobius function. Alias for mobius.",
     {"mu(6) = 1", NULL, NULL, NULL, NULL, NULL},
     "mobius"},
    
    {"omega", FH_OMEGA, 1, 1, FF_SCALAR, "Number Theory", "omega(n)",
     "Number of distinct prime factors.",
     {"omega(12) = 2", "omega(30) = 3", "omega(64) = 1", NULL, NULL, NULL},
     "bigomega, factor"},
    
    {"Omega", FH_BIGOMEGA, 1, 1, FF_SCALAR, "Number Theory", "Omega(n)",
     "Number of prime factors with multiplicity.",
     {"Omega(12) = 3", "Omega(8) = 3", "Omega(30) = 3", NULL, NULL, NULL},
     "omega, factor"},
    
    {"bigomega", FH_BIGOMEGA, 1, 1, FF_SCALAR, "Number Theory", "bigomega(n)",
     "Number of prime factors with multiplicity. Alias for Omega.",
     {"bigomega(12) = 3", NULL, NULL, NULL, NULL, NULL},
     "Omega"},
    
    {"ndigits", FH_NDIGITS, 1, 1, FF_SCALAR, "Number Theory", "ndigits(n)",
     "Number of digits in n.",
     {"ndigits(12345) = 5", "ndigits(1000) = 4", "ndigits(0) = 1", NULL, NULL, NULL},
     "digitsum"},
    
    {"numdigits", FH_NDIGITS, 1, 1, FF_SCALAR, "Number Theory", "numdigits(n)",
     "Number of digits. Alias for ndigits.",
     {"numdigits(12345) = 5", NULL, NULL, NULL, NULL, NULL},
     "ndigits"},
    
    {"digitsum", FH_DIGITSUM, 1, 1, FF_SCALAR, "Number Theory", "digitsum(n)",
     "Sum of digits.",
     {"digitsum(12345) = 15", "digitsum(999) = 27", NULL, NULL, NULL, NULL},
     "digitroot, ndigits"},
    
    {"digsum", FH_DIGITSUM, 1, 1, FF_SCALAR, "Number Theory", "digsum(n)",
     "Sum of digits. Alias for digitsum.",
     {"digsum(12345) = 15", NULL, NULL, NULL, NULL, NULL},
     "digitsum"},
    
    {"digitroot", FH_DIGITROOT, 1, 1, FF_SCALAR, "Number Theory", "digitroot(n)",
     "Digital root (repeated digit sum until single digit).",
     {"digitroot(12345) = 6", "digitroot(999) = 9", "digitroot(123456789) = 9", NULL, NULL, NULL},
     "digitsum"},
    
    {"digroot", FH_DIGITROOT, 1, 1, FF_SCALAR, "Number Theory", "digroot(n)",
     "Digital root. Alias for digitroot.",
     {"digroot(12345) = 6", NULL, NULL, NULL, NULL, NULL},
     "digitroot"},
    
    {"digitalroot", FH_DIGITROOT, 1, 1, FF_SCALAR, "Number Theory", "digitalroot(n)",
     "Digital root. Alias for digitroot.",
     {"digitalroot(12345) = 6", NULL, NULL, NULL, NULL, NULL},
     "digitroot"},
    
    {"isperfect", FH_ISPERFECT, 1, 1, FF_SCALAR, "Number Theory", "isperfect(n)",
     "Test if n is a perfect number (sum of proper divisors equals n).",
     {"isperfect(6) = true", "isperfect(28) = true", "isperfect(12) = false", NULL, NULL, NULL},
     "divisorsum, isabundant"},
    
    {"isabundant", FH_ISABUNDANT, 1, 1, FF_SCALAR, "Number Theory", "isabundant(n)",
     "Test if n is abundant (sum of proper divisors > n).",
     {"isabundant(12) = true", "isabundant(6) = false", NULL, NULL, NULL, NULL},
     "isperfect, isdeficient"},
    
    {"isdeficient", FH_ISDEFICIENT, 1, 1, FF_SCALAR, "Number Theory", "isdeficient(n)",
     "Test if n is deficient (sum of proper divisors < n).",
     {"isdeficient(8) = true", "isdeficient(12) = false", NULL, NULL, NULL, NULL},
     "isperfect, isabundant"},
    
    {"issquarefree", FH_ISSQUAREFREE, 1, 1, FF_SCALAR, "Number Theory", "issquarefree(n)",
     "Test if n has no repeated prime factors.",
     {"issquarefree(6) = true", "issquarefree(12) = false", "issquarefree(30) = true", NULL, NULL, NULL},
     "factor, mobius"},
    
    {"fib", FH_FIBONACCI, 1, 1, FF_SCALAR, "Number Theory", "fib(n)",
     "Fibonacci number. Alias for fibonacci.",
     {"fib(10) = 55", NULL, NULL, NULL, NULL, NULL},
     "fibonacci"},
    
    {"fact2", FH_FACTORIAL2, 1, 1, FF_SCALAR, "Number Theory", "fact2(n)",
     "Double factorial. n!! = n*(n-2)*(n-4)*...",
     {"fact2(5) = 15", "fact2(6) = 48", "fact2(7) = 105", NULL, NULL, NULL},
     "factorial2, fact"},
    
    {"factorial2", FH_FACTORIAL2, 1, 1, FF_SCALAR, "Number Theory", "factorial2(n)",
     "Double factorial. n!! = n*(n-2)*(n-4)*...",
     {"factorial2(5) = 15", "factorial2(6) = 48", NULL, NULL, NULL, NULL},
     "fact, fact2"},
    
    {"ffactorial", FH_FACTORIAL2, 1, 1, FF_SCALAR, "Number Theory", "ffactorial(n)",
     "Double factorial. Alias for factorial2.",
     {"ffactorial(5) = 15", NULL, NULL, NULL, NULL, NULL},
     "factorial2"},
    
    {"choose", FH_NCR, 2, 2, FF_SCALAR, "Number Theory", "choose(n, r)",
     "Binomial coefficient. Alias for ncr.",
     {"choose(5, 2) = 10", NULL, NULL, NULL, NULL, NULL},
     "ncr, comb"},
    
    {"pochhammer", FH_POCHHAMMER, 2, 2, FF_SCALAR, "Number Theory", "pochhammer(x, n)",
     "Pochhammer symbol (rising factorial). (x)_n = x*(x+1)*...*(x+n-1).",
     {"pochhammer(1, 5) = 120", "pochhammer(3, 4) = 360", NULL, NULL, NULL, NULL},
     "falling, fact"},
    
    {"rising", FH_POCHHAMMER, 2, 2, FF_SCALAR, "Number Theory", "rising(x, n)",
     "Rising factorial. Alias for pochhammer.",
     {"rising(1, 5) = 120", NULL, NULL, NULL, NULL, NULL},
     "pochhammer"},
    
    {"falling", FH_FALLING, 2, 2, FF_SCALAR, "Number Theory", "falling(x, n)",
     "Falling factorial. x*(x-1)*...*(x-n+1).",
     {"falling(5, 3) = 60", "falling(10, 4) = 5040", NULL, NULL, NULL, NULL},
     "pochhammer, npr"},
    
    {"primepi", FH_PRIMEPI, 1, 1, FF_SCALAR, "Number Theory", "primepi(n)",
     "Prime counting function. Number of primes <= n.",
     {"primepi(10) = 4", "primepi(100) = 25", "primepi(1000) = 168", NULL, NULL, NULL},
     "isprime, nthprime"},
    
    {"nthprime", FH_NTHPRIME, 1, 1, FF_SCALAR, "Number Theory", "nthprime(n)",
     "The nth prime number.",
     {"nthprime(1) = 2", "nthprime(10) = 29", "nthprime(100) = 541", NULL, NULL, NULL},
     "primepi, isprime"},
    
    {"partition", FH_PARTITION, 1, 1, FF_SCALAR, "Number Theory", "partition(n)",
     "Number of integer partitions of n.",
     {"partition(5) = 7", "partition(10) = 42", "partition(50) = 204226", NULL, NULL, NULL},
     "partitions"},
    
    {"partitions", FH_PARTITION, 1, 1, FF_SCALAR, "Number Theory", "partitions(n)",
     "Number of integer partitions. Alias for partition.",
     {"partitions(5) = 7", NULL, NULL, NULL, NULL, NULL},
     "partition"},
    
    /* ===== MORE SPECIAL FUNCTIONS ===== */
    {"digamma", FH_DIGAMMA, 1, 1, FF_SCALAR, "Special Functions", "digamma(x)",
     "Digamma function. Logarithmic derivative of gamma: psi(x) = gamma'(x)/gamma(x).",
     {"digamma(1) = -0.5772", "digamma(2) = 0.4228", "digamma(0.5) = -1.9635", NULL, NULL, NULL},
     "gamma, psi"},
    
    {"psi", FH_DIGAMMA, 1, 1, FF_SCALAR, "Special Functions", "psi(x)",
     "Digamma function. Alias for digamma.",
     {"psi(1) = -0.5772", NULL, NULL, NULL, NULL, NULL},
     "digamma"},
    
    {"harmonic", FH_HARMONIC, 1, 1, FF_SCALAR, "Special Functions", "harmonic(n)",
     "nth harmonic number. H(n) = 1 + 1/2 + 1/3 + ... + 1/n.",
     {"harmonic(1) = 1", "harmonic(10) = 2.9290", "harmonic(100) = 5.1874", NULL, NULL, NULL},
     "zeta"},
    
    {"genharmonic", FH_GENHARMONIC, 2, 2, FF_SCALAR, "Special Functions", "genharmonic(n, m)",
     "Generalized harmonic number. H(n,m) = sum(1/k^m, k=1..n).",
     {"genharmonic(10, 2) = 1.5498", "genharmonic(100, 2) = 1.6350", NULL, NULL, NULL, NULL},
     "harmonic, zeta"},
    
    {"harmonic2", FH_GENHARMONIC, 2, 2, FF_SCALAR, "Special Functions", "harmonic2(n, m)",
     "Generalized harmonic. Alias for genharmonic.",
     {"harmonic2(10, 2) = 1.5498", NULL, NULL, NULL, NULL, NULL},
     "genharmonic"},
    
    {"tgamma", FH_GAMMA, 1, 1, FF_SCALAR, "Special Functions", "tgamma(x)",
     "Gamma function. Alias for gamma.",
     {"tgamma(5) = 24", NULL, NULL, NULL, NULL, NULL},
     "gamma"},
    
    {"softplus", FH_SOFTPLUS, 1, 1, FF_SCALAR, "Special Functions", "softplus(x)",
     "Softplus function. softplus(x) = log(1 + exp(x)).",
     {"softplus(0) = 0.6931", "softplus(1) = 1.3133", "softplus(-10) = 0.00005", NULL, NULL, NULL},
     "sigmoid"},
    
    {"step", FH_HEAVISIDE, 1, 1, FF_SCALAR, "Special Functions", "step(x)",
     "Step function. Alias for heaviside.",
     {"step(-1) = 0", "step(0) = 0.5", "step(1) = 1", NULL, NULL, NULL},
     "heaviside"},
    
    {"rect", FH_RECT, 1, 1, FF_SCALAR, "Special Functions", "rect(x)",
     "Rectangle function. 1 for |x| < 0.5, 0.5 at |x| = 0.5, 0 otherwise.",
     {"rect(0) = 1", "rect(0.3) = 1", "rect(0.5) = 0.5", "rect(1) = 0", NULL, NULL},
     "tri, heaviside"},
    
    {"rectangle", FH_RECT, 1, 1, FF_SCALAR, "Special Functions", "rectangle(x)",
     "Rectangle function. Alias for rect.",
     {"rectangle(0) = 1", NULL, NULL, NULL, NULL, NULL},
     "rect"},
    
    {"tri", FH_TRI, 1, 1, FF_SCALAR, "Special Functions", "tri(x)",
     "Triangle function. 1-|x| for |x| < 1, 0 otherwise.",
     {"tri(0) = 1", "tri(0.5) = 0.5", "tri(1) = 0", "tri(-0.5) = 0.5", NULL, NULL},
     "rect"},
    
    {"triangle", FH_TRI, 1, 1, FF_SCALAR, "Special Functions", "triangle(x)",
     "Triangle function. Alias for tri.",
     {"triangle(0) = 1", NULL, NULL, NULL, NULL, NULL},
     "tri"},
    
    {"logistic", FH_SIGMOID, 1, 1, FF_SCALAR, "Special Functions", "logistic(x)",
     "Logistic function. Alias for sigmoid.",
     {"logistic(0) = 0.5", NULL, NULL, NULL, NULL, NULL},
     "sigmoid"},
    
    {"erfinv", FH_ERFINV, 1, 1, FF_SCALAR, "Special Functions", "erfinv(x)",
     "Inverse error function.",
     {"erfinv(0) = 0", "erfinv(0.5) = 0.4769", "erfinv(0.9) = 1.1631", NULL, NULL, NULL},
     "erf, erfc"},
    
    {"besselJ", FH_BESSELJ, 2, 2, FF_SCALAR, "Special Functions", "besselJ(n, x)",
     "Bessel function of first kind. Alias for besselj.",
     {"besselJ(0, 1) = 0.7652", NULL, NULL, NULL, NULL, NULL},
     "besselj"},
    
    {"besselY", FH_BESSELY, 2, 2, FF_SCALAR, "Special Functions", "besselY(n, x)",
     "Bessel function of second kind. Alias for bessely.",
     {"besselY(0, 1) = 0.0883", NULL, NULL, NULL, NULL, NULL},
     "bessely"},
    
    {"riemann_zeta", FH_ZETA, 1, 1, FF_SCALAR, "Special Functions", "riemann_zeta(s)",
     "Riemann zeta function. Alias for zeta.",
     {"riemann_zeta(2) = 1.6449", NULL, NULL, NULL, NULL, NULL},
     "zeta"},
    
    /* ===== MORE PROBABILITY ===== */
    {"chi2cdf", FH_CHI2CDF, 2, 2, FF_SCALAR, "Probability", "chi2cdf(x, df)",
     "Chi-squared cumulative distribution function.",
     {"chi2cdf(3.84, 1) = 0.95", "chi2cdf(5.99, 2) = 0.95", NULL, NULL, NULL, NULL},
     "chi2inv, normcdf"},
    
    {"chicdf", FH_CHI2CDF, 2, 2, FF_SCALAR, "Probability", "chicdf(x, df)",
     "Chi-squared CDF. Alias for chi2cdf.",
     {"chicdf(3.84, 1) = 0.95", NULL, NULL, NULL, NULL, NULL},
     "chi2cdf"},
    
    {"tcdf", FH_TCDF, 2, 2, FF_SCALAR, "Probability", "tcdf(x, df)",
     "Student's t cumulative distribution function.",
     {"tcdf(0, 10) = 0.5", "tcdf(2.228, 10) = 0.975", NULL, NULL, NULL, NULL},
     "tinv, normcdf"},
    
    {"student_cdf", FH_TCDF, 2, 2, FF_SCALAR, "Probability", "student_cdf(x, df)",
     "Student's t CDF. Alias for tcdf.",
     {"student_cdf(0, 10) = 0.5", NULL, NULL, NULL, NULL, NULL},
     "tcdf"},
    
    {"normalpdf", FH_NORMPDF, 1, 3, FF_SCALAR, "Probability", "normalpdf(x)",
     "Normal PDF. Alias for normpdf.",
     {"normalpdf(0) = 0.3989", NULL, NULL, NULL, NULL, NULL},
     "normpdf"},
    
    {"normalcdf", FH_NORMCDF, 1, 3, FF_SCALAR, "Probability", "normalcdf(x)",
     "Normal CDF. Alias for normcdf.",
     {"normalcdf(0) = 0.5", NULL, NULL, NULL, NULL, NULL},
     "normcdf"},
    
    {"normalinv", FH_NORMINV, 1, 3, FF_SCALAR, "Probability", "normalinv(p)",
     "Inverse normal CDF. Alias for norminv.",
     {"normalinv(0.5) = 0", NULL, NULL, NULL, NULL, NULL},
     "norminv"},
    
    {"probit", FH_NORMINV, 1, 1, FF_SCALAR, "Probability", "probit(p)",
     "Probit function (inverse normal CDF). Alias for norminv.",
     {"probit(0.975) = 1.96", NULL, NULL, NULL, NULL, NULL},
     "norminv"},
    
    {"randperm", FH_RANDPERM, 1, 2, FF_MATRIX, "Probability", "randperm(n) or randperm(n, k)",
     "Random permutation.",
     {"randperm(5)", "randperm(10, 3)", NULL, NULL, NULL, NULL},
     "rand, randi"},
    
    /* ===== MORE STATISTICS ===== */
    {"geomean", FH_GEOMEAN, 1, -1, FF_MATRIX | FF_VARIADIC, "Statistics", "geomean(v)",
     "Geometric mean.",
     {"geomean([1,2,4,8]) = 2.8284", "geomean([4,9]) = 6", NULL, NULL, NULL, NULL},
     "mean, harmmean"},
    
    {"harmmean", FH_HARMMEAN, 1, -1, FF_MATRIX | FF_VARIADIC, "Statistics", "harmmean(v)",
     "Harmonic mean.",
     {"harmmean([1,2,4]) = 1.7143", "harmmean([2,3,6]) = 3", NULL, NULL, NULL, NULL},
     "mean, geomean"},
    
    {"skewness", FH_SKEWNESS, 1, -1, FF_MATRIX | FF_VARIADIC, "Statistics", "skewness(v)",
     "Sample skewness (asymmetry measure).",
     {"skewness([1,2,3,4,5]) = 0", "skewness([1,2,2,3,10]) = 1.4626", NULL, NULL, NULL, NULL},
     "kurtosis, std"},
    
    {"kurtosis", FH_KURTOSIS, 1, -1, FF_MATRIX | FF_VARIADIC, "Statistics", "kurtosis(v)",
     "Sample kurtosis (tail heaviness measure).",
     {"kurtosis([1,2,3,4,5])", "kurtosis([1,1,1,1,10])", NULL, NULL, NULL, NULL},
     "skewness, std"},
    
    {"mad", FH_MAD, 1, -1, FF_MATRIX | FF_VARIADIC, "Statistics", "mad(v)",
     "Mean absolute deviation.",
     {"mad([1,2,3,4,5]) = 1.2", "mad([1,1,1,1,10]) = 2.88", NULL, NULL, NULL, NULL},
     "std, iqr"},
    
    {"iqr", FH_IQR, 1, -1, FF_MATRIX | FF_VARIADIC, "Statistics", "iqr(v)",
     "Interquartile range (75th percentile - 25th percentile).",
     {"iqr([1,2,3,4,5,6,7,8,9,10]) = 5", NULL, NULL, NULL, NULL, NULL},
     "prctile, mad"},
    
    {"rms", FH_RMS, 1, -1, FF_MATRIX | FF_VARIADIC, "Statistics", "rms(v)",
     "Root mean square.",
     {"rms([3,4]) = 3.5355", "rms([1,2,3,4,5]) = 3.3166", NULL, NULL, NULL, NULL},
     "mean, sumsq"},
    
    {"sumsq", FH_SUMSQ, 1, -1, FF_MATRIX | FF_VARIADIC, "Statistics", "sumsq(v)",
     "Sum of squares.",
     {"sumsq([1,2,3]) = 14", "sumsq([3,4]) = 25", NULL, NULL, NULL, NULL},
     "sum, meansq"},
    
    {"meansq", FH_MEANSQ, 1, -1, FF_MATRIX | FF_VARIADIC, "Statistics", "meansq(v)",
     "Mean of squares.",
     {"meansq([1,2,3]) = 4.6667", "meansq([3,4]) = 12.5", NULL, NULL, NULL, NULL},
     "sumsq, rms"},
    
    {"prctile", FH_PRCTILE, 2, 2, FF_MATRIX, "Statistics", "prctile(v, p)",
     "Percentile. p is percentage (0-100).",
     {"prctile([1,2,3,4,5], 50) = 3", "prctile([1,2,3,4,5], 25) = 1.5", NULL, NULL, NULL, NULL},
     "iqr, median"},
    
    {"pct", FH_PRCTILE, 2, 2, FF_MATRIX, "Statistics", "pct(v, p)",
     "Percentile. Alias for prctile.",
     {"pct([1,2,3,4,5], 50) = 3", NULL, NULL, NULL, NULL, NULL},
     "prctile"},
    
    {"zscore", FH_ZSCORE, 1, 1, FF_MATRIX, "Statistics", "zscore(v)",
     "Standardized scores (subtract mean, divide by std).",
     {"zscore([1,2,3,4,5])", "mean(zscore([1,2,3,4,5])) = 0", NULL, NULL, NULL, NULL},
     "std, mean"},
    
    {"cov", FH_COV, 2, 2, FF_MATRIX, "Statistics", "cov(x, y)",
     "Covariance between x and y.",
     {"cov([1,2,3], [4,5,6]) = 1", "cov([1,2,3], [3,2,1]) = -1", NULL, NULL, NULL, NULL},
     "corrcoef, var"},
    
    {"corrcoef", FH_CORRCOEF, 2, 2, FF_MATRIX, "Statistics", "corrcoef(x, y)",
     "Correlation coefficient between x and y.",
     {"corrcoef([1,2,3], [4,5,6]) = 1", "corrcoef([1,2,3], [3,2,1]) = -1", NULL, NULL, NULL, NULL},
     "cov"},
    
    {"sd", FH_STD, 1, -1, FF_MATRIX | FF_VARIADIC, "Statistics", "sd(v)",
     "Standard deviation. Alias for std.",
     {"sd([1,2,3,4,5]) = 1.5811", NULL, NULL, NULL, NULL, NULL},
     "std"},
    
    /* ===== MORE LINEAR ALGEBRA ===== */
    {"tr", FH_TRACE, 1, 1, FF_MATRIX, "Linear Algebra", "tr(A)",
     "Trace. Alias for trace.",
     {"tr(eye(5)) = 5", NULL, NULL, NULL, NULL, NULL},
     "trace"},
    
    {"cond", FH_COND, 1, 2, FF_MATRIX, "Linear Algebra", "cond(A) or cond(A, p)",
     "Condition number.",
     {"cond(eye(3)) = 1", "cond([1,2;3,4])", NULL, NULL, NULL, NULL},
     "svd, norm"},
    
    {"inverse", FH_INV, 1, 1, FF_MATRIX, "Linear Algebra", "inverse(A)",
     "Matrix inverse. Alias for inv.",
     {"inverse([1,2;3,4])", NULL, NULL, NULL, NULL, NULL},
     "inv"},
    
    {"null", FH_NULL, 1, 1, FF_MATRIX, "Linear Algebra", "null(A)",
     "Null space basis.",
     {"null([1,2;2,4])", "null([1,0,0;0,1,0])", NULL, NULL, NULL, NULL},
     "rank"},
    
    {"nullspace", FH_NULL, 1, 1, FF_MATRIX, "Linear Algebra", "nullspace(A)",
     "Null space. Alias for null.",
     {"nullspace([1,2;2,4])", NULL, NULL, NULL, NULL, NULL},
     "null"},
    
    {"cholesky", FH_CHOL, 1, 1, FF_MATRIX, "Linear Algebra", "cholesky(A)",
     "Cholesky decomposition. Alias for chol.",
     {"cholesky([4,2;2,2])", NULL, NULL, NULL, NULL, NULL},
     "chol"},
    
    {"eigenvalues", FH_EIG, 1, 1, FF_MATRIX, "Linear Algebra", "eigenvalues(A)",
     "Eigenvalues. Alias for eig.",
     {"eigenvalues([1,2;2,1])", NULL, NULL, NULL, NULL, NULL},
     "eig"},
    
    {"schur", FH_SCHUR, 1, 1, FF_MATRIX, "Linear Algebra", "schur(A)",
     "Schur decomposition.",
     {"schur([1,2;0,3])", "schur(eye(3))", NULL, NULL, NULL, NULL},
     "eig, qr"},
    
    {"mldivide", FH_MLDIVIDE, 2, 2, FF_MATRIX, "Linear Algebra", "mldivide(A, b)",
     "Matrix left division (solve Ax=b).",
     {"mldivide([1,0;0,1], [3;4])", NULL, NULL, NULL, NULL, NULL},
     "linsolve, ldivide"},
    
    {"ldivide", FH_MLDIVIDE, 2, 2, FF_MATRIX, "Linear Algebra", "ldivide(A, b)",
     "Left divide. Alias for mldivide.",
     {"ldivide([1,0;0,1], [3;4])", NULL, NULL, NULL, NULL, NULL},
     "mldivide"},
    
    /* ===== MORE MATRIX FUNCTIONS ===== */
    {"identity", FH_EYE, 1, 2, FF_MATRIX, "Matrix", "identity(n)",
     "Identity matrix. Alias for eye.",
     {"identity(3)", NULL, NULL, NULL, NULL, NULL},
     "eye"},
    
    {"magic", FH_MAGIC, 1, 1, FF_MATRIX, "Matrix", "magic(n)",
     "Magic square matrix (rows, cols, diagonals sum to same value).",
     {"magic(3)", "sum(magic(3))", "sum(magic(5))", NULL, NULL, NULL},
     "pascal, hilb"},
    
    {"pascal", FH_PASCAL, 1, 1, FF_MATRIX, "Matrix", "pascal(n)",
     "Pascal matrix (binomial coefficients).",
     {"pascal(4)", "pascal(5)", NULL, NULL, NULL, NULL},
     "magic, hilb"},
    
    {"hilb", FH_HILB, 1, 1, FF_MATRIX, "Matrix", "hilb(n)",
     "Hilbert matrix. H(i,j) = 1/(i+j-1).",
     {"hilb(3)", "hilb(4)", "cond(hilb(5))", NULL, NULL, NULL},
     "invhilb, pascal"},
    
    {"invhilb", FH_INVHILB, 1, 1, FF_MATRIX, "Matrix", "invhilb(n)",
     "Inverse Hilbert matrix.",
     {"invhilb(3)", "hilb(3) * invhilb(3)", NULL, NULL, NULL, NULL},
     "hilb"},
    
    {"vander", FH_VANDER, 1, 2, FF_MATRIX, "Matrix", "vander(v) or vander(v, n)",
     "Vandermonde matrix.",
     {"vander([1,2,3])", "vander([1,2,3,4])", NULL, NULL, NULL, NULL},
     "polyval"},
    
    {"toeplitz", FH_TOEPLITZ, 1, 2, FF_MATRIX, "Matrix", "toeplitz(c) or toeplitz(c, r)",
     "Toeplitz matrix (constant diagonals).",
     {"toeplitz([1,2,3])", "toeplitz([1,2,3], [1,4,5])", NULL, NULL, NULL, NULL},
     "hankel"},
    
    {"hankel", FH_HANKEL, 1, 2, FF_MATRIX, "Matrix", "hankel(c) or hankel(c, r)",
     "Hankel matrix (constant anti-diagonals).",
     {"hankel([1,2,3])", "hankel([1,2,3], [3,4,5])", NULL, NULL, NULL, NULL},
     "toeplitz"},
    
    {"compan", FH_COMPAN, 1, 1, FF_MATRIX, "Matrix", "compan(p)",
     "Companion matrix for polynomial coefficients.",
     {"compan([1,0,-1])", "eig(compan([1,-6,11,-6]))", NULL, NULL, NULL, NULL},
     "roots, polyval"},
    
    {"blkdiag", FH_BLKDIAG, 1, -1, FF_MATRIX | FF_VARIADIC, "Matrix", "blkdiag(A, B, ...)",
     "Block diagonal matrix.",
     {"blkdiag([1,2;3,4], [5,6;7,8])", "blkdiag(eye(2), eye(3))", NULL, NULL, NULL, NULL},
     "diag"},
    
    {"cat", FH_CAT, 2, -1, FF_MATRIX | FF_VARIADIC, "Matrix", "cat(dim, A, B, ...)",
     "Concatenate arrays along dimension dim.",
     {"cat(1, [1,2], [3,4])", "cat(2, [1;2], [3;4])", NULL, NULL, NULL, NULL},
     "horzcat, vertcat"},
    
    {"horzcat", FH_HORZCAT, 2, -1, FF_MATRIX | FF_VARIADIC, "Matrix", "horzcat(A, B, ...)",
     "Horizontal concatenation.",
     {"horzcat([1,2], [3,4])", "horzcat([1;2], [3;4])", NULL, NULL, NULL, NULL},
     "vertcat, cat"},
    
    {"vertcat", FH_VERTCAT, 2, -1, FF_MATRIX | FF_VARIADIC, "Matrix", "vertcat(A, B, ...)",
     "Vertical concatenation.",
     {"vertcat([1,2], [3,4])", "vertcat([1;2], [3;4])", NULL, NULL, NULL, NULL},
     "horzcat, cat"},
    
    {"circshift", FH_CIRCSHIFT, 2, 3, FF_MATRIX, "Matrix", "circshift(A, k) or circshift(A, k, dim)",
     "Circular shift of array elements.",
     {"circshift([1,2,3,4], 1)", "circshift([1,2,3,4], -1)", NULL, NULL, NULL, NULL},
     "flip, rot90"},
    
    {"repmat", FH_REPMAT, 2, 3, FF_MATRIX, "Matrix", "repmat(A, m, n)",
     "Replicate and tile array.",
     {"repmat([1,2], 2, 3)", "repmat([1;2], 3, 2)", NULL, NULL, NULL, NULL},
     "kron"},
    
    {"triu", FH_TRIU, 1, 2, FF_MATRIX, "Matrix", "triu(A) or triu(A, k)",
     "Upper triangular part.",
     {"triu([1,2,3;4,5,6;7,8,9])", "triu([1,2,3;4,5,6;7,8,9], 1)", NULL, NULL, NULL, NULL},
     "tril"},
    
    {"tril", FH_TRIL, 1, 2, FF_MATRIX, "Matrix", "tril(A) or tril(A, k)",
     "Lower triangular part.",
     {"tril([1,2,3;4,5,6;7,8,9])", "tril([1,2,3;4,5,6;7,8,9], -1)", NULL, NULL, NULL, NULL},
     "triu"},
    
    {"sortrows", FH_SORTROWS, 1, 2, FF_MATRIX, "Matrix", "sortrows(A) or sortrows(A, col)",
     "Sort rows of matrix.",
     {"sortrows([3,1;1,2;2,3])", "sortrows([3,1;1,2;2,3], 2)", NULL, NULL, NULL, NULL},
     "sort"},
    
    {"find", FH_FIND, 1, 2, FF_MATRIX, "Matrix", "find(A) or find(A, n)",
     "Find indices of nonzero elements.",
     {"find([0,1,0,2,0])", "find([1,0;0,1])", NULL, NULL, NULL, NULL},
     "nnz"},
    
    /* ===== MORE MATRIX QUERY ===== */
    {"ndims", FH_NDIMS, 1, 1, FF_MATRIX, "Matrix Query", "ndims(A)",
     "Number of dimensions (always 2 for matrices).",
     {"ndims([1,2,3]) = 2", "ndims([1;2;3]) = 2", NULL, NULL, NULL, NULL},
     "size"},
    
    {"isrow", FH_ISROW, 1, 1, FF_MATRIX, "Matrix Query", "isrow(A)",
     "Test if row vector.",
     {"isrow([1,2,3]) = true", "isrow([1;2;3]) = false", NULL, NULL, NULL, NULL},
     "iscolumn, isvector"},
    
    {"iscolumn", FH_ISCOLUMN, 1, 1, FF_MATRIX, "Matrix Query", "iscolumn(A)",
     "Test if column vector.",
     {"iscolumn([1;2;3]) = true", "iscolumn([1,2,3]) = false", NULL, NULL, NULL, NULL},
     "isrow, isvector"},
    
    {"ismatrix", FH_ISMATRIX, 1, 1, FF_MATRIX, "Matrix Query", "ismatrix(A)",
     "Test if 2D matrix.",
     {"ismatrix([1,2;3,4]) = true", "ismatrix([1,2,3]) = true", NULL, NULL, NULL, NULL},
     "isvector, isscalar"},
    
    {"issorted", FH_ISSORTED, 1, 1, FF_MATRIX, "Matrix Query", "issorted(v)",
     "Test if vector is sorted.",
     {"issorted([1,2,3,4]) = true", "issorted([1,3,2,4]) = false", NULL, NULL, NULL, NULL},
     "sort"},
    
    {"issymmetric", FH_ISSYMMETRIC, 1, 1, FF_MATRIX, "Matrix Query", "issymmetric(A)",
     "Test if matrix is symmetric.",
     {"issymmetric([1,2;2,1]) = true", "issymmetric([1,2;3,4]) = false", NULL, NULL, NULL, NULL},
     "issquare"},
    
    /* ===== MORE SIGNAL PROCESSING ===== */
    {"deconv", FH_DECONV, 2, 2, FF_MATRIX, "Signal Processing", "deconv(u, v)",
     "Deconvolution (polynomial division).",
     {"deconv([1,3,3,1], [1,1])", "deconv(conv([1,2],[1,1]), [1,1])", NULL, NULL, NULL, NULL},
     "conv"},
    
    {"gradient", FH_GRADIENT, 1, 2, FF_MATRIX, "Signal Processing", "gradient(v) or gradient(v, h)",
     "Numerical gradient.",
     {"gradient([1,4,9,16])", "gradient([1,4,9,16], 0.5)", NULL, NULL, NULL, NULL},
     "diff"},
    
    {"cummin", FH_CUMMIN, 1, 1, FF_MATRIX, "Signal Processing", "cummin(v)",
     "Cumulative minimum.",
     {"cummin([3,1,4,1,5])", "cummin([5,4,3,2,1])", NULL, NULL, NULL, NULL},
     "cummax, min"},
    
    {"cummax", FH_CUMMAX, 1, 1, FF_MATRIX, "Signal Processing", "cummax(v)",
     "Cumulative maximum.",
     {"cummax([1,3,2,5,4])", "cummax([1,2,3,4,5])", NULL, NULL, NULL, NULL},
     "cummin, max"},
    
    {"movmean", FH_MOVMEAN, 2, 2, FF_MATRIX, "Signal Processing", "movmean(v, k)",
     "Moving average with window size k.",
     {"movmean([1,2,3,4,5], 3)", "movmean([1,2,10,2,1], 3)", NULL, NULL, NULL, NULL},
     "movsum, mean"},
    
    {"movsum", FH_MOVSUM, 2, 2, FF_MATRIX, "Signal Processing", "movsum(v, k)",
     "Moving sum with window size k.",
     {"movsum([1,2,3,4,5], 3)", "movsum([1,1,1,1,1], 2)", NULL, NULL, NULL, NULL},
     "movmean, movmax, sum"},
    
    {"movmax", FH_MOVMAX, 2, 2, FF_MATRIX, "Signal Processing", "movmax(v, k)",
     "Moving maximum with window size k.",
     {"movmax([1,3,2,5,4], 3)", "movmax([5,4,3,2,1], 2)", NULL, NULL, NULL, NULL},
     "movmin, movmean, max"},
    
    {"movmin", FH_MOVMIN, 2, 2, FF_MATRIX, "Signal Processing", "movmin(v, k)",
     "Moving minimum with window size k.",
     {"movmin([1,3,2,5,4], 3)", "movmin([5,4,3,2,1], 2)", NULL, NULL, NULL, NULL},
     "movmax, movmean, min"},
    
    {"movstd", FH_MOVSTD, 2, 2, FF_MATRIX, "Signal Processing", "movstd(v, k)",
     "Moving standard deviation with window size k.",
     {"movstd([1,2,3,4,5], 3)", "movstd(randn(1,20), 5)", NULL, NULL, NULL, NULL},
     "movmean, std"},
    
    {"interp1", FH_INTERP1, 3, 4, FF_MATRIX, "Signal Processing", "interp1(x, y, xq)",
     "1-D interpolation.",
     {"interp1([1,2,3], [1,4,9], 1.5)", "interp1([0,1,2], [0,1,4], [0.5,1.5])", NULL, NULL, NULL, NULL},
     "lerp"},
    
    /* ===== POLYNOMIALS ===== */
    {"polyval", FH_POLYVAL, 2, 2, FF_MATRIX, "Polynomials", "polyval(p, x)",
     "Evaluate polynomial at x. p = [a_n, ..., a_1, a_0].",
     {"polyval([1,0,-1], 2) = 3", "polyval([1,2,1], 3) = 16", NULL, NULL, NULL, NULL},
     "polyder, roots"},
    
    {"polyder", FH_POLYDER, 1, 1, FF_MATRIX, "Polynomials", "polyder(p)",
     "Polynomial derivative.",
     {"polyder([1,0,-1])", "polyder([3,2,1])", NULL, NULL, NULL, NULL},
     "polyval, polyint"},
    
    {"polyint", FH_POLYINT, 1, 2, FF_MATRIX, "Polynomials", "polyint(p) or polyint(p, k)",
     "Polynomial integration.",
     {"polyint([3,2,1])", "polyint([1,0], 5)", NULL, NULL, NULL, NULL},
     "polyder, polyval"},
    
    {"roots2", FH_ROOTS2, 3, 3, FF_MATRIX, "Polynomials", "roots2(a, b, c)",
     "Roots of quadratic ax^2 + bx + c = 0.",
     {"roots2(1, 0, -1)", "roots2(1, -5, 6)", "roots2(1, 0, 1)", NULL, NULL, NULL},
     "quadroots"},
    
    {"quadroots", FH_ROOTS2, 3, 3, FF_MATRIX, "Polynomials", "quadroots(a, b, c)",
     "Quadratic roots. Alias for roots2.",
     {"quadroots(1, 0, -1)", NULL, NULL, NULL, NULL, NULL},
     "roots2"},
    
    {"quadratic", FH_ROOTS2, 3, 3, FF_MATRIX, "Polynomials", "quadratic(a, b, c)",
     "Quadratic roots. Alias for roots2.",
     {"quadratic(1, 0, -1)", NULL, NULL, NULL, NULL, NULL},
     "roots2"},
    
    /* ===== ANGLE/TIME CONVERSION ===== */
    {"wrapToPi", FH_WRAPTOPI, 1, 1, FF_SCALAR, "Angle Conversion", "wrapToPi(x)",
     "Wrap angle to [-pi, pi].",
     {"wrapToPi(4) = -2.2832", "wrapToPi(-4) = 2.2832", NULL, NULL, NULL, NULL},
     "wrapTo2Pi, wrap180"},
    
    {"wrap", FH_WRAPTOPI, 1, 1, FF_SCALAR, "Angle Conversion", "wrap(x)",
     "Wrap angle to [-pi, pi]. Alias for wrapToPi.",
     {"wrap(4) = -2.2832", NULL, NULL, NULL, NULL, NULL},
     "wrapToPi"},
    
    {"wrapTo180", FH_WRAPTO180, 1, 1, FF_SCALAR, "Angle Conversion", "wrapTo180(x)",
     "Wrap angle to [-180, 180] degrees.",
     {"wrapTo180(270) = -90", "wrapTo180(-270) = 90", NULL, NULL, NULL, NULL},
     "wrapTo360, wrapToPi"},
    
    {"wrap180", FH_WRAPTO180, 1, 1, FF_SCALAR, "Angle Conversion", "wrap180(x)",
     "Wrap angle to [-180, 180]. Alias for wrapTo180.",
     {"wrap180(270) = -90", NULL, NULL, NULL, NULL, NULL},
     "wrapTo180"},
    
    {"wrapTo360", FH_WRAPTO360, 1, 1, FF_SCALAR, "Angle Conversion", "wrapTo360(x)",
     "Wrap angle to [0, 360] degrees.",
     {"wrapTo360(-90) = 270", "wrapTo360(450) = 90", NULL, NULL, NULL, NULL},
     "wrapTo180"},
    
    {"wrap360", FH_WRAPTO360, 1, 1, FF_SCALAR, "Angle Conversion", "wrap360(x)",
     "Wrap angle to [0, 360]. Alias for wrapTo360.",
     {"wrap360(-90) = 270", NULL, NULL, NULL, NULL, NULL},
     "wrapTo360"},
    
    /* ===== SET OPERATIONS ===== */
    {"union", FH_UNION, 2, 2, FF_MATRIX, "Set Operations", "union(a, b)",
     "Set union of vectors.",
     {"union([1,2,3], [2,3,4])", "union([1,1,2], [2,3,3])", NULL, NULL, NULL, NULL},
     "intersect, setdiff"},
    
    {"intersect", FH_INTERSECT, 2, 2, FF_MATRIX, "Set Operations", "intersect(a, b)",
     "Set intersection of vectors.",
     {"intersect([1,2,3], [2,3,4])", "intersect([1,2,3], [4,5,6])", NULL, NULL, NULL, NULL},
     "union, setdiff"},
    
    {"setdiff", FH_SETDIFF, 2, 2, FF_MATRIX, "Set Operations", "setdiff(a, b)",
     "Set difference (elements in a but not b).",
     {"setdiff([1,2,3,4], [2,4])", "setdiff([1,2,3], [1,2,3])", NULL, NULL, NULL, NULL},
     "union, intersect"},
    
    {"setxor", FH_SETXOR, 2, 2, FF_MATRIX, "Set Operations", "setxor(a, b)",
     "Set exclusive or (elements in a or b but not both).",
     {"setxor([1,2,3], [2,3,4])", "setxor([1,2], [3,4])", NULL, NULL, NULL, NULL},
     "union, setdiff"},
    
    /* ===== LOGICAL OPERATORS ===== */
    {"and", FH_AND, 2, 2, FF_SCALAR, "Logical", "and(a, b)",
     "Logical AND.",
     {"and(1, 1) = true", "and(1, 0) = false", NULL, NULL, NULL, NULL},
     "or, not, xor"},
    
    {"band", FH_AND, 2, 2, FF_SCALAR, "Logical", "band(a, b)",
     "Logical AND. Alias for and.",
     {"band(1, 1) = true", NULL, NULL, NULL, NULL, NULL},
     "and"},
    
    {"or", FH_OR, 2, 2, FF_SCALAR, "Logical", "or(a, b)",
     "Logical OR.",
     {"or(1, 0) = true", "or(0, 0) = false", NULL, NULL, NULL, NULL},
     "and, not, xor"},
    
    {"bor", FH_OR, 2, 2, FF_SCALAR, "Logical", "bor(a, b)",
     "Logical OR. Alias for or.",
     {"bor(1, 0) = true", NULL, NULL, NULL, NULL, NULL},
     "or"},
    
    {"not", FH_NOT, 1, 1, FF_SCALAR, "Logical", "not(a)",
     "Logical NOT.",
     {"not(0) = true", "not(1) = false", NULL, NULL, NULL, NULL},
     "and, or"},
    
    {"bnot", FH_NOT, 1, 1, FF_SCALAR, "Logical", "bnot(a)",
     "Logical NOT. Alias for not.",
     {"bnot(0) = true", NULL, NULL, NULL, NULL, NULL},
     "not"},
    
    {"xor", FH_XOR, 2, 2, FF_SCALAR, "Logical", "xor(a, b)",
     "Logical XOR.",
     {"xor(1, 0) = true", "xor(1, 1) = false", NULL, NULL, NULL, NULL},
     "and, or"},
    
    {"bxor", FH_XOR, 2, 2, FF_SCALAR, "Logical", "bxor(a, b)",
     "Logical XOR. Alias for xor.",
     {"bxor(1, 0) = true", NULL, NULL, NULL, NULL, NULL},
     "xor, xnor"},
    
    {"xnor", FH_XNOR, 2, 2, FF_SCALAR, "Logical", "xnor(a, b)",
     "Logical XNOR (not xor). True if both equal.",
     {"xnor(1, 1) = true", "xnor(0, 0) = true", "xnor(1, 0) = false", NULL, NULL, NULL},
     "xor, nand, nor"},
    
    {"nand", FH_NAND, 2, 2, FF_SCALAR, "Logical", "nand(a, b)",
     "Logical NAND (not and).",
     {"nand(1, 1) = false", "nand(1, 0) = true", NULL, NULL, NULL, NULL},
     "and, nor"},
    
    {"nor", FH_NOR, 2, 2, FF_SCALAR, "Logical", "nor(a, b)",
     "Logical NOR (not or).",
     {"nor(0, 0) = true", "nor(1, 0) = false", NULL, NULL, NULL, NULL},
     "or, nand"},
    
    {"implies", FH_IMPLIES, 2, 2, FF_SCALAR, "Logical", "implies(a, b)",
     "Logical implication (not a or b).",
     {"implies(0, 0) = true", "implies(1, 0) = false", "implies(0, 1) = true", NULL, NULL, NULL},
     "and, or"},
    
    {"imply", FH_IMPLIES, 2, 2, FF_SCALAR, "Logical", "imply(a, b)",
     "Logical implication. Alias for implies.",
     {"imply(0, 1) = true", NULL, NULL, NULL, NULL, NULL},
     "implies"},
    
    /* ===== COMPARISON ===== */
    {"eq", FH_EQ, 2, 2, FF_SCALAR | FF_MATRIX, "Logical", "eq(a, b)",
     "Equal (a == b).",
     {"eq(3, 3) = true", "eq(3, 4) = false", NULL, NULL, NULL, NULL},
     "ne, lt, le, gt, ge"},
    
    {"ne", FH_NE, 2, 2, FF_SCALAR | FF_MATRIX, "Logical", "ne(a, b)",
     "Not equal (a != b).",
     {"ne(3, 4) = true", "ne(3, 3) = false", NULL, NULL, NULL, NULL},
     "eq"},
    
    {"lt", FH_LT, 2, 2, FF_SCALAR | FF_MATRIX, "Logical", "lt(a, b)",
     "Less than (a < b).",
     {"lt(3, 4) = true", "lt(4, 3) = false", NULL, NULL, NULL, NULL},
     "le, gt, ge"},
    
    {"le", FH_LE, 2, 2, FF_SCALAR | FF_MATRIX, "Logical", "le(a, b)",
     "Less than or equal (a <= b).",
     {"le(3, 3) = true", "le(4, 3) = false", NULL, NULL, NULL, NULL},
     "lt, ge"},
    
    {"gt", FH_GT, 2, 2, FF_SCALAR | FF_MATRIX, "Logical", "gt(a, b)",
     "Greater than (a > b).",
     {"gt(4, 3) = true", "gt(3, 4) = false", NULL, NULL, NULL, NULL},
     "ge, lt, le"},
    
    {"ge", FH_GE, 2, 2, FF_SCALAR | FF_MATRIX, "Logical", "ge(a, b)",
     "Greater than or equal (a >= b).",
     {"ge(3, 3) = true", "ge(3, 4) = false", NULL, NULL, NULL, NULL},
     "gt, le"},
    
    {"approxeq", FH_APPROXEQ, 2, 3, FF_SCALAR, "Logical", "approxeq(a, b) or approxeq(a, b, tol)",
     "Approximately equal within tolerance.",
     {"approxeq(0.1+0.2, 0.3) = true", "approxeq(1, 1.01, 0.1) = true", NULL, NULL, NULL, NULL},
     "eq"},
    
    {"isapprox", FH_APPROXEQ, 2, 3, FF_SCALAR, "Logical", "isapprox(a, b)",
     "Approximately equal. Alias for approxeq.",
     {"isapprox(0.1+0.2, 0.3) = true", NULL, NULL, NULL, NULL, NULL},
     "approxeq"},
    
    {"isequal", FH_ISEQUAL, 2, 2, FF_SCALAR | FF_MATRIX, "Logical", "isequal(a, b)",
     "Test equality of arrays.",
     {"isequal([1,2,3], [1,2,3]) = true", "isequal([1,2], [1,3]) = false", NULL, NULL, NULL, NULL},
     "eq, approxeq"},
    
    /* ===== COORDINATE TRANSFORMS ===== */
    {"cart2sph", FH_CART2SPH, 3, 3, FF_SCALAR, "Coordinates", "cart2sph(x, y, z)",
     "Cartesian to spherical coordinates. Returns [r, theta, phi].",
     {"cart2sph(1, 0, 0)", "cart2sph(0, 0, 1)", NULL, NULL, NULL, NULL},
     "sph2cart, cart2pol"},
    
    {"sph2cart", FH_SPH2CART, 3, 3, FF_SCALAR, "Coordinates", "sph2cart(r, theta, phi)",
     "Spherical to Cartesian coordinates. Returns [x, y, z].",
     {"sph2cart(1, 0, 0)", "sph2cart(1, pi/2, 0)", NULL, NULL, NULL, NULL},
     "cart2sph, pol2cart"},
    
    /* ===== DATA SCIENCE ===== */
    {"pca", FH_PCA, 1, 1, FF_MATRIX, "Data Science", "pca(X)",
     "Principal component analysis. Returns principal components.",
     {"pca([1,2;3,4;5,6])", "pca([1,0;0,1;1,1])", NULL, NULL, NULL, NULL},
     "svd, cov"},
    
    {"pcareduce", FH_PCAREDUCE, 2, 2, FF_MATRIX, "Data Science", "pcareduce(X, k)",
     "Reduce to k principal components.",
     {"pcareduce([1,2,3;4,5,6;7,8,9], 2)", NULL, NULL, NULL, NULL, NULL},
     "pca"},
    
    {"kmeans", FH_KMEANS, 2, 2, FF_MATRIX, "Data Science", "kmeans(X, k)",
     "K-means clustering.",
     {"kmeans([1,1;1,2;8,8;8,9], 2)", NULL, NULL, NULL, NULL, NULL},
     "pdist, silhouette, scatter, lineplot"},
    
    {"lineplot", FH_LINEPLOT, 1, 2, FF_NONE, "Data Science", "lineplot(y) or lineplot(x, y)",
     "ASCII line plot of vector y, optionally with x values.",
     {"lineplot(sin(linspace(0,2*pi,50)'))", "lineplot(1:10, rand(10,1))", NULL, NULL, NULL, NULL},
     "scatter, plot"},
    
    {"scatter", FH_SCATTER, 2, 4, FF_NONE, "Data Science", "scatter(x, y, [groups], [chars])",
     "ASCII scatter plot. Groups color points (1,2,3...). Chars sets markers.",
     {"scatter(rand(20,1), rand(20,1))", NULL, NULL, NULL, NULL, NULL},
     "kmeans, pca"},
    
    {"silhouette", FH_SILHOUETTE, 2, 2, FF_MATRIX, "Data Science", "silhouette(X, idx)",
     "Silhouette values for clustering quality.",
     {"silhouette([1,1;2,2;8,8;9,9], [1;1;2;2])", NULL, NULL, NULL, NULL, NULL},
     "kmeans"},
    
    {"pdist", FH_PDIST, 1, 1, FF_MATRIX, "Data Science", "pdist(X)",
     "Pairwise distances between rows.",
     {"pdist([0,0;3,4;0,5])", "pdist([1,0;0,1;1,1])", NULL, NULL, NULL, NULL},
     "dist, kmeans"},
    
    {"dist", FH_PDIST, 2, 2, FF_MATRIX, "Data Science", "dist(a, b)",
     "Euclidean distance between vectors.",
     {"dist([0,0], [3,4]) = 5", "dist([1,2,3], [4,5,6])", NULL, NULL, NULL, NULL},
     "pdist, norm"},
    
    {"distance", FH_PDIST, 2, 2, FF_MATRIX, "Data Science", "distance(a, b)",
     "Distance. Alias for dist.",
     {"distance([0,0], [3,4]) = 5", NULL, NULL, NULL, NULL, NULL},
     "dist"},
    
    {"manhattan", FH_MANHATTAN, 2, 2, FF_MATRIX, "Data Science", "manhattan(a, b)",
     "Manhattan (L1) distance.",
     {"manhattan([0,0], [3,4]) = 7", "manhattan([1,2], [4,6]) = 7", NULL, NULL, NULL, NULL},
     "dist"},
    
    /* ===== MISC ===== */
    {"int", FH_TRUNC, 1, 1, FF_SCALAR, "Rounding", "int(x)",
     "Integer part. Alias for trunc.",
     {"int(3.7) = 3", "int(-3.7) = -3", NULL, NULL, NULL, NULL},
     "trunc, floor"},
    
    {"isinteger", FH_ISINTEGER, 1, 1, FF_SCALAR, "Logical", "isinteger(x)",
     "Test if x is an integer.",
     {"isinteger(5) = true", "isinteger(5.5) = false", NULL, NULL, NULL, NULL},
     "floor"},
    
    {"isnumeric", FH_ISNUMERIC, 1, 1, FF_SCALAR | FF_MATRIX, "Logical", "isnumeric(x)",
     "Test if numeric (always true for numbers).",
     {"isnumeric(5) = true", "isnumeric([1,2,3]) = true", NULL, NULL, NULL, NULL},
     "isfinite"},
    
    {"islogical", FH_ISLOGICAL, 1, 1, FF_SCALAR, "Logical", "islogical(x)",
     "Test if x is 0 or 1.",
     {"islogical(0) = true", "islogical(1) = true", "islogical(2) = false", NULL, NULL, NULL},
     "isinteger"},
    
    {"isnegative", FH_ISNEGATIVE, 1, 1, FF_SCALAR, "Logical", "isnegative(x)",
     "Test if x < 0.",
     {"isnegative(-5) = true", "isnegative(5) = false", NULL, NULL, NULL, NULL},
     "ispositive, sign"},
    
    {"ispositive", FH_ISPOSITIVE, 1, 1, FF_SCALAR, "Logical", "ispositive(x)",
     "Test if x > 0.",
     {"ispositive(5) = true", "ispositive(-5) = false", NULL, NULL, NULL, NULL},
     "isnegative, sign"},
    
    /* ===== DATETIME FUNCTIONS ===== */
    {"datenum", FH_DATENUM, 1, 6, FF_SCALAR, "Datetime", "datenum(y,m,d,h,mi,s)",
     "Convert date to Unix timestamp (seconds since 1970-01-01).",
     {"datenum(2024,1,1,0,0,0)", "datenum(2024,6,15,12,30,0)", NULL, NULL, NULL, NULL},
     "datetime, year, month, day"},
    
    {"datetime", FH_DATETIME, 1, 1, FF_SCALAR, "Datetime", "datetime(timestamps)",
     "Create datetime array from cell array of date strings.",
     {"datetime({'2024-01-15'; '2024-02-20'})", NULL, NULL, NULL, NULL, NULL},
     "datenum, dateshift"},
    
    {"dateshift", FH_DATESHIFT, 1, 3, FF_MATRIX, "Datetime", "dateshift(t, 'start', unit)",
     "Shift datetime to start of period (month/year/day). Vectorized.",
     {"dateshift(t, 'start', 'month')", "months = dateshift(TT.times)", NULL, NULL, NULL, NULL},
     "startofmonth, startofyear, datetime"},
    
    {"year", FH_YEAR, 1, 1, FF_SCALAR, "Datetime", "year(timestamp)",
     "Extract year from Unix timestamp.",
     {"year(datenum(2024,6,15,0,0,0)) = 2024", NULL, NULL, NULL, NULL, NULL},
     "month, day, hour"},
    
    {"month", FH_MONTH, 1, 1, FF_SCALAR, "Datetime", "month(timestamp)",
     "Extract month (1-12) from Unix timestamp.",
     {"month(datenum(2024,6,15,0,0,0)) = 6", NULL, NULL, NULL, NULL, NULL},
     "year, day, hour"},
    
    {"day", FH_DAY, 1, 1, FF_SCALAR, "Datetime", "day(timestamp)",
     "Extract day of month (1-31) from Unix timestamp.",
     {"day(datenum(2024,6,15,0,0,0)) = 15", NULL, NULL, NULL, NULL, NULL},
     "year, month, hour"},
    
    {"hour", FH_HOUR, 1, 1, FF_SCALAR, "Datetime", "hour(timestamp)",
     "Extract hour (0-23) from Unix timestamp.",
     {"hour(datenum(2024,6,15,14,30,0)) = 14", NULL, NULL, NULL, NULL, NULL},
     "minute, second, day"},
    
    {"minute", FH_MINUTE, 1, 1, FF_SCALAR, "Datetime", "minute(timestamp)",
     "Extract minute (0-59) from Unix timestamp.",
     {"minute(datenum(2024,6,15,14,30,0)) = 30", NULL, NULL, NULL, NULL, NULL},
     "hour, second"},
    
    {"second", FH_SECOND, 1, 1, FF_SCALAR, "Datetime", "second(timestamp)",
     "Extract second (0-59) from Unix timestamp.",
     {"second(datenum(2024,6,15,14,30,45)) = 45", NULL, NULL, NULL, NULL, NULL},
     "hour, minute"},
    
    {"days", FH_DAYS, 1, 1, FF_SCALAR, "Datetime", "days(n)",
     "Convert days to seconds.",
     {"days(1) = 86400", "days(7) = 604800", NULL, NULL, NULL, NULL},
     "weeks, hours, minutes, seconds"},
    
    {"weeks", FH_WEEKS, 1, 1, FF_SCALAR, "Datetime", "weeks(n)",
     "Convert weeks to seconds.",
     {"weeks(1) = 604800", "weeks(2) = 1209600", NULL, NULL, NULL, NULL},
     "days, hours, minutes"},
    
    {"hours", FH_HOURS, 1, 1, FF_SCALAR, "Datetime", "hours(n)",
     "Convert hours to seconds.",
     {"hours(1) = 3600", "hours(24) = 86400", NULL, NULL, NULL, NULL},
     "days, minutes, seconds"},
    
    {"minutes", FH_MINUTES, 1, 1, FF_SCALAR, "Datetime", "minutes(n)",
     "Convert minutes to seconds.",
     {"minutes(1) = 60", "minutes(60) = 3600", NULL, NULL, NULL, NULL},
     "hours, seconds, days"},
    
    {"seconds", FH_SECONDS, 1, 1, FF_SCALAR, "Datetime", "seconds(n)",
     "Return seconds (identity function for consistency).",
     {"seconds(60) = 60", NULL, NULL, NULL, NULL, NULL},
     "minutes, hours, days"},
    
    {"startofmonth", FH_STARTOFMONTH, 1, 1, FF_SCALAR, "Datetime", "startofmonth(t)",
     "Get timestamp for start of month containing t. Alias: som.",
     {"startofmonth(datenum(2024,6,15,0,0,0))", NULL, NULL, NULL, NULL, NULL},
     "dateshift, startofyear, startofday"},
    
    {"som", FH_STARTOFMONTH, 1, 1, FF_SCALAR, "Datetime", "som(t)",
     "Alias for startofmonth.",
     {"som(t)", NULL, NULL, NULL, NULL, NULL},
     "startofmonth"},
    
    {"startofyear", FH_STARTOFYEAR, 1, 1, FF_SCALAR, "Datetime", "startofyear(t)",
     "Get timestamp for start of year containing t. Alias: soy.",
     {"startofyear(datenum(2024,6,15,0,0,0))", NULL, NULL, NULL, NULL, NULL},
     "startofmonth, startofday"},
    
    {"soy", FH_STARTOFYEAR, 1, 1, FF_SCALAR, "Datetime", "soy(t)",
     "Alias for startofyear.",
     {"soy(t)", NULL, NULL, NULL, NULL, NULL},
     "startofyear"},
    
    {"startofday", FH_STARTOFDAY, 1, 1, FF_SCALAR, "Datetime", "startofday(t)",
     "Get timestamp for start of day containing t. Alias: sod.",
     {"startofday(datenum(2024,6,15,14,30,0))", NULL, NULL, NULL, NULL, NULL},
     "startofmonth, startofyear"},
    
    {"sod", FH_STARTOFDAY, 1, 1, FF_SCALAR, "Datetime", "sod(t)",
     "Alias for startofday.",
     {"sod(t)", NULL, NULL, NULL, NULL, NULL},
     "startofday"},
    
    /* ===== TABLE FUNCTIONS ===== */
    {"timetable", FH_TIMETABLE, 0, -1, FF_MATRIX, "Tables", "timetable name col1 col2 ...",
     "Create timetable from column variables. First column is row times.",
     {"timetable TT times custid usage", "timetable Events t values", NULL, NULL, NULL, NULL},
     "table, groupsummary, sortrows, height"},
    
    {"table", FH_TABLE, 0, -1, FF_MATRIX, "Tables", "table name col1 col2 ...",
     "Create table from column variables.",
     {"table T x y z", NULL, NULL, NULL, NULL, NULL},
     "timetable, groupsummary"},
    
    {"height", FH_HEIGHT, 1, 1, FF_SCALAR, "Tables", "height(T)",
     "Number of rows in table T.",
     {"height(TT)", "n = height(MyTable)", NULL, NULL, NULL, NULL},
     "width, size, numel"},
    
    {"width", FH_WIDTH, 1, 1, FF_SCALAR, "Tables", "width(T)",
     "Number of columns in table T.",
     {"width(TT)", NULL, NULL, NULL, NULL, NULL},
     "height, size"},
    
    {"groupsummary", FH_GROUPSUMMARY, 0, -1, FF_MATRIX, "Tables", "groupsummary Result T groupcol method datacol",
     "Group table by column and apply aggregation (sum/mean/min/max/numel).",
     {"groupsummary MonthlySum TT month sum usage", "groupsummary Counts TT custid numel custid", NULL, NULL, NULL, NULL},
     "timetable, sortrows, unique"},
    
    {"groupcounts", FH_GROUPCOUNTS, 0, -1, FF_MATRIX, "Tables", "groupcounts Result T groupcol",
     "Count rows in each group.",
     {"groupcounts CustCounts TT custid", NULL, NULL, NULL, NULL, NULL},
     "groupsummary, unique, numel"},
    
    {"sortrows", FH_SORTROWS_TBL, 1, 2, FF_MATRIX, "Tables", "sortrows tablename",
     "Sort table by first column (ascending).",
     {"sortrows MonthlyUsage", NULL, NULL, NULL, NULL, NULL},
     "sort, timetable"},
    
    {"head", FH_HEAD, 1, 2, FF_MATRIX, "Tables", "head(T, n)",
     "First n rows of table (default n=5).",
     {"head(TT)", "head(TT, 10)", NULL, NULL, NULL, NULL},
     "tail, height"},
    
    {"tail", FH_TAIL, 1, 2, FF_MATRIX, "Tables", "tail(T, n)",
     "Last n rows of table (default n=5).",
     {"tail(TT)", "tail(TT, 10)", NULL, NULL, NULL, NULL},
     "head, height"},
    
    /* ===== FORECASTING FUNCTIONS ===== */
    {"movavg", FH_MOVAVG, 2, 2, FF_MATRIX, "Forecasting", "movavg(x, window)",
     "Simple moving average with specified window size.",
     {"movavg([1;2;3;4;5], 3)", "movavg([10;20;30;40;50], 2)", NULL, NULL, NULL, NULL},
     "ewma, movmean, movsum"},
    
    {"ewma", FH_EWMA, 2, 2, FF_MATRIX, "Forecasting", "ewma(x, alpha)",
     "Exponentially weighted moving average. Alpha in (0,1).",
     {"ewma([1;2;3;4;5], 0.3)", "ewma([10;20;30;40;50], 0.5)", NULL, NULL, NULL, NULL},
     "movavg, expsmooth"},
    
    {"trend", FH_TREND, 1, 1, FF_MATRIX, "Forecasting", "trend(x)",
     "Extract linear trend component from time series.",
     {"trend([1;3;5;7;9])", "trend([10;12;15;18;20])", NULL, NULL, NULL, NULL},
     "detrend, season"},
    
    {"detrend", FH_DETREND, 1, 2, FF_MATRIX, "Forecasting", "detrend(x, type)",
     "Remove trend from data. Type: 'linear' (default) or 'constant'.",
     {"detrend([1;3;5;7;9])", "detrend([10;12;14;16;18])", NULL, NULL, NULL, NULL},
     "trend, diff"},
    
    {"lag", FH_LAG, 1, 2, FF_MATRIX, "Forecasting", "lag(x, k)",
     "Lag time series by k periods (default k=1). Prepends NaN.",
     {"lag([1;2;3;4;5])", "lag([10;20;30;40;50], 2)", NULL, NULL, NULL, NULL},
     "lead, diff"},
    
    {"lead", FH_LEAD, 1, 2, FF_MATRIX, "Forecasting", "lead(x, k)",
     "Lead time series by k periods (default k=1). Appends NaN.",
     {"lead([1;2;3;4;5])", "lead([10;20;30;40;50], 2)", NULL, NULL, NULL, NULL},
     "lag, diff"},
    
    {"autocorr", FH_AUTOCORR, 1, 2, FF_MATRIX, "Forecasting", "autocorr(x, maxlag)",
     "Autocorrelation function up to maxlag.",
     {"autocorr([1;2;3;4;5])", "autocorr([1;2;1;2;1;2], 3)", NULL, NULL, NULL, NULL},
     "crosscorr, corrcoef"},
    
    {"crosscorr", FH_CROSSCORR, 2, 3, FF_MATRIX, "Forecasting", "crosscorr(x, y, maxlag)",
     "Cross-correlation between two series.",
     {"crosscorr([1;2;3], [1;2;3])", NULL, NULL, NULL, NULL, NULL},
     "autocorr, corrcoef"},
    
    {"expsmooth", FH_EXPSMOOTH, 2, 2, FF_MATRIX, "Forecasting", "expsmooth(x, alpha)",
     "Simple exponential smoothing for forecasting.",
     {"expsmooth([1;2;3;4;5], 0.3)", NULL, NULL, NULL, NULL, NULL},
     "ewma, movavg, forecast"},
    
    {"forecast", FH_FORECAST, 2, 3, FF_MATRIX, "Forecasting", "forecast(x, periods, method)",
     "Forecast future values. Methods: 'linear', 'exp', 'naive'.",
     {"forecast([1;2;3;4;5], 3)", NULL, NULL, NULL, NULL, NULL},
     "trend, expsmooth"},
    
    {"pctchange", FH_DIFF_TS, 1, 1, FF_MATRIX, "Forecasting", "pctchange(x)",
     "Percentage change between consecutive elements.",
     {"pctchange([100;110;121;133])", NULL, NULL, NULL, NULL, NULL},
     "diff, lag"},
    
    {"percentchange", FH_DIFF_TS, 1, 1, FF_MATRIX, "Forecasting", "percentchange(x)",
     "Alias for pctchange.",
     {"percentchange([100;110;121])", NULL, NULL, NULL, NULL, NULL},
     "pctchange"},
    
    /* ===== NUMBER THEORY - ADDITIONAL ===== */
    {"triangular", FH_TRIANGULAR, 1, 1, FF_SCALAR, "Number Theory", "triangular(n)",
     "nth triangular number: n*(n+1)/2.",
     {"triangular(1) = 1", "triangular(5) = 15", "triangular(10) = 55", NULL, NULL, NULL},
     "pentagonal, hexagonal, istriangular"},
    
    {"pentagonal", FH_PENTAGONAL, 1, 1, FF_SCALAR, "Number Theory", "pentagonal(n)",
     "nth pentagonal number: n*(3n-1)/2.",
     {"pentagonal(1) = 1", "pentagonal(5) = 35", NULL, NULL, NULL, NULL},
     "triangular, hexagonal"},
    
    {"hexagonal", FH_HEXAGONAL, 1, 1, FF_SCALAR, "Number Theory", "hexagonal(n)",
     "nth hexagonal number: n*(2n-1).",
     {"hexagonal(1) = 1", "hexagonal(5) = 45", NULL, NULL, NULL, NULL},
     "triangular, pentagonal, square_num"},
    
    {"square_num", FH_SQUARE_NUM, 1, 1, FF_SCALAR, "Number Theory", "square_num(n)",
     "nth square number: n^2.",
     {"square_num(1) = 1", "square_num(5) = 25", "square_num(10) = 100", NULL, NULL, NULL},
     "cube_num, triangular"},
    
    {"cube_num", FH_CUBE_NUM, 1, 1, FF_SCALAR, "Number Theory", "cube_num(n)",
     "nth cube number: n^3.",
     {"cube_num(1) = 1", "cube_num(3) = 27", "cube_num(5) = 125", NULL, NULL, NULL},
     "square_num, triangular"},
    
    {"istriangular", FH_ISTRIANGULAR, 1, 1, FF_SCALAR, "Number Theory", "istriangular(n)",
     "Test if n is a triangular number.",
     {"istriangular(10) = true", "istriangular(11) = false", NULL, NULL, NULL, NULL},
     "triangular, isprime"},
    
    {"isperfectsquare", FH_ISPERFECTSQUARE, 1, 1, FF_SCALAR, "Number Theory", "isperfectsquare(n)",
     "Test if n is a perfect square.",
     {"isperfectsquare(16) = true", "isperfectsquare(17) = false", NULL, NULL, NULL, NULL},
     "sqrt, isprime"},
    
    {"ispower", FH_ISPOWER, 1, 1, FF_SCALAR, "Number Theory", "ispower(n)",
     "Test if n is a perfect power (n = m^k for some m, k >= 2).",
     {"ispower(8) = true", "ispower(7) = false", "ispower(27) = true", NULL, NULL, NULL},
     "isperfectsquare, nthroot"},
    
    {"isperfectpower", FH_ISPERFECTPOWER, 1, 1, FF_SCALAR, "Number Theory", "isperfectpower(n)",
     "Alias for ispower.",
     {"isperfectpower(8) = true", NULL, NULL, NULL, NULL, NULL},
     "ispower"},
    
    {"reversedigits", FH_REVERSEDIGITS, 1, 1, FF_SCALAR, "Number Theory", "reversedigits(n)",
     "Reverse the digits of n.",
     {"reversedigits(123) = 321", "reversedigits(1000) = 1", NULL, NULL, NULL, NULL},
     "digitsum, ndigits"},
    
    {"revdigits", FH_REVERSEDIGITS, 1, 1, FF_SCALAR, "Number Theory", "revdigits(n)",
     "Alias for reversedigits.",
     {"revdigits(123) = 321", NULL, NULL, NULL, NULL, NULL},
     "reversedigits"},
    
    {"ispalindrome", FH_ISPALINDROME, 1, 1, FF_SCALAR, "Number Theory", "ispalindrome(n)",
     "Test if n reads same forwards and backwards.",
     {"ispalindrome(121) = true", "ispalindrome(123) = false", NULL, NULL, NULL, NULL},
     "reversedigits"},
    
    {"subfactorial", FH_SUBFACTORIAL, 1, 1, FF_SCALAR, "Number Theory", "subfactorial(n)",
     "Subfactorial (derangements): !n = n! * sum((-1)^k/k!, k=0..n).",
     {"subfactorial(0) = 1", "subfactorial(3) = 2", "subfactorial(5) = 44", NULL, NULL, NULL},
     "factorial, derangements"},
    
    {"derangements", FH_DERANGEMENTS, 1, 1, FF_SCALAR, "Number Theory", "derangements(n)",
     "Alias for subfactorial. Number of permutations with no fixed points.",
     {"derangements(5) = 44", NULL, NULL, NULL, NULL, NULL},
     "subfactorial"},
    
    {"bell", FH_BELL, 1, 1, FF_SCALAR, "Number Theory", "bell(n)",
     "nth Bell number: number of partitions of a set of n elements.",
     {"bell(0) = 1", "bell(3) = 5", "bell(5) = 52", NULL, NULL, NULL},
     "stirling2, partition"},
    
    {"stirling2", FH_STIRLING2, 2, 2, FF_SCALAR, "Number Theory", "stirling2(n, k)",
     "Stirling number of second kind: partitions of n into k non-empty subsets.",
     {"stirling2(4, 2) = 7", "stirling2(5, 3) = 25", NULL, NULL, NULL, NULL},
     "bell, ncr"},
    
    {"partition", FH_PARTITION, 1, 1, FF_SCALAR, "Number Theory", "partition(n)",
     "Number of ways to write n as sum of positive integers.",
     {"partition(5) = 7", "partition(10) = 42", NULL, NULL, NULL, NULL},
     "bell, stirling2"},
    
    /* ===== UTILITY FUNCTIONS ===== */
    {"rng", FH_RAND, 1, 1, FF_SCALAR, "Random", "rng(seed)",
     "Set random number generator seed for reproducibility.",
     {"rng(42)", "rng(0)", NULL, NULL, NULL, NULL},
     "rand, randn, randi"},
    
    {"categorical", FH_NONE, 1, 1, FF_MATRIX, "Data Types", "categorical(x)",
     "Create categorical array from numeric or cell array.",
     {"c = categorical([1;2;1;3])", NULL, NULL, NULL, NULL, NULL},
     "unique, groupsummary"},
    
    /* ===== ANGLE/TIME CONVERSIONS ===== */
    {"deg2dms", FH_NONE, 1, 1, FF_SCALAR, "Angle Conversion", "deg2dms(deg)",
     "Convert decimal degrees to [degrees, minutes, seconds] vector.",
     {"deg2dms(45.5125) = [45, 30, 45]", NULL, NULL, NULL, NULL, NULL},
     "dms2deg, deg2rad"},
    
    {"dms2deg", FH_NONE, 3, 3, FF_SCALAR, "Angle Conversion", "dms2deg(d, m, s)",
     "Convert degrees, minutes, seconds to decimal degrees.",
     {"dms2deg(45, 30, 45) = 45.5125", NULL, NULL, NULL, NULL, NULL},
     "deg2dms, deg2rad"},
    
    {"dms", FH_NONE, 1, 1, FF_SCALAR, "Angle Conversion", "dms(deg)",
     "Alias for deg2dms. Convert decimal degrees to DMS.",
     {"dms(45.5125)", NULL, NULL, NULL, NULL, NULL},
     "deg2dms"},
    
    {"hms2hr", FH_NONE, 3, 3, FF_SCALAR, "Time Conversion", "hms2hr(h, m, s)",
     "Convert hours, minutes, seconds to decimal hours.",
     {"hms2hr(2, 30, 0) = 2.5", NULL, NULL, NULL, NULL, NULL},
     "hr2hms"},
    
    {"hr2hms", FH_NONE, 1, 1, FF_SCALAR, "Time Conversion", "hr2hms(hours)",
     "Convert decimal hours to [hours, minutes, seconds] vector.",
     {"hr2hms(2.5) = [2, 30, 0]", NULL, NULL, NULL, NULL, NULL},
     "hms2hr"},
    
    {"hms", FH_NONE, 1, 1, FF_SCALAR, "Time Conversion", "hms(hours)",
     "Alias for hr2hms.",
     {"hms(2.5)", NULL, NULL, NULL, NULL, NULL},
     "hr2hms"},
    
    {"hr", FH_NONE, 3, 3, FF_SCALAR, "Time Conversion", "hr(h, m, s)",
     "Alias for hms2hr.",
     {"hr(2, 30, 0)", NULL, NULL, NULL, NULL, NULL},
     "hms2hr"},
    
    {"ms", FH_NONE, 1, 1, FF_SCALAR, "Datetime", "ms(n)",
     "Alias for milliseconds.",
     {"ms(1000)", NULL, NULL, NULL, NULL, NULL},
     "milliseconds"},
    
    {"milliseconds", FH_NONE, 1, 1, FF_SCALAR, "Datetime", "milliseconds(n)",
     "Convert milliseconds to seconds.",
     {"milliseconds(1000) = 1", NULL, NULL, NULL, NULL, NULL},
     "seconds, minutes"},
    
    /* ===== MATRIX PREDICATES ===== */
    {"isdiag", FH_NONE, 1, 1, FF_MATRIX, "Matrix", "isdiag(A)",
     "Test if A is a diagonal matrix.",
     {"isdiag(eye(3)) = true", "isdiag([1 2; 3 4]) = false", NULL, NULL, NULL, NULL},
     "diag, istriu, istril"},
    
    {"istriu", FH_NONE, 1, 1, FF_MATRIX, "Matrix", "istriu(A)",
     "Test if A is upper triangular.",
     {"istriu(triu(A))", NULL, NULL, NULL, NULL, NULL},
     "triu, istril, isdiag"},
    
    {"istril", FH_NONE, 1, 1, FF_MATRIX, "Matrix", "istril(A)",
     "Test if A is lower triangular.",
     {"istril(tril(A))", NULL, NULL, NULL, NULL, NULL},
     "tril, istriu, isdiag"},
    
    {"squeeze", FH_NONE, 1, 1, FF_MATRIX, "Matrix", "squeeze(A)",
     "Remove singleton dimensions from matrix.",
     {"squeeze(A)", NULL, NULL, NULL, NULL, NULL},
     "reshape, size"},
    
    /* ===== INDEX CONVERSION ===== */
    {"ind2sub", FH_NONE, 2, 2, FF_SCALAR, "Matrix", "ind2sub(size, idx)",
     "Convert linear index to row, column subscripts.",
     {"ind2sub([3,3], 5) = [2, 2]", NULL, NULL, NULL, NULL, NULL},
     "sub2ind"},
    
    {"sub2ind", FH_NONE, 3, 3, FF_SCALAR, "Matrix", "sub2ind(size, row, col)",
     "Convert row, column subscripts to linear index.",
     {"sub2ind([3,3], 2, 2) = 5", NULL, NULL, NULL, NULL, NULL},
     "ind2sub"},
    
    /* ===== BITWISE OPERATIONS ===== */
    {"lsl", FH_NONE, 2, 2, FF_SCALAR, "Bitwise", "lsl(x, n)",
     "Logical shift left. Alias for bitshift(x, n).",
     {"lsl(1, 4) = 16", NULL, NULL, NULL, NULL, NULL},
     "lsr, bitshift"},
    
    {"lsr", FH_NONE, 2, 2, FF_SCALAR, "Bitwise", "lsr(x, n)",
     "Logical shift right. Alias for bitshift(x, -n).",
     {"lsr(16, 4) = 1", NULL, NULL, NULL, NULL, NULL},
     "lsl, bitshift"},
    
    /* ===== NUMBER THEORY ADDITIONS ===== */
    {"numdivisors", FH_NONE, 1, 1, FF_SCALAR, "Number Theory", "numdivisors(n)",
     "Count of divisors of n. Alias for sigma0.",
     {"numdivisors(12) = 6", NULL, NULL, NULL, NULL, NULL},
     "divisors, divisorsum"},
    
    {"sigma0", FH_NONE, 1, 1, FF_SCALAR, "Number Theory", "sigma0(n)",
     "Divisor count function (number of divisors).",
     {"sigma0(12) = 6", NULL, NULL, NULL, NULL, NULL},
     "numdivisors, sigma"},
    
    {"sigmak", FH_NONE, 2, 2, FF_SCALAR, "Number Theory", "sigmak(n, k)",
     "Sum of kth powers of divisors.",
     {"sigmak(12, 1) = 28", "sigmak(12, 2) = 210", NULL, NULL, NULL, NULL},
     "sigma, divisorsum"},
    
    {"divisorsum_k", FH_NONE, 2, 2, FF_SCALAR, "Number Theory", "divisorsum_k(n, k)",
     "Alias for sigmak.",
     {"divisorsum_k(12, 2)", NULL, NULL, NULL, NULL, NULL},
     "sigmak"},
    
    {"stirling", FH_NONE, 2, 2, FF_SCALAR, "Number Theory", "stirling(n, k)",
     "Stirling number of first kind (unsigned).",
     {"stirling(4, 2) = 11", NULL, NULL, NULL, NULL, NULL},
     "stirling2, bell"},
    
    {"radical", FH_NONE, 1, 1, FF_SCALAR, "Number Theory", "radical(n)",
     "Product of distinct prime factors of n.",
     {"radical(12) = 6", "radical(100) = 10", NULL, NULL, NULL, NULL},
     "factor, isprime"},
    
    {"silver", FH_SILVER, 0, 0, FF_CONST, "Constants", "silver",
     "Silver ratio: 1 + sqrt(2) ≈ 2.414.",
     {"silver = 2.4142", NULL, NULL, NULL, NULL, NULL},
     "phi, sqrt2"},
    
    /* ===== FIGURATE NUMBERS ===== */
    {"tri_num", FH_TRIANGULAR, 1, 1, FF_SCALAR, "Number Theory", "tri_num(n)",
     "Alias for triangular.",
     {"tri_num(5) = 15", NULL, NULL, NULL, NULL, NULL},
     "triangular"},
    
    {"pent_num", FH_PENTAGONAL, 1, 1, FF_SCALAR, "Number Theory", "pent_num(n)",
     "Alias for pentagonal.",
     {"pent_num(5) = 35", NULL, NULL, NULL, NULL, NULL},
     "pentagonal"},
    
    {"hex_num", FH_HEXAGONAL, 1, 1, FF_SCALAR, "Number Theory", "hex_num(n)",
     "Alias for hexagonal.",
     {"hex_num(5) = 45", NULL, NULL, NULL, NULL, NULL},
     "hexagonal"},
    
    /* ===== ARRAY OPERATIONS ===== */
    {"maxk", FH_NONE, 2, 2, FF_MATRIX, "Statistics", "maxk(x, k)",
     "Return k largest elements.",
     {"maxk([3,1,4,1,5,9], 3) = [9,5,4]", NULL, NULL, NULL, NULL, NULL},
     "mink, max, sort"},
    
    {"mink", FH_NONE, 2, 2, FF_MATRIX, "Statistics", "mink(x, k)",
     "Return k smallest elements.",
     {"mink([3,1,4,1,5,9], 3) = [1,1,3]", NULL, NULL, NULL, NULL, NULL},
     "maxk, min, sort"},
    
    {"midpoint", FH_NONE, 2, 2, FF_SCALAR, "Arithmetic", "midpoint(a, b)",
     "Midpoint of a and b: (a+b)/2.",
     {"midpoint(2, 8) = 5", NULL, NULL, NULL, NULL, NULL},
     "mean, lerp"},
    
    {"mix", FH_NONE, 3, 3, FF_SCALAR, "Arithmetic", "mix(a, b, t)",
     "Linear interpolation: a*(1-t) + b*t. Alias for lerp.",
     {"mix(0, 10, 0.3) = 3", NULL, NULL, NULL, NULL, NULL},
     "lerp, midpoint"},
    
    {"remap", FH_NONE, 5, 5, FF_SCALAR, "Arithmetic", "remap(x, a1, b1, a2, b2)",
     "Map x from [a1,b1] to [a2,b2].",
     {"remap(5, 0, 10, 0, 100) = 50", NULL, NULL, NULL, NULL, NULL},
     "lerp, clamp"},
    
    {"map", FH_NONE, 5, 5, FF_SCALAR, "Arithmetic", "map(x, a1, b1, a2, b2)",
     "Alias for remap.",
     {"map(5, 0, 10, 0, 100) = 50", NULL, NULL, NULL, NULL, NULL},
     "remap"},
    
    /* ===== DISPLAY/OUTPUT ===== */
    {"print", FH_NONE, 1, 1, FF_SCALAR, "Display", "print(x)",
     "Display value of x.",
     {"print(42)", NULL, NULL, NULL, NULL, NULL},
     "disp, format"},
    
    {"printf", FH_NONE, 1, -1, FF_SCALAR, "Display", "printf(fmt, ...)",
     "Formatted output (C-style).",
     {"printf('x = %g', 3.14)", NULL, NULL, NULL, NULL, NULL},
     "fprintf, sprintf"},
    
    {"fprintf", FH_NONE, 1, -1, FF_SCALAR, "Display", "fprintf(fmt, ...)",
     "Formatted output to stdout.",
     {"fprintf('Result: %g\\n', x)", NULL, NULL, NULL, NULL, NULL},
     "printf"},
    
    {"printbin", FH_NONE, 1, 1, FF_SCALAR, "Display", "printbin(n)",
     "Print integer in binary.",
     {"printbin(42) = 101010", NULL, NULL, NULL, NULL, NULL},
     "printhex, printdec"},
    
    {"printhex", FH_NONE, 1, 1, FF_SCALAR, "Display", "printhex(n)",
     "Print integer in hexadecimal.",
     {"printhex(255) = FF", NULL, NULL, NULL, NULL, NULL},
     "printbin, printdec"},
    
    {"printdec", FH_NONE, 1, 1, FF_SCALAR, "Display", "printdec(n)",
     "Print integer in decimal (default).",
     {"printdec(42) = 42", NULL, NULL, NULL, NULL, NULL},
     "printbin, printhex"},
    
    {"ieee", FH_NONE, 1, 1, FF_SCALAR, "Display", "ieee(x)",
     "Show IEEE 754 representation of x.",
     {"ieee(1.5)", NULL, NULL, NULL, NULL, NULL},
     "format, hex"},
    
    {"accuracy", FH_NONE, 0, 0, FF_SCALAR, "Display", "accuracy",
     "Show current precision settings.",
     {"accuracy", NULL, NULL, NULL, NULL, NULL},
     "format, precision"},
    
    /* ===== ML FUNCTIONS ===== */
    {"fitcknn", FH_NONE, 3, 3, FF_MATRIX, "Machine Learning", "fitcknn(X, y, k)",
     "Train k-nearest neighbors classifier.",
     {"model = fitcknn(X, y, 5)", NULL, NULL, NULL, NULL, NULL},
     "predict, fitcsvm"},
    
    {"fitcsvm", FH_NONE, 2, 2, FF_MATRIX, "Machine Learning", "fitcsvm(X, y)",
     "Train support vector machine classifier.",
     {"model = fitcsvm(X, y)", NULL, NULL, NULL, NULL, NULL},
     "predict, fitcknn"},
    
    {"fitcnb", FH_NONE, 2, 2, FF_MATRIX, "Machine Learning", "fitcnb(X, y)",
     "Train naive Bayes classifier.",
     {"model = fitcnb(X, y)", NULL, NULL, NULL, NULL, NULL},
     "predict, fitcknn"},
    
    {"fitctree", FH_NONE, 2, 2, FF_MATRIX, "Machine Learning", "fitctree(X, y)",
     "Train decision tree classifier.",
     {"model = fitctree(X, y)", NULL, NULL, NULL, NULL, NULL},
     "predict, fitcknn"},
    
    {"predict", FH_NONE, 2, 2, FF_MATRIX, "Machine Learning", "predict(model, X)",
     "Predict labels using trained model.",
     {"y_pred = predict(knn_model, X_test)", NULL, NULL, NULL, NULL, NULL},
     "fitcknn, fitcsvm"},
    
    {"cvpartition", FH_NONE, 2, 2, FF_MATRIX, "Machine Learning", "cvpartition(n, k)",
     "Create k-fold cross-validation partition.",
     {"cv = cvpartition(100, 5)", NULL, NULL, NULL, NULL, NULL},
     "crossval"},
    
    {"confusionmat", FH_NONE, 2, 2, FF_MATRIX, "Machine Learning", "confusionmat(y_true, y_pred)",
     "Compute confusion matrix.",
     {"C = confusionmat(y, y_pred)", NULL, NULL, NULL, NULL, NULL},
     "accuracy, precision"},
    
    {"randindex", FH_NONE, 2, 2, FF_SCALAR, "Machine Learning", "randindex(labels1, labels2)",
     "Rand index for clustering comparison.",
     {"ri = randindex(cluster1, cluster2)", NULL, NULL, NULL, NULL, NULL},
     "kmeans, silhouette"},
    
    /* ===== STATISTICS ADDITIONS ===== */
    {"corr", FH_NONE, 2, 2, FF_SCALAR, "Statistics", "corr(x, y)",
     "Pearson correlation coefficient. Alias for corrcoef.",
     {"corr(x, y)", NULL, NULL, NULL, NULL, NULL},
     "corrcoef, cov"},
    
    {"crosstab", FH_NONE, 2, 2, FF_MATRIX, "Statistics", "crosstab(a, b)",
     "Cross-tabulation of two categorical vectors.",
     {"crosstab(gender, department)", NULL, NULL, NULL, NULL, NULL},
     "tabulate, groupcounts"},
    
    {"percent", FH_NONE, 2, 2, FF_SCALAR, "Arithmetic", "percent(part, whole)",
     "Calculate percentage: 100 * part / whole.",
     {"percent(25, 100) = 25", NULL, NULL, NULL, NULL, NULL},
     "pctchange"},
    
    /* ===== DECOMPOSITION OUTPUTS ===== */
    {"qrq", FH_NONE, 1, 1, FF_MATRIX, "Linear Algebra", "qrq(A)",
     "Q factor from QR decomposition.",
     {"Q = qrq(A)", NULL, NULL, NULL, NULL, NULL},
     "qr, qrr"},
    
    {"schurq", FH_NONE, 1, 1, FF_MATRIX, "Linear Algebra", "schurq(A)",
     "Unitary factor from Schur decomposition.",
     {"U = schurq(A)", NULL, NULL, NULL, NULL, NULL},
     "schur"},
    
    /* ===== MISC UTILITIES ===== */
    {"colon", FH_NONE, 2, 3, FF_MATRIX, "Matrix", "colon(a, b) or colon(a, step, b)",
     "Create range vector. Same as a:b or a:step:b.",
     {"colon(1, 5) = [1,2,3,4,5]", "colon(0, 0.5, 2)", NULL, NULL, NULL, NULL},
     "linspace, range"},
    
    {"linspace2", FH_NONE, 3, 3, FF_MATRIX, "Matrix", "linspace2(a, b, n)",
     "Alias for linspace.",
     {"linspace2(0, 1, 5)", NULL, NULL, NULL, NULL, NULL},
     "linspace"},
    
    {"breakdown", FH_NONE, 1, 1, FF_SCALAR, "Number Theory", "breakdown(n)",
     "Show prime factorization breakdown.",
     {"breakdown(360)", NULL, NULL, NULL, NULL, NULL},
     "factor, isprime"},
    
    {"explain", FH_NONE, 1, 1, FF_SCALAR, "Help", "explain(func)",
     "Explain function usage in detail.",
     {"explain(sin)", NULL, NULL, NULL, NULL, NULL},
     "help, demo"},
    
    {"plus", FH_NONE, 2, 2, FF_SCALAR, "Arithmetic", "plus(a, b)",
     "Addition: a + b.",
     {"plus(2, 3) = 5", NULL, NULL, NULL, NULL, NULL},
     "minus, times"},
    
    {"minus", FH_NONE, 2, 2, FF_SCALAR, "Arithmetic", "minus(a, b)",
     "Subtraction: a - b.",
     {"minus(5, 3) = 2", NULL, NULL, NULL, NULL, NULL},
     "plus, times"},
    
    {"times", FH_NONE, 2, 2, FF_SCALAR, "Arithmetic", "times(a, b)",
     "Element-wise multiplication: a .* b.",
     {"times(2, 3) = 6", NULL, NULL, NULL, NULL, NULL},
     "plus, rdivide"},
    
    {"rdivide", FH_NONE, 2, 2, FF_SCALAR, "Arithmetic", "rdivide(a, b)",
     "Element-wise right division: a ./ b.",
     {"rdivide(6, 3) = 2", NULL, NULL, NULL, NULL, NULL},
     "times, ldivide"},
    
    {"log_", FH_NONE, 2, 2, FF_SCALAR, "Logarithm", "log_(x, base)",
     "Logarithm of x with specified base.",
     {"log_(8, 2) = 3", NULL, NULL, NULL, NULL, NULL},
     "log, log10, log2"},
    
    /* ===== DATA CLEANING FUNCTIONS ===== */
    {"fillmissing", FH_FILLMISSING, 1, 2, FF_MATRIX, "Data Cleaning", "fillmissing(x, method)",
     "Fill NaN values. Methods: 'linear' (default), 'mean', 'previous'.",
     {"fillmissing([1;NaN;3], 'linear')", "fillmissing(data, 'mean')", NULL, NULL, NULL, NULL},
     "rmmissing, isnan, smoothdata"},
    
    {"rmmissing", FH_RMMISSING, 1, 1, FF_MATRIX, "Data Cleaning", "rmmissing(x)",
     "Remove NaN values from vector.",
     {"rmmissing([1;NaN;3;NaN;5])", NULL, NULL, NULL, NULL, NULL},
     "fillmissing, isnan"},
    
    {"smoothdata", FH_SMOOTHDATA, 1, 3, FF_MATRIX, "Data Cleaning", "smoothdata(x, method, window)",
     "Smooth data. Methods: 'movmean' (default), 'gaussian'. Window default=5.",
     {"smoothdata(noisy_data)", "smoothdata(x, 'gaussian', 7)", NULL, NULL, NULL, NULL},
     "movmean, fillmissing, movavg"},
    
    /* ===== CURVE FITTING FUNCTIONS ===== */
    {"polyfit", FH_POLYFIT, 3, 3, FF_MATRIX, "Curve Fitting", "polyfit(x, y, degree)",
     "Fit polynomial of specified degree. Returns coefficients [highest...lowest].",
     {"polyfit([1,2,3,4], [1,4,9,16], 2) = [1,0,0]", "coeffs = polyfit(x, y, 3)", NULL, NULL, NULL, NULL},
     "polyval, linreg"},
    
    {"linreg", FH_LINREG, 2, 2, FF_MATRIX, "Curve Fitting", "linreg(x, y)",
     "Simple linear regression. Returns [slope, intercept, R-squared].",
     {"linreg([1,2,3,4,5], [2,4,6,8,10]) = [2,0,1]", NULL, NULL, NULL, NULL, NULL},
     "polyfit, corrcoef, trend"},
    
    {"regress", FH_LINREG, 2, 2, FF_MATRIX, "Curve Fitting", "regress(x, y)",
     "Alias for linreg. Returns [slope, intercept, R-squared].",
     {"regress(x, y)", NULL, NULL, NULL, NULL, NULL},
     "linreg"},
    
    {"cumtrapz", FH_CUMTRAPZ, 1, 2, FF_MATRIX, "Integration", "cumtrapz(y) or cumtrapz(x, y)",
     "Cumulative trapezoidal numerical integration.",
     {"cumtrapz([1,2,3,4,5]) = [0,1.5,4,7.5,12]", NULL, NULL, NULL, NULL, NULL},
     "trapz, diff, cumsum"},
    
    {"quantile", FH_QUANTILE, 2, 2, FF_SCALAR, "Statistics", "quantile(x, p)",
     "Quantile at probability p (0 to 1). quantile(x, 0.5) = median.",
     {"quantile([1:10], 0.5) = 5", "quantile(data, 0.25)", NULL, NULL, NULL, NULL},
     "prctile, median, iqr"},
    
    /* ===== SAAS METRICS ===== */
    {"mrr", FH_MRR, 1, 1, FF_MATRIX, "SaaS", "mrr(data)",
     "Monthly Recurring Revenue from billing data [timestamp, custid, ..., billing].",
     {"mrr(rand(100,4))", "size(mrr(rand(1000,4)))", NULL, NULL, NULL, NULL},
     "arpu, mrrbridge, nrr"},
    
    {"arpu", FH_ARPU, 1, 1, FF_MATRIX, "SaaS", "arpu(data)",
     "Average Revenue Per User = MRR / active_customers per month.",
     {"arpu(rand(100,4))", NULL, NULL, NULL, NULL, NULL},
     "mrr, customercount"},
    
    {"mrrbridge", FH_MRRBRIDGE, 1, 1, FF_MATRIX, "SaaS", "mrrbridge(data)",
     "MRR movement: [month, new, expansion, contraction, churned, net].",
     {"mrrbridge(rand(100,4))", NULL, NULL, NULL, NULL, NULL},
     "mrr, nrr, churn"},
    
    {"customercount", FH_CUSTOMERCOUNT, 1, 1, FF_MATRIX, "SaaS", "customercount(data)",
     "Count unique active customers per month.",
     {"customercount(rand(100,4))", NULL, NULL, NULL, NULL, NULL},
     "newcustomers, churn"},
    
    {"newcustomers", FH_NEWCUSTOMERS, 1, 1, FF_MATRIX, "SaaS", "newcustomers(data)",
     "Count customers with first appearance each month.",
     {"newcustomers(rand(100,4))", NULL, NULL, NULL, NULL, NULL},
     "churn, reactivated"},
    
    {"churn", FH_CHURN, 1, 1, FF_MATRIX, "SaaS", "churn(data)",
     "Customers who churned (active last month, not this month).",
     {"churn(rand(100,4))", NULL, NULL, NULL, NULL, NULL},
     "churnrate, nrr"},
    
    {"reactivated", FH_REACTIVATED, 1, 1, FF_MATRIX, "SaaS", "reactivated(data)",
     "Returning customers (not active last month, active now, not new).",
     {"reactivated(rand(100,4))", NULL, NULL, NULL, NULL, NULL},
     "churn, newcustomers"},
    
    {"churnrate", FH_CHURNRATE, 1, 1, FF_MATRIX, "SaaS", "churnrate(data)",
     "Logo churn rate: churned / previous_customers * 100 (percent).",
     {"churnrate(rand(100,4))", NULL, NULL, NULL, NULL, NULL},
     "churn, nrr, grr"},
    
    {"nrr", FH_NRR, 1, 1, FF_MATRIX, "SaaS", "nrr(data)",
     "Net Revenue Retention %. >100% means existing customers growing.",
     {"nrr(rand(100,4))", NULL, NULL, NULL, NULL, NULL},
     "grr, churnrate"},
    
    {"grr", FH_GRR, 1, 1, FF_MATRIX, "SaaS", "grr(data)",
     "Gross Revenue Retention %. Like NRR but ignores expansion (max 100%).",
     {"grr(rand(100,4))", NULL, NULL, NULL, NULL, NULL},
     "nrr, churnrate"},
    
    {"retention", FH_RETENTION, 1, 1, FF_MATRIX, "SaaS", "retention(data)",
     "Cohort retention matrix. Row=cohort, Col=months since start, Value=%.",
     {"retention(rand(100,4))", NULL, NULL, NULL, NULL, NULL},
     "churn, tenure"},
    
    {"tenure", FH_TENURE, 1, 1, FF_MATRIX, "SaaS", "tenure(data)",
     "Customer tenure: [customerid, months, first_month, last_month].",
     {"tenure(rand(100,4))", NULL, NULL, NULL, NULL, NULL},
     "retention, churn"},
    
    {"toprevenue", FH_TOPREVENUE, 1, 2, FF_MATRIX, "SaaS", "toprevenue(data, n)",
     "Top N customers by revenue: [customerid, total_revenue, pct_of_total].",
     {"toprevenue(rand(100,4), 10)", NULL, NULL, NULL, NULL, NULL},
     "concentration, ltv"},
    
    {"concentration", FH_CONCENTRATION, 1, 1, FF_MATRIX, "SaaS", "concentration(data)",
     "Revenue concentration: [top1%, top5%, top10%, top20%, herfindahl_index].",
     {"concentration(rand(100,4))", NULL, NULL, NULL, NULL, NULL},
     "toprevenue"},
    
    {"revchurn", FH_REVCHURN, 1, 1, FF_MATRIX, "SaaS", "revchurn(data)",
     "Revenue churn rate: [month, churned_mrr, start_mrr, rev_churn%].",
     {"revchurn(rand(100,4))", NULL, NULL, NULL, NULL, NULL},
     "churnrate, netchurn"},
    
    {"netchurn", FH_NETCHURN, 1, 1, FF_MATRIX, "SaaS", "netchurn(data)",
     "Net churn rate: (churn + contraction - expansion) / start_mrr.",
     {"netchurn(rand(100,4))", NULL, NULL, NULL, NULL, NULL},
     "revchurn, quickratio"},
    
    {"quickratio", FH_QUICKRATIO, 1, 1, FF_MATRIX, "SaaS", "quickratio(data)",
     "SaaS Quick Ratio: (new + expansion) / (churn + contraction). >4 excellent.",
     {"quickratio(rand(100,4))", NULL, NULL, NULL, NULL, NULL},
     "mrrbridge, netchurn"},
    
    {"ltv", FH_LTV, 1, 1, FF_MATRIX, "SaaS", "ltv(data)",
     "Customer LTV: [id, total_rev, tenure, monthly_avg, estimated_ltv].",
     {"ltv(rand(100,4))", NULL, NULL, NULL, NULL, NULL},
     "tenure, arpu"},
    
    /* ===== FINANCIAL FUNCTIONS ===== */
    {"npv", FH_NONE, 2, 2, FF_SCALAR, "Financial", "npv(rate, cashflows)",
     "Net Present Value of cashflows at discount rate.",
     {"npv(0.1, [-1000;300;300;300;300])", NULL, NULL, NULL, NULL, NULL},
     "irr, payback, pv"},
    
    {"irr", FH_NONE, 1, 1, FF_SCALAR, "Financial", "irr(cashflows)",
     "Internal Rate of Return using Newton-Raphson.",
     {"irr([-1000;300;300;300;300]) = 0.077", NULL, NULL, NULL, NULL, NULL},
     "npv, payback"},
    
    {"payback", FH_NONE, 1, 1, FF_SCALAR, "Financial", "payback(cashflows)",
     "Payback period (when cumulative cashflow becomes positive).",
     {"payback([-1000;300;300;300;300]) = 4.33", NULL, NULL, NULL, NULL, NULL},
     "npv, irr"},
    
    {"cagr", FH_NONE, 3, 3, FF_SCALAR, "Financial", "cagr(start, end, years)",
     "Compound Annual Growth Rate.",
     {"cagr(100, 200, 5) = 0.1487", NULL, NULL, NULL, NULL, NULL},
     "compound, roi"},
    
    {"compound", FH_NONE, 4, 4, FF_SCALAR, "Financial", "compound(principal, rate, n, t)",
     "Compound interest: P*(1+r/n)^(n*t).",
     {"compound(1000, 0.05, 12, 10) = 1647", NULL, NULL, NULL, NULL, NULL},
     "cagr, pv, fv"},
    
    {"roi", FH_NONE, 2, 2, FF_SCALAR, "Financial", "roi(gain, cost)",
     "Return on Investment as percentage.",
     {"roi(1500, 1000) = 50", NULL, NULL, NULL, NULL, NULL},
     "cagr, irr"},
    
    {"sharpe", FH_SHARPE, 1, 2, FF_MATRIX, "Financial", "sharpe(returns) or sharpe(returns, rf)",
     "Sharpe ratio. Measures risk-adjusted return. rf is risk-free rate (default 0).",
     {"sharpe([0.1,0.05,-0.02,0.08])", "sharpe([0.1,0.05,-0.02,0.08], 0.01)", NULL, NULL, NULL, NULL},
     "std, mean, drawdown"},
    
    {"sharperatio", FH_SHARPE, 1, 2, FF_MATRIX, "Financial", "sharperatio(returns)",
     "Sharpe ratio. Alias for sharpe.",
     {"sharperatio([0.1,0.05,-0.02,0.08])", NULL, NULL, NULL, NULL, NULL},
     "sharpe"},
    
    {"drawdown", FH_DRAWDOWN, 1, 1, FF_MATRIX, "Financial", "drawdown(prices)",
     "Maximum drawdown. Largest peak-to-trough decline as percentage.",
     {"drawdown([100,110,105,120,90,95])", "drawdown([1,1.1,1.05,0.95])", NULL, NULL, NULL, NULL},
     "cumret, max, min"},
    
    {"cumret", FH_CUMRET, 1, 1, FF_MATRIX, "Financial", "cumret(returns)",
     "Cumulative returns. Converts period returns to cumulative wealth.",
     {"cumret([0.1,0.05,-0.02])", "cumret([0.01,0.02,-0.01,0.03])", NULL, NULL, NULL, NULL},
     "cumprod, drawdown"},
    
    {"cumreturn", FH_CUMRET, 1, 1, FF_MATRIX, "Financial", "cumreturn(returns)",
     "Cumulative returns. Alias for cumret.",
     {"cumreturn([0.1,0.05,-0.02])", NULL, NULL, NULL, NULL, NULL},
     "cumret"},
    
    {"winsorize", FH_WINSORIZE, 2, 3, FF_MATRIX, "Statistics", "winsorize(x, p) or winsorize(x, lo, hi)",
     "Winsorize data by limiting extreme values. p is percentile (e.g., 0.05 for 5%).",
     {"winsorize([1,2,3,100], 0.1)", "winsorize([1,2,3,4,5], 0.2, 0.8)", NULL, NULL, NULL, NULL},
     "clip, clamp, isoutlier"},
    
    {"bounds", FH_BOUNDS, 1, 1, FF_MATRIX, "Statistics", "bounds(x)",
     "Return [min, max] bounds of data.",
     {"bounds([3,1,4,1,5,9])", "bounds([1,2;3,4])", NULL, NULL, NULL, NULL},
     "min, max, range"},
    
    {"minmax", FH_MINMAX, 1, 1, FF_MATRIX, "Statistics", "minmax(x)",
     "Return [min, max] bounds. Alias for bounds.",
     {"minmax([3,1,4,1,5])", NULL, NULL, NULL, NULL, NULL},
     "bounds, min, max"},
    
    /* ===== SIGNAL ANALYSIS ===== */
    {"peaks", FH_NONE, 1, 1, FF_MATRIX, "Signal", "peaks(x)",
     "Find indices of local maxima.",
     {"peaks([1;3;2;5;4]) = [2;4]", NULL, NULL, NULL, NULL, NULL},
     "valleys, findpeaks, max"},
    
    {"findpeaks", FH_NONE, 1, 1, FF_MATRIX, "Signal", "findpeaks(x)",
     "Alias for peaks. Find indices of local maxima.",
     {"findpeaks([1;3;2;5;4])", NULL, NULL, NULL, NULL, NULL},
     "peaks, valleys"},
    
    {"valleys", FH_NONE, 1, 1, FF_MATRIX, "Signal", "valleys(x)",
     "Find indices of local minima.",
     {"valleys([3;1;4;2;5]) = [2;4]", NULL, NULL, NULL, NULL, NULL},
     "peaks, min"},
    
    {"findvalleys", FH_NONE, 1, 1, FF_MATRIX, "Signal", "findvalleys(x)",
     "Alias for valleys. Find indices of local minima.",
     {"findvalleys([3;1;4;2;5])", NULL, NULL, NULL, NULL, NULL},
     "valleys, peaks"},
    
    {"crossover", FH_NONE, 2, 2, FF_MATRIX, "Signal", "crossover(x, y)",
     "Find indices where x crosses y.",
     {"crossover(price, sma)", NULL, NULL, NULL, NULL, NULL},
     "peaks, valleys"},
    
    /* ===== DATA NORMALIZATION ===== */
    {"standardize", FH_NONE, 1, 1, FF_MATRIX, "Statistics", "standardize(x)",
     "Z-score normalization: (x - mean) / std.",
     {"standardize([1;2;3;4;5])", NULL, NULL, NULL, NULL, NULL},
     "zscore, minmax_scale, rescale"},
    
    {"zscore_vec", FH_NONE, 1, 1, FF_MATRIX, "Statistics", "zscore_vec(x)",
     "Alias for standardize. Z-score normalization for vectors.",
     {"zscore_vec([1;2;3;4;5])", NULL, NULL, NULL, NULL, NULL},
     "standardize"},
    
    {"minmax_scale", FH_NONE, 1, 3, FF_MATRIX, "Statistics", "minmax_scale(x, a, b)",
     "Scale data to [a, b] range (default [0, 1]).",
     {"minmax_scale([1;2;3;4;5])", "minmax_scale(x, 0, 100)", NULL, NULL, NULL, NULL},
     "standardize, rescale"},
    
    {"minmaxscale", FH_NONE, 1, 3, FF_MATRIX, "Statistics", "minmaxscale(x, a, b)",
     "Alias for minmax_scale.",
     {"minmaxscale([1;2;3;4;5], 0, 100)", NULL, NULL, NULL, NULL, NULL},
     "minmax_scale"},
    
    /* ===== CURVE FITTING ===== */
    {"polyfit", FH_POLYFIT, 3, 3, FF_MATRIX, "Polynomials", "polyfit(x, y, n)",
     "Fit polynomial of degree n to data points (x, y). Returns coefficients.",
     {"polyfit([1;2;3;4;5], [1;4;9;16;25], 2)", "p = polyfit(x, y, 3)", NULL, NULL, NULL, NULL},
     "polyval, polyder, interp1"},
    
    {"cumtrapz", FH_CUMTRAPZ, 1, 2, FF_MATRIX, "Integration", "cumtrapz(y) or cumtrapz(x, y)",
     "Cumulative trapezoidal integration.",
     {"cumtrapz([1;2;3;4;5])", "cumtrapz(x, y)", NULL, NULL, NULL, NULL},
     "trapz, diff, cumsum"},
    
    {"diff2", FH_DIFF2, 1, 1, FF_MATRIX, "Calculus", "diff2(v)",
     "Second difference (discrete second derivative).",
     {"diff2([1;4;9;16;25])", NULL, NULL, NULL, NULL, NULL},
     "diff, gradient"},
    
    {"standardize", FH_STANDARDIZE, 1, 1, FF_MATRIX, "Statistics", "standardize(v)",
     "Z-score standardization: (x - mean) / std. Alias: zscore_vec.",
     {"standardize([1;2;3;4;5])", NULL, NULL, NULL, NULL, NULL},
     "center, rescale, normalize"},
    
    {"zscore_vec", FH_STANDARDIZE, 1, 1, FF_MATRIX, "Statistics", "zscore_vec(v)",
     "Z-score standardization (alias for standardize).",
     {"zscore_vec([1;2;3;4;5])", NULL, NULL, NULL, NULL, NULL},
     "standardize"},
    
    {"center", FH_CENTER, 1, 1, FF_MATRIX, "Statistics", "center(v)",
     "Center data by subtracting mean: x - mean(x).",
     {"center([1;2;3;4;5])", NULL, NULL, NULL, NULL, NULL},
     "standardize, detrend"},
    
    /* Sentinel - must be last */
    {NULL, FH_NONE, 0, 0, 0, NULL, NULL, NULL, {NULL}, NULL}
};

/* Registry size (computed once) */
static int registry_size = -1;

/* ========== API Implementation ========== */

void func_registry_init(void)
{
    int i;
    if (registry_size < 0) {
        for (i = 0; func_registry[i].name != NULL; i++) {
            /* Just count */
        }
        registry_size = i;
    }
}

const FuncEntry *func_lookup(const char *name)
{
    int i;
    if (name == NULL) return NULL;
    
    func_registry_init();
    
    for (i = 0; i < registry_size; i++) {
        if (str_eq(func_registry[i].name, name)) {
            return &func_registry[i];
        }
    }
    return NULL;
}

int func_exists(const char *name)
{
    return func_lookup(name) != NULL;
}

FuncHandler func_get_handler(const char *name)
{
    const FuncEntry *entry = func_lookup(name);
    return entry ? entry->handler : FH_NONE;
}

int func_count(void)
{
    func_registry_init();
    return registry_size;
}

const FuncEntry *func_get_all(int *count)
{
    func_registry_init();
    if (count) *count = registry_size;
    return func_registry;
}

void func_help(const char *name)
{
    const FuncEntry *f = func_lookup(name);
    int i;
    
    if (!f) {
        printf("Unknown function: %s\n", name);
        return;
    }
    
    printf("\n%s - %s\n", f->name, f->category);
    printf("\342\224\201\342\224\201\342\224\201\342\224\201\342\224\201\342\224\201\342\224\201\342\224\201\342\224\201\342\224\201\342\224\201\342\224\201\342\224\201\342\224\201\342\224\201\342\224\201\342\224\201\342\224\201\342\224\201\342\224\201\342\224\201\342\224\201\342\224\201\342\224\201\342\224\201\342\224\201\342\224\201\342\224\201\342\224\201\342\224\201\342\224\201\342\224\201\342\224\201\342\224\201\342\224\201\342\224\201\342\224\201\342\224\201\342\224\201\342\224\201\n");
    printf("Syntax:      %s\n", f->syntax);
    printf("Description: %s\n", f->description);
    
    /* Examples */
    printf("\nExamples:\n");
    for (i = 0; i < 6 && f->examples[i] != NULL; i++) {
        printf("  %s\n", f->examples[i]);
    }
    
    if (f->see_also && f->see_also[0]) {
        printf("\nSee also: %s\n", f->see_also);
    }
    printf("\n");
}

void func_demo(const char *name)
{
    const FuncEntry *f = func_lookup(name);
    int i;
    char expr[256];
    const char *ex;
    size_t len;
    extern void mat_arena_reset(void);
    extern int cmd_load_dataset(const char *name);
    
    if (!f) {
        printf("Unknown function: %s\n", name);
        return;
    }
    
    printf("\n=== Demo: %s ===\n", f->name);
    printf("Description: %s\n\n", f->description);
    
    /* Special handling for SaaS category - load billing data first */
    if (str_eq(f->category, "SaaS")) {
        printf(">>> load hourlybilling  %% Loading SaaS demo data...\n");
        mat_arena_reset();
        if (!cmd_load_dataset("hourlybilling")) {
            printf("Error: Could not load hourlybilling dataset.\n");
            printf("Run 'generate' command first to create demo data.\n\n");
            return;
        }
        printf("\n");
        
        /* Now run the function on data */
        sprintf(expr, "%s(data)", f->name);
        printf(">>> %s\n", expr);
        eval_expr_line(expr, 0);
        printf("\n");
        return;
    }
    
    for (i = 0; i < 6 && f->examples[i] != NULL; i++) {
        ex = f->examples[i];
        
        /* Reset arena between examples to prevent memory exhaustion */
        mat_arena_reset();
        
        /* Extract expression (before '=') */
        len = 0;
        while (ex[len] && ex[len] != '=' && len < 254) {
            len++;
        }
        
        if (len > 0) {
            memcpy(expr, ex, len);
            expr[len] = '\0';
            
            /* Trim trailing whitespace */
            while (len > 0 && (expr[len-1] == ' ' || expr[len-1] == '\t')) {
                expr[--len] = '\0';
            }
            
            printf(">>> %s\n", expr);
            eval_expr_line(expr, 0);
            printf("\n");
        }
    }
}

void func_list(void)
{
    const char *printed_cats[30];
    int num_printed = 0;
    int i, j, c;
    int count;
    
    func_registry_init();
    
    printf("\nBuilt-in Functions (%d total):\n", registry_size);
    printf("=================================\n\n");
    
    /* Print by category */
    for (i = 0; i < registry_size; i++) {
        const char *cat = func_registry[i].category;
        int already = 0;
        
        /* Skip if deprecated */
        if (func_registry[i].flags & FF_DEPRECATED) continue;
        
        /* Check if we've printed this category */
        for (c = 0; c < num_printed; c++) {
            if (str_eq(printed_cats[c], cat)) {
                already = 1;
                break;
            }
        }
        
        if (!already && num_printed < 30) {
            /* New category - print header and all functions in it */
            printf("%s:\n  ", cat);
            printed_cats[num_printed++] = cat;
            count = 0;
            
            for (j = 0; j < registry_size; j++) {
                if (str_eq(func_registry[j].category, cat) &&
                    !(func_registry[j].flags & FF_DEPRECATED)) {
                    printf("%s ", func_registry[j].name);
                    count++;
                    if (count % 10 == 0) {
                        printf("\n  ");
                    }
                }
            }
            printf("\n\n");
        }
    }
    
    printf("Use 'help <func>' for details, 'demo <func>' for examples.\n\n");
}
