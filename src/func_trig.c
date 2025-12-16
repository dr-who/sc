/*
 * func_trig.c - Trigonometric function documentation
 * 
 * All trig functions with doxygen-style documentation for
 * automatic help, demo, and test generation.
 *
 * C89 compliant
 */
#include "func_defs.h"

/**
 * @func sin
 * @category Trigonometry
 * @syntax sin(x)
 * @desc Sine of x. Input angle is in radians by default.
 *       For complex x, sin(x) = (e^(ix) - e^(-ix)) / (2i).
 *       Use sind(x) for degrees input.
 * @example sin(0) = 0
 * @example sin(pi/6) = 0.5
 * @example sin(pi/2) = 1
 * @example sin(pi) = 0
 * @example sin(1+2i)
 * @see_also cos, tan, asin, sind
 */

/**
 * @func cos
 * @category Trigonometry
 * @syntax cos(x)
 * @desc Cosine of x. Input angle is in radians by default.
 *       For complex x, cos(x) = (e^(ix) + e^(-ix)) / 2.
 *       Use cosd(x) for degrees input.
 * @example cos(0) = 1
 * @example cos(pi/3) = 0.5
 * @example cos(pi/2) = 0
 * @example cos(pi) = -1
 * @example cos(2*pi) = 1
 * @see_also sin, tan, acos, cosd
 */

/**
 * @func tan
 * @category Trigonometry
 * @syntax tan(x)
 * @desc Tangent of x. Equal to sin(x)/cos(x).
 *       Undefined at x = pi/2 + n*pi.
 * @example tan(0) = 0
 * @example tan(pi/4) = 1
 * @example tan(pi/6) = 0.5774
 * @example tan(-pi/4) = -1
 * @see_also sin, cos, atan, tand
 */

/**
 * @func asin
 * @category Trigonometry
 * @syntax asin(x)
 * @desc Inverse sine (arcsine). Returns angle in radians.
 *       Domain: [-1, 1]. Range: [-pi/2, pi/2].
 *       For |x| > 1, returns complex result.
 * @example asin(0) = 0
 * @example asin(0.5) = 0.5236
 * @example asin(1) = 1.5708
 * @example asin(-1) = -1.5708
 * @example 6*asin(0.5)/pi = 1
 * @see_also sin, acos, atan, asind
 */

/**
 * @func acos
 * @category Trigonometry
 * @syntax acos(x)
 * @desc Inverse cosine (arccosine). Returns angle in radians.
 *       Domain: [-1, 1]. Range: [0, pi].
 *       For |x| > 1, returns complex result.
 * @example acos(1) = 0
 * @example acos(0.5) = 1.0472
 * @example acos(0) = 1.5708
 * @example acos(-1) = 3.1416
 * @example 3*acos(0.5)/pi = 1
 * @see_also cos, asin, atan, acosd
 */

/**
 * @func atan
 * @category Trigonometry
 * @syntax atan(x)
 * @desc Inverse tangent (arctangent). Returns angle in radians.
 *       Domain: all real numbers. Range: (-pi/2, pi/2).
 * @example atan(0) = 0
 * @example atan(1) = 0.7854
 * @example atan(-1) = -0.7854
 * @example 4*atan(1) = 3.1416
 * @example atan(1e10) = 1.5708
 * @see_also tan, asin, acos, atan2, atand
 */

/**
 * @func atan2
 * @category Trigonometry
 * @syntax atan2(y, x)
 * @desc Two-argument arctangent. Returns angle of point (x,y) from origin.
 *       Handles all quadrants correctly. Range: (-pi, pi].
 *       More robust than atan(y/x) for x near zero.
 * @example atan2(0, 1) = 0
 * @example atan2(1, 1) = 0.7854
 * @example atan2(1, 0) = 1.5708
 * @example atan2(1, -1) = 2.3562
 * @example atan2(-1, -1) = -2.3562
 * @see_also atan, atan2d
 */

/**
 * @func sinh
 * @category Trigonometry
 * @syntax sinh(x)
 * @desc Hyperbolic sine. sinh(x) = (e^x - e^(-x)) / 2.
 *       Inverse of asinh.
 * @example sinh(0) = 0
 * @example sinh(1) = 1.1752
 * @example sinh(-1) = -1.1752
 * @example sinh(2) = 3.6269
 * @see_also cosh, tanh, asinh
 */

/**
 * @func cosh
 * @category Trigonometry
 * @syntax cosh(x)
 * @desc Hyperbolic cosine. cosh(x) = (e^x + e^(-x)) / 2.
 *       Always >= 1 for real x. Inverse of acosh.
 * @example cosh(0) = 1
 * @example cosh(1) = 1.5431
 * @example cosh(-1) = 1.5431
 * @example cosh(2) = 3.7622
 * @see_also sinh, tanh, acosh
 */

/**
 * @func tanh
 * @category Trigonometry
 * @syntax tanh(x)
 * @desc Hyperbolic tangent. tanh(x) = sinh(x)/cosh(x).
 *       Range: (-1, 1). Approaches +/-1 as x -> +/-inf.
 * @example tanh(0) = 0
 * @example tanh(1) = 0.7616
 * @example tanh(-1) = -0.7616
 * @example tanh(10) = 1
 * @see_also sinh, cosh, atanh
 */

/**
 * @func asinh
 * @category Trigonometry
 * @syntax asinh(x)
 * @desc Inverse hyperbolic sine. asinh(x) = ln(x + sqrt(x^2 + 1)).
 *       Domain: all real numbers.
 * @example asinh(0) = 0
 * @example asinh(1) = 0.8814
 * @example asinh(-1) = -0.8814
 * @example asinh(sinh(2)) = 2
 * @see_also sinh, acosh, atanh
 */

/**
 * @func acosh
 * @category Trigonometry
 * @syntax acosh(x)
 * @desc Inverse hyperbolic cosine. acosh(x) = ln(x + sqrt(x^2 - 1)).
 *       Domain: [1, inf). Returns complex for x < 1.
 * @example acosh(1) = 0
 * @example acosh(2) = 1.3170
 * @example acosh(10) = 2.9932
 * @example acosh(cosh(3)) = 3
 * @see_also cosh, asinh, atanh
 */

/**
 * @func atanh
 * @category Trigonometry
 * @syntax atanh(x)
 * @desc Inverse hyperbolic tangent. atanh(x) = 0.5 * ln((1+x)/(1-x)).
 *       Domain: (-1, 1). Returns complex for |x| > 1.
 * @example atanh(0) = 0
 * @example atanh(0.5) = 0.5493
 * @example atanh(-0.5) = -0.5493
 * @example atanh(tanh(1)) = 1
 * @see_also tanh, asinh, acosh
 */

/**
 * @func sec
 * @category Trigonometry
 * @syntax sec(x)
 * @desc Secant. sec(x) = 1/cos(x).
 *       Undefined where cos(x) = 0.
 * @example sec(0) = 1
 * @example sec(pi/3) = 2
 * @example sec(pi/4) = 1.4142
 * @see_also cos, csc, cot, asec
 */

/**
 * @func csc
 * @category Trigonometry
 * @syntax csc(x)
 * @desc Cosecant. csc(x) = 1/sin(x).
 *       Undefined where sin(x) = 0.
 * @example csc(pi/2) = 1
 * @example csc(pi/6) = 2
 * @example csc(pi/4) = 1.4142
 * @see_also sin, sec, cot, acsc
 */

/**
 * @func cot
 * @category Trigonometry
 * @syntax cot(x)
 * @desc Cotangent. cot(x) = 1/tan(x) = cos(x)/sin(x).
 *       Undefined where sin(x) = 0.
 * @example cot(pi/4) = 1
 * @example cot(pi/2) = 0
 * @example cot(pi/6) = 1.7321
 * @see_also tan, sec, csc, acot
 */

/**
 * @func asec
 * @category Trigonometry
 * @syntax asec(x)
 * @desc Inverse secant. asec(x) = acos(1/x).
 *       Domain: |x| >= 1.
 * @example asec(1) = 0
 * @example asec(2) = 1.0472
 * @example asec(-1) = 3.1416
 * @see_also sec, acsc, acot
 */

/**
 * @func acsc
 * @category Trigonometry
 * @syntax acsc(x)
 * @desc Inverse cosecant. acsc(x) = asin(1/x).
 *       Domain: |x| >= 1.
 * @example acsc(1) = 1.5708
 * @example acsc(2) = 0.5236
 * @example acsc(-1) = -1.5708
 * @see_also csc, asec, acot
 */

/**
 * @func acot
 * @category Trigonometry
 * @syntax acot(x)
 * @desc Inverse cotangent. acot(x) = atan(1/x).
 *       For x = 0, returns pi/2.
 * @example acot(1) = 0.7854
 * @example acot(0) = 1.5708
 * @example acot(-1) = -0.7854
 * @see_also cot, asec, acsc
 */

/**
 * @func sech
 * @category Trigonometry
 * @syntax sech(x)
 * @desc Hyperbolic secant. sech(x) = 1/cosh(x).
 *       Maximum value is 1 at x = 0.
 * @example sech(0) = 1
 * @example sech(1) = 0.6481
 * @example sech(2) = 0.2658
 * @see_also cosh, csch, coth
 */

/**
 * @func csch
 * @category Trigonometry
 * @syntax csch(x)
 * @desc Hyperbolic cosecant. csch(x) = 1/sinh(x).
 *       Undefined at x = 0.
 * @example csch(1) = 0.8509
 * @example csch(2) = 0.2757
 * @example csch(-1) = -0.8509
 * @see_also sinh, sech, coth
 */

/**
 * @func coth
 * @category Trigonometry
 * @syntax coth(x)
 * @desc Hyperbolic cotangent. coth(x) = 1/tanh(x) = cosh(x)/sinh(x).
 *       Undefined at x = 0.
 * @example coth(1) = 1.3130
 * @example coth(2) = 1.0373
 * @example coth(-1) = -1.3130
 * @see_also tanh, sech, csch
 */

/**
 * @func asech
 * @category Trigonometry
 * @syntax asech(x)
 * @desc Inverse hyperbolic secant. asech(x) = acosh(1/x).
 *       Domain: (0, 1].
 * @example asech(1) = 0
 * @example asech(0.5) = 1.3170
 * @example asech(0.1) = 2.9932
 * @see_also sech, acsch, acoth
 */

/**
 * @func acsch
 * @category Trigonometry
 * @syntax acsch(x)
 * @desc Inverse hyperbolic cosecant. acsch(x) = asinh(1/x).
 *       Domain: x != 0.
 * @example acsch(1) = 0.8814
 * @example acsch(2) = 0.4812
 * @example acsch(-1) = -0.8814
 * @see_also csch, asech, acoth
 */

/**
 * @func acoth
 * @category Trigonometry
 * @syntax acoth(x)
 * @desc Inverse hyperbolic cotangent. acoth(x) = atanh(1/x).
 *       Domain: |x| > 1.
 * @example acoth(2) = 0.5493
 * @example acoth(10) = 0.1003
 * @example acoth(-2) = -0.5493
 * @see_also coth, asech, acsch
 */

/**
 * @func sind
 * @category Trigonometry
 * @syntax sind(x)
 * @desc Sine of x degrees. More accurate than sin(x*pi/180)
 *       for multiples of 90 degrees.
 * @example sind(0) = 0
 * @example sind(30) = 0.5
 * @example sind(90) = 1
 * @example sind(180) = 0
 * @example sind(270) = -1
 * @see_also sin, cosd, tand
 */

/**
 * @func cosd
 * @category Trigonometry
 * @syntax cosd(x)
 * @desc Cosine of x degrees. More accurate than cos(x*pi/180)
 *       for multiples of 90 degrees.
 * @example cosd(0) = 1
 * @example cosd(60) = 0.5
 * @example cosd(90) = 0
 * @example cosd(180) = -1
 * @example cosd(360) = 1
 * @see_also cos, sind, tand
 */

/**
 * @func tand
 * @category Trigonometry
 * @syntax tand(x)
 * @desc Tangent of x degrees.
 * @example tand(0) = 0
 * @example tand(45) = 1
 * @example tand(60) = 1.7321
 * @example tand(-45) = -1
 * @see_also tan, sind, cosd
 */

/**
 * @func asind
 * @category Trigonometry
 * @syntax asind(x)
 * @desc Inverse sine, result in degrees.
 *       Domain: [-1, 1]. Range: [-90, 90].
 * @example asind(0) = 0
 * @example asind(0.5) = 30
 * @example asind(1) = 90
 * @example asind(-1) = -90
 * @see_also asin, acosd, atand
 */

/**
 * @func acosd
 * @category Trigonometry
 * @syntax acosd(x)
 * @desc Inverse cosine, result in degrees.
 *       Domain: [-1, 1]. Range: [0, 180].
 * @example acosd(1) = 0
 * @example acosd(0.5) = 60
 * @example acosd(0) = 90
 * @example acosd(-1) = 180
 * @see_also acos, asind, atand
 */

/**
 * @func atand
 * @category Trigonometry
 * @syntax atand(x)
 * @desc Inverse tangent, result in degrees.
 *       Range: (-90, 90).
 * @example atand(0) = 0
 * @example atand(1) = 45
 * @example atand(-1) = -45
 * @example atand(1e10) = 90
 * @see_also atan, asind, acosd, atan2d
 */

/**
 * @func atan2d
 * @category Trigonometry
 * @syntax atan2d(y, x)
 * @desc Two-argument arctangent, result in degrees.
 *       Range: (-180, 180].
 * @example atan2d(0, 1) = 0
 * @example atan2d(1, 1) = 45
 * @example atan2d(1, 0) = 90
 * @example atan2d(0, -1) = 180
 * @example atan2d(-1, 0) = -90
 * @see_also atan2, atand
 */

/**
 * @func sinc
 * @category Trigonometry
 * @syntax sinc(x)
 * @desc Normalized sinc function. sinc(x) = sin(pi*x) / (pi*x).
 *       sinc(0) = 1 by definition (limit).
 *       Used in signal processing and Fourier analysis.
 * @example sinc(0) = 1
 * @example sinc(1) = 0
 * @example sinc(0.5) = 0.6366
 * @example sinc(-1) = 0
 * @see_also sin
 */

/**
 * @func sinpi
 * @category Trigonometry
 * @syntax sinpi(x)
 * @desc sin(pi * x). More accurate than sin(pi*x) for integer x.
 *       sinpi(n) = 0 for integer n.
 * @example sinpi(0) = 0
 * @example sinpi(0.5) = 1
 * @example sinpi(1) = 0
 * @example sinpi(1.5) = -1
 * @see_also cospi, sin
 */

/**
 * @func cospi
 * @category Trigonometry
 * @syntax cospi(x)
 * @desc cos(pi * x). More accurate than cos(pi*x) for half-integer x.
 *       cospi(n + 0.5) = 0 for integer n.
 * @example cospi(0) = 1
 * @example cospi(0.5) = 0
 * @example cospi(1) = -1
 * @example cospi(2) = 1
 * @see_also sinpi, cos
 */

/* Function documentation table */
const FuncDoc func_trig_docs[] = {
    {"sin", CAT_TRIG, "sin(x)", 
     "Sine of x (radians). For complex x, sin(x) = (e^(ix) - e^(-ix)) / (2i).",
     {"sin(0) = 0", "sin(pi/6) = 0.5", "sin(pi/2) = 1", "sin(pi) = 0", NULL, NULL},
     "cos, tan, asin, sind"},
    
    {"cos", CAT_TRIG, "cos(x)",
     "Cosine of x (radians). For complex x, cos(x) = (e^(ix) + e^(-ix)) / 2.",
     {"cos(0) = 1", "cos(pi/3) = 0.5", "cos(pi/2) = 0", "cos(pi) = -1", NULL, NULL},
     "sin, tan, acos, cosd"},
    
    {"tan", CAT_TRIG, "tan(x)",
     "Tangent of x. Equal to sin(x)/cos(x). Undefined at x = pi/2 + n*pi.",
     {"tan(0) = 0", "tan(pi/4) = 1", "tan(pi/6) = 0.5774", NULL, NULL, NULL},
     "sin, cos, atan, tand"},
    
    {"asin", CAT_TRIG, "asin(x)",
     "Inverse sine. Domain: [-1,1]. Range: [-pi/2, pi/2]. Complex for |x|>1.",
     {"asin(0) = 0", "asin(0.5) = 0.5236", "asin(1) = 1.5708", "asin(-1) = -1.5708", NULL, NULL},
     "sin, acos, atan, asind"},
    
    {"acos", CAT_TRIG, "acos(x)",
     "Inverse cosine. Domain: [-1,1]. Range: [0, pi]. Complex for |x|>1.",
     {"acos(1) = 0", "acos(0.5) = 1.0472", "acos(0) = 1.5708", "acos(-1) = 3.1416", NULL, NULL},
     "cos, asin, atan, acosd"},
    
    {"atan", CAT_TRIG, "atan(x)",
     "Inverse tangent. Domain: all reals. Range: (-pi/2, pi/2).",
     {"atan(0) = 0", "atan(1) = 0.7854", "atan(-1) = -0.7854", "4*atan(1) = 3.1416", NULL, NULL},
     "tan, asin, acos, atan2, atand"},
    
    {"atan2", CAT_TRIG, "atan2(y, x)",
     "Two-argument arctangent. Returns angle of (x,y) from origin. Range: (-pi, pi].",
     {"atan2(0, 1) = 0", "atan2(1, 1) = 0.7854", "atan2(1, 0) = 1.5708", "atan2(-1, -1) = -2.3562", NULL, NULL},
     "atan, atan2d"},
    
    {"sinh", CAT_TRIG, "sinh(x)",
     "Hyperbolic sine. sinh(x) = (e^x - e^(-x)) / 2.",
     {"sinh(0) = 0", "sinh(1) = 1.1752", "sinh(-1) = -1.1752", "sinh(2) = 3.6269", NULL, NULL},
     "cosh, tanh, asinh"},
    
    {"cosh", CAT_TRIG, "cosh(x)",
     "Hyperbolic cosine. cosh(x) = (e^x + e^(-x)) / 2. Always >= 1 for real x.",
     {"cosh(0) = 1", "cosh(1) = 1.5431", "cosh(-1) = 1.5431", "cosh(2) = 3.7622", NULL, NULL},
     "sinh, tanh, acosh"},
    
    {"tanh", CAT_TRIG, "tanh(x)",
     "Hyperbolic tangent. tanh(x) = sinh(x)/cosh(x). Range: (-1, 1).",
     {"tanh(0) = 0", "tanh(1) = 0.7616", "tanh(-1) = -0.7616", "tanh(10) = 1", NULL, NULL},
     "sinh, cosh, atanh"},
    
    {"asinh", CAT_TRIG, "asinh(x)",
     "Inverse hyperbolic sine. asinh(x) = ln(x + sqrt(x^2 + 1)).",
     {"asinh(0) = 0", "asinh(1) = 0.8814", "asinh(-1) = -0.8814", NULL, NULL, NULL},
     "sinh, acosh, atanh"},
    
    {"acosh", CAT_TRIG, "acosh(x)",
     "Inverse hyperbolic cosine. Domain: [1, inf). Complex for x < 1.",
     {"acosh(1) = 0", "acosh(2) = 1.3170", "acosh(10) = 2.9932", NULL, NULL, NULL},
     "cosh, asinh, atanh"},
    
    {"atanh", CAT_TRIG, "atanh(x)",
     "Inverse hyperbolic tangent. Domain: (-1, 1). Complex for |x| > 1.",
     {"atanh(0) = 0", "atanh(0.5) = 0.5493", "atanh(-0.5) = -0.5493", NULL, NULL, NULL},
     "tanh, asinh, acosh"},
    
    {"sec", CAT_TRIG, "sec(x)",
     "Secant. sec(x) = 1/cos(x).",
     {"sec(0) = 1", "sec(pi/3) = 2", "sec(pi/4) = 1.4142", NULL, NULL, NULL},
     "cos, csc, cot, asec"},
    
    {"csc", CAT_TRIG, "csc(x)",
     "Cosecant. csc(x) = 1/sin(x).",
     {"csc(pi/2) = 1", "csc(pi/6) = 2", "csc(pi/4) = 1.4142", NULL, NULL, NULL},
     "sin, sec, cot, acsc"},
    
    {"cot", CAT_TRIG, "cot(x)",
     "Cotangent. cot(x) = 1/tan(x) = cos(x)/sin(x).",
     {"cot(pi/4) = 1", "cot(pi/2) = 0", "cot(pi/6) = 1.7321", NULL, NULL, NULL},
     "tan, sec, csc, acot"},
    
    {"asec", CAT_TRIG, "asec(x)",
     "Inverse secant. asec(x) = acos(1/x). Domain: |x| >= 1.",
     {"asec(1) = 0", "asec(2) = 1.0472", "asec(-1) = 3.1416", NULL, NULL, NULL},
     "sec, acsc, acot"},
    
    {"acsc", CAT_TRIG, "acsc(x)",
     "Inverse cosecant. acsc(x) = asin(1/x). Domain: |x| >= 1.",
     {"acsc(1) = 1.5708", "acsc(2) = 0.5236", "acsc(-1) = -1.5708", NULL, NULL, NULL},
     "csc, asec, acot"},
    
    {"acot", CAT_TRIG, "acot(x)",
     "Inverse cotangent. acot(x) = atan(1/x). acot(0) = pi/2.",
     {"acot(1) = 0.7854", "acot(0) = 1.5708", "acot(-1) = -0.7854", NULL, NULL, NULL},
     "cot, asec, acsc"},
    
    {"sech", CAT_TRIG, "sech(x)",
     "Hyperbolic secant. sech(x) = 1/cosh(x). Maximum 1 at x=0.",
     {"sech(0) = 1", "sech(1) = 0.6481", "sech(2) = 0.2658", NULL, NULL, NULL},
     "cosh, csch, coth"},
    
    {"csch", CAT_TRIG, "csch(x)",
     "Hyperbolic cosecant. csch(x) = 1/sinh(x). Undefined at x=0.",
     {"csch(1) = 0.8509", "csch(2) = 0.2757", "csch(-1) = -0.8509", NULL, NULL, NULL},
     "sinh, sech, coth"},
    
    {"coth", CAT_TRIG, "coth(x)",
     "Hyperbolic cotangent. coth(x) = 1/tanh(x). Undefined at x=0.",
     {"coth(1) = 1.3130", "coth(2) = 1.0373", "coth(-1) = -1.3130", NULL, NULL, NULL},
     "tanh, sech, csch"},
    
    {"asech", CAT_TRIG, "asech(x)",
     "Inverse hyperbolic secant. Domain: (0, 1].",
     {"asech(1) = 0", "asech(0.5) = 1.3170", "asech(0.1) = 2.9932", NULL, NULL, NULL},
     "sech, acsch, acoth"},
    
    {"acsch", CAT_TRIG, "acsch(x)",
     "Inverse hyperbolic cosecant. Domain: x != 0.",
     {"acsch(1) = 0.8814", "acsch(2) = 0.4812", "acsch(-1) = -0.8814", NULL, NULL, NULL},
     "csch, asech, acoth"},
    
    {"acoth", CAT_TRIG, "acoth(x)",
     "Inverse hyperbolic cotangent. Domain: |x| > 1.",
     {"acoth(2) = 0.5493", "acoth(10) = 0.1003", "acoth(-2) = -0.5493", NULL, NULL, NULL},
     "coth, asech, acsch"},
    
    {"sind", CAT_TRIG, "sind(x)",
     "Sine of x degrees. Exact for multiples of 90.",
     {"sind(0) = 0", "sind(30) = 0.5", "sind(90) = 1", "sind(180) = 0", "sind(270) = -1", NULL},
     "sin, cosd, tand"},
    
    {"cosd", CAT_TRIG, "cosd(x)",
     "Cosine of x degrees. Exact for multiples of 90.",
     {"cosd(0) = 1", "cosd(60) = 0.5", "cosd(90) = 0", "cosd(180) = -1", NULL, NULL},
     "cos, sind, tand"},
    
    {"tand", CAT_TRIG, "tand(x)",
     "Tangent of x degrees.",
     {"tand(0) = 0", "tand(45) = 1", "tand(60) = 1.7321", "tand(-45) = -1", NULL, NULL},
     "tan, sind, cosd"},
    
    {"asind", CAT_TRIG, "asind(x)",
     "Inverse sine in degrees. Range: [-90, 90].",
     {"asind(0) = 0", "asind(0.5) = 30", "asind(1) = 90", "asind(-1) = -90", NULL, NULL},
     "asin, acosd, atand"},
    
    {"acosd", CAT_TRIG, "acosd(x)",
     "Inverse cosine in degrees. Range: [0, 180].",
     {"acosd(1) = 0", "acosd(0.5) = 60", "acosd(0) = 90", "acosd(-1) = 180", NULL, NULL},
     "acos, asind, atand"},
    
    {"atand", CAT_TRIG, "atand(x)",
     "Inverse tangent in degrees. Range: (-90, 90).",
     {"atand(0) = 0", "atand(1) = 45", "atand(-1) = -45", NULL, NULL, NULL},
     "atan, asind, acosd, atan2d"},
    
    {"atan2d", CAT_TRIG, "atan2d(y, x)",
     "Two-argument arctangent in degrees. Range: (-180, 180].",
     {"atan2d(0, 1) = 0", "atan2d(1, 1) = 45", "atan2d(1, 0) = 90", "atan2d(-1, 0) = -90", NULL, NULL},
     "atan2, atand"},
    
    {"sinc", CAT_TRIG, "sinc(x)",
     "Normalized sinc: sin(pi*x)/(pi*x). sinc(0) = 1.",
     {"sinc(0) = 1", "sinc(1) = 0", "sinc(0.5) = 0.6366", "sinc(-1) = 0", NULL, NULL},
     "sin"},
    
    {"sinpi", CAT_TRIG, "sinpi(x)",
     "sin(pi*x). Exact for integers.",
     {"sinpi(0) = 0", "sinpi(0.5) = 1", "sinpi(1) = 0", "sinpi(1.5) = -1", NULL, NULL},
     "cospi, sin"},
    
    {"cospi", CAT_TRIG, "cospi(x)",
     "cos(pi*x). Exact for half-integers.",
     {"cospi(0) = 1", "cospi(0.5) = 0", "cospi(1) = -1", "cospi(2) = 1", NULL, NULL},
     "sinpi, cos"},
    
    /* Sentinel */
    {NULL, NULL, NULL, NULL, {NULL, NULL, NULL, NULL, NULL, NULL}, NULL}
};
