/*
 * func_complex.c - Complex number function documentation
 * C89 compliant
 */
#include "func_defs.h"

const FuncDoc func_complex_docs[] = {
    {"re", CAT_COMPLEX, "re(z)",
     "Real part of complex number z.",
     {"re(3+4i) = 3", "re(5) = 5", "re(2i) = 0", "re(-3-4i) = -3", NULL, NULL},
     "im, conj, abs"},
    
    {"real", CAT_COMPLEX, "real(z)",
     "Real part of complex number. Alias for re.",
     {"real(3+4i) = 3", "real(5) = 5", "real(-2i) = 0", NULL, NULL, NULL},
     "re, imag"},
    
    {"im", CAT_COMPLEX, "im(z)",
     "Imaginary part of complex number z.",
     {"im(3+4i) = 4", "im(5) = 0", "im(2i) = 2", "im(-3-4i) = -4", NULL, NULL},
     "re, conj, abs"},
    
    {"imag", CAT_COMPLEX, "imag(z)",
     "Imaginary part of complex number. Alias for im.",
     {"imag(3+4i) = 4", "imag(5) = 0", "imag(-2i) = -2", NULL, NULL, NULL},
     "im, real"},
    
    {"conj", CAT_COMPLEX, "conj(z)",
     "Complex conjugate. conj(a + bi) = a - bi.",
     {"conj(3+4i) = 3-4i", "conj(5) = 5", "conj(2i) = -2i", "conj(-3-4i) = -3+4i", NULL, NULL},
     "re, im, abs"},
    
    {"abs", CAT_COMPLEX, "abs(z)",
     "Absolute value (magnitude). For complex: sqrt(re^2 + im^2).",
     {"abs(-5) = 5", "abs(3+4i) = 5", "abs(-3.14) = 3.14", "abs(5-12i) = 13", "abs(0) = 0", NULL},
     "arg, conj"},
    
    {"arg", CAT_COMPLEX, "arg(z)",
     "Phase angle (argument) in radians. Range: (-pi, pi].",
     {"arg(1) = 0", "arg(i) = 1.5708", "arg(-1) = 3.1416", "arg(1+i) = 0.7854", "arg(-i) = -1.5708", NULL},
     "abs, angle, phase"},
    
    {"angle", CAT_COMPLEX, "angle(z)",
     "Phase angle in radians. Alias for arg.",
     {"angle(1) = 0", "angle(i) = 1.5708", "angle(1+i) = 0.7854", NULL, NULL, NULL},
     "arg, phase"},
    
    {"phase", CAT_COMPLEX, "phase(z)",
     "Phase angle in radians. Alias for arg.",
     {"phase(1) = 0", "phase(-1) = 3.1416", "phase(3+4i) = 0.9273", NULL, NULL, NULL},
     "arg, angle"},
    
    {"complex", CAT_COMPLEX, "complex(re, im)",
     "Create complex number from real and imaginary parts.",
     {"complex(3, 4) = 3+4i", "complex(1, 0) = 1", "complex(0, 1) = i", "complex(-2, 5) = -2+5i", NULL, NULL},
     "re, im"},
    
    {"cart2pol", CAT_COMPLEX, "cart2pol(x, y)",
     "Convert Cartesian (x,y) to polar (r, theta). Returns [r, theta].",
     {"cart2pol(1, 0)", "cart2pol(0, 1)", "cart2pol(1, 1)", "cart2pol(3, 4)", NULL, NULL},
     "pol2cart, abs, arg"},
    
    {"pol2cart", CAT_COMPLEX, "pol2cart(r, theta)",
     "Convert polar (r, theta) to Cartesian (x, y). Returns [x, y].",
     {"pol2cart(1, 0)", "pol2cart(1, pi/2)", "pol2cart(5, 0.9273)", NULL, NULL, NULL},
     "cart2pol"},
    
    /* Sentinel */
    {NULL, NULL, NULL, NULL, {NULL, NULL, NULL, NULL, NULL, NULL}, NULL}
};
