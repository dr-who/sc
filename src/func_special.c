/*
 * func_special.c - Special mathematical function documentation
 * C89 compliant
 */
#include "func_defs.h"

const FuncDoc func_special_docs[] = {
    {"gamma", CAT_SPECIAL, "gamma(x)",
     "Gamma function. gamma(n) = (n-1)! for positive integers. Extends factorial to reals.",
     {"gamma(5) = 24", "gamma(1) = 1", "gamma(0.5) = 1.7725", "gamma(6) = 120", "gamma(1.5) = 0.8862", NULL},
     "lgamma, fact, beta"},
    
    {"lgamma", CAT_SPECIAL, "lgamma(x)",
     "Log-gamma function. ln(|gamma(x)|). More stable for large x.",
     {"lgamma(5) = 3.1781", "lgamma(100) = 359.13", "lgamma(1) = 0", "lgamma(10) = 12.802", NULL, NULL},
     "gamma"},
    
    {"tgamma", CAT_SPECIAL, "tgamma(x)",
     "True gamma function. Alias for gamma.",
     {"tgamma(5) = 24", "tgamma(0.5) = 1.7725", "tgamma(3) = 2", NULL, NULL, NULL},
     "gamma"},
    
    {"beta", CAT_SPECIAL, "beta(a, b)",
     "Beta function. beta(a,b) = gamma(a)*gamma(b)/gamma(a+b).",
     {"beta(2, 3) = 0.0833", "beta(1, 1) = 1", "beta(0.5, 0.5) = 3.1416", "beta(3, 4) = 0.0167", NULL, NULL},
     "gamma, betainc"},
    
    {"digamma", CAT_SPECIAL, "digamma(x)",
     "Digamma function. psi(x) = d/dx ln(gamma(x)).",
     {"digamma(1) = -0.5772", "digamma(2) = 0.4228", "digamma(5) = 1.5061", "digamma(10) = 2.2518", NULL, NULL},
     "gamma, psi"},
    
    {"psi", CAT_SPECIAL, "psi(x)",
     "Digamma function. Alias for digamma.",
     {"psi(1) = -0.5772", "psi(2) = 0.4228", "psi(3) = 0.9228", NULL, NULL, NULL},
     "digamma"},
    
    {"erf", CAT_SPECIAL, "erf(x)",
     "Error function. erf(x) = (2/sqrt(pi)) * integral from 0 to x of exp(-t^2) dt.",
     {"erf(0) = 0", "erf(1) = 0.8427", "erf(2) = 0.9953", "erf(-1) = -0.8427", "erf(0.5) = 0.5205", NULL},
     "erfc, normcdf"},
    
    {"erfc", CAT_SPECIAL, "erfc(x)",
     "Complementary error function. erfc(x) = 1 - erf(x).",
     {"erfc(0) = 1", "erfc(1) = 0.1573", "erfc(2) = 0.0047", "erfc(-1) = 1.8427", NULL, NULL},
     "erf"},
    
    {"besselj", CAT_SPECIAL, "besselj(nu, x)",
     "Bessel function of first kind J_nu(x).",
     {"besselj(0, 0) = 1", "besselj(0, 1) = 0.7652", "besselj(1, 1) = 0.4401", "besselj(0, 2.4) = 0", NULL, NULL},
     "bessely"},
    
    {"bessely", CAT_SPECIAL, "bessely(nu, x)",
     "Bessel function of second kind Y_nu(x).",
     {"bessely(0, 1) = 0.0883", "bessely(1, 1) = -0.7812", "bessely(0, 2) = 0.5104", NULL, NULL, NULL},
     "besselj"},
    
    {"zeta", CAT_SPECIAL, "zeta(s)",
     "Riemann zeta function. zeta(s) = sum of n^(-s) for n=1 to infinity.",
     {"zeta(2) = 1.6449", "zeta(3) = 1.2021", "zeta(4) = 1.0823", "zeta(6) = 1.0173", NULL, NULL},
     "harmonic"},
    
    {"harmonic", CAT_SPECIAL, "harmonic(n)",
     "Harmonic number H_n = 1 + 1/2 + 1/3 + ... + 1/n.",
     {"harmonic(1) = 1", "harmonic(2) = 1.5", "harmonic(10) = 2.9290", "harmonic(100) = 5.1874", NULL, NULL},
     "zeta"},
    
    {"sigmoid", CAT_SPECIAL, "sigmoid(x)",
     "Logistic sigmoid. sigmoid(x) = 1 / (1 + exp(-x)).",
     {"sigmoid(0) = 0.5", "sigmoid(1) = 0.7311", "sigmoid(-1) = 0.2689", "sigmoid(10) = 1", "sigmoid(-10) = 0", NULL},
     "softplus"},
    
    {"softplus", CAT_SPECIAL, "softplus(x)",
     "Softplus function. softplus(x) = ln(1 + exp(x)). Smooth approximation to ReLU.",
     {"softplus(0) = 0.6931", "softplus(1) = 1.3133", "softplus(-1) = 0.3133", "softplus(10) = 10", NULL, NULL},
     "sigmoid"},
    
    {"heaviside", CAT_SPECIAL, "heaviside(x)",
     "Heaviside step function. Returns 0 for x<0, 0.5 for x=0, 1 for x>0.",
     {"heaviside(-1) = 0", "heaviside(0) = 0.5", "heaviside(1) = 1", "heaviside(5) = 1", NULL, NULL},
     "step, rect"},
    
    {"step", CAT_SPECIAL, "step(x)",
     "Unit step function. Returns 0 for x<0, 1 for x>=0.",
     {"step(-1) = 0", "step(0) = 1", "step(1) = 1", "step(0.001) = 1", NULL, NULL},
     "heaviside, rect"},
    
    {"rect", CAT_SPECIAL, "rect(x)",
     "Rectangle function. Returns 1 for |x| < 0.5, 0.5 for |x| = 0.5, 0 otherwise.",
     {"rect(0) = 1", "rect(0.3) = 1", "rect(0.5) = 0.5", "rect(1) = 0", "rect(-0.25) = 1", NULL},
     "tri, step"},
    
    {"tri", CAT_SPECIAL, "tri(x)",
     "Triangle function. Returns 1-|x| for |x| < 1, 0 otherwise.",
     {"tri(0) = 1", "tri(0.5) = 0.5", "tri(1) = 0", "tri(-0.5) = 0.5", "tri(2) = 0", NULL},
     "rect"},
    
    /* Sentinel */
    {NULL, NULL, NULL, NULL, {NULL, NULL, NULL, NULL, NULL, NULL}, NULL}
};
