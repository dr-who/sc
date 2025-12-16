/* apfx.h - Extended math functions for APF
 * Pure C89, no external dependencies
 */
#ifndef APFX_H
#define APFX_H

#include "apf.h"

/* Initialization - MUST be called before using any apfx functions on DOS */
void apfx_init(void);

/* Constants */
void apfx_pi(apf *r);
void apfx_e(apf *r);

/* Transcendental functions (Phase 1: Core) */
void apfx_exp(apf *r, const apf *x);
void apfx_log(apf *r, const apf *x);
void apfx_sin(apf *r, const apf *x);
void apfx_cos(apf *r, const apf *x);
void apfx_sincos(apf *sin_r, apf *cos_r, const apf *x);
void apfx_tan(apf *r, const apf *x);
void apfx_atan(apf *r, const apf *x);
void apfx_atan2(apf *r, const apf *y, const apf *x);
void apfx_asin(apf *r, const apf *x);
void apfx_acos(apf *r, const apf *x);
void apfx_sinh(apf *r, const apf *x);
void apfx_cosh(apf *r, const apf *x);
void apfx_tanh(apf *r, const apf *x);
void apfx_sinhcosh(apf *sinh_r, apf *cosh_r, const apf *x);
void apfx_asinh(apf *r, const apf *x);
void apfx_acosh(apf *r, const apf *x);
void apfx_atanh(apf *r, const apf *x);

/* Power functions */
void apfx_pow(apf *r, const apf *base, const apf *exp);

/* Factorial */
void apfx_fact(apf *r, long n);

/* Phase 2: Scientific Essentials */

/* Gamma functions */
void apfx_lgamma(apf *r, const apf *x);   /* Log-gamma: ln(|Gamma(x)|) */
void apfx_tgamma(apf *r, const apf *x);   /* Gamma function */

/* Error functions */
void apfx_erf(apf *r, const apf *x);      /* Error function */
void apfx_erfc(apf *r, const apf *x);     /* Complementary error function: 1-erf(x) */

/* Normal distribution */
void apfx_normcdf(apf *r, const apf *x);  /* Standard normal CDF: Phi(x) */
void apfx_norminv(apf *r, const apf *p);  /* Inverse normal CDF (probit) */

/* Student's t distribution */
void apfx_tcdf(apf *r, const apf *t, long df);    /* t-distribution CDF */
void apfx_tinv(apf *r, const apf *p, long df);    /* Inverse t-distribution */

/* Chi-square distribution */
void apfx_chi2cdf(apf *r, const apf *x, long df); /* Chi-square CDF */
void apfx_chi2inv(apf *r, const apf *p, long df); /* Inverse chi-square */

/* F distribution */
void apfx_fcdf(apf *r, const apf *x, long df1, long df2); /* F-distribution CDF */

/* Binomial distribution */
void apfx_binompdf(apf *r, long k, long n, const apf *p); /* Binomial PMF */
void apfx_binomcdf(apf *r, long k, long n, const apf *p); /* Binomial CDF */

/* Poisson distribution */
void apfx_poisspdf(apf *r, long k, const apf *lambda);    /* Poisson PMF */
void apfx_poisscdf(apf *r, long k, const apf *lambda);    /* Poisson CDF */

/* Bessel functions - First kind */
void apfx_j0(apf *r, const apf *x);       /* Bessel J_0(x) */
void apfx_j1(apf *r, const apf *x);       /* Bessel J_1(x) */

/* Bessel functions - Second kind (Neumann functions) */
void apfx_y0(apf *r, const apf *x);       /* Bessel Y_0(x) */
void apfx_y1(apf *r, const apf *x);       /* Bessel Y_1(x) */

/* Modified Bessel functions - First kind */
void apfx_i0(apf *r, const apf *x);       /* Modified Bessel I_0(x) */
void apfx_i1(apf *r, const apf *x);       /* Modified Bessel I_1(x) */

/* Modified Bessel functions - Second kind */
void apfx_k0(apf *r, const apf *x);       /* Modified Bessel K_0(x) */
void apfx_k1(apf *r, const apf *x);       /* Modified Bessel K_1(x) */

/* Elliptic integrals (AGM method) */
void apfx_ellipk(apf *r, const apf *m);   /* Complete elliptic K(m) */
void apfx_ellipe(apf *r, const apf *m);   /* Complete elliptic E(m) */

/* Phase 4: Advanced */

/* Lambert W function */
void apfx_lambertw(apf *r, const apf *x); /* Principal branch W_0(x) */

/* Gauss-Legendre quadrature helpers */
void apfx_gl_node(apf *r, int n, int i);  /* Get node i for n-point quadrature */
void apfx_gl_weight(apf *r, int n, int i); /* Get weight i for n-point quadrature */

/* Beta function and incomplete beta (for distributions) */
void apfx_beta(apf *r, const apf *a, const apf *b);     /* Beta(a,b) = Gamma(a)*Gamma(b)/Gamma(a+b) */
void apfx_betainc(apf *r, const apf *x, const apf *a, const apf *b); /* Regularized incomplete beta I_x(a,b) */

/* Incomplete gamma (for distributions) */
void apfx_gammainc(apf *r, const apf *a, const apf *x); /* Regularized lower incomplete gamma P(a,x) */

#endif /* APFX_H */
