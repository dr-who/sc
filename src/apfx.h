/* apfx.h - Extended math functions for APF */
#ifndef APFX_H
#define APFX_H

#include "apf.h"

/* Constants */
void apfx_pi(apf *r);
void apfx_e(apf *r);

/* Transcendental functions */
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

/* Power functions */
void apfx_pow(apf *r, const apf *base, const apf *exp);

/* Factorial */
void apfx_fact(apf *r, long n);

#endif /* APFX_H */
