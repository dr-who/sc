/* apfc.h - Arbitrary Precision Complex numbers */
#ifndef APFC_H
#define APFC_H

#include "apf.h"
#include "apfx.h"

typedef struct {
    apf re;
    apf im;
} apfc;

/* Basic operations */
void apfc_zero(apfc *z);
void apfc_from_real(apfc *z, const apf *r);
void apfc_copy(apfc *dest, const apfc *src);
int apfc_is_zero(const apfc *z);
int apfc_is_real(const apfc *z);

/* Arithmetic */
void apfc_add(apfc *r, const apfc *a, const apfc *b);
void apfc_sub(apfc *r, const apfc *a, const apfc *b);
void apfc_mul(apfc *r, const apfc *a, const apfc *b);
void apfc_div(apfc *r, const apfc *a, const apfc *b);
void apfc_neg(apfc *r, const apfc *a);
void apfc_conj(apfc *r, const apfc *a);
void apfc_abs(apf *r, const apfc *z);
void apfc_arg(apf *r, const apfc *z);

/* Functions */
void apfc_sqrt(apfc *r, const apfc *z);
void apfc_exp(apfc *r, const apfc *z);
void apfc_log(apfc *r, const apfc *z);
void apfc_pow(apfc *r, const apfc *base, const apfc *exp);
void apfc_sin(apfc *r, const apfc *z);
void apfc_cos(apfc *r, const apfc *z);
void apfc_tan(apfc *r, const apfc *z);
void apfc_sinh(apfc *r, const apfc *z);
void apfc_cosh(apfc *r, const apfc *z);

/* String conversion */
void apfc_to_str(char *buf, int bufsize, const apfc *z, int max_frac);

#endif /* APFC_H */
