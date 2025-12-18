/* apf_native.c - Native float conversion functions
 * This file is NOT part of the pure soft-float library.
 * Use this in applications that need to interface with native doubles.
 */

#include "apf.h"
#include <stdio.h>
#include <stdlib.h>

/* Convert from double - useful for interfacing with native math */
void apf_from_double(apf *r, double d)
{
    char buf[64];
    
    /* Check for special values first - using != comparison to detect NaN */
    if (d != d) {
        /* NaN: d != d is only true for NaN */
        apf_set_nan(r);
        return;
    }
    if (d == 0.0) {
        apf_zero(r);
        return;
    }
    /* Check for infinity */
    if (d > 1e308) {
        apf_set_inf(r, 0);  /* +Inf */
        return;
    }
    if (d < -1e308) {
        apf_set_inf(r, 1);  /* -Inf */
        return;
    }
    /* Use sprintf to convert, then parse */
    sprintf(buf, "%.17g", d);
    apf_from_str(r, buf);
}

/* Convert to double - loses precision but useful for plotting */
double apf_to_double(const apf *a)
{
    char buf[64];
    if (a->cls == APF_CLASS_ZERO) return 0.0;
    if (a->cls == APF_CLASS_INF) return a->sign ? -1e308 : 1e308;
    if (a->cls == APF_CLASS_NAN) return 0.0;  /* Can't represent NaN portably */
    
    apf_to_str(buf, sizeof(buf), a, 15);
    return atof(buf);
}

/* Debug dump to FILE* - for debugging only */
void apf_dump(const apf *x, void *fp)
{
    FILE *f = fp ? (FILE *)fp : stdout;
    int i;

    if (x->cls == APF_CLASS_NAN) {
        fprintf(f, "NaN\n");
        return;
    }
    if (x->cls == APF_CLASS_INF) {
        fprintf(f, "%sInf\n", x->sign ? "-" : "+");
        return;
    }
    if (x->cls == APF_CLASS_ZERO) {
        fprintf(f, "%s0\n", x->sign ? "-" : "+");
        return;
    }

    fprintf(f, "NORMAL sign=%d exp=%ld mant=0x", x->sign, x->exp);
    for (i = AP_LIMBS - 1; i >= 0; --i)
        fprintf(f, "%04X", x->mant[i] & 0xFFFFu);
    fprintf(f, "\n");
}
