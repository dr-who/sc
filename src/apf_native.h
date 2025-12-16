/* apf_native.h - Native float conversion functions
 * This is NOT part of the pure soft-float library.
 * Include this in applications that need double/float interop.
 */
#ifndef APF_NATIVE_H
#define APF_NATIVE_H

#include "apf.h"
#include <stdio.h>

/* Convert between apf and native double */
void apf_from_double(apf *r, double d);
double apf_to_double(const apf *a);

/* Debug dump to FILE* */
void apf_dump(const apf *x, void *fp);

#endif /* APF_NATIVE_H */
