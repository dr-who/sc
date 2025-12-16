/* apf.h - Arbitrary Precision Float (soft float)
 * Pure C89, standalone library - no external dependencies
 * Uses only 16/32-bit integer arithmetic
 * No float, double, or stdio
 */
#ifndef APF_H
#define APF_H

/*
 * Portable integer types for 16-bit DOS compatibility.
 * On DOS (Watcom C): int=16-bit, long=32-bit
 * On modern systems: int=32-bit, long=32/64-bit
 */
typedef int  INT16;   /* 16-bit on DOS, may be larger elsewhere */
typedef long INT32;   /* 32-bit everywhere */
typedef unsigned int  UINT16;
typedef unsigned long UINT32;

/* Configuration - number of 16-bit limbs */
#ifndef AP_LIMBS
#define AP_LIMBS 8   /* Default: 8 x 16-bit = 128 bits */
#endif

#define AP_BITS (AP_LIMBS * 16)

/* Types: use only 16-bit and 32-bit integers */
typedef unsigned int limb_t;     /* 16-bit limb */
typedef unsigned long wide_t;    /* 32-bit for intermediate results */

/* Classification */
#define APF_CLASS_ZERO   0
#define APF_CLASS_NORMAL 1
#define APF_CLASS_INF    2
#define APF_CLASS_NAN    3

/* APF structure */
typedef struct {
    limb_t mant[AP_LIMBS];   /* mantissa, mant[0] = LSW */
    long   exp;              /* binary exponent */
    int    sign;             /* 0 = positive, 1 = negative */
    int    cls;              /* APF_CLASS_xxx */
} apf;

/* Basic operations */
void apf_zero(apf *x);
void apf_set_nan(apf *x);
void apf_set_inf(apf *x, int sign);
void apf_from_int(apf *x, long v);
long apf_to_long(const apf *x);
void apf_copy(apf *dst, const apf *src);

/* Predicates */
int apf_isnan(const apf *x);
int apf_isinf(const apf *x);
int apf_iszero(const apf *x);
int apf_is_zero(const apf *x);
int apf_is_neg(const apf *x);

/* Special values */
void apf_set_max(apf *x);      /* largest finite positive */
void apf_set_min(apf *x);      /* smallest positive non-zero */

/* Limits - exponent range for 32-bit long */
#define APF_EXP_MAX  2147483647L
#define APF_EXP_MIN  (-2147483647L - 1L)

/* Normalization */
void apf_norm(apf *x);

/* Comparison */
int apf_cmp(const apf *a, const apf *b);
int apf_cmp_int(const apf *a, long n);
int apf_eq(const apf *a, const apf *b);
int apf_ne(const apf *a, const apf *b);
int apf_lt(const apf *a, const apf *b);
int apf_le(const apf *a, const apf *b);
int apf_gt(const apf *a, const apf *b);
int apf_ge(const apf *a, const apf *b);

/* Arithmetic */
void apf_neg(apf *r, const apf *a);
void apf_abs(apf *r, const apf *a);
void apf_add(apf *r, const apf *a, const apf *b);
void apf_sub(apf *r, const apf *a, const apf *b);
void apf_mul(apf *r, const apf *a, const apf *b);
void apf_div(apf *r, const apf *a, const apf *b);
void apf_mod(apf *r, const apf *a, const apf *b);
void apf_sqrt(apf *r, const apf *a);

/* Bitwise operations */
void apf_lsl(apf *r, const apf *a, int shift);
void apf_lsr(apf *r, const apf *a, int shift);
void apf_band(apf *r, const apf *a, const apf *b);
void apf_bor(apf *r, const apf *a, const apf *b);
void apf_bxor(apf *r, const apf *a, const apf *b);
void apf_bnot(apf *r, const apf *a);

/* Rounding */
void apf_trunc(apf *r, const apf *a);
void apf_floor(apf *r, const apf *a);
void apf_ceil(apf *r, const apf *a);

/* Random */
void apf_random(apf *r);

/* Output rounding mode for apf_to_str */
#define APF_ROUND_NEAREST 0   /* Round to nearest, ties to even (default) */
#define APF_ROUND_CEIL    1   /* Round toward +infinity */
#define APF_ROUND_FLOOR   2   /* Round toward -infinity */
#define APF_ROUND_TRUNC   3   /* Round toward zero */
#define APF_ROUND_AWAY    4   /* Round to nearest, ties away from zero */

void apf_set_round_mode(int mode);
int apf_get_round_mode(void);

/* String conversion */
char *apf_to_str(char *buf, int bufsize, const apf *val, int max_digits);
void apf_from_str(apf *x, const char *s);

/* Debug - dump internal representation to string buffer */
char *apf_dump_str(char *buf, int bufsize, const apf *x);

#endif /* APF_H */
