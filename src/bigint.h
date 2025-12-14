/* bigint.h - Signed 128-bit integer arithmetic
 * Lean C89 for VIC-20 (cc65), DOS (Watcom), Linux (gcc)
 */
#ifndef BIGINT_H
#define BIGINT_H

#if defined(__CC65__) || (defined(__WATCOMC__) && defined(__I86__))
  typedef unsigned int  bi_u16;
  typedef unsigned long bi_u32;
#else
  typedef unsigned short bi_u16;
  typedef unsigned long  bi_u32;
#endif

/* Configurable: 8=128bit, 4=64bit, 2=32bit */
#ifndef BI_LIMBS
#define BI_LIMBS 8
#endif

typedef struct {
    bi_u16 m[BI_LIMBS];
    char sign;
} s128;

void bi_zero(s128 *r);
void bi_neg(s128 *r);
int  bi_is_zero(const s128 *a);
void bi_add(s128 *r, const s128 *a, const s128 *b);
void bi_sub(s128 *r, const s128 *a, const s128 *b);
void bi_mul(s128 *r, const s128 *a, const s128 *b);
void bi_div(s128 *q, s128 *rem, const s128 *a, const s128 *b);
void bi_parse(s128 *r, const char *s);
char *bi_fmt(char *end, s128 *a);

#endif
