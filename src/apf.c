/* apf.c - arbitrary-precision bigfloat core using only 16/32-bit arithmetic
 * C89 compliant for Watcom C / DOS
 */

#include "apf.h"
#include <stdio.h>

/* Global output rounding mode */
static int apf_round_mode = APF_ROUND_NEAREST;

void apf_set_round_mode(int mode)
{
    if (mode >= APF_ROUND_NEAREST && mode <= APF_ROUND_AWAY) {
        apf_round_mode = mode;
    }
}

int apf_get_round_mode(void)
{
    return apf_round_mode;
}

/* ---------- limb helpers (AP_LIMBS x 16-bit) ---------- */

static void limb_zero(limb_t *a)
{
    int i;
    for (i = 0; i < AP_LIMBS; ++i) a[i] = 0;
}

static void limb_copy(limb_t *dst, const limb_t *src)
{
    int i;
    for (i = 0; i < AP_LIMBS; ++i) dst[i] = src[i];
}

static int limb_is_zero(const limb_t *a)
{
    int i;
    for (i = 0; i < AP_LIMBS; ++i)
        if (a[i] != 0) return 0;
    return 1;
}

/* compare mantissa |a| ? |b|, return -1,0,1 */
static int limb_cmp(const limb_t *a, const limb_t *b)
{
    int i;
    for (i = AP_LIMBS - 1; i >= 0; --i) {
        if (a[i] < b[i]) return -1;
        if (a[i] > b[i]) return 1;
    }
    return 0;
}

/* r = a + b (AP_BITS-bit add), return carry (0/1) */
static limb_t limb_add(limb_t *r, const limb_t *a, const limb_t *b)
{
    wide_t carry = 0;
    int i;
    for (i = 0; i < AP_LIMBS; ++i) {
        wide_t s = (wide_t)a[i] + (wide_t)b[i] + carry;
        r[i] = (limb_t)(s & 0xFFFFu);
        carry = s >> 16;
    }
    return (limb_t)carry;
}

/* r = a - b, assuming a >= b; return borrow (0/1, normally 0) */
static limb_t limb_sub(limb_t *r, const limb_t *a, const limb_t *b)
{
    long borrow = 0;
    int i;
    for (i = 0; i < AP_LIMBS; ++i) {
        long s = (long)a[i] - (long)b[i] - borrow;
        if (s < 0) {
            s += (1L << 16);
            borrow = 1;
        } else {
            borrow = 0;
        }
        r[i] = (limb_t)(s & 0xFFFFu);
    }
    return (limb_t)borrow;
}

/* shift left by bits (0<bits<16) within limbs */
static void limb_shl_bits(limb_t *a, int bits)
{
    limb_t carry = 0;
    int i;
    if (bits <= 0 || bits >= 16) return;
    for (i = 0; i < AP_LIMBS; ++i) {
        wide_t w = (wide_t)a[i];
        wide_t nw = (w << bits) | carry;
        a[i] = (limb_t)(nw & 0xFFFFu);
        carry = (limb_t)(nw >> 16);
    }
}

/* shift right by bits (0<bits<16) within limbs */
static void limb_shr_bits(limb_t *a, int bits)
{
    limb_t carry = 0;
    int i;
    if (bits <= 0 || bits >= 16) return;
    for (i = AP_LIMBS - 1; i >= 0; --i) {
        wide_t w = (wide_t)a[i];
        wide_t nw = (w >> bits) | ((wide_t)carry << (16 - bits));
        carry = (limb_t)(w & ((1u << bits) - 1u));
        a[i] = (limb_t)(nw & 0xFFFFu);
    }
}

/* a <<= k bits (k>=0) */
static void limb_shl(limb_t *a, long k)
{
    long w;
    int b, i;
    
    if (k <= 0) return;
    w = k / 16;
    b = (int)(k % 16);

    if (w >= (long)AP_LIMBS) {
        limb_zero(a);
        return;
    }
    /* word shift */
    for (i = AP_LIMBS - 1; i >= 0; --i) {
        long src = (long)i - w;
        a[i] = (src >= 0) ? a[src] : 0;
    }
    /* bit shift */
    limb_shl_bits(a, b);
}

/* a >>= k bits (k>=0) */
static void limb_shr(limb_t *a, long k)
{
    long w;
    int b, i;
    
    if (k <= 0) return;
    w = k / 16;
    b = (int)(k % 16);

    if (w >= (long)AP_LIMBS) {
        limb_zero(a);
        return;
    }
    /* word shift */
    for (i = 0; i < AP_LIMBS; ++i) {
        long src = (long)i + w;
        a[i] = (src < (long)AP_LIMBS) ? a[src] : 0;
    }
    /* bit shift */
    limb_shr_bits(a, b);
}

/* highest set bit index (0..AP_BITS-1), or -1 if zero */
static long limb_highbit(const limb_t *a, int n_limbs)
{
    int i, bit;
    for (i = n_limbs - 1; i >= 0; --i) {
        limb_t w = a[i];
        if (w) {
            for (bit = 15; bit >= 0; --bit) {
                if (w & (limb_t)(1u << bit))
                    return (long)i * 16 + bit;
            }
        }
    }
    return -1;
}

/* ---------- basic apf ---------- */

void apf_zero(apf *x)
{
    x->cls  = APF_CLASS_ZERO;
    x->sign = 0;
    x->exp  = 0;
    limb_zero(x->mant);
}

void apf_set_nan(apf *x)
{
    x->cls  = APF_CLASS_NAN;
    x->sign = 0;
    x->exp  = 0;
    limb_zero(x->mant);
}

void apf_set_inf(apf *x, int sign)
{
    x->cls  = APF_CLASS_INF;
    x->sign = sign ? 1 : 0;
    x->exp  = 0;
    limb_zero(x->mant);
}

void apf_from_int(apf *x, long v)
{
    apf_zero(x);
    if (v == 0) return;

    x->cls  = APF_CLASS_NORMAL;
    x->sign = (v < 0) ? 1 : 0;
    if (v < 0) v = -v;

    limb_zero(x->mant);
    /* store v into low limbs (assuming v fits in 32 bits) */
    x->mant[0] = (limb_t)(v & 0xFFFFL);
    x->mant[1] = (limb_t)((v >> 16) & 0xFFFFL);
    x->exp     = 0;

    apf_norm(x);
}

/* Convert from double - useful for interfacing with native math */
void apf_from_double(apf *r, double d)
{
    char buf[64];
    if (d == 0.0) {
        apf_zero(r);
        return;
    }
    /* Use sprintf to convert, then parse - simple and portable */
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

int apf_isnan(const apf *x)  { return x->cls == APF_CLASS_NAN; }
int apf_isinf(const apf *x)  { return x->cls == APF_CLASS_INF; }
int apf_iszero(const apf *x) { return x->cls == APF_CLASS_ZERO; }

/* Set to largest finite positive value */
void apf_set_max(apf *x)
{
    int i;
    x->cls = APF_CLASS_NORMAL;
    x->sign = 0;
    x->exp = APF_EXP_MAX - (AP_BITS - 1);  /* Adjust for normalized mantissa */
    for (i = 0; i < AP_LIMBS; i++) {
        x->mant[i] = 0xFFFFu;
    }
}

/* Set to smallest positive non-zero value */
void apf_set_min(apf *x)
{
    int i;
    x->cls = APF_CLASS_NORMAL;
    x->sign = 0;
    x->exp = APF_EXP_MIN;
    for (i = 0; i < AP_LIMBS - 1; i++) {
        x->mant[i] = 0;
    }
    x->mant[AP_LIMBS - 1] = 0x8000u;  /* MSB set only */
}

/* ---------- normalization ---------- */

void apf_norm(apf *x)
{
    long hb, shift;

    if (x->cls != APF_CLASS_NORMAL)
        return;

    if (limb_is_zero(x->mant)) {
        apf_zero(x);
        return;
    }

    hb = limb_highbit(x->mant, AP_LIMBS);
    if (hb < 0) {
        apf_zero(x);
        return;
    }

    shift = (AP_BITS - 1) - hb; /* want MSB at AP_BITS-1 */
    if (shift > 0) {
        limb_shl(x->mant, shift);
        x->exp -= shift;
    } else if (shift < 0) {
        limb_shr(x->mant, -shift);
        x->exp += -shift;
    }
}

/* ---------- comparisons ---------- */

int apf_cmp(const apf *a, const apf *b)
{
    apf A = *a, B = *b;
    int c;

    if (A.cls == APF_CLASS_NAN || B.cls == APF_CLASS_NAN)
        return 0;

    if (A.cls == APF_CLASS_ZERO && B.cls == APF_CLASS_ZERO)
        return 0;

    if (A.cls == APF_CLASS_INF || B.cls == APF_CLASS_INF) {
        if (A.cls == APF_CLASS_INF && B.cls == APF_CLASS_INF) {
            if (A.sign == B.sign) return 0;
            return A.sign ? -1 : 1;
        }
        if (A.cls == APF_CLASS_INF)
            return A.sign ? -1 : 1;
        else
            return B.sign ? 1 : -1;
    }

    if (A.cls == APF_CLASS_ZERO && B.cls == APF_CLASS_NORMAL)
        return B.sign ? 1 : -1;
    if (B.cls == APF_CLASS_ZERO && A.cls == APF_CLASS_NORMAL)
        return A.sign ? -1 : 1;

    if (A.cls == APF_CLASS_NORMAL && B.cls == APF_CLASS_NORMAL) {
        apf_norm(&A);
        apf_norm(&B);
        if (A.sign != B.sign)
            return A.sign ? -1 : 1;
        if (A.exp < B.exp) return A.sign ? 1 : -1;
        if (A.exp > B.exp) return A.sign ? -1 : 1;
        c = limb_cmp(A.mant, B.mant);
        if (c == 0) return 0;
        return A.sign ? -c : c;
    }

    if (A.cls == APF_CLASS_ZERO) return B.sign ? 1 : -1;
    if (B.cls == APF_CLASS_ZERO) return A.sign ? -1 : 1;

    return 0;
}

int apf_eq(const apf *a, const apf *b) { return !apf_isnan(a) && !apf_isnan(b) && apf_cmp(a,b) == 0; }
int apf_ne(const apf *a, const apf *b) { return  apf_isnan(a) ||  apf_isnan(b) || apf_cmp(a,b) != 0; }
int apf_lt(const apf *a, const apf *b) { return !apf_isnan(a) && !apf_isnan(b) && apf_cmp(a,b) < 0; }
int apf_le(const apf *a, const apf *b) { return !apf_isnan(a) && !apf_isnan(b) && apf_cmp(a,b) <= 0; }
int apf_gt(const apf *a, const apf *b) { return !apf_isnan(a) && !apf_isnan(b) && apf_cmp(a,b) > 0; }
int apf_ge(const apf *a, const apf *b) { return !apf_isnan(a) && !apf_isnan(b) && apf_cmp(a,b) >= 0; }

/* ---------- arithmetic ---------- */

void apf_neg(apf *r, const apf *a)
{
    apf_copy(r, a);
    if (r->cls != APF_CLASS_NAN) {
        r->sign ^= 1;
    }
}

void apf_add(apf *r, const apf *a, const apf *b)
{
    apf A = *a, B = *b;
    long diff;
    int i;
    static limb_t sum[AP_LIMBS];      /* static to save stack */
    static limb_t tmp[AP_LIMBS + 1];  /* static to save stack */
    limb_t carry;

    /* handle specials */
    if (A.cls == APF_CLASS_NAN || B.cls == APF_CLASS_NAN) {
        apf_set_nan(r);
        return;
    }
    if (A.cls == APF_CLASS_INF) {
        if (B.cls == APF_CLASS_INF && A.sign != B.sign) {
            apf_set_nan(r);
        } else {
            apf_set_inf(r, A.sign);
        }
        return;
    }
    if (B.cls == APF_CLASS_INF) {
        apf_set_inf(r, B.sign);
        return;
    }
    if (A.cls == APF_CLASS_ZERO) {
        *r = B;
        return;
    }
    if (B.cls == APF_CLASS_ZERO) {
        *r = A;
        return;
    }

    /* both NORMAL */
    apf_norm(&A);
    apf_norm(&B);

    /* align exponents: shift the smaller towards the larger */
    if (A.exp < B.exp) {
        diff = B.exp - A.exp;
        if (diff >= AP_BITS) {
            *r = B;
            return;
        }
        limb_shr(A.mant, diff);
        A.exp = B.exp;
    } else if (B.exp < A.exp) {
        diff = A.exp - B.exp;
        if (diff >= AP_BITS) {
            *r = A;
            return;
        }
        limb_shr(B.mant, diff);
        B.exp = A.exp;
    }

    r->cls = APF_CLASS_NORMAL;

    if (A.sign == B.sign) {
        /* addition of magnitudes */
        carry = limb_add(sum, A.mant, B.mant);
        if (carry) {
            for (i = 0; i < AP_LIMBS; ++i) tmp[i] = sum[i];
            tmp[AP_LIMBS] = carry;
            /* shift (AP_LIMBS+1) limbs right by 1 bit */
            {
                limb_t c = 0;
                for (i = AP_LIMBS; i >= 0; --i) {
                    wide_t w = (wide_t)tmp[i];
                    wide_t nw = (w >> 1) | ((wide_t)c << 15);
                    c = (limb_t)(w & 1u);
                    tmp[i] = (limb_t)(nw & 0xFFFFu);
                }
            }
            for (i = 0; i < AP_LIMBS; ++i) r->mant[i] = tmp[i];
            r->exp = A.exp + 1;
            /* Check for exponent overflow */
            if (r->exp > APF_EXP_MAX - (AP_BITS - 1)) {
                apf_set_inf(r, A.sign);
                return;
            }
        } else {
            limb_copy(r->mant, sum);
            r->exp = A.exp;
        }
        r->sign = A.sign;
    } else {
        /* subtraction of magnitudes */
        int cmp = limb_cmp(A.mant, B.mant);
        if (cmp == 0) {
            apf_zero(r);
            return;
        } else if (cmp > 0) {
            limb_sub(r->mant, A.mant, B.mant);
            r->sign = A.sign;
            r->exp  = A.exp;
        } else {
            limb_sub(r->mant, B.mant, A.mant);
            r->sign = B.sign;
            r->exp  = A.exp;
        }
    }

    apf_norm(r);
}

void apf_sub(apf *r, const apf *a, const apf *b)
{
    apf nb = *b;
    nb.sign ^= 1;
    apf_add(r, a, &nb);
}

/* ---------- helpers ---------- */

long apf_to_long(const apf *x)
{
    apf tmp = *x;
    int neg;
    long hb, val_exp;
    wide_t top32;

    if (tmp.cls == APF_CLASS_NAN) return 0;
    if (tmp.cls == APF_CLASS_INF) return tmp.sign ? (long)0x80000000L : (long)0x7FFFFFFFL;
    if (tmp.cls == APF_CLASS_ZERO) return 0;

    apf_norm(&tmp);
    neg = tmp.sign;

    /* After normalization, MSB is at bit AP_BITS-1.
     * Value = mant * 2^exp
     * The highest bit of the value is at position (AP_BITS-1) + exp.
     */
    hb = limb_highbit(tmp.mant, AP_LIMBS);
    if (hb < 0) return 0;  /* zero mantissa */

    val_exp = hb + tmp.exp;  /* position of MSB in actual value */

    if (val_exp < 0) {
        /* Value < 1, integer part is 0 */
        return 0;
    }
    if (val_exp >= 31) {
        /* Overflow */
        return neg ? (long)0x80000000L : (long)0x7FFFFFFFL;
    }

    /* Extract top 32 bits of mantissa */
    top32 = ((wide_t)tmp.mant[AP_LIMBS - 1] << 16) | (wide_t)tmp.mant[AP_LIMBS - 2];

    /* top32 represents bits [AP_BITS-1 .. AP_BITS-32] of mantissa.
     * We need to shift to get the integer value.
     * 
     * The value is: top32 * 2^(exp + AP_BITS - 32 - (AP_BITS-1))
     *             = top32 * 2^(exp - 31)
     *
     * For integer conversion:
     * - If exp >= 31, shift left (but we'd overflow, handled above)
     * - If exp < 31, shift right by (31 - exp) = (31 - exp) - (hb - (AP_BITS-1))
     *
     * Actually simpler: the MSB of top32 is at bit 31.
     * Value's MSB is at position val_exp.
     * So we need to shift right by (31 - val_exp) to get integer.
     */
    if (val_exp < 31) {
        top32 >>= (31 - val_exp);
    }

    /* Clamp to 31 bits for signed long */
    if (top32 > 0x7FFFFFFFUL) {
        top32 = 0x7FFFFFFFUL;
    }

    return neg ? -(long)top32 : (long)top32;
}

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

/* full product: out[0..2*n-1] = a*b, where a,b are n-limb values */
static void limb_mul_full(limb_t *out, const limb_t *a, const limb_t *b, int n)
{
    int i, j, k;
    wide_t carry, acc;
    
    for (i = 0; i < 2*n; ++i) out[i] = 0;

    for (i = 0; i < n; ++i) {
        carry = 0;
        for (j = 0; j < n; ++j) {
            k = i + j;
            acc = (wide_t)out[k]
                + (wide_t)a[i] * (wide_t)b[j]
                + carry;
            out[k] = (limb_t)(acc & 0xFFFFu);
            carry  = acc >> 16;
        }
        out[i + n] = (limb_t)((wide_t)out[i + n] + carry);
    }
}

/* shift a 2*n-limb array right by k bits (k>=0) */
static void limb2n_shr(limb_t *a, long k, int n)
{
    int total = 2 * n;
    long w;
    int b, i;
    static limb_t tmp[2*AP_LIMBS];  /* static to save stack */

    if (k <= 0) return;
    w = k / 16;
    b = (int)(k % 16);

    if (w >= (long)total) {
        for (i = 0; i < total; ++i) a[i] = 0;
        return;
    }

    /* word shift */
    for (i = 0; i < total; ++i) {
        long src = (long)i + w;
        tmp[i] = (src < (long)total) ? a[src] : 0;
    }
    for (i = 0; i < total; ++i) a[i] = tmp[i];

    /* bit shift */
    if (b > 0 && b < 16) {
        limb_t carry = 0;
        for (i = total - 1; i >= 0; --i) {
            wide_t wv = (wide_t)a[i];
            wide_t nv = (wv >> b) | ((wide_t)carry << (16 - b));
            carry = (limb_t)(wv & ((1u << b) - 1u));
            a[i] = (limb_t)(nv & 0xFFFFu);
        }
    }
}

/* Multiply: r = a * b */
void apf_mul(apf *r, const apf *a, const apf *b)
{
    apf A = *a, B = *b;
    static limb_t prod[2*AP_LIMBS];  /* static to save stack */
    long hb, shift;
    int i;

    /* specials */
    if (A.cls == APF_CLASS_NAN || B.cls == APF_CLASS_NAN) {
        apf_set_nan(r);
        return;
    }
    if ((A.cls == APF_CLASS_INF && B.cls == APF_CLASS_ZERO) ||
        (B.cls == APF_CLASS_INF && A.cls == APF_CLASS_ZERO)) {
        apf_set_nan(r);
        return;
    }
    if (A.cls == APF_CLASS_INF || B.cls == APF_CLASS_INF) {
        apf_set_inf(r, A.sign ^ B.sign);
        return;
    }
    if (A.cls == APF_CLASS_ZERO || B.cls == APF_CLASS_ZERO) {
        apf_zero(r);
        r->sign = A.sign ^ B.sign;
        return;
    }

    /* NORMAL x NORMAL */
    apf_norm(&A);
    apf_norm(&B);

    limb_mul_full(prod, A.mant, B.mant, AP_LIMBS);

    /* find highest bit in the 2*AP_LIMBS-limb product */
    hb = limb_highbit(prod, 2*AP_LIMBS);
    if (hb < 0) {
        apf_zero(r);
        return;
    }

    /* we want top bit at AP_BITS-1 => shift right by hb-(AP_BITS-1) */
    shift = hb - (AP_BITS - 1);
    if (shift < 0) shift = 0;
    limb2n_shr(prod, shift, AP_LIMBS);

    r->cls  = APF_CLASS_NORMAL;
    r->sign = A.sign ^ B.sign;
    r->exp  = A.exp + B.exp + shift;

    /* Check for exponent overflow -> infinity */
    if (r->exp > APF_EXP_MAX - (AP_BITS - 1)) {
        apf_set_inf(r, r->sign);
        return;
    }
    /* Check for exponent underflow -> zero */
    if (r->exp < APF_EXP_MIN) {
        apf_zero(r);
        r->sign = A.sign ^ B.sign;
        return;
    }

    /* copy low AP_LIMBS limbs into mant */
    for (i = 0; i < AP_LIMBS; ++i) {
        r->mant[i] = prod[i];
    }

    apf_norm(r);
}

/* ---------- Division helpers ---------- */

/* Compare n-limb numbers: return -1, 0, or 1 */
static int limbn_cmp(const limb_t *a, const limb_t *b, int n)
{
    int i;
    for (i = n - 1; i >= 0; --i) {
        if (a[i] < b[i]) return -1;
        if (a[i] > b[i]) return 1;
    }
    return 0;
}

/* Subtract n-limb: r = a - b. Assumes a >= b. */
static void limbn_sub(limb_t *r, const limb_t *a, const limb_t *b, int n)
{
    long borrow = 0;
    int i;
    for (i = 0; i < n; ++i) {
        long s = (long)a[i] - (long)b[i] - borrow;
        if (s < 0) {
            s += (1L << 16);
            borrow = 1;
        } else {
            borrow = 0;
        }
        r[i] = (limb_t)(s & 0xFFFFu);
    }
}

/* Shift n-limb array left by 1 bit, return carry out */
static limb_t limbn_shl1(limb_t *a, int n)
{
    int i;
    limb_t carry = 0;
    limb_t new_carry;
    for (i = 0; i < n; ++i) {
        new_carry = (a[i] >> 15) & 1u;
        a[i] = (limb_t)(((a[i] << 1) | carry) & 0xFFFFu);
        carry = new_carry;
    }
    return carry;
}

/* Long division of normalized mantissas.
 * Both num and den have MSB set (bit AP_BITS-1).
 * Produces AP_BITS quotient bits in quot.
 * Returns exp_adjust: 0 if num >= den, -1 if num < den
 */
static int limb_div(limb_t *quot, const limb_t *num, const limb_t *den)
{
    static limb_t remainder[AP_LIMBS];  /* static to save stack */
    int i, bit;
    int exp_adjust;
    limb_t carry;
    int limb_idx, bit_idx;

    /* Copy num to remainder */
    for (i = 0; i < AP_LIMBS; ++i) {
        remainder[i] = num[i];
    }

    /* Clear quotient */
    for (i = 0; i < AP_LIMBS; ++i) {
        quot[i] = 0;
    }

    carry = 0;

    /* Check if num >= den */
    if (limbn_cmp(remainder, den, AP_LIMBS) >= 0) {
        /* num >= den: quotient will be in [1, 2) */
        exp_adjust = 0;
        /* Subtract den from remainder, set MSB of quotient */
        limbn_sub(remainder, remainder, den, AP_LIMBS);
        quot[AP_LIMBS - 1] = 0x8000u;
        
        /* Now produce remaining AP_BITS-1 quotient bits */
        for (bit = AP_BITS - 2; bit >= 0; --bit) {
            limb_idx = bit / 16;
            bit_idx = bit % 16;

            /* Shift remainder left by 1 */
            carry = limbn_shl1(remainder, AP_LIMBS);

            /* If carry set OR remainder >= den, subtract and set quotient bit */
            if (carry || limbn_cmp(remainder, den, AP_LIMBS) >= 0) {
                limbn_sub(remainder, remainder, den, AP_LIMBS);
                quot[limb_idx] |= (limb_t)(1u << bit_idx);
                carry = 0;
            }
        }
    } else {
        /* num < den: quotient will be in [0.5, 1) */
        exp_adjust = -1;
        
        /* Produce AP_BITS quotient bits, starting with MSB */
        for (bit = AP_BITS - 1; bit >= 0; --bit) {
            limb_idx = bit / 16;
            bit_idx = bit % 16;

            /* Shift remainder left by 1 */
            carry = limbn_shl1(remainder, AP_LIMBS);

            /* If carry set OR remainder >= den, subtract and set quotient bit */
            if (carry || limbn_cmp(remainder, den, AP_LIMBS) >= 0) {
                limbn_sub(remainder, remainder, den, AP_LIMBS);
                quot[limb_idx] |= (limb_t)(1u << bit_idx);
                carry = 0;
            }
        }
    }

    return exp_adjust;
}

/* Division: r = a / b */
void apf_div(apf *r, const apf *a, const apf *b)
{
    apf A = *a, B = *b;
    static limb_t quot[AP_LIMBS];  /* static to save stack */
    int exp_adjust;

    /* Specials: follow IEEE-ish behaviour */

    if (A.cls == APF_CLASS_NAN || B.cls == APF_CLASS_NAN) {
        apf_set_nan(r);
        return;
    }

    /* 0/0, Inf/Inf -> NaN */
    if ((A.cls == APF_CLASS_ZERO && B.cls == APF_CLASS_ZERO) ||
        (A.cls == APF_CLASS_INF  && B.cls == APF_CLASS_INF)) {
        apf_set_nan(r);
        return;
    }

    /* finite / 0 -> Inf */
    if (B.cls == APF_CLASS_ZERO) {
        if (A.cls == APF_CLASS_ZERO) {
            apf_set_nan(r);
        } else {
            apf_set_inf(r, A.sign ^ B.sign);
        }
        return;
    }

    /* 0 / finite -> 0 */
    if (A.cls == APF_CLASS_ZERO) {
        apf_zero(r);
        r->sign = A.sign ^ B.sign;
        return;
    }

    /* Inf / finite -> Inf */
    if (A.cls == APF_CLASS_INF && B.cls == APF_CLASS_NORMAL) {
        apf_set_inf(r, A.sign ^ B.sign);
        return;
    }

    /* finite / Inf -> 0 */
    if (A.cls == APF_CLASS_NORMAL && B.cls == APF_CLASS_INF) {
        apf_zero(r);
        r->sign = A.sign ^ B.sign;
        return;
    }

    /* NORMAL / NORMAL: use long division */
    apf_norm(&A);
    apf_norm(&B);

    /* Perform mantissa division */
    exp_adjust = limb_div(quot, A.mant, B.mant);

    /* Set up result */
    r->cls = APF_CLASS_NORMAL;
    r->sign = A.sign ^ B.sign;
    limb_copy(r->mant, quot);

    /* Compute exponent */
    r->exp = A.exp - B.exp - (AP_BITS - 1) + exp_adjust;

    apf_norm(r);
}

/* Square root via Newton-Raphson: x_{n+1} = (x_n + S/x_n) / 2 */
void apf_sqrt(apf *r, const apf *a)
{
    apf S;
    apf x, x_new;
    apf quotient, sum;
    apf two;
    int i, iterations;

    /* Handle special cases */
    if (a->cls == APF_CLASS_NAN || a->sign) {
        /* sqrt(NaN) = NaN, sqrt(negative) = NaN */
        apf_set_nan(r);
        return;
    }
    if (a->cls == APF_CLASS_ZERO) {
        apf_zero(r);
        return;
    }
    if (a->cls == APF_CLASS_INF) {
        apf_set_inf(r, 0);
        return;
    }

    S = *a;
    apf_norm(&S);

    /* Initial guess: x0 = 2^(exp/2) 
     * For normalized S: S = mant * 2^exp where mant in [0.5, 1)
     * sqrt(S) ~ sqrt(mant) * 2^(exp/2) ~ 2^(exp/2) (rough)
     *
     * We set x0.mant = 0x8000...0 (normalized 1.0 equivalent)
     * and x0.exp = S.exp / 2
     */
    apf_from_int(&x, 1);
    /* Adjust exponent for initial guess.
     * After from_int(1), x represents 1.0 with exp = -(AP_BITS-1)
     * We want x ~ 2^(S.exp/2), so x.exp should give us that value.
     * Current value = 2^(AP_BITS-1) * 2^x.exp = 1 (so x.exp = -(AP_BITS-1))
     * We want value = 2^((S.exp + (AP_BITS-1))/2)
     * Hmm, let's think differently...
     *
     * S = S.mant * 2^S.exp, with S.mant having MSB at bit AP_BITS-1
     * So S ~ 2^(AP_BITS-1) * 2^S.exp = 2^(AP_BITS-1 + S.exp)
     * sqrt(S) ~ 2^((AP_BITS-1 + S.exp)/2)
     *
     * For x0: x0 = x0.mant * 2^x0.exp ~ 2^(AP_BITS-1 + x0.exp)
     * We want: AP_BITS-1 + x0.exp = (AP_BITS-1 + S.exp)/2
     * So: x0.exp = (AP_BITS-1 + S.exp)/2 - (AP_BITS-1)
     *            = (S.exp - AP_BITS + 1) / 2
     *            = (S.exp - (AP_BITS-1)) / 2
     */
    x.exp = (S.exp - (AP_BITS - 1)) / 2;

    apf_from_int(&two, 2);

    /* Newton iterations: x = (x + S/x) / 2
     * Convergence is quadratic, so log2(AP_BITS) + 2 iterations suffice.
     * 10 iterations handles up to 512+ bits safely.
     */
    iterations = 10;
    if (AP_BITS > 256) iterations = 12;
    if (AP_BITS > 512) iterations = 14;

    for (i = 0; i < iterations; i++) {
        apf_div(&quotient, &S, &x);      /* S / x */
        apf_add(&sum, &x, &quotient);    /* x + S/x */
        apf_div(&x_new, &sum, &two);     /* (x + S/x) / 2 */
        x = x_new;
    }

    *r = x;
}

/* Convert apf to decimal string
 * Uses scientific notation for very large or very small numbers
 */
char *apf_to_str(char *buf, int bufsize, const apf *val, int max_digits)
{
    apf x, ten, digit_apf, temp;
    char *p = buf;
    char *end = buf + bufsize - 1;
    static char digits[128];  /* static to save stack */
    int num_digits;
    int i, digit;
    int max_frac;
    int has_fraction;
    int decimal_pos;
    long dec_exp;          /* decimal exponent for scientific notation */
    int use_sci;           /* use scientific notation? */
    int sig_digits;        /* significant digits to show */

    if (bufsize < 2) {
        buf[0] = '\0';
        return buf;
    }

    if (val->cls == APF_CLASS_NAN) {
        if (p + 3 <= end) { *p++ = 'N'; *p++ = 'a'; *p++ = 'N'; }
        *p = '\0';
        return buf;
    }
    if (val->cls == APF_CLASS_INF) {
        if (val->sign && p < end) *p++ = '-';
        if (p + 3 <= end) { *p++ = 'I'; *p++ = 'n'; *p++ = 'f'; }
        *p = '\0';
        return buf;
    }
    if (val->cls == APF_CLASS_ZERO) {
        *p++ = '0';
        *p = '\0';
        return buf;
    }

    x = *val;

    /* Handle sign */
    if (x.sign) {
        if (p < end) *p++ = '-';
        x.sign = 0;
    }

    apf_from_int(&ten, 10);
    
    /* Estimate decimal exponent from binary exponent
     * dec_exp ~ (binary_exp + AP_BITS - 1) * log10(2)
     * log10(2) = 0.30102999566... ≈ 30103 / 100000
     * For normalized mantissa, effective binary exp is exp + (AP_BITS - 1)
     * 
     * Split calculation carefully to avoid INT32 overflow:
     * bin_exp * 30103 / 100000 can overflow, so we split 30103 = 301*100 + 3
     * and compute in parts that each fit in INT32.
     */
    {
        INT32 bin_exp = x.exp + (AP_BITS - 1);
        INT32 rem;
        /* Split: 30103/100000 = 301/1000 + 3/100000 */
        /* Part 1: bin_exp * 301 / 1000 */
        dec_exp = (bin_exp / 1000L) * 301L;
        rem = bin_exp % 1000L;
        dec_exp += (rem * 301L) / 1000L;
        /* Part 2: bin_exp * 3 / 100000 (small correction term) */
        dec_exp += (bin_exp / 100000L) * 3L;
        rem = bin_exp % 100000L;
        dec_exp += (rem * 3L) / 100000L;
    }
    
    /* Determine significant digits based on mantissa precision */
    sig_digits = (max_digits > 0) ? max_digits : (int)((long)AP_BITS * 301L / 1000L);
    if (sig_digits > 100) sig_digits = 100;
    if (sig_digits < 1) sig_digits = 1;
    
    /* Use scientific notation if exponent is large */
    use_sci = (dec_exp > sig_digits + 5) || (dec_exp < -5);
    
    if (use_sci) {
        /* Scientific notation: extract sig_digits significant figures */
        apf one;
        
        apf_from_int(&one, 1);
        
        /* Scale x to be in range [1, 10) by multiplying/dividing by powers of 10 */
        /* We want x_scaled = x / 10^dec_exp to be in [1, 10) */
        
        /* Build 10^|dec_exp| using binary exponentiation */
        {
            apf pow10, base10;
            long n = dec_exp;
            int neg_scale = 0;
            
            if (n < 0) {
                neg_scale = 1;
                n = -n;
            }
            
            apf_from_int(&pow10, 1);
            apf_from_int(&base10, 10);
            
            while (n > 0) {
                if (n & 1) {
                    apf_mul(&temp, &pow10, &base10);
                    pow10 = temp;
                }
                apf_mul(&temp, &base10, &base10);
                base10 = temp;
                n >>= 1;
            }
            
            /* Scale x */
            if (neg_scale) {
                apf_mul(&x, &x, &pow10);
            } else {
                apf_div(&x, &x, &pow10);
            }
        }
        
        /* Adjust if x is not in [1, 10) */
        while (apf_ge(&x, &ten)) {
            apf_div(&temp, &x, &ten);
            x = temp;
            dec_exp++;
        }
        while (apf_lt(&x, &one) && !apf_is_zero(&x)) {
            apf_mul(&temp, &x, &ten);
            x = temp;
            dec_exp--;
        }
        
        /* Extract significant digits */
        num_digits = 0;
        for (i = 0; i < sig_digits + 1 && i < 100; i++) {
            digit = apf_to_long(&x);
            if (digit < 0) digit = 0;
            if (digit > 9) digit = 9;
            digits[num_digits++] = '0' + digit;
            
            apf_from_int(&digit_apf, digit);
            apf_sub(&temp, &x, &digit_apf);
            apf_mul(&x, &temp, &ten);
        }
        
        /* Round based on the extra digit and rounding mode */
        if (num_digits > sig_digits) {
            int round_digit = digits[sig_digits] - '0';
            int do_round_up = 0;
            int is_negative = val->sign;
            
            num_digits = sig_digits;
            
            switch (apf_round_mode) {
                case APF_ROUND_CEIL:
                    /* Round toward +infinity: round up if positive and any extra digits */
                    do_round_up = !is_negative && round_digit > 0;
                    break;
                case APF_ROUND_FLOOR:
                    /* Round toward -infinity: round up (in magnitude) if negative and any extra */
                    do_round_up = is_negative && round_digit > 0;
                    break;
                case APF_ROUND_TRUNC:
                    /* Round toward zero: never round up, just truncate */
                    do_round_up = 0;
                    break;
                case APF_ROUND_AWAY:
                    /* Round to nearest, ties away from zero */
                    do_round_up = round_digit >= 5;
                    break;
                case APF_ROUND_NEAREST:
                default:
                    /* Round to nearest, ties to even */
                    if (round_digit > 5) {
                        do_round_up = 1;
                    } else if (round_digit == 5) {
                        /* Tie: round to even (round up if last kept digit is odd) */
                        int last_digit = (num_digits > 0) ? (digits[num_digits - 1] - '0') : 0;
                        do_round_up = (last_digit % 2) == 1;
                    } else {
                        do_round_up = 0;
                    }
                    break;
            }
            
            if (do_round_up) {
                int carry = 1;
                for (i = num_digits - 1; i >= 0 && carry; i--) {
                    int d = digits[i] - '0' + carry;
                    if (d > 9) {
                        digits[i] = '0';
                        carry = 1;
                    } else {
                        digits[i] = '0' + d;
                        carry = 0;
                    }
                }
                if (carry) {
                    /* Overflow: shift and adjust exponent */
                    for (i = num_digits; i > 0; i--) {
                        digits[i] = digits[i - 1];
                    }
                    digits[0] = '1';
                    dec_exp++;
                }
            }
        }
        
        /* Trim trailing zeros but keep at least one digit */
        while (num_digits > 1 && digits[num_digits - 1] == '0') {
            num_digits--;
        }
        
        /* Output: d.dddddde+/-ddd */
        if (num_digits > 0 && p < end) *p++ = digits[0];
        if (num_digits > 1 && p < end) {
            *p++ = '.';
            for (i = 1; i < num_digits && p < end; i++) {
                *p++ = digits[i];
            }
        }
        
        /* Exponent */
        if (p < end) *p++ = 'e';
        if (dec_exp >= 0) {
            if (p < end) *p++ = '+';
        }
        /* Convert exponent to string */
        {
            char exp_buf[16];
            int exp_len = 0;
            long e = dec_exp;
            
            if (e < 0) {
                e = -e;
                if (p < end) *p++ = '-';
            }
            
            if (e == 0) {
                exp_buf[exp_len++] = '0';
            } else {
                while (e > 0) {
                    exp_buf[exp_len++] = '0' + (int)(e % 10);
                    e /= 10;
                }
            }
            /* Reverse and copy */
            for (i = exp_len - 1; i >= 0 && p < end; i--) {
                *p++ = exp_buf[i];
            }
        }
        
        *p = '\0';
        return buf;
    }
    
    /* Non-scientific notation for "normal" sized numbers */
    num_digits = 0;
    {
        apf one, count_x;
        apf_from_int(&one, 1);
        count_x = x;
        
        if (apf_lt(&count_x, &one)) {
            digits[num_digits++] = '0';
            decimal_pos = 1;
        } else {
            apf magnitude, divisor, q;
            int int_digits = 0;
            
            /* Count integer digits - now limited since we only get here for small exponents */
            magnitude = x;
            while (apf_ge(&magnitude, &ten) && int_digits < 50) {
                apf_div(&temp, &magnitude, &ten);
                magnitude = temp;
                int_digits++;
            }
            int_digits++;
            
            for (i = 0; i < int_digits; i++) {
                int j;
                apf_from_int(&divisor, 1);
                for (j = 0; j < int_digits - 1 - i; j++) {
                    apf_mul(&temp, &divisor, &ten);
                    divisor = temp;
                }
                
                apf_div(&q, &x, &divisor);
                digit = apf_to_long(&q);
                if (digit < 0) digit = 0;
                if (digit > 9) digit = 9;
                digits[num_digits++] = '0' + digit;
                
                apf_from_int(&digit_apf, digit);
                apf_mul(&temp, &digit_apf, &divisor);
                apf_sub(&x, &x, &temp);
            }
            decimal_pos = num_digits;
        }
    }

    has_fraction = !apf_iszero(&x);

    if (has_fraction) {
        int total_sig;  /* Total significant digits so far */
        int frac_sig;   /* Significant figures found in fraction */
        int round_digit;
        
        /* Count significant digits already extracted (skip leading zero) */
        total_sig = decimal_pos;
        if (total_sig == 1 && digits[0] == '0') total_sig = 0;
        
        /* For fractional part with leading zeros, we need to keep extracting
         * until we have max_digits significant figures plus one for rounding */
        max_frac = (max_digits > 0) ? 100 : (AP_BITS * 3) / 10;
        if (max_frac > 100) max_frac = 100;

        /* Extract fractional digits until we have enough significant figures */
        frac_sig = 0;
        for (i = 0; i < max_frac && !apf_iszero(&x); i++) {
            apf_mul(&temp, &x, &ten);
            x = temp;
            digit = apf_to_long(&x);
            if (digit < 0) digit = 0;
            if (digit > 9) digit = 9;
            digits[num_digits++] = '0' + digit;

            apf_from_int(&digit_apf, digit);
            apf_sub(&temp, &x, &digit_apf);
            x = temp;
            
            /* Track significant figures in fraction */
            if (digit != 0 || frac_sig > 0) frac_sig++;
            
            /* Stop when we have enough for rounding (max_digits + 1 total) */
            if (max_digits > 0 && total_sig + frac_sig >= max_digits + 1) {
                break;
            }
        }

        /* Apply rounding if we have max_digits set */
        if (max_digits > 0 && num_digits > 0) {
            int sig_count = 0;
            int round_pos = -1;
            
            /* Find position after max_digits significant figures */
            for (i = 0; i < num_digits; i++) {
                /* Start counting sig figs from first non-zero, or from decimal */
                if (i < decimal_pos) {
                    /* Integer part: always count unless leading zero */
                    if (digits[i] != '0' || sig_count > 0) sig_count++;
                } else {
                    /* Fraction part: count after first non-zero */
                    if (digits[i] != '0' || sig_count > 0) sig_count++;
                }
                if (sig_count == max_digits + 1 && round_pos < 0) {
                    round_pos = i;
                    round_digit = digits[i] - '0';
                    break;
                }
            }
            
            /* Round if needed - respect rounding mode */
            if (round_pos > 0) {
                int do_round_up = 0;
                int is_negative = val->sign;
                
                switch (apf_round_mode) {
                    case APF_ROUND_CEIL:
                        /* Round toward +infinity: round up if positive and any extra digits */
                        do_round_up = !is_negative && round_digit > 0;
                        break;
                    case APF_ROUND_FLOOR:
                        /* Round toward -infinity: round up (in magnitude) if negative and any extra */
                        do_round_up = is_negative && round_digit > 0;
                        break;
                    case APF_ROUND_TRUNC:
                        /* Round toward zero: never round up, just truncate */
                        do_round_up = 0;
                        break;
                    case APF_ROUND_AWAY:
                        /* Round to nearest, ties away from zero */
                        do_round_up = round_digit >= 5;
                        break;
                    case APF_ROUND_NEAREST:
                    default:
                        /* Round to nearest, ties to even */
                        if (round_digit > 5) {
                            do_round_up = 1;
                        } else if (round_digit == 5) {
                            /* Tie: round to even (round up if last kept digit is odd) */
                            int last_digit = (round_pos > 0) ? (digits[round_pos - 1] - '0') : 0;
                            do_round_up = (last_digit % 2) == 1;
                        } else {
                            do_round_up = 0;
                        }
                        break;
                }
                
                if (do_round_up) {
                    int carry = 1;
                    for (i = round_pos - 1; i >= 0 && carry; i--) {
                        int d = digits[i] - '0' + carry;
                        if (d > 9) {
                            digits[i] = '0';
                            carry = 1;
                        } else {
                            digits[i] = '0' + d;
                            carry = 0;
                        }
                    }
                    /* Handle carry overflow (e.g., 0.9999 -> 1.0000) */
                    if (carry && decimal_pos > 0) {
                        /* Shift everything right and insert 1 */
                        for (i = num_digits; i > 0; i--) {
                            digits[i] = digits[i - 1];
                        }
                        digits[0] = '1';
                        decimal_pos++;
                        round_pos++;
                    }
                }
            }
            
            /* Truncate after max_digits significant figures */
            if (round_pos > 0) {
                num_digits = round_pos;
            }
        }

        /* Trim trailing zeros */
        while (num_digits > decimal_pos && digits[num_digits - 1] == '0') {
            num_digits--;
        }
    }

    /* Copy to output buffer */
    for (i = 0; i < num_digits && p < end; i++) {
        if (i == decimal_pos && has_fraction) {
            *p++ = '.';
            if (p >= end) break;
        }
        *p++ = digits[i];
    }
    *p = '\0';

    return buf;
}

/* ---------- additional functions for scalc compatibility ---------- */

void apf_copy(apf *dst, const apf *src)
{
    int i;
    dst->cls = src->cls;
    dst->sign = src->sign;
    dst->exp = src->exp;
    for (i = 0; i < AP_LIMBS; i++) {
        dst->mant[i] = src->mant[i];
    }
}

int apf_is_zero(const apf *x)
{
    return x->cls == APF_CLASS_ZERO;
}

int apf_is_neg(const apf *x)
{
    return x->sign && x->cls == APF_CLASS_NORMAL;
}

void apf_abs(apf *r, const apf *a)
{
    apf_copy(r, a);
    r->sign = 0;
}

/* Parse a decimal string into an apf
 * Handles: integers, decimals, negative numbers
 * Uses arbitrary precision for integer part (no overflow)
 */
void apf_from_str(apf *x, const char *s)
{
    apf ten, digit_apf, temp;
    int neg = 0;
    
    /* Skip leading whitespace */
    while (*s == ' ' || *s == '\t') s++;
    
    /* Handle sign */
    if (*s == '-') {
        neg = 1;
        s++;
    } else if (*s == '+') {
        s++;
    }
    
    /* Parse integer part using arbitrary precision */
    apf_from_int(x, 0);
    apf_from_int(&ten, 10);
    
    while (*s >= '0' && *s <= '9') {
        int d = *s - '0';
        apf_mul(&temp, x, &ten);
        apf_from_int(&digit_apf, d);
        apf_add(x, &temp, &digit_apf);
        s++;
    }
    
    /* Parse fractional part using arbitrary precision */
    if (*s == '.') {
        apf frac, divisor;
        apf_from_int(&frac, 0);
        apf_from_int(&divisor, 1);
        s++;
        
        while (*s >= '0' && *s <= '9') {
            int d = *s - '0';
            /* frac = frac * 10 + d */
            apf_mul(&temp, &frac, &ten);
            apf_from_int(&digit_apf, d);
            apf_add(&frac, &temp, &digit_apf);
            /* divisor *= 10 */
            apf_mul(&temp, &divisor, &ten);
            apf_copy(&divisor, &temp);
            s++;
        }
        
        /* Add frac/divisor to result */
        if (!apf_is_zero(&frac)) {
            apf_div(&temp, &frac, &divisor);
            apf_add(&frac, x, &temp);
            apf_copy(x, &frac);
        }
    }
    
    /* Apply sign */
    if (neg && x->cls == APF_CLASS_NORMAL) {
        x->sign = 1;
    }
}

