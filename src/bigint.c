/* bigint.c - Signed 128-bit integer arithmetic (lean version)
 * Portable C89 for VIC-20 (cc65), DOS (Watcom), Linux (gcc)
 */
#include "bigint.h"

void bi_zero(s128 *r) {
    char i = BI_LIMBS;
    while (i--) r->m[(int)i] = 0;
    r->sign = 0;
}

int bi_is_zero(const s128 *a) {
    char i = BI_LIMBS;
    while (i--) if (a->m[(int)i]) return 0;
    return 1;
}

void bi_neg(s128 *r) {
    if (!bi_is_zero(r)) r->sign = !r->sign;
}

/* Compare magnitudes: -1 if a<b, 0 if ==, 1 if a>b */
static int mcmp(const s128 *a, const s128 *b) {
    char i = BI_LIMBS;
    while (i--) {
        if (a->m[(int)i] > b->m[(int)i]) return 1;
        if (a->m[(int)i] < b->m[(int)i]) return -1;
    }
    return 0;
}

/* Add magnitudes: r = a + b */
static void madd(s128 *r, const s128 *a, const s128 *b) {
    bi_u32 c = 0;
    char i;
    for (i = 0; i < BI_LIMBS; i++) {
        c += (bi_u32)a->m[(int)i] + b->m[(int)i];
        r->m[(int)i] = (bi_u16)c;
        c >>= 16;
    }
}

/* Sub magnitudes: r = a - b (assumes a >= b) */
static void msub(s128 *r, const s128 *a, const s128 *b) {
    bi_u32 c = 0;
    char i;
    for (i = 0; i < BI_LIMBS; i++) {
        c = (bi_u32)a->m[(int)i] - b->m[(int)i] - c;
        r->m[(int)i] = (bi_u16)c;
        c = (c >> 16) & 1;
    }
}

void bi_add(s128 *r, const s128 *a, const s128 *b) {
    int c;
    if (a->sign == b->sign) {
        madd(r, a, b);
        r->sign = a->sign;
    } else {
        c = mcmp(a, b);
        if (c == 0) bi_zero(r);
        else if (c > 0) { msub(r, a, b); r->sign = a->sign; }
        else { msub(r, b, a); r->sign = b->sign; }
    }
}

void bi_sub(s128 *r, const s128 *a, const s128 *b) {
    s128 t;
    char i = BI_LIMBS;
    while (i--) t.m[(int)i] = b->m[(int)i];
    t.sign = !b->sign;
    if (bi_is_zero(&t)) t.sign = 0;
    bi_add(r, a, &t);
}

void bi_mul(s128 *r, const s128 *a, const s128 *b) {
    char i, j, rs;
    bi_zero(r);
    rs = a->sign != b->sign;
    for (i = 0; i < BI_LIMBS; i++) {
        bi_u32 c = 0;
        bi_u16 m = b->m[(int)i];
        if (m) {
            for (j = 0; (int)i + j < BI_LIMBS; j++) {
                c += (bi_u32)a->m[(int)j] * m + r->m[(int)i + j];
                r->m[(int)i + j] = (bi_u16)c;
                c >>= 16;
            }
        }
    }
    r->sign = rs && !bi_is_zero(r);
}

void bi_div(s128 *q, s128 *rem, const s128 *a, const s128 *b) {
    int i;
    char j, qs, rs;
    bi_zero(q);
    bi_zero(rem);
    if (bi_is_zero(b)) return;
    qs = a->sign != b->sign;
    rs = a->sign;
    for (i = BI_LIMBS * 16 - 1; i >= 0; i--) {
        bi_u32 c = 0;
        for (j = 0; j < BI_LIMBS; j++) {
            c |= (bi_u32)rem->m[(int)j] << 1;
            rem->m[(int)j] = (bi_u16)c;
            c >>= 16;
        }
        rem->m[0] |= (a->m[i >> 4] >> (i & 15)) & 1;
        if (mcmp(rem, b) >= 0) {
            msub(rem, rem, b);
            q->m[i >> 4] |= (bi_u16)1 << (i & 15);
        }
    }
    q->sign = qs && !bi_is_zero(q);
    rem->sign = rs && !bi_is_zero(rem);
}

void bi_parse(s128 *r, const char *s) {
    char neg = 0;
    bi_u32 c;
    char i;
    bi_zero(r);
    while (*s == ' ') s++;
    if (*s == '-') { neg = 1; s++; }
    else if (*s == '+') s++;
    while (*s >= '0' && *s <= '9') {
        c = 0;
        for (i = 0; i < BI_LIMBS; i++) {
            c += (bi_u32)r->m[(int)i] * 10;
            r->m[(int)i] = (bi_u16)c;
            c >>= 16;
        }
        r->m[0] += (bi_u16)(*s++ - '0');
    }
    r->sign = neg && !bi_is_zero(r);
}

/* Format: writes backwards from 'end', returns pointer to start */
char *bi_fmt(char *end, s128 *a) {
    bi_u32 r;
    char i, neg;
    *--end = '\0';
    if (bi_is_zero(a)) { *--end = '0'; return end; }
    neg = a->sign;
    a->sign = 0;
    while (!bi_is_zero(a)) {
        r = 0;
        i = BI_LIMBS;
        while (i--) {
            r = (r << 16) | a->m[(int)i];
            a->m[(int)i] = (bi_u16)(r / 10);
            r %= 10;
        }
        *--end = '0' + (char)r;
    }
    if (neg) *--end = '-';
    return end;
}
