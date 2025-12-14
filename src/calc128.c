/* calc128.c - 128-bit integer calculator (lean version)
 * Portable C89 for VIC-20 (cc65), DOS (Watcom), Linux (gcc)
 * Supports: + - * / () unary-
 */
#include <stdio.h>
#include "bigint.h"

#ifdef __CC65__
#define LINESZ 40
#else
#define LINESZ 80
#endif

static char line[LINESZ];
static char buf[44];
static const char *p;

static int expr(s128 *r);

static void skipws(void) { while (*p == ' ') p++; }

static int atom(s128 *r) {
    skipws();
    if (*p == '(') {
        p++;
        if (!expr(r)) return 0;
        skipws();
        if (*p == ')') p++;
        return 1;
    }
    if (*p == '-' || *p == '+' || (*p >= '0' && *p <= '9')) {
        bi_parse(r, p);
        if (*p == '-' || *p == '+') p++;
        while (*p >= '0' && *p <= '9') p++;
        return 1;
    }
    return 0;
}

static int unary(s128 *r) {
    skipws();
    if (*p == '-') { p++; if (!unary(r)) return 0; bi_neg(r); return 1; }
    if (*p == '+') { p++; return unary(r); }
    return atom(r);
}

static int term(s128 *r) {
    s128 b, t, rem;
    char op;
    if (!unary(r)) return 0;
    for (;;) {
        skipws(); op = *p;
        if (op != '*' && op != '/') break;
        p++;
        if (!unary(&b)) return 0;
        if (op == '*') { bi_mul(&t, r, &b); *r = t; }
        else { bi_div(&t, &rem, r, &b); *r = t; }
    }
    return 1;
}

static int expr(s128 *r) {
    s128 b, t;
    char op;
    if (!term(r)) return 0;
    for (;;) {
        skipws(); op = *p;
        if (op != '+' && op != '-') break;
        p++;
        if (!term(&b)) return 0;
        if (op == '+') bi_add(&t, r, &b); else bi_sub(&t, r, &b);
        *r = t;
    }
    return 1;
}

int main(void) {
    s128 r;
    int i;
#ifdef __CC65__
    puts("128-BIT CALC +-*/()");
#else
    puts("128-bit calc +-*/() q=quit");
#endif
    for (;;) {
#ifndef __CC65__
        putchar('>'); putchar(' ');
#endif
        for (i = 0; i < LINESZ - 1; i++) {
            int c = getchar();
            if (c < 0 || c == '\n' || c == '\r') break;
            line[i] = (char)c;
        }
        line[i] = 0;
        if (!i) continue;
        if (line[0] == 'q') break;
        p = line;
        if (expr(&r)) puts(bi_fmt(buf + sizeof(buf), &r));
        else puts("?");
    }
    return 0;
}
