/* scalc20.c - Scientific Calculator for VIC-20
 * C89/cc65 compatible - 128-bit integers
 * Supports: Infix (algebraic) and RPN modes
 * 
 * Commands:
 *   rpn    - switch to RPN mode
 *   alg    - switch to algebraic mode
 *   .s     - show stack (RPN)
 *   clr    - clear stack
 *   q      - quit
 *
 * RPN: 5 3 + 2 *  -> 16
 * Infix: (5+3)*2  -> 16
 */

#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include "bigint.h"

#ifdef __CC65__
#define LINESZ 40
#define STACKSZ 4
#else
#define LINESZ 80
#define STACKSZ 8
#endif

/* ========== Global State ========== */
static char line[LINESZ];
static char buf[44];
static const char *p;
static int rpn_mode = 0;

/* RPN Stack */
static s128 stack[STACKSZ];
static int sp = 0;
static s128 lastx;

/* ========== RPN Stack Operations ========== */

static void push(const s128 *v) {
    if (sp < STACKSZ) {
        stack[sp++] = *v;
    } else {
        /* Stack full - shift down */
        int i;
        for (i = 0; i < STACKSZ - 1; i++) {
            stack[i] = stack[i + 1];
        }
        stack[STACKSZ - 1] = *v;
    }
}

static int pop(s128 *v) {
    if (sp > 0) {
        *v = stack[--sp];
        return 1;
    }
    return 0;
}

static void show_stack(void) {
    int i;
    s128 tmp;
    if (sp == 0) {
        puts("(empty)");
        return;
    }
    for (i = 0; i < sp; i++) {
        tmp = stack[i];
        printf("%d: %s\n", i, bi_fmt(buf + sizeof(buf), &tmp));
    }
}

static void save_lastx(void) {
    if (sp > 0) lastx = stack[sp - 1];
}

/* ========== Infix Parser ========== */

static int expr(s128 *r);

static void skipws(void) { 
    while (*p == ' ') p++; 
}

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
    if (*p == '-') { 
        p++; 
        if (!unary(r)) return 0; 
        bi_neg(r); 
        return 1; 
    }
    if (*p == '+') { 
        p++; 
        return unary(r); 
    }
    return atom(r);
}

static int term(s128 *r) {
    s128 b, t, rem;
    char op;
    if (!unary(r)) return 0;
    for (;;) {
        skipws(); 
        op = *p;
        if (op != '*' && op != '/') break;
        p++;
        if (!unary(&b)) return 0;
        if (op == '*') { 
            bi_mul(&t, r, &b); 
            *r = t; 
        } else { 
            bi_div(&t, &rem, r, &b); 
            *r = t; 
        }
    }
    return 1;
}

static int expr(s128 *r) {
    s128 b, t;
    char op;
    if (!term(r)) return 0;
    for (;;) {
        skipws(); 
        op = *p;
        if (op != '+' && op != '-') break;
        p++;
        if (!term(&b)) return 0;
        if (op == '+') 
            bi_add(&t, r, &b); 
        else 
            bi_sub(&t, r, &b);
        *r = t;
    }
    return 1;
}

/* ========== RPN Token Parser ========== */

static int parse_num(const char *s, s128 *r) {
    if (*s == '-' || *s == '+' || isdigit((unsigned char)*s)) {
        bi_parse(r, s);
        return 1;
    }
    return 0;
}

static int rpn_binop(char op) {
    s128 a, b, r, rem;
    if (sp < 2) {
        puts("need 2");
        return 0;
    }
    save_lastx();
    pop(&b);
    pop(&a);
    switch (op) {
        case '+': bi_add(&r, &a, &b); break;
        case '-': bi_sub(&r, &a, &b); break;
        case '*': bi_mul(&r, &a, &b); break;
        case '/': bi_div(&r, &rem, &a, &b); break;
        default: return 0;
    }
    push(&r);
    return 1;
}

static int rpn_token(const char *tok) {
    s128 val, t;
    
    /* Stack commands */
    if (strcmp(tok, "drop") == 0) {
        if (sp > 0) sp--;
        return 1;
    }
    if (strcmp(tok, "dup") == 0 || strcmp(tok, "enter") == 0) {
        if (sp > 0 && sp < STACKSZ) {
            stack[sp] = stack[sp - 1];
            sp++;
        }
        return 1;
    }
    if (strcmp(tok, "swap") == 0) {
        if (sp >= 2) {
            t = stack[sp - 1];
            stack[sp - 1] = stack[sp - 2];
            stack[sp - 2] = t;
        }
        return 1;
    }
    if (strcmp(tok, "clr") == 0 || strcmp(tok, "clear") == 0) {
        sp = 0;
        return 1;
    }
    if (strcmp(tok, ".s") == 0 || strcmp(tok, "stack") == 0) {
        show_stack();
        return 1;
    }
    if (strcmp(tok, "lastx") == 0) {
        push(&lastx);
        return 1;
    }
    if (strcmp(tok, "neg") == 0 || strcmp(tok, "chs") == 0) {
        if (sp > 0) {
            save_lastx();
            bi_neg(&stack[sp - 1]);
        }
        return 1;
    }
    
    /* Modulo */
    if (strcmp(tok, "mod") == 0) {
        s128 a, b, q, rem;
        if (sp < 2) {
            puts("need 2");
            return 1;
        }
        save_lastx();
        pop(&b);
        pop(&a);
        bi_div(&q, &rem, &a, &b);
        push(&rem);
        return 1;
    }
    
    /* Single-char operators */
    if (strlen(tok) == 1) {
        char c = tok[0];
        if (c == '+' || c == '-' || c == '*' || c == '/') {
            rpn_binop(c);
            return 1;
        }
    }
    
    /* Try as number */
    if (parse_num(tok, &val)) {
        push(&val);
        return 1;
    }
    
    return 0;
}

static void rpn_eval(const char *ln) {
    char tokbuf[LINESZ];
    const char *s = ln;
    char *d;
    s128 tmp;
    
    while (*s) {
        /* Skip whitespace */
        while (*s && *s == ' ') s++;
        if (!*s) break;
        
        /* Copy token */
        d = tokbuf;
        while (*s && *s != ' ' && (d - tokbuf) < LINESZ - 1) {
            *d++ = *s++;
        }
        *d = '\0';
        
        if (!rpn_token(tokbuf)) {
            printf("? %s\n", tokbuf);
        }
    }
    
    /* Show X register */
    if (sp > 0) {
        tmp = stack[sp - 1];
        printf("X: %s\n", bi_fmt(buf + sizeof(buf), &tmp));
    }
}

/* ========== Main ========== */

int main(void) {
    s128 r;
    int i;
    
    bi_zero(&lastx);
    
#ifdef __CC65__
    puts("SCALC20 128-BIT");
    puts("RPN/ALG .S CLR Q");
#else
    puts("scalc20 - 128-bit calculator");
    puts("Commands: rpn alg .s clr q");
    puts("RPN: 5 3 +   Infix: 5+3");
#endif

    for (;;) {
#ifndef __CC65__
        printf("%s> ", rpn_mode ? "rpn" : "alg");
#endif
        /* Read line */
        for (i = 0; i < LINESZ - 1; i++) {
            int c = getchar();
            if (c < 0 || c == '\n' || c == '\r') break;
            line[i] = (char)c;
        }
        line[i] = '\0';
        
        if (i == 0) continue;
        
        /* Check commands */
        if (strcmp(line, "q") == 0 || strcmp(line, "quit") == 0) break;
        
        if (strcmp(line, "rpn") == 0) {
            rpn_mode = 1;
            puts("RPN mode");
            continue;
        }
        if (strcmp(line, "alg") == 0) {
            rpn_mode = 0;
            puts("Algebraic mode");
            continue;
        }
        
        /* Evaluate */
        if (rpn_mode) {
            rpn_eval(line);
        } else {
            p = line;
            if (expr(&r)) {
                puts(bi_fmt(buf + sizeof(buf), &r));
            } else {
                puts("?");
            }
        }
    }
    
    return 0;
}
