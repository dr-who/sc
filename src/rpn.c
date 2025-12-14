/* rpn.c - Reverse Polish Notation stack engine
 * C89 portable for DOS, Linux, VIC-20
 */

#include "config.h"

#ifdef HAVE_RPN

#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include "rpn.h"
#include "apf.h"
#include "apfc.h"
#include "apfx.h"

/* External display precision setting */
extern int display_digits;

/* Stack storage */
apfc rpn_stack[RPN_STACK_SIZE];
int rpn_sp = 0;
int rpn_mode = 0;
apfc rpn_lastx;  /* LASTX register */

/* Forward declaration for formatting */
extern void format_value_rpn(const apfc *val);

/* Save X to LASTX before operations */
static void save_lastx(void) {
    if (rpn_sp > 0) {
        apfc_copy(&rpn_lastx, &rpn_stack[rpn_sp - 1]);
    }
}

void rpn_init(void) {
    int i;
    rpn_sp = 0;
    for (i = 0; i < RPN_STACK_SIZE; i++) {
        apfc_zero(&rpn_stack[i]);
    }
    apfc_zero(&rpn_lastx);
}

void rpn_push(const apfc *val) {
    if (rpn_sp < RPN_STACK_SIZE) {
        apfc_copy(&rpn_stack[rpn_sp], val);
        rpn_sp++;
    } else {
        /* Stack full - shift down, lose bottom */
        int i;
        for (i = 0; i < RPN_STACK_SIZE - 1; i++) {
            apfc_copy(&rpn_stack[i], &rpn_stack[i + 1]);
        }
        apfc_copy(&rpn_stack[RPN_STACK_SIZE - 1], val);
    }
}

int rpn_pop(apfc *val) {
    if (rpn_sp > 0) {
        rpn_sp--;
        apfc_copy(val, &rpn_stack[rpn_sp]);
        return 1;
    }
    return 0;
}

int rpn_peek(apfc *val) {
    if (rpn_sp > 0) {
        apfc_copy(val, &rpn_stack[rpn_sp - 1]);
        return 1;
    }
    return 0;
}

void rpn_drop(void) {
    if (rpn_sp > 0) rpn_sp--;
}

void rpn_dup(void) {
    if (rpn_sp > 0 && rpn_sp < RPN_STACK_SIZE) {
        apfc_copy(&rpn_stack[rpn_sp], &rpn_stack[rpn_sp - 1]);
        rpn_sp++;
    }
}

void rpn_swap(void) {
    if (rpn_sp >= 2) {
        apfc tmp;
        apfc_copy(&tmp, &rpn_stack[rpn_sp - 1]);
        apfc_copy(&rpn_stack[rpn_sp - 1], &rpn_stack[rpn_sp - 2]);
        apfc_copy(&rpn_stack[rpn_sp - 2], &tmp);
    }
}

void rpn_clear(void) {
    rpn_sp = 0;
}

void rpn_roll_down(void) {
    if (rpn_sp >= 2) {
        apfc tmp;
        int i;
        apfc_copy(&tmp, &rpn_stack[rpn_sp - 1]);
        for (i = rpn_sp - 1; i > 0; i--) {
            apfc_copy(&rpn_stack[i], &rpn_stack[i - 1]);
        }
        apfc_copy(&rpn_stack[0], &tmp);
    }
}

void rpn_roll_up(void) {
    if (rpn_sp >= 2) {
        apfc tmp;
        int i;
        apfc_copy(&tmp, &rpn_stack[0]);
        for (i = 0; i < rpn_sp - 1; i++) {
            apfc_copy(&rpn_stack[i], &rpn_stack[i + 1]);
        }
        apfc_copy(&rpn_stack[rpn_sp - 1], &tmp);
    }
}

void rpn_show_stack(void) {
    int i;
    char labels[] = "TZYX";
    int label_start;
    
    if (rpn_sp == 0) {
        printf("  (empty)\n");
        return;
    }
    
    /* Show stack with T,Z,Y,X labels for top 4 */
    label_start = (rpn_sp > 4) ? rpn_sp - 4 : 0;
    
    for (i = 0; i < rpn_sp; i++) {
        if (i >= label_start && rpn_sp <= 4) {
            printf("  %c: ", labels[4 - rpn_sp + i]);
        } else if (i >= rpn_sp - 4) {
            printf("  %c: ", labels[i - (rpn_sp - 4)]);
        } else {
            printf(" %2d: ", i);
        }
        format_value_rpn(&rpn_stack[i]);
        printf("\n");
    }
}

/* Parse a number from string */
static int parse_number(const char *s, apfc *val) {
    const char *p = s;
    
    apfc_zero(val);
    
    /* Check for leading minus/plus - apf_from_str handles it */
    if (*p == '-' || *p == '+') {
        p++;
    }
    
    /* Must start with digit or decimal point */
    if (!isdigit((unsigned char)*p) && *p != '.') {
        return 0;
    }
    
    /* Parse real part */
    apf_from_str(&val->re, s);
    
    /* Check for trailing 'i' for imaginary */
    while (*p && (isdigit((unsigned char)*p) || *p == '.' || *p == 'e' || *p == 'E' || *p == '-' || *p == '+')) {
        p++;
    }
    if (*p == 'i' && *(p+1) == '\0') {
        /* It's a pure imaginary number */
        apf_copy(&val->im, &val->re);
        apf_zero(&val->re);
    }
    
    return 1;
}

/* Binary operators */
static int do_binary_op(char op) {
    apfc a, b, r;
    
    if (rpn_sp < 2) {
        printf("  Error: need 2 operands\n");
        return 0;
    }
    
    save_lastx();  /* Save X before operation */
    rpn_pop(&b);
    rpn_pop(&a);
    
    switch (op) {
        case '+': apfc_add(&r, &a, &b); break;
        case '-': apfc_sub(&r, &a, &b); break;
        case '*': apfc_mul(&r, &a, &b); break;
        case '/': apfc_div(&r, &a, &b); break;
#ifdef HAVE_POW
        case '^': apfc_pow(&r, &a, &b); break;
#endif
        default:
            rpn_push(&a);
            rpn_push(&b);
            return 0;
    }
    
    rpn_push(&r);
    return 1;
}

/* Unary functions */
static int do_unary_func(const char *name) {
    apfc a, r;
    
    if (rpn_sp < 1) {
        printf("  Error: need 1 operand\n");
        return 0;
    }
    
    save_lastx();  /* Save X before operation */
    rpn_pop(&a);
    
    if (strcmp(name, "neg") == 0 || strcmp(name, "chs") == 0) {
        apfc_neg(&r, &a);
    }
    /* HP-41C style functions */
    else if (strcmp(name, "inv") == 0 || strcmp(name, "1/x") == 0) {
        apfc one;
        apf_from_int(&one.re, 1);
        apf_zero(&one.im);
        apfc_div(&r, &one, &a);
    }
    else if (strcmp(name, "sq") == 0 || strcmp(name, "x2") == 0) {
        apfc_mul(&r, &a, &a);
    }
    else if (strcmp(name, "int") == 0) {
        /* Integer part - truncate toward zero */
        long n = apf_to_long(&a.re);
        apf_from_int(&r.re, n);
        apf_zero(&r.im);
    }
    else if (strcmp(name, "frac") == 0) {
        /* Fractional part */
        apf intpart;
        long n = apf_to_long(&a.re);
        apf_from_int(&intpart, n);
        apf_sub(&r.re, &a.re, &intpart);
        apf_zero(&r.im);
    }
#ifdef HAVE_CONST
    else if (strcmp(name, "deg") == 0 || strcmp(name, ">deg") == 0) {
        /* Radians to degrees: x * 180 / pi */
        apf pi, d180, tmp;
        apfx_pi(&pi);
        apf_from_int(&d180, 180);
        apf_mul(&tmp, &a.re, &d180);
        apf_div(&r.re, &tmp, &pi);
        apf_zero(&r.im);
    }
    else if (strcmp(name, "rad") == 0 || strcmp(name, ">rad") == 0) {
        /* Degrees to radians: x * pi / 180 */
        apf pi, d180, tmp;
        apfx_pi(&pi);
        apf_from_int(&d180, 180);
        apf_mul(&tmp, &a.re, &pi);
        apf_div(&r.re, &tmp, &d180);
        apf_zero(&r.im);
    }
#endif
#ifdef HAVE_SQRT
    else if (strcmp(name, "sqrt") == 0) {
        apfc_sqrt(&r, &a);
    }
#endif
#ifdef HAVE_EXP
    else if (strcmp(name, "exp") == 0) {
        apfc_exp(&r, &a);
    }
    else if (strcmp(name, "ln") == 0 || strcmp(name, "log") == 0) {
        apfc_log(&r, &a);
    }
    else if (strcmp(name, "log10") == 0) {
        /* log10(x) = ln(x) / ln(10) */
        apf ten, ln10;
        apf_from_int(&ten, 10);
        apfx_log(&ln10, &ten);
        apfc_log(&r, &a);
        apf_div(&r.re, &r.re, &ln10);
    }
    else if (strcmp(name, "10x") == 0 || strcmp(name, "10^x") == 0) {
        /* 10^x = exp(x * ln(10)) */
        apf ten, ln10, tmp;
        apf_from_int(&ten, 10);
        apfx_log(&ln10, &ten);
        apf_mul(&tmp, &a.re, &ln10);
        apfx_exp(&r.re, &tmp);
        apf_zero(&r.im);
    }
#endif
#ifdef HAVE_TRIG
    else if (strcmp(name, "sin") == 0) {
        apfc_sin(&r, &a);
    }
    else if (strcmp(name, "cos") == 0) {
        apfc_cos(&r, &a);
    }
    else if (strcmp(name, "tan") == 0) {
        apfc_tan(&r, &a);
    }
    else if (strcmp(name, "atan") == 0) {
        apfx_atan(&r.re, &a.re);
        apf_zero(&r.im);
    }
#endif
#ifdef HAVE_HYPER
    else if (strcmp(name, "sinh") == 0) {
        apf_copy(&r.re, &a.re);
        apf_zero(&r.im);
        apfx_sinh(&r.re, &r.re);
    }
    else if (strcmp(name, "cosh") == 0) {
        apf_copy(&r.re, &a.re);
        apf_zero(&r.im);
        apfx_cosh(&r.re, &r.re);
    }
    else if (strcmp(name, "tanh") == 0) {
        apf_copy(&r.re, &a.re);
        apf_zero(&r.im);
        apfx_tanh(&r.re, &r.re);
    }
    else if (strcmp(name, "asin") == 0) {
        apf_copy(&r.re, &a.re);
        apf_zero(&r.im);
        apfx_asin(&r.re, &r.re);
    }
    else if (strcmp(name, "acos") == 0) {
        apf_copy(&r.re, &a.re);
        apf_zero(&r.im);
        apfx_acos(&r.re, &r.re);
    }
#endif
#ifdef HAVE_COMPLEX
    else if (strcmp(name, "conj") == 0) {
        apfc_conj(&r, &a);
    }
    else if (strcmp(name, "re") == 0) {
        apf_copy(&r.re, &a.re);
        apf_zero(&r.im);
    }
    else if (strcmp(name, "im") == 0) {
        apf_copy(&r.re, &a.im);
        apf_zero(&r.im);
    }
    else if (strcmp(name, "arg") == 0) {
        apfc_arg(&r.re, &a);
        apf_zero(&r.im);
    }
    else if (strcmp(name, "abs") == 0) {
        apfc_abs(&r.re, &a);
        apf_zero(&r.im);
    }
#endif
#ifdef HAVE_FACTORIAL
    else if (strcmp(name, "fact") == 0) {
        long n = apf_to_long(&a.re);
        apfx_fact(&r.re, n);
        apf_zero(&r.im);
    }
#endif
    else {
        /* Unknown function */
        rpn_push(&a);
        return 0;
    }
    
    rpn_push(&r);
    return 1;
}

int rpn_process(const char *token) {
    apfc val;
    
    /* Skip empty tokens */
    if (!token || !*token) return 0;
    
    /* Stack commands */
    if (strcmp(token, "drop") == 0) { rpn_drop(); return 1; }
    if (strcmp(token, "dup") == 0 || strcmp(token, "enter") == 0) { 
        rpn_dup(); 
        return 1; 
    }
    if (strcmp(token, "swap") == 0 || strcmp(token, "x<>y") == 0) { 
        rpn_swap(); 
        return 1; 
    }
    if (strcmp(token, "clear") == 0 || strcmp(token, "clst") == 0) { 
        rpn_clear(); 
        return 1; 
    }
    if (strcmp(token, "clx") == 0) {
        /* Clear X only */
        if (rpn_sp > 0) {
            apfc_zero(&rpn_stack[rpn_sp - 1]);
        }
        return 1;
    }
    if (strcmp(token, "lastx") == 0) {
        rpn_push(&rpn_lastx);
        return 1;
    }
    if (strcmp(token, ".s") == 0 || strcmp(token, "stack") == 0) { 
        rpn_show_stack(); 
        return 1; 
    }
    if (strcmp(token, "roll") == 0 || strcmp(token, "r>") == 0) { 
        rpn_roll_down(); 
        return 1; 
    }
    if (strcmp(token, "rollu") == 0 || strcmp(token, "r<") == 0) { 
        rpn_roll_up(); 
        return 1; 
    }
    
    /* Percent: Y * X / 100 */
    if (strcmp(token, "%") == 0) {
        apfc a, b, r, hundred;
        if (rpn_sp < 2) {
            printf("  Error: need 2 operands for %%\n");
            return 1;
        }
        save_lastx();
        rpn_pop(&b);  /* X = percent */
        rpn_pop(&a);  /* Y = base */
        apf_from_int(&hundred.re, 100);
        apf_zero(&hundred.im);
        apfc_mul(&r, &a, &b);
        apfc_div(&r, &r, &hundred);
        rpn_push(&a);  /* Keep Y on stack */
        rpn_push(&r);
        return 1;
    }
    
    /* Delta percent: (X - Y) / Y * 100 */
    if (strcmp(token, "d%") == 0 || strcmp(token, "delta%") == 0) {
        apfc a, b, r, hundred;
        if (rpn_sp < 2) {
            printf("  Error: need 2 operands for d%%\n");
            return 1;
        }
        save_lastx();
        rpn_pop(&b);  /* X = new */
        rpn_pop(&a);  /* Y = old */
        apf_from_int(&hundred.re, 100);
        apf_zero(&hundred.im);
        apfc_sub(&r, &b, &a);
        apfc_div(&r, &r, &a);
        apfc_mul(&r, &r, &hundred);
        rpn_push(&r);
        return 1;
    }
    
    /* Modulo */
    if (strcmp(token, "mod") == 0) {
        apfc a, b, r;
        if (rpn_sp < 2) {
            printf("  Error: need 2 operands for mod\n");
            return 1;
        }
        save_lastx();
        rpn_pop(&b);
        rpn_pop(&a);
        /* r = a mod b (integer) */
        {
            long ai = apf_to_long(&a.re);
            long bi = apf_to_long(&b.re);
            if (bi != 0) {
                apf_from_int(&r.re, ai % bi);
            } else {
                apf_zero(&r.re);
            }
            apf_zero(&r.im);
        }
        rpn_push(&r);
        return 1;
    }
    
#ifdef HAVE_TRIG
    /* Rectangular to Polar: (x,y) -> (r,θ) */
    if (strcmp(token, ">pol") == 0 || strcmp(token, "pol") == 0) {
        apfc x, y;
        apf r, theta;
        if (rpn_sp < 2) {
            printf("  Error: need X and Y for >pol\n");
            return 1;
        }
        save_lastx();
        rpn_pop(&y);  /* Y coordinate in Y register */
        rpn_pop(&x);  /* X coordinate in X register */
        /* r = sqrt(x² + y²) */
        {
            apf x2, y2, sum;
            apf_mul(&x2, &x.re, &x.re);
            apf_mul(&y2, &y.re, &y.re);
            apf_add(&sum, &x2, &y2);
            apf_sqrt(&r, &sum);
        }
        /* θ = atan2(y, x) */
        apfx_atan2(&theta, &y.re, &x.re);
        /* Push θ then r (so r is in X, θ in Y) */
        apf_copy(&y.re, &theta);
        apf_zero(&y.im);
        rpn_push(&y);
        apf_copy(&x.re, &r);
        apf_zero(&x.im);
        rpn_push(&x);
        return 1;
    }
    
    /* Polar to Rectangular: (r,θ) -> (x,y) */
    if (strcmp(token, ">rec") == 0 || strcmp(token, "rec") == 0) {
        apfc r_val, theta;
        apf x, y, c, s;
        if (rpn_sp < 2) {
            printf("  Error: need r and theta for >rec\n");
            return 1;
        }
        save_lastx();
        rpn_pop(&r_val);  /* r in X */
        rpn_pop(&theta);   /* θ in Y */
        /* x = r * cos(θ), y = r * sin(θ) */
        apfx_cos(&c, &theta.re);
        apfx_sin(&s, &theta.re);
        apf_mul(&x, &r_val.re, &c);
        apf_mul(&y, &r_val.re, &s);
        /* Push y then x (so x is in X, y in Y) */
        {
            apfc tmp;
            apf_copy(&tmp.re, &y);
            apf_zero(&tmp.im);
            rpn_push(&tmp);
            apf_copy(&tmp.re, &x);
            rpn_push(&tmp);
        }
        return 1;
    }
#endif
    
    /* Constants */
#ifdef HAVE_CONST
    if (strcmp(token, "pi") == 0) {
        apfx_pi(&val.re);
        apf_zero(&val.im);
        rpn_push(&val);
        return 1;
    }
    if (strcmp(token, "e") == 0) {
        apfx_e(&val.re);
        apf_zero(&val.im);
        rpn_push(&val);
        return 1;
    }
#endif
    if (strcmp(token, "i") == 0) {
        apf_zero(&val.re);
        apf_from_int(&val.im, 1);
        rpn_push(&val);
        return 1;
    }
    
    /* Single-char operators */
    if (strlen(token) == 1) {
        char c = token[0];
        if (c == '+' || c == '-' || c == '*' || c == '/' || c == '^') {
            return do_binary_op(c);
        }
    }
    
    /* Try as number */
    if (parse_number(token, &val)) {
        rpn_push(&val);
        return 1;
    }
    
    /* Try as unary function */
    if (do_unary_func(token)) {
        return 1;
    }
    
    return 0;
}

void rpn_eval_line(const char *line) {
    char buf[256];
    char *token;
    char *p;
    int len;
    
    /* Copy line to buffer for tokenizing */
    len = (int)strlen(line);
    if (len > 255) len = 255;
    memcpy(buf, line, len);
    buf[len] = '\0';
    
    /* Tokenize by spaces */
    p = buf;
    while (*p) {
        /* Skip whitespace */
        while (*p && isspace((unsigned char)*p)) p++;
        if (!*p) break;
        
        /* Find token end */
        token = p;
        while (*p && !isspace((unsigned char)*p)) p++;
        if (*p) *p++ = '\0';
        
        /* Process token */
        if (!rpn_process(token)) {
            printf("  Unknown: %s\n", token);
        }
    }
    
    /* Show top of stack as result */
    if (rpn_sp > 0) {
        printf("  X: ");
        format_value_rpn(&rpn_stack[rpn_sp - 1]);
        printf("\n");
    }
}

/* Format a value for RPN display (simpler than algebraic) */
void format_value_rpn(const apfc *val) {
    char buf[64];
    int digits = display_digits ? display_digits : 10;  /* Default 10 for stack */
    
    if (apf_is_zero(&val->im)) {
        /* Real number */
        apf_to_str(buf, sizeof(buf), &val->re, digits);
        printf("%s", buf);
    } else if (apf_is_zero(&val->re)) {
        /* Pure imaginary */
        apf_to_str(buf, sizeof(buf), &val->im, digits);
        printf("%si", buf);
    } else {
        /* Complex */
        apf tmp;
        apf_to_str(buf, sizeof(buf), &val->re, digits);
        printf("%s", buf);
        apf_copy(&tmp, &val->im);
        if (apf_is_neg(&tmp)) {
            apf_neg(&tmp, &tmp);
            apf_to_str(buf, sizeof(buf), &tmp, digits);
            printf(" - %si", buf);
        } else {
            apf_to_str(buf, sizeof(buf), &val->im, digits);
            printf(" + %si", buf);
        }
    }
}

#endif /* HAVE_RPN */
