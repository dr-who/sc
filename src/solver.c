/* solver.c - Equation solvers
 * C89 portable for DOS, Linux
 * 
 * Quadratic formula: x = (-b ± √(b²-4ac)) / 2a
 */

#include "config.h"

#ifdef HAVE_SOLVER

#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include "solver.h"
#include "apf.h"
#include "apfc.h"
#include "apfx.h"

/* Forward declaration for formatting */
extern void format_complex_result(const apfc *val);

int solve_quadratic(const apfc *a, const apfc *b, const apfc *c,
                    apfc *x1, apfc *x2) {
    apfc disc, sqrt_disc, neg_b, two_a, tmp1, tmp2;
    apfc four, two;
    
    /* Check for a = 0 (linear equation, not quadratic) */
    if (apfc_is_zero(a)) {
        /* bx + c = 0  =>  x = -c/b */
        if (apfc_is_zero(b)) {
            return 0;  /* No solution or infinite solutions */
        }
        apfc_neg(&tmp1, c);
        apfc_div(x1, &tmp1, b);
        apfc_copy(x2, x1);
        return 1;
    }
    
    /* Calculate discriminant: b² - 4ac */
    apf_from_int(&four.re, 4);
    apf_zero(&four.im);
    
    apfc_mul(&disc, b, b);           /* b² */
    apfc_mul(&tmp1, &four, a);       /* 4a */
    apfc_mul(&tmp2, &tmp1, c);       /* 4ac */
    apfc_sub(&disc, &disc, &tmp2);   /* b² - 4ac */
    
    /* Calculate √(discriminant) - works for complex too */
    apfc_sqrt(&sqrt_disc, &disc);
    
    /* Calculate -b */
    apfc_neg(&neg_b, b);
    
    /* Calculate 2a */
    apf_from_int(&two.re, 2);
    apf_zero(&two.im);
    apfc_mul(&two_a, &two, a);
    
    /* x1 = (-b + √disc) / 2a */
    apfc_add(&tmp1, &neg_b, &sqrt_disc);
    apfc_div(x1, &tmp1, &two_a);
    
    /* x2 = (-b - √disc) / 2a */
    apfc_sub(&tmp1, &neg_b, &sqrt_disc);
    apfc_div(x2, &tmp1, &two_a);
    
    /* Check if roots are equal (discriminant ≈ 0) */
    if (apfc_is_zero(&disc)) {
        return 1;  /* Double root */
    }
    
    return 2;  /* Two distinct roots */
}

/* Parse a number or complex from string */
static int parse_coeff(const char **pp, apfc *val) {
    const char *p = *pp;
    char buf[64];
    int i = 0;
    int has_i = 0;
    
    /* Skip whitespace */
    while (*p && isspace((unsigned char)*p)) p++;
    
    if (!*p) return 0;
    
    /* Handle sign */
    if (*p == '-' || *p == '+') {
        buf[i++] = *p++;
    }
    
    /* Read number part */
    while (*p && (isdigit((unsigned char)*p) || *p == '.')) {
        if (i < 62) buf[i++] = *p;
        p++;
    }
    
    /* Check for 'i' (imaginary) */
    if (*p == 'i') {
        has_i = 1;
        p++;
    }
    
    buf[i] = '\0';
    *pp = p;
    
    if (i == 0 || (i == 1 && (buf[0] == '-' || buf[0] == '+'))) {
        /* Just sign, means 1 or -1 */
        if (has_i) {
            apf_zero(&val->re);
            apf_from_int(&val->im, (buf[0] == '-') ? -1 : 1);
        }
        return has_i;
    }
    
    if (has_i) {
        apf_zero(&val->re);
        apf_from_str(&val->im, buf);
    } else {
        apf_from_str(&val->re, buf);
        apf_zero(&val->im);
    }
    
    return 1;
}

void cmd_quadratic(const char *args) {
    apfc a, b, c, x1, x2;
    const char *p = args;
    int n;
    
    /* Initialize to zero */
    apfc_zero(&a);
    apfc_zero(&b);
    apfc_zero(&c);
    
    /* Parse three coefficients */
    if (!parse_coeff(&p, &a)) {
        printf("Usage: quad a b c\n");
        printf("  Solves ax^2 + bx + c = 0\n");
        printf("Examples:\n");
        printf("  quad 1 -5 6    -> x = 2, 3\n");
        printf("  quad 1 0 1     -> x = i, -i\n");
        printf("  quad 1 2 5     -> x = -1+2i, -1-2i\n");
        return;
    }
    if (!parse_coeff(&p, &b)) {
        printf("Error: need 3 coefficients (a b c)\n");
        return;
    }
    if (!parse_coeff(&p, &c)) {
        printf("Error: need 3 coefficients (a b c)\n");
        return;
    }
    
    /* Solve */
    n = solve_quadratic(&a, &b, &c, &x1, &x2);
    
    if (n == 0) {
        printf("Error: not a valid equation (a=0, b=0)\n");
        return;
    }
    
    /* Display equation */
    {
        char buf[32];
        printf("Solving: ");
        if (!apfc_is_zero(&a)) {
            apf_to_str(buf, sizeof(buf), &a.re, 6);
            printf("%sx^2 ", buf);
        }
        if (!apfc_is_zero(&b)) {
            if (apf_is_neg(&b.re)) {
                apf tmp;
                apf_neg(&tmp, &b.re);
                apf_to_str(buf, sizeof(buf), &tmp, 6);
                printf("- %sx ", buf);
            } else {
                apf_to_str(buf, sizeof(buf), &b.re, 6);
                printf("+ %sx ", buf);
            }
        }
        if (!apfc_is_zero(&c)) {
            if (apf_is_neg(&c.re)) {
                apf tmp;
                apf_neg(&tmp, &c.re);
                apf_to_str(buf, sizeof(buf), &tmp, 6);
                printf("- %s ", buf);
            } else {
                apf_to_str(buf, sizeof(buf), &c.re, 6);
                printf("+ %s ", buf);
            }
        }
        printf("= 0\n");
    }
    
    /* Display roots */
    if (n == 1) {
        printf("Double root:\n");
        printf("  x = ");
        format_complex_result(&x1);
        printf("\n");
    } else {
        printf("Roots:\n");
        printf("  x1 = ");
        format_complex_result(&x1);
        printf("\n");
        printf("  x2 = ");
        format_complex_result(&x2);
        printf("\n");
    }
    
    /* Also show discriminant info */
    {
        apfc disc, four, tmp1, tmp2;
        apf_from_int(&four.re, 4);
        apf_zero(&four.im);
        apfc_mul(&disc, &b, &b);
        apfc_mul(&tmp1, &four, &a);
        apfc_mul(&tmp2, &tmp1, &c);
        apfc_sub(&disc, &disc, &tmp2);
        
        printf("Discriminant (b^2-4ac): ");
        format_complex_result(&disc);
        if (apf_is_zero(&disc.im)) {
            if (apf_is_neg(&disc.re)) {
                printf(" (negative -> complex roots)");
            } else if (apf_is_zero(&disc.re)) {
                printf(" (zero -> double root)");
            } else {
                printf(" (positive -> real roots)");
            }
        }
        printf("\n");
    }
}

void cmd_solve(const char *args) {
    /* Skip whitespace */
    while (*args && isspace((unsigned char)*args)) args++;
    
    /* For now, only quadratic solver */
    /* Could add: linear, cubic, newton-raphson, etc. */
    cmd_quadratic(args);
}

/* Format complex result nicely */
void format_complex_result(const apfc *val) {
    char buf[64];
    
    if (apf_is_zero(&val->im)) {
        /* Pure real */
        apf_to_str(buf, sizeof(buf), &val->re, 10);
        printf("%s", buf);
    } else if (apf_is_zero(&val->re)) {
        /* Pure imaginary */
        apf one;
        apf tmp;
        apf_from_int(&one, 1);
        apf_copy(&tmp, &val->im);
        
        if (apf_is_neg(&tmp)) {
            apf_neg(&tmp, &tmp);
            if (apf_cmp(&tmp, &one) == 0) {
                printf("-i");
            } else {
                apf_to_str(buf, sizeof(buf), &tmp, 10);
                printf("-%si", buf);
            }
        } else {
            if (apf_cmp(&tmp, &one) == 0) {
                printf("i");
            } else {
                apf_to_str(buf, sizeof(buf), &val->im, 10);
                printf("%si", buf);
            }
        }
    } else {
        /* Complex a + bi */
        apf_to_str(buf, sizeof(buf), &val->re, 10);
        printf("%s", buf);
        
        if (apf_is_neg(&val->im)) {
            apf tmp;
            apf_neg(&tmp, &val->im);
            apf_to_str(buf, sizeof(buf), &tmp, 10);
            printf(" - %si", buf);
        } else {
            apf_to_str(buf, sizeof(buf), &val->im, 10);
            printf(" + %si", buf);
        }
    }
}

#endif /* HAVE_SOLVER */
