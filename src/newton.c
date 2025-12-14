/* newton.c - Newton-Raphson root finder
 * C89 portable for DOS, Linux
 * 
 * Solves f(x) = 0 using Newton-Raphson iteration:
 *   x_{n+1} = x_n - f(x_n) / f'(x_n)
 * 
 * Derivative is computed numerically.
 */

#include "config.h"

#ifdef HAVE_NEWTON

#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include "newton.h"
#include "apf.h"
#include "apfc.h"
#include "sc.h"

/* External display digits */
extern int display_digits;

/* External parser interface */
extern const char *input_ptr;
extern token_t current_token;
extern void next_token(void);

/* Evaluate user function at x */
static int eval_func(int func_idx, const apf *x, apf *result) {
    apfc cx;
    int param_idx;
    apfc save_val;
    int save_defined;
    const char *save_input;
    
    if (func_idx < 0 || func_idx >= MAX_FUNCTIONS) return 0;
    if (!user_funcs[func_idx].defined) return 0;
    
    /* Save the parameter variable */
    param_idx = user_funcs[func_idx].param - 'a';
    if (param_idx < 0 || param_idx >= MAX_SCALAR_VARS) return 0;
    
    save_defined = scalar_vars[param_idx].defined;
    if (save_defined) {
        apfc_copy(&save_val, &scalar_vars[param_idx].val);
    }
    
    /* Set parameter to x */
    apf_copy(&cx.re, x);
    apf_zero(&cx.im);
    scalar_vars[param_idx].defined = 1;
    apfc_copy(&scalar_vars[param_idx].val, &cx);
    
    /* Save and set up parser state */
    save_input = input_ptr;
    input_ptr = user_funcs[func_idx].body;
    next_token();
    
    /* Evaluate */
    {
        apfc expr_result;
        if (parse_expr(&expr_result)) {
            apf_copy(result, &expr_result.re);
        } else {
            /* Restore and fail */
            input_ptr = save_input;
            if (save_defined) {
                apfc_copy(&scalar_vars[param_idx].val, &save_val);
            } else {
                scalar_vars[param_idx].defined = 0;
            }
            return 0;
        }
    }
    
    /* Restore parser state */
    input_ptr = save_input;
    
    /* Restore parameter */
    if (save_defined) {
        apfc_copy(&scalar_vars[param_idx].val, &save_val);
    } else {
        scalar_vars[param_idx].defined = 0;
    }
    
    return 1;
}

int newton_solve(int func_idx, const apf *x0, apf *result) {
    apf x, x_new, fx, fx_plus, df, delta;
    apf h, tolerance, abs_fx, abs_delta;
    int iter;
    
    apf_copy(&x, x0);
    
    /* h for numerical derivative */
    apf_from_str(&h, "1e-10");
    
    /* Convergence tolerance */
    apf_from_str(&tolerance, "1e-14");
    
    for (iter = 0; iter < NEWTON_MAX_ITER; iter++) {
        /* Evaluate f(x) */
        if (!eval_func(func_idx, &x, &fx)) {
            return NEWTON_NO_CONVERGE;
        }
        
        /* Check if we're close enough */
        apf_abs(&abs_fx, &fx);
        if (apf_cmp(&abs_fx, &tolerance) < 0) {
            apf_copy(result, &x);
            return NEWTON_OK;
        }
        
        /* Numerical derivative: f'(x) ≈ (f(x+h) - f(x)) / h */
        {
            apf x_plus;
            apf_add(&x_plus, &x, &h);
            if (!eval_func(func_idx, &x_plus, &fx_plus)) {
                return NEWTON_NO_CONVERGE;
            }
            apf_sub(&df, &fx_plus, &fx);
            apf_div(&df, &df, &h);
        }
        
        /* Check for zero derivative */
        {
            apf abs_df, small;
            apf_abs(&abs_df, &df);
            apf_from_str(&small, "1e-30");
            if (apf_cmp(&abs_df, &small) < 0) {
                return NEWTON_ZERO_DERIV;
            }
        }
        
        /* x_new = x - f(x) / f'(x) */
        apf_div(&delta, &fx, &df);
        apf_sub(&x_new, &x, &delta);
        
        /* Check for convergence */
        apf_abs(&abs_delta, &delta);
        if (apf_cmp(&abs_delta, &tolerance) < 0) {
            apf_copy(result, &x_new);
            return NEWTON_OK;
        }
        
        apf_copy(&x, &x_new);
    }
    
    /* Didn't converge but return best guess */
    apf_copy(result, &x);
    return NEWTON_NO_CONVERGE;
}

void cmd_newton(const char *args) {
    char buf[64];
    apf x0, result;
    int func_idx;
    int ret;
    const char *p = args;
    int digits = display_digits ? display_digits : 10;
    
    /* Skip whitespace */
    while (*p && isspace((unsigned char)*p)) p++;
    
    if (*p == '\0') {
        printf("Usage: solve f [x0]\n");
        printf("  Solve f(x) = 0 using Newton-Raphson\n");
        printf("  f   - function name (a-z)\n");
        printf("  x0  - initial guess (default: 1)\n");
        printf("\nExample:\n");
        printf("  f(x) = x^2 - 2\n");
        printf("  solve f 1      -> x = 1.414... (sqrt(2))\n");
        return;
    }
    
    /* Get function name */
    if (*p >= 'a' && *p <= 'z') {
        func_idx = *p - 'a';
    } else if (*p >= 'A' && *p <= 'Z') {
        func_idx = *p - 'A';
    } else {
        printf("Error: expected function name (a-z)\n");
        return;
    }
    p++;
    
    if (!user_funcs[func_idx].defined) {
        printf("Error: function '%c' not defined\n", 'a' + func_idx);
        return;
    }
    
    /* Skip whitespace */
    while (*p && isspace((unsigned char)*p)) p++;
    
    /* Get initial guess (default 1) */
    if (*p) {
        apf_from_str(&x0, p);
    } else {
        apf_from_int(&x0, 1);
    }
    
    /* Solve */
    ret = newton_solve(func_idx, &x0, &result);
    
    if (ret == NEWTON_OK) {
        apf_to_str(buf, sizeof(buf), &result, digits);
        printf("x = %s\n", buf);
        
        /* Verify */
        {
            apf check;
            if (eval_func(func_idx, &result, &check)) {
                apf_to_str(buf, sizeof(buf), &check, digits);
                printf("%c(%s) = %s\n", 'a' + func_idx, 
                       apf_to_str(buf, sizeof(buf), &result, digits), buf);
            }
        }
    } else if (ret == NEWTON_ZERO_DERIV) {
        printf("Error: zero derivative encountered\n");
    } else {
        apf_to_str(buf, sizeof(buf), &result, digits);
        printf("Warning: may not have converged. Best guess: x = %s\n", buf);
    }
}

#endif /* HAVE_NEWTON */
