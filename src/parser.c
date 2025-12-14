/* parser.c - Expression parser with built-in functions
 * Recursive descent parser with precedence climbing
 * C89 compliant for Watcom C / DOS
 * 16-bit clean: int is 16-bit, long is 32-bit
 */
#include "sc.h"

/* External runtime accessors */
extern int is_in_user_func(void);
extern char get_param_name(void);
extern apfc *get_param_value(void);

/*
 * Shared static temporaries for parse_value_* functions.
 * These are static to avoid stack overflow on DOS (value_t is ~1.5KB).
 * Using a shared pool keeps BSS under the 64KB DGROUP limit.
 */
static value_t pv_arg;        /* Shared function argument */
static matrix_t pv_tmp_mat;   /* Shared temporary matrix */
static matrix_t pv_tmp_mat2;  /* Second shared temporary matrix */

/* Forward declarations for scalar parsing */
static int parse_term(apfc *result);
static int parse_unary(apfc *result);
static int parse_power(apfc *result);
static int parse_postfix(apfc *result);
static int parse_factor(apfc *result);

/* Forward declarations for value parsing */
static int parse_value_term(value_t *result);
static int parse_value_unary(value_t *result);
static int parse_value_power(value_t *result);
static int parse_value_postfix(value_t *result);
static int parse_value_factor(value_t *result);
static int parse_value_expr(value_t *result);

/* ========== Scalar Expression Parser ========== */

int parse_expr(apfc *result)
{
    apfc right, temp;
    int op;

    if (!parse_term(result)) return 0;

    while (current_token.type == TOK_PLUS || current_token.type == TOK_MINUS) {
        op = current_token.type;
        next_token();
        if (!parse_term(&right)) return 0;

        if (op == TOK_PLUS) {
            apfc_add(&temp, result, &right);
        } else {
            apfc_sub(&temp, result, &right);
        }
        *result = temp;
    }
    return 1;
}

static int parse_term(apfc *result)
{
    apfc right, temp;
    int op;

    if (!parse_unary(result)) return 0;

    while (current_token.type == TOK_MUL || current_token.type == TOK_DIV ||
           current_token.type == TOK_BACKSLASH) {
        op = current_token.type;
        next_token();
        if (!parse_unary(&right)) return 0;

        if (op == TOK_MUL) {
            apfc_mul(&temp, result, &right);
        } else if (op == TOK_DIV) {
            apfc_div(&temp, result, &right);
        } else {
            /* Backslash for scalars: a\b = b/a */
            apfc_div(&temp, &right, result);
        }
        *result = temp;
    }
    return 1;
}

/* Unary minus/plus - applied AFTER exponentiation */
static int parse_unary(apfc *result)
{
    int neg = 0;

    while (current_token.type == TOK_MINUS) {
        neg = !neg;
        next_token();
    }
    while (current_token.type == TOK_PLUS) {
        next_token();
    }

    if (!parse_power(result)) return 0;

    if (neg) {
        apfc_neg(result, result);
    }
    return 1;
}

/* Power operator - right associative */
static int parse_power(apfc *result)
{
    apfc right, temp;

    if (!parse_postfix(result)) return 0;

    if (current_token.type == TOK_POW) {
        next_token();
        if (!parse_unary(&right)) return 0;
        apfc_pow(&temp, result, &right);
        *result = temp;
    }
    return 1;
}

/* Postfix factorial */
static int parse_postfix(apfc *result)
{
    long n;
    apfc temp;

    if (!parse_factor(result)) return 0;

    while (current_token.type == TOK_FACT) {
        next_token();
        if (!apfc_is_real(result)) {
            printf("Error: factorial requires real integer\n");
            return 0;
        }
        n = apf_to_long(&result->re);
        if (n < 0 || n > 1000) {
            printf("Error: factorial argument out of range\n");
            return 0;
        }
        apfx_fact(&temp.re, n);
        apf_from_int(&temp.im, 0);
        *result = temp;
    }
    return 1;
}

/* Parse random functions with various argument patterns */
static int parse_random_func(apfc *result, const char *name)
{
    if (current_token.type != TOK_LPAREN) {
        /* rand without parens - single scalar */
        if (str_eq(name, "rand")) {
            rand_uniform(result);
            return 1;
        } else if (str_eq(name, "randn")) {
            rand_normal(result);
            return 1;
        }
        printf("Error: expected '(' after '%s'\n", name);
        return 0;
    }
    
    next_token();
    
    /* Check for empty parens: rand() */
    if (current_token.type == TOK_RPAREN) {
        next_token();
        if (str_eq(name, "rand")) {
            rand_uniform(result);
        } else if (str_eq(name, "randn")) {
            rand_normal(result);
        } else if (str_eq(name, "randi")) {
            rand_int(result, 1);
        }
        return 1;
    }
    
    /* Parse first argument */
    {
        apfc arg1;
        if (!parse_expr(&arg1)) return 0;
        
        /* Check for comma (multiple args) */
        if (current_token.type == TOK_COMMA) {
            apfc arg2;
            next_token();
            if (!parse_expr(&arg2)) return 0;
            
            if (str_eq(name, "rand")) {
                if (current_token.type != TOK_RPAREN) {
                    printf("Error: expected ')'\n");
                    return 0;
                }
                next_token();
                rand_uniform_range(result, &arg1, &arg2);
            } else if (str_eq(name, "randi")) {
                if (current_token.type != TOK_RPAREN) {
                    printf("Error: expected ')'\n");
                    return 0;
                }
                next_token();
                rand_int_range(result, apf_to_long(&arg1.re), apf_to_long(&arg2.re));
            } else {
                printf("Error: %s doesn't support two arguments\n", name);
                return 0;
            }
        } else if (current_token.type == TOK_RPAREN) {
            next_token();
            if (str_eq(name, "randi")) {
                rand_int(result, apf_to_long(&arg1.re));
            } else {
                printf("Error: %s requires 0 or 2 arguments\n", name);
                return 0;
            }
        } else {
            printf("Error: expected ')' or ','\n");
            return 0;
        }
    }
    
    return 1;
}

/* Factor: numbers, constants, functions, parentheses */
static int parse_factor(apfc *result)
{
    apfc arg;

    /* Parentheses */
    if (current_token.type == TOK_LPAREN) {
        next_token();
        if (!parse_expr(result)) return 0;
        if (current_token.type != TOK_RPAREN) {
            printf("Error: expected ')'\n");
            return 0;
        }
        next_token();
        return 1;
    }

    /* Pure imaginary number (e.g., 3i) */
    if (current_token.type == TOK_IMAG) {
        int i_idx = 'i' - 'a';
        if (scalar_vars[i_idx].defined) {
            *result = scalar_vars[i_idx].val;
            next_token();
            return 1;
        }
        apf_zero(&result->re);
        apf_copy(&result->im, &current_token.value);
        next_token();
        return 1;
    }

    /* Real number */
    if (current_token.type == TOK_NUM) {
        apf_copy(&result->re, &current_token.value);
        apf_zero(&result->im);
        next_token();
        return 1;
    }

    /* Previous answer (ans) */
    if (current_token.type == TOK_ANS) {
        if (!last_ans_valid) {
            printf("Error: no previous answer\n");
            return 0;
        }
        if (last_ans.type == VAL_SCALAR) {
            *result = last_ans.v.scalar;
        } else {
            /* Matrix ans - check if 1x1 */
            if (last_ans.v.matrix.rows == 1 && last_ans.v.matrix.cols == 1) {
                *result = MAT_AT(&last_ans.v.matrix, 0, 0);
            } else {
                printf("Error: ans is a matrix, use in matrix context\n");
                return 0;
            }
        }
        next_token();
        return 1;
    }

    /* Function or constant */
    if (current_token.type == TOK_FUNC) {
        char name[16];
        int func_idx, var_idx;
        strcpy(name, current_token.func_name);
        next_token();

        /* Parameter in user function */
        if (is_in_user_func() && name[0] == get_param_name() && name[1] == '\0') {
            *result = *get_param_value();
            return 1;
        }

        /* Variable lookup - scalar only in scalar context */
        var_idx = get_var_index(name);
        if (var_idx >= 0 && current_token.type != TOK_LPAREN) {
            if (scalar_vars[var_idx].defined) {
                *result = scalar_vars[var_idx].val;
                return 1;
            } else if (is_var_matrix(var_idx)) {
                /* Check if it's a 1x1 matrix */
                int mi;
                for (mi = 0; mi < MAX_MATRIX_VARS; mi++) {
                    if (matrix_vars[mi].defined && matrix_vars[mi].name == 'a' + var_idx) {
                        if (matrix_vars[mi].val.rows == 1 && matrix_vars[mi].val.cols == 1) {
                            *result = MAT_AT(&matrix_vars[mi].val, 0, 0);
                            return 1;
                        }
                        printf("Error: matrix '%s' used in scalar context\n", name);
                        return 0;
                    }
                }
            }
        }

        /* Constants */
        if (str_eq(name, "pi")) {
            apfx_pi(&result->re);
            apf_zero(&result->im);
            return 1;
        }
        if (str_eq(name, "e")) {
            apfx_e(&result->re);
            apf_zero(&result->im);
            return 1;
        }
        if (str_eq(name, "i")) {
            apf_zero(&result->re);
            apf_from_int(&result->im, 1);
            return 1;
        }
        if (str_eq(name, "Inf") || str_eq(name, "inf") || str_eq(name, "INF")) {
            apf_set_inf(&result->re, 0);
            apf_zero(&result->im);
            return 1;
        }
        if (str_eq(name, "NaN") || str_eq(name, "nan") || str_eq(name, "NAN")) {
            apf_set_nan(&result->re);
            apf_zero(&result->im);
            return 1;
        }

        /* Random functions */
        if (str_eq(name, "rand") || str_eq(name, "randn") || str_eq(name, "randi")) {
            return parse_random_func(result, name);
        }

        /* User function */
        func_idx = get_func_index(name);
        if (func_idx >= 0 && user_funcs[func_idx].defined) {
            if (current_token.type != TOK_LPAREN) {
                printf("Error: expected '(' after function '%s'\n", name);
                return 0;
            }
            next_token();
            if (!parse_expr(&arg)) return 0;
            if (current_token.type != TOK_RPAREN) {
                printf("Error: expected ')'\n");
                return 0;
            }
            next_token();
            return eval_user_func(result, func_idx, &arg);
        }

        /* Matrix functions that return scalar */
        if (str_eq(name, "det") || str_eq(name, "trace") || str_eq(name, "tr") ||
            str_eq(name, "norm") || str_eq(name, "size") || str_eq(name, "length") ||
            str_eq(name, "rows") || str_eq(name, "cols") ||
            str_eq(name, "sum") || str_eq(name, "mean") ||
            str_eq(name, "min") || str_eq(name, "max") ||
            str_eq(name, "median") || str_eq(name, "sd") || str_eq(name, "std")) {
            
            /* Uses shared static: pv_arg */
            
            if (current_token.type != TOK_LPAREN) {
                printf("Error: expected '(' after function '%s'\n", name);
                return 0;
            }
            next_token();
            if (!parse_value(&pv_arg)) return 0;
            if (current_token.type != TOK_RPAREN) {
                printf("Error: expected ')'\n");
                return 0;
            }
            next_token();
            
            if (pv_arg.type == VAL_SCALAR) {
                if (str_eq(name, "det") || str_eq(name, "trace") || str_eq(name, "tr") ||
                    str_eq(name, "sum") || str_eq(name, "mean") ||
                    str_eq(name, "min") || str_eq(name, "max") || str_eq(name, "median")) {
                    *result = pv_arg.v.scalar;
                } else if (str_eq(name, "sd") || str_eq(name, "std")) {
                    /* SD of single value is 0 */
                    apf_zero(&result->re);
                    apf_zero(&result->im);
                } else if (str_eq(name, "norm")) {
                    apfc_abs(&result->re, &pv_arg.v.scalar);
                    apf_zero(&result->im);
                } else {
                    apf_from_int(&result->re, 1);
                    apf_zero(&result->im);
                }
            } else {
                matrix_t *m = &pv_arg.v.matrix;
                
                if (str_eq(name, "det")) {
                    mat_det(result, m);
                } else if (str_eq(name, "trace") || str_eq(name, "tr")) {
                    mat_trace(result, m);
                } else if (str_eq(name, "norm")) {
                    mat_norm_frobenius(result, m);
                } else if (str_eq(name, "size")) {
                    apf_from_int(&result->re, m->rows * m->cols);
                    apf_zero(&result->im);
                } else if (str_eq(name, "length")) {
                    apf_from_int(&result->re, mat_length(m));
                    apf_zero(&result->im);
                } else if (str_eq(name, "rows")) {
                    apf_from_int(&result->re, m->rows);
                    apf_zero(&result->im);
                } else if (str_eq(name, "cols")) {
                    apf_from_int(&result->re, m->cols);
                    apf_zero(&result->im);
                } else if (str_eq(name, "sum")) {
                    mat_sum(result, m);
                } else if (str_eq(name, "mean")) {
                    mat_mean(result, m);
                } else if (str_eq(name, "min")) {
                    mat_min(result, m);
                } else if (str_eq(name, "max")) {
                    mat_max(result, m);
                } else if (str_eq(name, "median")) {
                    mat_median(result, m);
                } else if (str_eq(name, "sd") || str_eq(name, "std")) {
                    mat_std(result, m);
                }
            }
            return 1;
        }

        /* eig function */
        if (str_eq(name, "eig") || str_eq(name, "eigenvalues")) {
            /* Uses shared statics: pv_arg, pv_tmp_mat */
            
            if (current_token.type != TOK_LPAREN) {
                printf("Error: expected '(' after '%s'\n", name);
                return 0;
            }
            next_token();
            if (!parse_value(&pv_arg)) return 0;
            if (current_token.type != TOK_RPAREN) {
                printf("Error: expected ')'\n");
                return 0;
            }
            next_token();
            
            if (pv_arg.type == VAL_SCALAR) {
                *result = pv_arg.v.scalar;
            } else {
                matrix_t *m = &pv_arg.v.matrix;
                int i;
                int n;
                static char buf[128];
                
                n = mat_eig(&pv_tmp_mat, m);
                if (n == 0) {
                    return 0;
                }
                
                /* Print all eigenvalues */
                printf("eigenvalues: ");
                for (i = 0; i < n; i++) {
                    apfc_to_str(buf, sizeof(buf), &MAT_AT(&pv_tmp_mat, i, 0), 0);
                    printf("%s", buf);
                    if (i < n - 1) {
                        printf(", ");
                    }
                }
                printf("\n");
                
                /* Return the first eigenvalue as scalar result */
                *result = MAT_AT(&pv_tmp_mat, 0, 0);
            }
            return 1;
        }

        /* Built-in functions require '(' */
        if (current_token.type != TOK_LPAREN) {
            printf("Error: expected '(' after function '%s'\n", name);
            return 0;
        }
        next_token();
        if (!parse_expr(&arg)) return 0;
        
        /* Check for two-argument functions BEFORE consuming ')' */
#ifdef HAVE_GCD
        if (str_eq(name, "gcd")) {
            apfc arg2;
            if (current_token.type != TOK_COMMA) {
                printf("Error: gcd requires two arguments\n");
                return 0;
            }
            next_token();
            if (!parse_expr(&arg2)) return 0;
            if (current_token.type != TOK_RPAREN) {
                printf("Error: expected ')'\n");
                return 0;
            }
            next_token();
            apf_gcd(&result->re, &arg.re, &arg2.re);
            apf_zero(&result->im);
            return 1;
        }
        if (str_eq(name, "lcm")) {
            apfc arg2;
            if (current_token.type != TOK_COMMA) {
                printf("Error: lcm requires two arguments\n");
                return 0;
            }
            next_token();
            if (!parse_expr(&arg2)) return 0;
            if (current_token.type != TOK_RPAREN) {
                printf("Error: expected ')'\n");
                return 0;
            }
            next_token();
            apf_lcm(&result->re, &arg.re, &arg2.re);
            apf_zero(&result->im);
            return 1;
        }
#endif
#ifdef HAVE_COMB
        if (str_eq(name, "npr") || str_eq(name, "perm")) {
            apfc arg2;
            if (current_token.type != TOK_COMMA) {
                printf("Error: nPr requires two arguments\n");
                return 0;
            }
            next_token();
            if (!parse_expr(&arg2)) return 0;
            if (current_token.type != TOK_RPAREN) {
                printf("Error: expected ')'\n");
                return 0;
            }
            next_token();
            apf_npr(&result->re, apf_to_long(&arg.re), apf_to_long(&arg2.re));
            apf_zero(&result->im);
            return 1;
        }
        if (str_eq(name, "ncr") || str_eq(name, "comb") || str_eq(name, "choose")) {
            apfc arg2;
            if (current_token.type != TOK_COMMA) {
                printf("Error: nCr requires two arguments\n");
                return 0;
            }
            next_token();
            if (!parse_expr(&arg2)) return 0;
            if (current_token.type != TOK_RPAREN) {
                printf("Error: expected ')'\n");
                return 0;
            }
            next_token();
            apf_ncr(&result->re, apf_to_long(&arg.re), apf_to_long(&arg2.re));
            apf_zero(&result->im);
            return 1;
        }
#endif
#ifdef HAVE_BITWISE
        if (str_eq(name, "and") || str_eq(name, "band")) {
            apfc arg2;
            if (current_token.type != TOK_COMMA) {
                printf("Error: and requires two arguments\n");
                return 0;
            }
            next_token();
            if (!parse_expr(&arg2)) return 0;
            if (current_token.type != TOK_RPAREN) {
                printf("Error: expected ')'\n");
                return 0;
            }
            next_token();
            apf_and(&result->re, &arg.re, &arg2.re);
            apf_zero(&result->im);
            return 1;
        }
        if (str_eq(name, "or") || str_eq(name, "bor")) {
            apfc arg2;
            if (current_token.type != TOK_COMMA) {
                printf("Error: or requires two arguments\n");
                return 0;
            }
            next_token();
            if (!parse_expr(&arg2)) return 0;
            if (current_token.type != TOK_RPAREN) {
                printf("Error: expected ')'\n");
                return 0;
            }
            next_token();
            apf_or(&result->re, &arg.re, &arg2.re);
            apf_zero(&result->im);
            return 1;
        }
        if (str_eq(name, "xor") || str_eq(name, "bxor")) {
            apfc arg2;
            if (current_token.type != TOK_COMMA) {
                printf("Error: xor requires two arguments\n");
                return 0;
            }
            next_token();
            if (!parse_expr(&arg2)) return 0;
            if (current_token.type != TOK_RPAREN) {
                printf("Error: expected ')'\n");
                return 0;
            }
            next_token();
            apf_xor(&result->re, &arg.re, &arg2.re);
            apf_zero(&result->im);
            return 1;
        }
        if (str_eq(name, "lsl") || str_eq(name, "shl")) {
            apfc arg2;
            if (current_token.type != TOK_COMMA) {
                printf("Error: lsl requires two arguments\n");
                return 0;
            }
            next_token();
            if (!parse_expr(&arg2)) return 0;
            if (current_token.type != TOK_RPAREN) {
                printf("Error: expected ')'\n");
                return 0;
            }
            next_token();
            apf_lsl(&result->re, &arg.re, (int)apf_to_long(&arg2.re));
            apf_zero(&result->im);
            return 1;
        }
        if (str_eq(name, "lsr") || str_eq(name, "shr")) {
            apfc arg2;
            if (current_token.type != TOK_COMMA) {
                printf("Error: lsr requires two arguments\n");
                return 0;
            }
            next_token();
            if (!parse_expr(&arg2)) return 0;
            if (current_token.type != TOK_RPAREN) {
                printf("Error: expected ')'\n");
                return 0;
            }
            next_token();
            apf_lsr(&result->re, &arg.re, (int)apf_to_long(&arg2.re));
            apf_zero(&result->im);
            return 1;
        }
#endif

        /* Single-argument functions - require ')' */
        if (current_token.type != TOK_RPAREN) {
            printf("Error: expected ')'\n");
            return 0;
        }
        next_token();

        /* Built-in functions */
        if (str_eq(name, "sin")) {
            apfc_sin(result, &arg);
        } else if (str_eq(name, "cos")) {
            apfc_cos(result, &arg);
        } else if (str_eq(name, "tan")) {
            apfc_tan(result, &arg);
        } else if (str_eq(name, "exp")) {
            apfc_exp(result, &arg);
        } else if (str_eq(name, "log") || str_eq(name, "ln")) {
            apfc_log(result, &arg);
        } else if (str_eq(name, "log10")) {
            apfc log_arg, ln10;
            apf ten;
            apfc_log(&log_arg, &arg);
            apf_from_int(&ten, 10);
            apfc_from_real(&ln10, &ten);
            apfc_log(&ln10, &ln10);
            apfc_div(result, &log_arg, &ln10);
        } else if (str_eq(name, "sqrt")) {
            apfc_sqrt(result, &arg);
        } else if (str_eq(name, "roots2")) {
            static char buf1[256], buf2[256];
            apfc r1, r2;
            apfc_sqrt(&r1, &arg);
            apfc_neg(&r2, &r1);
            apfc_to_str(buf1, sizeof(buf1), &r1, 0);
            apfc_to_str(buf2, sizeof(buf2), &r2, 0);
            printf("roots: %s, %s\n", buf1, buf2);
            *result = r1;
        } else if (str_eq(name, "abs")) {
            apfc_abs(&result->re, &arg);
            apf_zero(&result->im);
        } else if (str_eq(name, "arg")) {
            apfc_arg(&result->re, &arg);
            apf_zero(&result->im);
        } else if (str_eq(name, "conj")) {
            apfc_conj(result, &arg);
        } else if (str_eq(name, "re") || str_eq(name, "real")) {
            result->re = arg.re;
            apf_zero(&result->im);
        } else if (str_eq(name, "im") || str_eq(name, "imag")) {
            result->re = arg.im;
            apf_zero(&result->im);
        } else if (str_eq(name, "sinh")) {
            apfc_sinh(result, &arg);
        } else if (str_eq(name, "cosh")) {
            apfc_cosh(result, &arg);
        } else if (str_eq(name, "atan") || str_eq(name, "arctan")) {
            /* atan for real numbers only (complex atan not implemented)
             * But allow NaN through - should return NaN */
            if (!apfc_is_real(&arg) && arg.re.cls != APF_CLASS_NAN) {
                printf("Error: atan requires real argument\n");
                return 0;
            }
            apfx_atan(&result->re, &arg.re);
            apf_zero(&result->im);
        } else if (str_eq(name, "asin") || str_eq(name, "arcsin")) {
            /* asin for real numbers only, but allow NaN */
            if (!apfc_is_real(&arg) && arg.re.cls != APF_CLASS_NAN) {
                printf("Error: asin requires real argument\n");
                return 0;
            }
            apfx_asin(&result->re, &arg.re);
            apf_zero(&result->im);
        } else if (str_eq(name, "acos") || str_eq(name, "arccos")) {
            /* acos for real numbers only, but allow NaN */
            if (!apfc_is_real(&arg) && arg.re.cls != APF_CLASS_NAN) {
                printf("Error: acos requires real argument\n");
                return 0;
            }
            apfx_acos(&result->re, &arg.re);
            apf_zero(&result->im);
        } else if (str_eq(name, "tanh")) {
            /* tanh for real numbers only, but allow NaN */
            if (!apfc_is_real(&arg) && arg.re.cls != APF_CLASS_NAN) {
                printf("Error: tanh requires real argument\n");
                return 0;
            }
            apfx_tanh(&result->re, &arg.re);
            apf_zero(&result->im);
        } else if (str_eq(name, "fact")) {
            long n;
            if (!apfc_is_real(&arg)) {
                printf("Error: factorial requires real integer\n");
                return 0;
            }
            n = apf_to_long(&arg.re);
            if (n < 0 || n > 1000) {
                printf("Error: factorial out of range\n");
                return 0;
            }
            apfx_fact(&result->re, n);
            apf_zero(&result->im);
        } else if (str_eq(name, "print") || str_eq(name, "printdec")) {
            static char pbuf[512];
            apfc_to_str(pbuf, sizeof(pbuf), &arg, 0);
            printf("%s\n", pbuf);
            *result = arg;
        } else if (str_eq(name, "printhex")) {
            apf_to_hex_str(&arg.re);
            *result = arg;
        } else if (str_eq(name, "printbin")) {
            apf_to_bin_str(&arg.re);
            *result = arg;
#ifdef HAVE_GCD
        } else if (str_eq(name, "isprime") || str_eq(name, "prime")) {
            long n = apf_to_long(&arg.re);
            apf_from_int(&result->re, is_prime_long(n));
            apf_zero(&result->im);
#endif
#ifdef HAVE_BITWISE
        } else if (str_eq(name, "not") || str_eq(name, "bnot")) {
            apf_not(&result->re, &arg.re);
            apf_zero(&result->im);
#endif
        } else {
            printf("Error: unknown function '%s'\n", name);
            return 0;
        }
        return 1;
    }

    /* Matrix literal */
    if (current_token.type == TOK_LBRACKET) {
        matrix_t mat;
        if (parse_matrix(&mat)) {
            if (mat.rows == 1 && mat.cols == 1) {
                *result = MAT_AT(&mat, 0, 0);
            } else {
                printf("Error: matrix in scalar expression\n");
                return 0;
            }
            return 1;
        }
        return 0;
    }

    if (current_token.type == TOK_END) {
        printf("Error: unexpected end of expression\n");
    } else {
        printf("Error: unexpected token\n");
    }
    return 0;
}

/* ========== Matrix Parsing ========== */

int parse_matrix(matrix_t *result)
{
    int row = 0, col = 0;
    int cols_per_row = -1;
    apfc elem;
    apfc temp_data[MAT_MAX_ELEM];
    int temp_idx = 0;
    
    if (current_token.type != TOK_LBRACKET) {
        printf("Error: expected '['\n");
        return 0;
    }
    next_token();
    
    while (current_token.type != TOK_RBRACKET && current_token.type != TOK_END) {
        if (!parse_expr(&elem)) {
            return 0;
        }
        
        if (temp_idx >= MAT_MAX_ELEM) {
            printf("Error: matrix too large\n");
            return 0;
        }
        
        temp_data[temp_idx++] = elem;
        col++;
        
        if (current_token.type == TOK_SEMI) {
            if (cols_per_row < 0) {
                cols_per_row = col;
            } else if (col != cols_per_row) {
                printf("Error: inconsistent row lengths\n");
                return 0;
            }
            row++;
            col = 0;
            next_token();
        } else if (current_token.type == TOK_COMMA) {
            next_token();
        }
    }
    
    if (current_token.type != TOK_RBRACKET) {
        printf("Error: expected ']'\n");
        return 0;
    }
    next_token();
    
    if (col > 0) {
        if (cols_per_row < 0) {
            cols_per_row = col;
        } else if (col != cols_per_row) {
            printf("Error: inconsistent row lengths\n");
            return 0;
        }
        row++;
    }
    
    result->rows = row;
    result->cols = cols_per_row > 0 ? cols_per_row : 1;
    
    {
        int r, c, src_idx = 0;
        for (r = 0; r < result->rows; r++) {
            for (c = 0; c < result->cols; c++) {
                MAT_AT(result, r, c) = temp_data[src_idx++];
            }
        }
    }
    
    return 1;
}

/* ========== Value Parser (scalar or matrix) ========== */

/* Parse matrix-creating functions */
static int parse_matrix_func(value_t *result, const char *name)
{
    apfc arg1, arg2, arg3;
    int rows, cols;
    long imax, imin;
    
    if (current_token.type != TOK_LPAREN) {
        printf("Error: expected '(' after '%s'\n", name);
        return 0;
    }
    next_token();
    
    if (!parse_expr(&arg1)) return 0;
    rows = apf_to_long(&arg1.re);
    
    if (current_token.type == TOK_COMMA) {
        next_token();
        if (!parse_expr(&arg2)) return 0;
        cols = apf_to_long(&arg2.re);
        
        if (str_eq(name, "rand")) {
            if (current_token.type != TOK_RPAREN) {
                printf("Error: expected ')'\n");
                return 0;
            }
            next_token();
            result->type = VAL_MATRIX;
            mat_rand(&result->v.matrix, rows, cols);
            return 1;
        } else if (str_eq(name, "randn")) {
            if (current_token.type != TOK_RPAREN) {
                printf("Error: expected ')'\n");
                return 0;
            }
            next_token();
            result->type = VAL_MATRIX;
            mat_randn(&result->v.matrix, rows, cols);
            return 1;
        } else if (str_eq(name, "randi")) {
            imax = rows;
            rows = cols;
            
            if (current_token.type == TOK_COMMA) {
                next_token();
                if (!parse_expr(&arg3)) return 0;
                
                if (current_token.type == TOK_COMMA) {
                    apfc arg4;
                    imin = imax;
                    imax = rows;
                    rows = apf_to_long(&arg3.re);
                    next_token();
                    if (!parse_expr(&arg4)) return 0;
                    cols = apf_to_long(&arg4.re);
                } else {
                    imin = 1;
                    cols = apf_to_long(&arg3.re);
                }
            } else {
                imin = 1;
                cols = rows;
            }
            
            if (current_token.type != TOK_RPAREN) {
                printf("Error: expected ')'\n");
                return 0;
            }
            next_token();
            result->type = VAL_MATRIX;
            mat_randi_range(&result->v.matrix, rows, cols, imin, imax);
            return 1;
        } else if (str_eq(name, "zeros")) {
            if (current_token.type != TOK_RPAREN) {
                printf("Error: expected ')'\n");
                return 0;
            }
            next_token();
            result->type = VAL_MATRIX;
            mat_zero(&result->v.matrix, rows, cols);
            return 1;
        } else if (str_eq(name, "ones")) {
            if (current_token.type != TOK_RPAREN) {
                printf("Error: expected ')'\n");
                return 0;
            }
            next_token();
            result->type = VAL_MATRIX;
            mat_ones(&result->v.matrix, rows, cols);
            return 1;
        }
    }
    
    if (current_token.type != TOK_RPAREN) {
        printf("Error: expected ')'\n");
        return 0;
    }
    next_token();
    
    result->type = VAL_MATRIX;
    if (str_eq(name, "zeros")) {
        mat_zero(&result->v.matrix, rows, rows);
    } else if (str_eq(name, "ones")) {
        mat_ones(&result->v.matrix, rows, rows);
    } else if (str_eq(name, "eye") || str_eq(name, "identity")) {
        mat_identity(&result->v.matrix, rows);
    } else if (str_eq(name, "rand")) {
        mat_rand(&result->v.matrix, rows, rows);
    } else if (str_eq(name, "randn")) {
        mat_randn(&result->v.matrix, rows, rows);
    } else if (str_eq(name, "randi")) {
        result->type = VAL_SCALAR;
        rand_int(&result->v.scalar, rows);
    }
    
    return 1;
}

/* Parse value factor (base element) */
static int parse_value_factor(value_t *result)
{
    /* Uses shared statics: pv_tmp_mat, pv_arg */
    
    if (current_token.type == TOK_LPAREN) {
        next_token();
        if (!parse_value_expr(result)) return 0;
        if (current_token.type != TOK_RPAREN) {
            printf("Error: expected ')'\n");
            return 0;
        }
        next_token();
        return 1;
    }
    
    if (current_token.type == TOK_LBRACKET) {
        if (!parse_matrix(&pv_tmp_mat)) return 0;
        result->type = VAL_MATRIX;
        mat_copy(&result->v.matrix, &pv_tmp_mat);
        return 1;
    }
    
    if (current_token.type == TOK_FUNC) {
        char name[16];
        int var_idx;
        const char *saved_pos = input_ptr;
        token_t saved_tok = current_token;
        
        strcpy(name, current_token.func_name);
        var_idx = get_var_index(name);
        
        /* Check for matrix variable */
        if (var_idx >= 0 && is_var_matrix(var_idx)) {
            next_token();
            if (current_token.type != TOK_LPAREN) {
                get_var_value(var_idx, result);
                return 1;
            }
            input_ptr = saved_pos;
            current_token = saved_tok;
        }
        
        /* Check for scalar variable (used as value_t) */
        if (var_idx >= 0 && scalar_vars[var_idx].defined) {
            next_token();
            if (current_token.type != TOK_LPAREN) {
                result->type = VAL_SCALAR;
                result->v.scalar = scalar_vars[var_idx].val;
                return 1;
            }
            input_ptr = saved_pos;
            current_token = saved_tok;
        }
        
        next_token();
        if (current_token.type == TOK_LPAREN) {
            if (str_eq(name, "zeros") || str_eq(name, "ones") || 
                str_eq(name, "eye") || str_eq(name, "identity") ||
                str_eq(name, "rand") || str_eq(name, "randn") || str_eq(name, "randi")) {
                input_ptr = saved_pos;
                current_token = saved_tok;
                next_token();
                return parse_matrix_func(result, name);
            }
            
            if (str_eq(name, "inv") || str_eq(name, "inverse") ||
                str_eq(name, "transpose") || str_eq(name, "trans") ||
                str_eq(name, "diag")) {
                
                next_token();
                
                if (!parse_value(&pv_arg)) return 0;
                if (current_token.type != TOK_RPAREN) {
                    printf("Error: expected ')'\n");
                    return 0;
                }
                next_token();
                
                if (pv_arg.type != VAL_MATRIX) {
                    mat_zero(&pv_arg.v.matrix, 1, 1);
                    MAT_AT(&pv_arg.v.matrix, 0, 0) = pv_arg.v.scalar;
                    pv_arg.type = VAL_MATRIX;
                }
                
                result->type = VAL_MATRIX;
                if (str_eq(name, "inv") || str_eq(name, "inverse")) {
                    if (!mat_inv(&result->v.matrix, &pv_arg.v.matrix)) {
                        return 0;
                    }
                } else if (str_eq(name, "transpose") || str_eq(name, "trans")) {
                    mat_transpose(&result->v.matrix, &pv_arg.v.matrix);
                } else if (str_eq(name, "diag")) {
                    if (mat_is_vector(&pv_arg.v.matrix)) {
                        mat_diag_create(&result->v.matrix, &pv_arg.v.matrix);
                    } else {
                        mat_diag_extract(&result->v.matrix, &pv_arg.v.matrix);
                    }
                }
                return 1;
            }
        }
        
        input_ptr = saved_pos;
        current_token = saved_tok;
    }
    
    /* Fall back to scalar factor parsing (not full expression!) */
    result->type = VAL_SCALAR;
    return parse_factor(&result->v.scalar);
}

/* Parse value postfix (transpose, factorial) */
static int parse_value_postfix(value_t *result)
{
    /* Uses shared static: pv_tmp_mat */
    
    if (!parse_value_factor(result)) return 0;
    
    while (current_token.type == TOK_TRANSPOSE || 
           current_token.type == TOK_DOT_TRANSPOSE ||
           current_token.type == TOK_FACT) {
        
        if (current_token.type == TOK_TRANSPOSE) {
            next_token();
            if (result->type == VAL_MATRIX) {
                mat_conj_transpose(&pv_tmp_mat, &result->v.matrix);
                mat_copy(&result->v.matrix, &pv_tmp_mat);
            } else {
                apfc_conj(&result->v.scalar, &result->v.scalar);
            }
        } else if (current_token.type == TOK_DOT_TRANSPOSE) {
            next_token();
            if (result->type == VAL_MATRIX) {
                mat_transpose(&pv_tmp_mat, &result->v.matrix);
                mat_copy(&result->v.matrix, &pv_tmp_mat);
            }
        } else if (current_token.type == TOK_FACT) {
            next_token();
            if (result->type != VAL_SCALAR) {
                printf("Error: factorial requires scalar\n");
                return 0;
            }
            if (!apfc_is_real(&result->v.scalar)) {
                printf("Error: factorial requires real integer\n");
                return 0;
            }
            {
                long n = apf_to_long(&result->v.scalar.re);
                if (n < 0 || n > 1000) {
                    printf("Error: factorial argument out of range\n");
                    return 0;
                }
                apfx_fact(&result->v.scalar.re, n);
                apf_zero(&result->v.scalar.im);
            }
        }
    }
    
    return 1;
}

/* Parse value power (^, .^) */
static int parse_value_power(value_t *result)
{
    /* Use local variable for exponent to avoid conflict with recursive calls */
    value_t exponent;
    
    if (!parse_value_postfix(result)) return 0;
    
    if (current_token.type == TOK_POW || current_token.type == TOK_DOT_POW) {
        int is_elemwise = (current_token.type == TOK_DOT_POW);
        
        next_token();
        /* Use parse_value_unary (which calls parse_value_power) for right-associativity */
        if (!parse_value_unary(&exponent)) return 0;
        
        if (is_elemwise) {
            if (result->type == VAL_MATRIX && exponent.type == VAL_MATRIX) {
                mat_elemwise_pow_mat(&pv_tmp_mat, &result->v.matrix, &exponent.v.matrix);
                mat_copy(&result->v.matrix, &pv_tmp_mat);
            } else if (result->type == VAL_MATRIX && exponent.type == VAL_SCALAR) {
                mat_elemwise_pow(&pv_tmp_mat, &result->v.matrix, &exponent.v.scalar);
                mat_copy(&result->v.matrix, &pv_tmp_mat);
            } else if (result->type == VAL_SCALAR && exponent.type == VAL_SCALAR) {
                apfc tmp;
                apfc_pow(&tmp, &result->v.scalar, &exponent.v.scalar);
                result->v.scalar = tmp;
            } else {
                printf("Error: invalid operand types for .^\n");
                return 0;
            }
        } else {
            if (result->type == VAL_MATRIX) {
                if (exponent.type != VAL_SCALAR) {
                    printf("Error: matrix power requires scalar exponent\n");
                    return 0;
                }
                {
                    long n = apf_to_long(&exponent.v.scalar.re);
                    if (!mat_pow(&pv_tmp_mat, &result->v.matrix, n)) return 0;
                    mat_copy(&result->v.matrix, &pv_tmp_mat);
                }
            } else {
                if (exponent.type != VAL_SCALAR) {
                    printf("Error: scalar power requires scalar exponent\n");
                    return 0;
                }
                {
                    apfc tmp;
                    apfc_pow(&tmp, &result->v.scalar, &exponent.v.scalar);
                    result->v.scalar = tmp;
                }
            }
        }
    }
    
    return 1;
}

/* Parse value unary (-, +) */
static int parse_value_unary(value_t *result)
{
    int neg = 0;
    
    while (current_token.type == TOK_MINUS) {
        neg = !neg;
        next_token();
    }
    while (current_token.type == TOK_PLUS) {
        next_token();
    }
    
    if (!parse_value_power(result)) return 0;
    
    if (neg) {
        if (result->type == VAL_MATRIX) {
            matrix_t tmp;
            mat_neg(&tmp, &result->v.matrix);
            mat_copy(&result->v.matrix, &tmp);
        } else {
            apfc_neg(&result->v.scalar, &result->v.scalar);
        }
    }
    
    return 1;
}

/* Parse value term (*, /, \, .*, ./) */
static int parse_value_term(value_t *result)
{
    /* Use local variable for RHS operand to avoid recursion conflicts */
    value_t rhs_operand;
    
    if (!parse_value_unary(result)) return 0;
    
    while (current_token.type == TOK_MUL || current_token.type == TOK_DIV ||
           current_token.type == TOK_BACKSLASH ||
           current_token.type == TOK_DOT_MUL || current_token.type == TOK_DOT_DIV) {
        
        int op = current_token.type;
        
        next_token();
        if (!parse_value_unary(&rhs_operand)) return 0;
        
        if (op == TOK_DOT_MUL) {
            if (result->type == VAL_MATRIX && rhs_operand.type == VAL_MATRIX) {
                mat_elemwise_mul(&pv_tmp_mat, &result->v.matrix, &rhs_operand.v.matrix);
                mat_copy(&result->v.matrix, &pv_tmp_mat);
            } else if (result->type == VAL_SCALAR && rhs_operand.type == VAL_SCALAR) {
                apfc tmp;
                apfc_mul(&tmp, &result->v.scalar, &rhs_operand.v.scalar);
                result->v.scalar = tmp;
            } else {
                printf("Error: .* requires matching types\n");
                return 0;
            }
        } else if (op == TOK_DOT_DIV) {
            if (result->type == VAL_MATRIX && rhs_operand.type == VAL_MATRIX) {
                mat_elemwise_div(&pv_tmp_mat, &result->v.matrix, &rhs_operand.v.matrix);
                mat_copy(&result->v.matrix, &pv_tmp_mat);
            } else if (result->type == VAL_SCALAR && rhs_operand.type == VAL_SCALAR) {
                apfc tmp;
                apfc_div(&tmp, &result->v.scalar, &rhs_operand.v.scalar);
                result->v.scalar = tmp;
            } else {
                printf("Error: ./ requires matching types\n");
                return 0;
            }
        } else if (op == TOK_MUL) {
            if (result->type == VAL_MATRIX && rhs_operand.type == VAL_MATRIX) {
                mat_mul(&pv_tmp_mat, &result->v.matrix, &rhs_operand.v.matrix);
                mat_copy(&result->v.matrix, &pv_tmp_mat);
            } else if (result->type == VAL_MATRIX && rhs_operand.type == VAL_SCALAR) {
                mat_scale(&pv_tmp_mat, &result->v.matrix, &rhs_operand.v.scalar);
                mat_copy(&result->v.matrix, &pv_tmp_mat);
            } else if (result->type == VAL_SCALAR && rhs_operand.type == VAL_MATRIX) {
                mat_scale(&pv_tmp_mat, &rhs_operand.v.matrix, &result->v.scalar);
                result->type = VAL_MATRIX;
                mat_copy(&result->v.matrix, &pv_tmp_mat);
            } else {
                apfc tmp;
                apfc_mul(&tmp, &result->v.scalar, &rhs_operand.v.scalar);
                result->v.scalar = tmp;
            }
        } else if (op == TOK_DIV) {
            if (result->type == VAL_MATRIX && rhs_operand.type == VAL_MATRIX) {
                if (!mat_mrdivide(&pv_tmp_mat, &result->v.matrix, &rhs_operand.v.matrix)) {
                    return 0;
                }
                mat_copy(&result->v.matrix, &pv_tmp_mat);
            } else if (result->type == VAL_MATRIX && rhs_operand.type == VAL_SCALAR) {
                apfc inv_s;
                apfc one;
                apf_from_int(&one.re, 1);
                apf_zero(&one.im);
                apfc_div(&inv_s, &one, &rhs_operand.v.scalar);
                mat_scale(&pv_tmp_mat, &result->v.matrix, &inv_s);
                mat_copy(&result->v.matrix, &pv_tmp_mat);
            } else if (result->type == VAL_SCALAR && rhs_operand.type == VAL_SCALAR) {
                apfc tmp;
                apfc_div(&tmp, &result->v.scalar, &rhs_operand.v.scalar);
                result->v.scalar = tmp;
            } else {
                printf("Error: invalid operand types for /\n");
                return 0;
            }
        } else if (op == TOK_BACKSLASH) {
            if (result->type == VAL_MATRIX && rhs_operand.type == VAL_MATRIX) {
                if (!mat_mldivide(&pv_tmp_mat, &result->v.matrix, &rhs_operand.v.matrix)) {
                    return 0;
                }
                result->type = VAL_MATRIX;
                mat_copy(&result->v.matrix, &pv_tmp_mat);
            } else if (result->type == VAL_SCALAR && rhs_operand.type == VAL_SCALAR) {
                apfc tmp;
                apfc_div(&tmp, &rhs_operand.v.scalar, &result->v.scalar);
                result->v.scalar = tmp;
            } else {
                printf("Error: invalid operand types for \\\n");
                return 0;
            }
        }
    }
    
    return 1;
}

/* Parse value expression (+, -) */
static int parse_value_expr(value_t *result)
{
    /* Use local variable for RHS operand to avoid recursion conflicts */
    value_t right_operand;
    
    if (!parse_value_term(result)) return 0;
    
    while (current_token.type == TOK_PLUS || current_token.type == TOK_MINUS) {
        int op = current_token.type;
        
        next_token();
        if (!parse_value_term(&right_operand)) return 0;
        
        if (op == TOK_PLUS) {
            if (result->type == VAL_MATRIX && right_operand.type == VAL_MATRIX) {
                mat_add(&pv_tmp_mat, &result->v.matrix, &right_operand.v.matrix);
                mat_copy(&result->v.matrix, &pv_tmp_mat);
            } else if (result->type == VAL_MATRIX && right_operand.type == VAL_SCALAR) {
                mat_add_scalar(&pv_tmp_mat, &result->v.matrix, &right_operand.v.scalar);
                mat_copy(&result->v.matrix, &pv_tmp_mat);
            } else if (result->type == VAL_SCALAR && right_operand.type == VAL_MATRIX) {
                mat_add_scalar(&pv_tmp_mat, &right_operand.v.matrix, &result->v.scalar);
                result->type = VAL_MATRIX;
                mat_copy(&result->v.matrix, &pv_tmp_mat);
            } else {
                apfc tmp;
                apfc_add(&tmp, &result->v.scalar, &right_operand.v.scalar);
                result->v.scalar = tmp;
            }
        } else {
            if (result->type == VAL_MATRIX && right_operand.type == VAL_MATRIX) {
                mat_sub(&pv_tmp_mat, &result->v.matrix, &right_operand.v.matrix);
                mat_copy(&result->v.matrix, &pv_tmp_mat);
            } else if (result->type == VAL_MATRIX && right_operand.type == VAL_SCALAR) {
                apfc neg_s;
                apfc_neg(&neg_s, &right_operand.v.scalar);
                mat_add_scalar(&pv_tmp_mat, &result->v.matrix, &neg_s);
                mat_copy(&result->v.matrix, &pv_tmp_mat);
            } else if (result->type == VAL_SCALAR && right_operand.type == VAL_MATRIX) {
                mat_neg(&pv_tmp_mat2, &right_operand.v.matrix);
                mat_add_scalar(&pv_tmp_mat, &pv_tmp_mat2, &result->v.scalar);
                result->type = VAL_MATRIX;
                mat_copy(&result->v.matrix, &pv_tmp_mat);
            } else {
                apfc tmp;
                apfc_sub(&tmp, &result->v.scalar, &right_operand.v.scalar);
                result->v.scalar = tmp;
            }
        }
    }
    
    return 1;
}

/* Parse a value (scalar or matrix) - main entry point */
int parse_value(value_t *result)
{
    return parse_value_expr(result);
}
