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
        if (n < 0 || n > 100000) {
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

    /* Parentheses - may contain boolean expressions */
    if (current_token.type == TOK_LPAREN) {
        int leading_not = 0;
        next_token();
        
        /* Handle leading 'not' inside parens */
        while (current_token.type == TOK_NOT) {
            leading_not = !leading_not;
            next_token();
        }
        
        if (!parse_expr(result)) return 0;
        
        /* Check for comparison operators inside parentheses */
        while (current_token.type == TOK_ASSIGN || current_token.type == TOK_EQUAL ||
               current_token.type == TOK_NE || current_token.type == TOK_LT ||
               current_token.type == TOK_LE || current_token.type == TOK_GT ||
               current_token.type == TOK_GE || current_token.type == TOK_APPROX ||
               current_token.type == TOK_AND || current_token.type == TOK_OR ||
               current_token.type == TOK_XOR) {
            
            token_type_t op = current_token.type;
            
            if (op == TOK_AND || op == TOK_OR || op == TOK_XOR) {
                /* Boolean operator */
                int lhs_bool = !apf_is_zero(&result->re);
                int rhs_bool, bool_result;
                apfc rhs;
                int rhs_not = 0;
                
                next_token();
                
                /* Handle 'not' before right side */
                while (current_token.type == TOK_NOT) {
                    rhs_not = !rhs_not;
                    next_token();
                }
                
                if (!parse_expr(&rhs)) return 0;
                
                /* Check for comparison on right side */
                if (current_token.type == TOK_ASSIGN || current_token.type == TOK_EQUAL ||
                    current_token.type == TOK_NE || current_token.type == TOK_LT ||
                    current_token.type == TOK_LE || current_token.type == TOK_GT ||
                    current_token.type == TOK_GE || current_token.type == TOK_APPROX) {
                    apfc rhs2;
                    int cmp_result = 0;
                    token_type_t cmp_op = current_token.type;
                    next_token();
                    if (!parse_expr(&rhs2)) return 0;
                    {
                        int cmp = apf_cmp(&rhs.re, &rhs2.re);
                        int eq = apf_eq(&rhs.re, &rhs2.re) && apf_eq(&rhs.im, &rhs2.im);
                        switch (cmp_op) {
                            case TOK_ASSIGN: case TOK_EQUAL: cmp_result = eq; break;
                            case TOK_NE: cmp_result = !eq; break;
                            case TOK_LT: cmp_result = (cmp < 0); break;
                            case TOK_LE: cmp_result = (cmp <= 0); break;
                            case TOK_GT: cmp_result = (cmp > 0); break;
                            case TOK_GE: cmp_result = (cmp >= 0); break;
                            case TOK_APPROX: {
                                apf abs_l, abs_r, up, lo, fh, fl;
                                apf_abs(&abs_l, &rhs.re);
                                apf_abs(&abs_r, &rhs2.re);
                                apf_from_str(&fh, "1.0500001");
                                apf_from_str(&fl, "0.9499999");
                                apf_mul(&up, &abs_l, &fh);
                                apf_mul(&lo, &abs_l, &fl);
                                if (apf_is_zero(&rhs.re)) cmp_result = apf_is_zero(&rhs2.re);
                                else if (rhs.re.sign != rhs2.re.sign) cmp_result = 0;
                                else cmp_result = (apf_cmp(&abs_r, &lo) >= 0 && apf_cmp(&abs_r, &up) <= 0);
                                break;
                            }
                            default: cmp_result = 0; break;
                        }
                    }
                    apf_from_int(&rhs.re, cmp_result);
                    apf_zero(&rhs.im);
                }
                
                rhs_bool = !apf_is_zero(&rhs.re);
                if (rhs_not) rhs_bool = !rhs_bool;
                
                switch (op) {
                    case TOK_AND: bool_result = lhs_bool && rhs_bool; break;
                    case TOK_OR:  bool_result = lhs_bool || rhs_bool; break;
                    case TOK_XOR: bool_result = lhs_bool != rhs_bool; break;
                    default: bool_result = 0; break;
                }
                apf_from_int(&result->re, bool_result);
                apf_zero(&result->im);
            } else {
                /* Comparison operator */
                apfc rhs;
                int cmp_result = 0;
                next_token();
                if (!parse_expr(&rhs)) return 0;
                {
                    int cmp = apf_cmp(&result->re, &rhs.re);
                    int eq = apf_eq(&result->re, &rhs.re) && apf_eq(&result->im, &rhs.im);
                    switch (op) {
                        case TOK_ASSIGN: case TOK_EQUAL: cmp_result = eq; break;
                        case TOK_NE: cmp_result = !eq; break;
                        case TOK_LT: cmp_result = (cmp < 0); break;
                        case TOK_LE: cmp_result = (cmp <= 0); break;
                        case TOK_GT: cmp_result = (cmp > 0); break;
                        case TOK_GE: cmp_result = (cmp >= 0); break;
                        case TOK_APPROX: {
                            apf abs_l, abs_r, up, lo, fh, fl;
                            apf_abs(&abs_l, &result->re);
                            apf_abs(&abs_r, &rhs.re);
                            apf_from_str(&fh, "1.0500001");
                            apf_from_str(&fl, "0.9499999");
                            apf_mul(&up, &abs_l, &fh);
                            apf_mul(&lo, &abs_l, &fl);
                            if (apf_is_zero(&result->re)) cmp_result = apf_is_zero(&rhs.re);
                            else if (result->re.sign != rhs.re.sign) cmp_result = 0;
                            else cmp_result = (apf_cmp(&abs_r, &lo) >= 0 && apf_cmp(&abs_r, &up) <= 0);
                            break;
                        }
                        default: cmp_result = 0; break;
                    }
                }
                apf_from_int(&result->re, cmp_result);
                apf_zero(&result->im);
            }
        }
        
        /* Apply leading 'not' */
        if (leading_not) {
            int val = apf_is_zero(&result->re) ? 1 : 0;
            apf_from_int(&result->re, val);
            apf_zero(&result->im);
        }
        
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
        if (str_eq(name, "true")) {
            apf_from_int(&result->re, 1);
            apf_zero(&result->im);
            result_is_boolean = 1;
            return 1;
        }
        if (str_eq(name, "false")) {
            apf_from_int(&result->re, 0);
            apf_zero(&result->im);
            result_is_boolean = 1;
            return 1;
        }
        if (str_eq(name, "eps")) {
            /* Machine epsilon: smallest x such that 1+x != 1 */
            /* For our arbitrary precision, use 2^(-52) like IEEE double */
            apf two;
            apf_from_int(&two, 2);
            apf_from_int(&result->re, -52);
            apfx_pow(&result->re, &two, &result->re);
            apf_zero(&result->im);
            return 1;
        }
        if (str_eq(name, "realmax")) {
            /* Largest finite floating point: ~1.8e308 for IEEE double */
            apf_from_str(&result->re, "1.7976931348623157e+308");
            apf_zero(&result->im);
            return 1;
        }
        if (str_eq(name, "realmin")) {
            /* Smallest positive normalized: ~2.2e-308 for IEEE double */
            apf_from_str(&result->re, "2.2250738585072014e-308");
            apf_zero(&result->im);
            return 1;
        }

        /* Random functions */
        if (str_eq(name, "rand") || str_eq(name, "randn") || str_eq(name, "randi")) {
            return parse_random_func(result, name);
        }

        /* disp() - MATLAB display function */
        if (str_eq(name, "disp")) {
            if (current_token.type != TOK_LPAREN) {
                printf("Error: expected '(' after 'disp'\n");
                return 0;
            }
            /* Check if next char is a string literal */
            while (*input_ptr == ' ') input_ptr++;
            if (*input_ptr == '\'' || *input_ptr == '"') {
                char delim = *input_ptr++;
                while (*input_ptr && *input_ptr != delim) {
                    putchar(*input_ptr++);
                }
                if (*input_ptr == delim) input_ptr++;
                putchar('\n');
                while (*input_ptr == ' ') input_ptr++;
                if (*input_ptr == ')') input_ptr++;
                next_token();
                apf_zero(&result->re);
                apf_zero(&result->im);
                return 1;
            }
            /* Otherwise parse expression */
            next_token();
            if (!parse_expr(&arg)) return 0;
            if (current_token.type != TOK_RPAREN) {
                printf("Error: expected ')'\n");
                return 0;
            }
            next_token();
            {
                char buf[256];
                apfc_to_str(buf, sizeof(buf), &arg, display_digits);
                printf("%s\n", buf);
            }
            *result = arg;
            return 1;
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
            str_eq(name, "rows") || str_eq(name, "cols") || str_eq(name, "numel") ||
            str_eq(name, "sum") || str_eq(name, "mean") || str_eq(name, "prod") ||
            str_eq(name, "min") || str_eq(name, "max") ||
            str_eq(name, "median") || str_eq(name, "sd") || str_eq(name, "std") ||
            str_eq(name, "var") || str_eq(name, "range") || str_eq(name, "rms") ||
            str_eq(name, "sumsq") || str_eq(name, "meansq") ||
            str_eq(name, "geomean") || str_eq(name, "harmmean") ||
            str_eq(name, "skewness") || str_eq(name, "kurtosis") ||
            str_eq(name, "mad") || str_eq(name, "iqr") ||
            str_eq(name, "any") || str_eq(name, "all") ||
            str_eq(name, "isempty") || str_eq(name, "isscalar") ||
            str_eq(name, "isvector") || str_eq(name, "ismatrix") ||
            str_eq(name, "isrow") || str_eq(name, "iscolumn") ||
            str_eq(name, "issquare") || str_eq(name, "issymmetric") ||
            str_eq(name, "isdiag") || str_eq(name, "istriu") || str_eq(name, "istril") ||
            str_eq(name, "issorted") || str_eq(name, "ndims")) {
            
            /* Use local variable to support nested calls like sum(sum(M)) */
            value_t local_arg;
            
            if (current_token.type != TOK_LPAREN) {
                printf("Error: expected '(' after function '%s'\n", name);
                return 0;
            }
            next_token();
            if (!parse_value(&local_arg)) return 0;
            if (current_token.type != TOK_RPAREN) {
                printf("Error: expected ')'\n");
                return 0;
            }
            next_token();
            
            if (local_arg.type == VAL_SCALAR) {
                if (str_eq(name, "det") || str_eq(name, "trace") || str_eq(name, "tr") ||
                    str_eq(name, "sum") || str_eq(name, "mean") || str_eq(name, "prod") ||
                    str_eq(name, "min") || str_eq(name, "max") || str_eq(name, "median") ||
                    str_eq(name, "geomean") || str_eq(name, "harmmean")) {
                    *result = local_arg.v.scalar;
                } else if (str_eq(name, "sd") || str_eq(name, "std") || str_eq(name, "var") ||
                           str_eq(name, "range") || str_eq(name, "skewness") || 
                           str_eq(name, "kurtosis") || str_eq(name, "mad") || str_eq(name, "iqr")) {
                    /* SD/Var/Range/skewness/kurtosis/mad/iqr of single value is 0 */
                    apf_zero(&result->re);
                    apf_zero(&result->im);
                } else if (str_eq(name, "sumsq") || str_eq(name, "meansq")) {
                    /* sumsq/meansq of scalar is x^2 */
                    apf_mul(&result->re, &local_arg.v.scalar.re, &local_arg.v.scalar.re);
                    apf_zero(&result->im);
                } else if (str_eq(name, "rms")) {
                    /* rms of scalar is |x| */
                    apfc_abs(&result->re, &local_arg.v.scalar);
                    apf_zero(&result->im);
                } else if (str_eq(name, "norm")) {
                    apfc_abs(&result->re, &local_arg.v.scalar);
                    apf_zero(&result->im);
                } else if (str_eq(name, "numel") || str_eq(name, "size") || str_eq(name, "length")) {
                    apf_from_int(&result->re, 1);
                    apf_zero(&result->im);
                } else if (str_eq(name, "any") || str_eq(name, "all")) {
                    int val = !apf_is_zero(&local_arg.v.scalar.re);
                    apf_from_int(&result->re, val);
                    apf_zero(&result->im);
                    result_is_boolean = 1;
                } else if (str_eq(name, "isempty")) {
                    apf_from_int(&result->re, 0);  /* scalar is not empty */
                    apf_zero(&result->im);
                    result_is_boolean = 1;
                } else if (str_eq(name, "isscalar")) {
                    apf_from_int(&result->re, 1);  /* scalar is scalar */
                    apf_zero(&result->im);
                    result_is_boolean = 1;
                } else if (str_eq(name, "isvector") || str_eq(name, "ismatrix")) {
                    apf_from_int(&result->re, 0);  /* scalar is not vector/matrix */
                    apf_zero(&result->im);
                    result_is_boolean = 1;
                } else if (str_eq(name, "isrow") || str_eq(name, "iscolumn")) {
                    apf_from_int(&result->re, 0);  /* scalar is not row/column vector */
                    apf_zero(&result->im);
                    result_is_boolean = 1;
                } else if (str_eq(name, "issquare") || str_eq(name, "issymmetric") ||
                           str_eq(name, "isdiag") || str_eq(name, "istriu") ||
                           str_eq(name, "istril")) {
                    apf_from_int(&result->re, 1);  /* scalar is trivially square/symmetric/diagonal/triangular */
                    apf_zero(&result->im);
                    result_is_boolean = 1;
                } else if (str_eq(name, "issorted")) {
                    apf_from_int(&result->re, 1);  /* single element is trivially sorted */
                    apf_zero(&result->im);
                    result_is_boolean = 1;
                } else if (str_eq(name, "ndims")) {
                    apf_from_int(&result->re, 2);  /* MATLAB convention: scalars are 2D */
                    apf_zero(&result->im);
                } else {
                    apf_from_int(&result->re, 1);
                    apf_zero(&result->im);
                }
            } else {
                matrix_t *m = &local_arg.v.matrix;
                
                if (str_eq(name, "det")) {
                    mat_det(result, m);
                } else if (str_eq(name, "trace") || str_eq(name, "tr")) {
                    mat_trace(result, m);
                } else if (str_eq(name, "norm")) {
                    mat_norm_frobenius(result, m);
                } else if (str_eq(name, "size") || str_eq(name, "numel")) {
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
                } else if (str_eq(name, "prod")) {
                    mat_prod(result, m);
                } else if (str_eq(name, "min")) {
                    mat_min(result, m);
                } else if (str_eq(name, "max")) {
                    mat_max(result, m);
                } else if (str_eq(name, "median")) {
                    mat_median(result, m);
                } else if (str_eq(name, "sd") || str_eq(name, "std")) {
                    mat_std(result, m);
                } else if (str_eq(name, "var")) {
                    mat_var(result, m);
                } else if (str_eq(name, "range")) {
                    /* range: max - min */
                    apfc min_val, max_val;
                    mat_min(&min_val, m);
                    mat_max(&max_val, m);
                    apf_sub(&result->re, &max_val.re, &min_val.re);
                    apf_zero(&result->im);
                } else if (str_eq(name, "sumsq")) {
                    /* sumsq: sum of squares */
                    int i, j;
                    apf sq;
                    apf_zero(&result->re);
                    for (i = 0; i < m->rows; i++) {
                        for (j = 0; j < m->cols; j++) {
                            apf_mul(&sq, &MAT_AT(m, i, j).re, &MAT_AT(m, i, j).re);
                            apf_add(&result->re, &result->re, &sq);
                        }
                    }
                    apf_zero(&result->im);
                } else if (str_eq(name, "meansq")) {
                    /* meansq: mean of squares */
                    int i, j, n;
                    apf sq, count;
                    apf_zero(&result->re);
                    n = m->rows * m->cols;
                    for (i = 0; i < m->rows; i++) {
                        for (j = 0; j < m->cols; j++) {
                            apf_mul(&sq, &MAT_AT(m, i, j).re, &MAT_AT(m, i, j).re);
                            apf_add(&result->re, &result->re, &sq);
                        }
                    }
                    apf_from_int(&count, n);
                    apf_div(&result->re, &result->re, &count);
                    apf_zero(&result->im);
                } else if (str_eq(name, "rms")) {
                    /* rms: root mean square = sqrt(meansq) */
                    int i, j, n;
                    apf sq, count;
                    apf_zero(&result->re);
                    n = m->rows * m->cols;
                    for (i = 0; i < m->rows; i++) {
                        for (j = 0; j < m->cols; j++) {
                            apf_mul(&sq, &MAT_AT(m, i, j).re, &MAT_AT(m, i, j).re);
                            apf_add(&result->re, &result->re, &sq);
                        }
                    }
                    apf_from_int(&count, n);
                    apf_div(&result->re, &result->re, &count);
                    apf_sqrt(&result->re, &result->re);
                    apf_zero(&result->im);
                } else if (str_eq(name, "any")) {
                    /* any: true if any element is nonzero */
                    int i, j, found = 0;
                    for (i = 0; i < m->rows && !found; i++) {
                        for (j = 0; j < m->cols && !found; j++) {
                            if (!apf_is_zero(&MAT_AT(m, i, j).re) || 
                                !apf_is_zero(&MAT_AT(m, i, j).im)) {
                                found = 1;
                            }
                        }
                    }
                    apf_from_int(&result->re, found);
                    apf_zero(&result->im);
                    result_is_boolean = 1;
                } else if (str_eq(name, "all")) {
                    /* all: true if all elements are nonzero */
                    int i, j, all_nonzero = 1;
                    for (i = 0; i < m->rows && all_nonzero; i++) {
                        for (j = 0; j < m->cols && all_nonzero; j++) {
                            if (apf_is_zero(&MAT_AT(m, i, j).re) && 
                                apf_is_zero(&MAT_AT(m, i, j).im)) {
                                all_nonzero = 0;
                            }
                        }
                    }
                    apf_from_int(&result->re, all_nonzero);
                    apf_zero(&result->im);
                    result_is_boolean = 1;
                } else if (str_eq(name, "isempty")) {
                    int is_empty = (m->rows == 0 || m->cols == 0);
                    apf_from_int(&result->re, is_empty);
                    apf_zero(&result->im);
                    result_is_boolean = 1;
                } else if (str_eq(name, "isscalar")) {
                    int is_scalar = (m->rows == 1 && m->cols == 1);
                    apf_from_int(&result->re, is_scalar);
                    apf_zero(&result->im);
                    result_is_boolean = 1;
                } else if (str_eq(name, "isvector")) {
                    int is_vec = (m->rows == 1 || m->cols == 1);
                    apf_from_int(&result->re, is_vec);
                    apf_zero(&result->im);
                    result_is_boolean = 1;
                } else if (str_eq(name, "ismatrix")) {
                    apf_from_int(&result->re, 1);  /* matrix is always a matrix */
                    apf_zero(&result->im);
                    result_is_boolean = 1;
                } else if (str_eq(name, "isrow")) {
                    int is_row = (m->rows == 1);
                    apf_from_int(&result->re, is_row);
                    apf_zero(&result->im);
                    result_is_boolean = 1;
                } else if (str_eq(name, "iscolumn")) {
                    int is_col = (m->cols == 1);
                    apf_from_int(&result->re, is_col);
                    apf_zero(&result->im);
                    result_is_boolean = 1;
                } else if (str_eq(name, "issquare")) {
                    int is_sq = (m->rows == m->cols);
                    apf_from_int(&result->re, is_sq);
                    apf_zero(&result->im);
                    result_is_boolean = 1;
                } else if (str_eq(name, "issymmetric")) {
                    int is_sym = 1;
                    int i, j;
                    if (m->rows != m->cols) {
                        is_sym = 0;
                    } else {
                        for (i = 0; i < m->rows && is_sym; i++) {
                            for (j = i + 1; j < m->cols && is_sym; j++) {
                                if (!apf_eq(&MAT_AT(m, i, j).re, &MAT_AT(m, j, i).re) ||
                                    !apf_eq(&MAT_AT(m, i, j).im, &MAT_AT(m, j, i).im)) {
                                    is_sym = 0;
                                }
                            }
                        }
                    }
                    apf_from_int(&result->re, is_sym);
                    apf_zero(&result->im);
                    result_is_boolean = 1;
                } else if (str_eq(name, "isdiag")) {
                    int is_diag = 1;
                    int i, j;
                    for (i = 0; i < m->rows && is_diag; i++) {
                        for (j = 0; j < m->cols && is_diag; j++) {
                            if (i != j) {
                                if (!apf_is_zero(&MAT_AT(m, i, j).re) ||
                                    !apf_is_zero(&MAT_AT(m, i, j).im)) {
                                    is_diag = 0;
                                }
                            }
                        }
                    }
                    apf_from_int(&result->re, is_diag);
                    apf_zero(&result->im);
                    result_is_boolean = 1;
                } else if (str_eq(name, "istriu")) {
                    int is_triu = 1;
                    int i, j;
                    for (i = 1; i < m->rows && is_triu; i++) {
                        for (j = 0; j < i && j < m->cols && is_triu; j++) {
                            if (!apf_is_zero(&MAT_AT(m, i, j).re) ||
                                !apf_is_zero(&MAT_AT(m, i, j).im)) {
                                is_triu = 0;
                            }
                        }
                    }
                    apf_from_int(&result->re, is_triu);
                    apf_zero(&result->im);
                    result_is_boolean = 1;
                } else if (str_eq(name, "istril")) {
                    int is_tril = 1;
                    int i, j;
                    for (i = 0; i < m->rows && is_tril; i++) {
                        for (j = i + 1; j < m->cols && is_tril; j++) {
                            if (!apf_is_zero(&MAT_AT(m, i, j).re) ||
                                !apf_is_zero(&MAT_AT(m, i, j).im)) {
                                is_tril = 0;
                            }
                        }
                    }
                    apf_from_int(&result->re, is_tril);
                    apf_zero(&result->im);
                    result_is_boolean = 1;
                } else if (str_eq(name, "issorted")) {
                    int is_sorted = 1;
                    int n = m->rows * m->cols;
                    int i;
                    for (i = 0; i < n - 1 && is_sorted; i++) {
                        int r1 = i % m->rows, c1 = i / m->rows;
                        int r2 = (i+1) % m->rows, c2 = (i+1) / m->rows;
                        if (apf_cmp(&MAT_AT(m, r1, c1).re, &MAT_AT(m, r2, c2).re) > 0) {
                            is_sorted = 0;
                        }
                    }
                    apf_from_int(&result->re, is_sorted);
                    apf_zero(&result->im);
                    result_is_boolean = 1;
                } else if (str_eq(name, "ndims")) {
                    apf_from_int(&result->re, 2);  /* always 2D in this implementation */
                    apf_zero(&result->im);
                } else if (str_eq(name, "geomean")) {
                    /* Geometric mean: (x1 * x2 * ... * xn)^(1/n) = exp(mean(log(x))) */
                    int i, j, n;
                    apf sum_log, log_val, n_val;
                    n = m->rows * m->cols;
                    apf_zero(&sum_log);
                    for (i = 0; i < m->rows; i++) {
                        for (j = 0; j < m->cols; j++) {
                            apfx_log(&log_val, &MAT_AT(m, i, j).re);
                            apf_add(&sum_log, &sum_log, &log_val);
                        }
                    }
                    apf_from_int(&n_val, n);
                    apf_div(&sum_log, &sum_log, &n_val);
                    apfx_exp(&result->re, &sum_log);
                    apf_zero(&result->im);
                } else if (str_eq(name, "harmmean")) {
                    /* Harmonic mean: n / (1/x1 + 1/x2 + ... + 1/xn) */
                    int i, j, n;
                    apf sum_inv, inv_val, n_val, one;
                    n = m->rows * m->cols;
                    apf_zero(&sum_inv);
                    apf_from_int(&one, 1);
                    for (i = 0; i < m->rows; i++) {
                        for (j = 0; j < m->cols; j++) {
                            apf_div(&inv_val, &one, &MAT_AT(m, i, j).re);
                            apf_add(&sum_inv, &sum_inv, &inv_val);
                        }
                    }
                    apf_from_int(&n_val, n);
                    apf_div(&result->re, &n_val, &sum_inv);
                    apf_zero(&result->im);
                } else if (str_eq(name, "skewness")) {
                    /* Skewness: E[(X-μ)³] / σ³ */
                    int i, j, n;
                    apf mean_val, sum_val, std_val, m3, tmp, n_val, var_val;
                    n = m->rows * m->cols;
                    apf_zero(&sum_val);
                    for (i = 0; i < m->rows; i++) {
                        for (j = 0; j < m->cols; j++) {
                            apf_add(&sum_val, &sum_val, &MAT_AT(m, i, j).re);
                        }
                    }
                    apf_from_int(&n_val, n);
                    apf_div(&mean_val, &sum_val, &n_val);
                    apf_zero(&var_val);
                    apf_zero(&m3);
                    for (i = 0; i < m->rows; i++) {
                        for (j = 0; j < m->cols; j++) {
                            apf diff, diff2, diff3;
                            apf_sub(&diff, &MAT_AT(m, i, j).re, &mean_val);
                            apf_mul(&diff2, &diff, &diff);
                            apf_add(&var_val, &var_val, &diff2);
                            apf_mul(&diff3, &diff2, &diff);
                            apf_add(&m3, &m3, &diff3);
                        }
                    }
                    apf_div(&var_val, &var_val, &n_val);
                    apf_sqrt(&std_val, &var_val);
                    apf_div(&m3, &m3, &n_val);
                    apf_mul(&tmp, &std_val, &var_val);
                    if (!apf_is_zero(&tmp)) {
                        apf_div(&result->re, &m3, &tmp);
                    } else {
                        apf_zero(&result->re);
                    }
                    apf_zero(&result->im);
                } else if (str_eq(name, "kurtosis")) {
                    /* Excess kurtosis: E[(X-μ)⁴] / σ⁴ - 3 */
                    int i, j, n;
                    apf mean_val, sum_val, m4, tmp, n_val, var_val, three;
                    n = m->rows * m->cols;
                    apf_zero(&sum_val);
                    for (i = 0; i < m->rows; i++) {
                        for (j = 0; j < m->cols; j++) {
                            apf_add(&sum_val, &sum_val, &MAT_AT(m, i, j).re);
                        }
                    }
                    apf_from_int(&n_val, n);
                    apf_div(&mean_val, &sum_val, &n_val);
                    apf_zero(&var_val);
                    apf_zero(&m4);
                    for (i = 0; i < m->rows; i++) {
                        for (j = 0; j < m->cols; j++) {
                            apf diff, diff2, diff4;
                            apf_sub(&diff, &MAT_AT(m, i, j).re, &mean_val);
                            apf_mul(&diff2, &diff, &diff);
                            apf_add(&var_val, &var_val, &diff2);
                            apf_mul(&diff4, &diff2, &diff2);
                            apf_add(&m4, &m4, &diff4);
                        }
                    }
                    apf_div(&var_val, &var_val, &n_val);
                    apf_div(&m4, &m4, &n_val);
                    apf_mul(&tmp, &var_val, &var_val);
                    if (!apf_is_zero(&tmp)) {
                        apf_div(&result->re, &m4, &tmp);
                        apf_from_int(&three, 3);
                        apf_sub(&result->re, &result->re, &three);
                    } else {
                        apf_zero(&result->re);
                    }
                    apf_zero(&result->im);
                } else if (str_eq(name, "mad")) {
                    /* Mean absolute deviation: mean(|x - mean(x)|) */
                    int i, j, n;
                    apf mean_val, sum_val, mad_sum, n_val;
                    n = m->rows * m->cols;
                    apf_zero(&sum_val);
                    for (i = 0; i < m->rows; i++) {
                        for (j = 0; j < m->cols; j++) {
                            apf_add(&sum_val, &sum_val, &MAT_AT(m, i, j).re);
                        }
                    }
                    apf_from_int(&n_val, n);
                    apf_div(&mean_val, &sum_val, &n_val);
                    apf_zero(&mad_sum);
                    for (i = 0; i < m->rows; i++) {
                        for (j = 0; j < m->cols; j++) {
                            apf diff, abs_diff;
                            apf_sub(&diff, &MAT_AT(m, i, j).re, &mean_val);
                            apf_abs(&abs_diff, &diff);
                            apf_add(&mad_sum, &mad_sum, &abs_diff);
                        }
                    }
                    apf_div(&result->re, &mad_sum, &n_val);
                    apf_zero(&result->im);
                } else if (str_eq(name, "iqr")) {
                    /* Interquartile range: Q3 - Q1 */
                    /* Sort and find Q1, Q3 */
                    int i, j, n, idx1, idx3;
                    apf q1, q3;
                    n = m->rows * m->cols;
                    /* Sort using pv_tmp_mat */
                    mat_copy(&pv_tmp_mat, m);
                    for (i = 0; i < n - 1; i++) {
                        for (j = 0; j < n - 1 - i; j++) {
                            int r1 = j % pv_tmp_mat.rows;
                            int c1 = j / pv_tmp_mat.rows;
                            int r2 = (j+1) % pv_tmp_mat.rows;
                            int c2 = (j+1) / pv_tmp_mat.rows;
                            if (apf_cmp(&MAT_AT(&pv_tmp_mat, r1, c1).re,
                                       &MAT_AT(&pv_tmp_mat, r2, c2).re) > 0) {
                                apfc tmp = MAT_AT(&pv_tmp_mat, r1, c1);
                                MAT_AT(&pv_tmp_mat, r1, c1) = MAT_AT(&pv_tmp_mat, r2, c2);
                                MAT_AT(&pv_tmp_mat, r2, c2) = tmp;
                            }
                        }
                    }
                    /* Q1 at 25%, Q3 at 75% */
                    idx1 = (n - 1) / 4;
                    idx3 = (3 * (n - 1)) / 4;
                    {
                        int r1 = idx1 % pv_tmp_mat.rows;
                        int c1 = idx1 / pv_tmp_mat.rows;
                        int r3 = idx3 % pv_tmp_mat.rows;
                        int c3 = idx3 / pv_tmp_mat.rows;
                        apf_copy(&q1, &MAT_AT(&pv_tmp_mat, r1, c1).re);
                        apf_copy(&q3, &MAT_AT(&pv_tmp_mat, r3, c3).re);
                    }
                    apf_sub(&result->re, &q3, &q1);
                    apf_zero(&result->im);
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

        /* rank(M) - matrix rank */
        if (str_eq(name, "rank")) {
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
                /* Scalar: rank is 1 if nonzero, 0 if zero */
                if (apf_is_zero(&pv_arg.v.scalar.re) && apf_is_zero(&pv_arg.v.scalar.im)) {
                    apf_zero(&result->re);
                } else {
                    apf_from_int(&result->re, 1);
                }
                apf_zero(&result->im);
            } else {
                /* Matrix: count nonzero pivots in LU decomposition */
                matrix_t *m = &pv_arg.v.matrix;
                matrix_t L, U;
                int perm[MAT_MAX_ROWS];
                int i, rank_val = 0;
                apf tol, abs_val;
                
                if (!mat_lu(&L, &U, perm, m)) {
                    printf("Error: LU decomposition failed\n");
                    return 0;
                }
                
                /* Tolerance for zero detection */
                apf_from_str(&tol, "1e-30");
                
                /* Count nonzero diagonal elements in U */
                for (i = 0; i < U.rows && i < U.cols; i++) {
                    apf_abs(&abs_val, &MAT_AT(&U, i, i).re);
                    if (apf_cmp(&abs_val, &tol) > 0) {
                        rank_val++;
                    }
                }
                
                apf_from_int(&result->re, rank_val);
                apf_zero(&result->im);
            }
            return 1;
        }

        /* cond(M) - condition number (ratio of largest to smallest singular value) */
        if (str_eq(name, "cond")) {
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
                /* Scalar: condition is 1 */
                apf_from_int(&result->re, 1);
                apf_zero(&result->im);
            } else {
                /* For 2x2, use SVD; otherwise use norm estimate */
                matrix_t *m = &pv_arg.v.matrix;
                if (m->rows == 2 && m->cols == 2) {
                    apfc s1, s2;
                    apf abs1, abs2;
                    mat_svd2(&s1, &s2, m);
                    apfc_abs(&abs1, &s1);
                    apfc_abs(&abs2, &s2);
                    if (apf_is_zero(&abs2)) {
                        apf_set_inf(&result->re, 0);
                    } else {
                        apf_div(&result->re, &abs1, &abs2);
                    }
                } else {
                    /* Use norm(A) * norm(inv(A)) estimate */
                    apfc norm_a, norm_inv;
                    matrix_t inv_m;
                    mat_norm_frobenius(&norm_a, m);
                    if (mat_inv(&inv_m, m)) {
                        mat_norm_frobenius(&norm_inv, &inv_m);
                        apf_mul(&result->re, &norm_a.re, &norm_inv.re);
                    } else {
                        apf_set_inf(&result->re, 0);
                    }
                }
                apf_zero(&result->im);
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
        /* mod(a, b) - modulo operation */
        if (str_eq(name, "mod")) {
            apfc arg2;
            if (current_token.type != TOK_COMMA) {
                printf("Error: mod requires two arguments: mod(a, b)\n");
                return 0;
            }
            next_token();
            if (!parse_expr(&arg2)) return 0;
            if (current_token.type != TOK_RPAREN) {
                printf("Error: expected ')'\n");
                return 0;
            }
            next_token();
            apf_mod(&result->re, &arg.re, &arg2.re);
            apf_zero(&result->im);
            return 1;
        }
        /* rem(a, b) - remainder (sign of dividend, unlike mod which uses divisor sign) */
        if (str_eq(name, "rem")) {
            apfc arg2;
            apf quot, truncated;
            if (current_token.type != TOK_COMMA) {
                printf("Error: rem requires two arguments: rem(a, b)\n");
                return 0;
            }
            next_token();
            if (!parse_expr(&arg2)) return 0;
            if (current_token.type != TOK_RPAREN) {
                printf("Error: expected ')'\n");
                return 0;
            }
            next_token();
            /* rem = a - trunc(a/b) * b */
            apf_div(&quot, &arg.re, &arg2.re);
            apf_trunc(&truncated, &quot);
            apf_mul(&truncated, &truncated, &arg2.re);
            apf_sub(&result->re, &arg.re, &truncated);
            apf_zero(&result->im);
            return 1;
        }
        /* dot(a1,a2,a3,b1,b2,b3) - 3D vector dot product */
        if (str_eq(name, "dot")) {
            apfc a2, a3, b1, b2, b3;
            apf t1, t2, t3;
            if (current_token.type != TOK_COMMA) goto dot_error;
            next_token(); if (!parse_expr(&a2)) return 0;
            if (current_token.type != TOK_COMMA) goto dot_error;
            next_token(); if (!parse_expr(&a3)) return 0;
            if (current_token.type != TOK_COMMA) goto dot_error;
            next_token(); if (!parse_expr(&b1)) return 0;
            if (current_token.type != TOK_COMMA) goto dot_error;
            next_token(); if (!parse_expr(&b2)) return 0;
            if (current_token.type != TOK_COMMA) goto dot_error;
            next_token(); if (!parse_expr(&b3)) return 0;
            if (current_token.type != TOK_RPAREN) goto dot_error;
            next_token();
            /* dot = a1*b1 + a2*b2 + a3*b3 */
            apf_mul(&t1, &arg.re, &b1.re);
            apf_mul(&t2, &a2.re, &b2.re);
            apf_mul(&t3, &a3.re, &b3.re);
            apf_add(&result->re, &t1, &t2);
            apf_add(&result->re, &result->re, &t3);
            apf_zero(&result->im);
            return 1;
            dot_error:
            printf("Error: dot requires 6 arguments: dot(a1,a2,a3,b1,b2,b3)\n");
            return 0;
        }
        /* cross(a1,a2,a3,b1,b2,b3) - 3D vector cross product (returns first component) */
        if (str_eq(name, "cross")) {
            apfc a2, a3, b1, b2, b3;
            apf t1, t2;
            if (current_token.type != TOK_COMMA) goto cross_error;
            next_token(); if (!parse_expr(&a2)) return 0;
            if (current_token.type != TOK_COMMA) goto cross_error;
            next_token(); if (!parse_expr(&a3)) return 0;
            if (current_token.type != TOK_COMMA) goto cross_error;
            next_token(); if (!parse_expr(&b1)) return 0;
            if (current_token.type != TOK_COMMA) goto cross_error;
            next_token(); if (!parse_expr(&b2)) return 0;
            if (current_token.type != TOK_COMMA) goto cross_error;
            next_token(); if (!parse_expr(&b3)) return 0;
            if (current_token.type != TOK_RPAREN) goto cross_error;
            next_token();
            /* cross x = a2*b3 - a3*b2 */
            apf_mul(&t1, &a2.re, &b3.re);
            apf_mul(&t2, &a3.re, &b2.re);
            apf_sub(&result->re, &t1, &t2);
            apf_zero(&result->im);
            /* Print all 3 components */
            {
                apf cx, cy, cz;
                char buf[64];
                apf_mul(&t1, &a2.re, &b3.re);
                apf_mul(&t2, &a3.re, &b2.re);
                apf_sub(&cx, &t1, &t2);
                apf_mul(&t1, &a3.re, &b1.re);
                apf_mul(&t2, &arg.re, &b3.re);
                apf_sub(&cy, &t1, &t2);
                apf_mul(&t1, &arg.re, &b2.re);
                apf_mul(&t2, &a2.re, &b1.re);
                apf_sub(&cz, &t1, &t2);
                printf("[ ");
                apf_to_str(buf, sizeof(buf), &cx, display_digits); printf("%s ", buf);
                apf_to_str(buf, sizeof(buf), &cy, display_digits); printf("%s ", buf);
                apf_to_str(buf, sizeof(buf), &cz, display_digits); printf("%s ", buf);
                printf("]\n");
            }
            return 1;
            cross_error:
            printf("Error: cross requires 6 arguments: cross(a1,a2,a3,b1,b2,b3)\n");
            return 0;
        }
        /* atan2(y, x) - two-argument arctangent */
        if (str_eq(name, "atan2")) {
            apfc arg2;
            if (current_token.type != TOK_COMMA) {
                printf("Error: atan2 requires two arguments: atan2(y, x)\n");
                return 0;
            }
            next_token();
            if (!parse_expr(&arg2)) return 0;
            if (current_token.type != TOK_RPAREN) {
                printf("Error: expected ')'\n");
                return 0;
            }
            next_token();
            apfx_atan2(&result->re, &arg.re, &arg2.re);
            if (angle_mode != ANGLE_RAD) {
                apf tmp;
                rad_to_angle(&tmp, &result->re, angle_mode);
                result->re = tmp;
            }
            apf_zero(&result->im);
            return 1;
        }
        /* atan2d(y, x) - two-argument arctangent returning degrees */
        if (str_eq(name, "atan2d")) {
            apfc arg2;
            apf pi_val, d180;
            if (current_token.type != TOK_COMMA) {
                printf("Error: atan2d requires two arguments: atan2d(y, x)\n");
                return 0;
            }
            next_token();
            if (!parse_expr(&arg2)) return 0;
            if (current_token.type != TOK_RPAREN) {
                printf("Error: expected ')'\n");
                return 0;
            }
            next_token();
            apfx_atan2(&result->re, &arg.re, &arg2.re);
            /* Convert radians to degrees */
            apfx_pi(&pi_val);
            apf_from_int(&d180, 180);
            apf_mul(&result->re, &result->re, &d180);
            apf_div(&result->re, &result->re, &pi_val);
            apf_zero(&result->im);
            return 1;
        }
        /* logb(x, base) - logarithm with arbitrary base */
        if (str_eq(name, "logb") || str_eq(name, "log_")) {
            apfc arg2;
            if (current_token.type != TOK_COMMA) {
                printf("Error: logb requires two arguments: logb(x, base)\n");
                return 0;
            }
            next_token();
            if (!parse_expr(&arg2)) return 0;
            if (current_token.type != TOK_RPAREN) {
                printf("Error: expected ')'\n");
                return 0;
            }
            next_token();
            if (!apfc_is_real(&arg) || !apfc_is_real(&arg2)) {
                printf("Error: logb requires real arguments\n");
                return 0;
            }
            apfx_logb(&result->re, &arg.re, &arg2.re);
            apf_zero(&result->im);
            return 1;
        }
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
        /* hypot(x,y) = sqrt(x^2 + y^2) */
        if (str_eq(name, "hypot")) {
            apfc arg2;
            apf x2, y2, sum;
            if (current_token.type != TOK_COMMA) {
                printf("Error: hypot requires two arguments\n");
                return 0;
            }
            next_token();
            if (!parse_expr(&arg2)) return 0;
            if (current_token.type != TOK_RPAREN) {
                printf("Error: expected ')'\n");
                return 0;
            }
            next_token();
            apf_mul(&x2, &arg.re, &arg.re);
            apf_mul(&y2, &arg2.re, &arg2.re);
            apf_add(&sum, &x2, &y2);
            apf_sqrt(&result->re, &sum);
            apf_zero(&result->im);
            return 1;
        }
        /* beta(a, b) - Beta function: gamma(a)*gamma(b)/gamma(a+b) */
        if (str_eq(name, "beta")) {
            apfc arg2;
            apf ga, gb, gab, ab, num;
            if (current_token.type != TOK_COMMA) {
                printf("Error: beta requires two arguments: beta(a, b)\n");
                return 0;
            }
            next_token();
            if (!parse_expr(&arg2)) return 0;
            if (current_token.type != TOK_RPAREN) {
                printf("Error: expected ')'\n");
                return 0;
            }
            next_token();
            /* B(a,b) = gamma(a)*gamma(b)/gamma(a+b) */
            apfx_tgamma(&ga, &arg.re);
            apfx_tgamma(&gb, &arg2.re);
            apf_add(&ab, &arg.re, &arg2.re);
            apfx_tgamma(&gab, &ab);
            apf_mul(&num, &ga, &gb);
            apf_div(&result->re, &num, &gab);
            apf_zero(&result->im);
            return 1;
        }
        /* besselj(n, x) - Bessel function of first kind */
        if (str_eq(name, "besselj") || str_eq(name, "besselJ")) {
            apfc arg2;
            int n;
            if (current_token.type != TOK_COMMA) {
                printf("Error: besselj requires two arguments: besselj(n, x)\n");
                return 0;
            }
            next_token();
            if (!parse_expr(&arg2)) return 0;
            if (current_token.type != TOK_RPAREN) {
                printf("Error: expected ')'\n");
                return 0;
            }
            next_token();
            n = (int)apf_to_long(&arg.re);
            apfx_besselj(&result->re, n, &arg2.re);
            apf_zero(&result->im);
            return 1;
        }
        /* bessely(n, x) - Bessel function of second kind */
        if (str_eq(name, "bessely") || str_eq(name, "besselY")) {
            apfc arg2;
            int n;
            if (current_token.type != TOK_COMMA) {
                printf("Error: bessely requires two arguments: bessely(n, x)\n");
                return 0;
            }
            next_token();
            if (!parse_expr(&arg2)) return 0;
            if (current_token.type != TOK_RPAREN) {
                printf("Error: expected ')'\n");
                return 0;
            }
            next_token();
            n = (int)apf_to_long(&arg.re);
            apfx_bessely(&result->re, n, &arg2.re);
            apf_zero(&result->im);
            return 1;
        }
        /* tcdf(t, df) - Student's t CDF */
        if (str_eq(name, "tcdf") || str_eq(name, "student_cdf")) {
            apfc arg2;
            if (current_token.type != TOK_COMMA) {
                printf("Error: tcdf requires two arguments: tcdf(t, df)\n");
                return 0;
            }
            next_token();
            if (!parse_expr(&arg2)) return 0;
            if (current_token.type != TOK_RPAREN) {
                printf("Error: expected ')'\n");
                return 0;
            }
            next_token();
            apfx_tcdf(&result->re, &arg.re, apf_to_long(&arg2.re));
            apf_zero(&result->im);
            return 1;
        }
        /* chi2cdf(x, df) - Chi-squared CDF */
        if (str_eq(name, "chi2cdf") || str_eq(name, "chicdf")) {
            apfc arg2;
            if (current_token.type != TOK_COMMA) {
                printf("Error: chi2cdf requires two arguments: chi2cdf(x, df)\n");
                return 0;
            }
            next_token();
            if (!parse_expr(&arg2)) return 0;
            if (current_token.type != TOK_RPAREN) {
                printf("Error: expected ')'\n");
                return 0;
            }
            next_token();
            apfx_chi2cdf(&result->re, &arg.re, apf_to_long(&arg2.re));
            apf_zero(&result->im);
            return 1;
        }
        /* nthroot(x, n) - real nth root preserving sign */
        if (str_eq(name, "nthroot")) {
            apfc arg2;
            long n;
            int neg;
            apf abs_x, inv_n, res;
            if (current_token.type != TOK_COMMA) {
                printf("Error: nthroot requires two arguments: nthroot(x, n)\n");
                return 0;
            }
            next_token();
            if (!parse_expr(&arg2)) return 0;
            if (current_token.type != TOK_RPAREN) {
                printf("Error: expected ')'\n");
                return 0;
            }
            next_token();
            n = apf_to_long(&arg2.re);
            if (n == 0) {
                printf("Error: nthroot root cannot be zero\n");
                return 0;
            }
            neg = apf_cmp_int(&arg.re, 0) < 0;
            if (neg && (n % 2 == 0)) {
                printf("Error: even root of negative number\n");
                return 0;
            }
            if (neg) {
                apf_neg(&abs_x, &arg.re);
            } else {
                abs_x = arg.re;
            }
            apf_from_int(&inv_n, 1);
            apf_div(&inv_n, &inv_n, &arg2.re);
            apfx_pow(&res, &abs_x, &inv_n);
            if (neg) {
                apf_neg(&result->re, &res);
            } else {
                result->re = res;
            }
            apf_zero(&result->im);
            return 1;
        }
        /* realpow(x, y) - real power, error if result would be complex */
        if (str_eq(name, "realpow")) {
            apfc arg2;
            if (current_token.type != TOK_COMMA) {
                printf("Error: realpow requires two arguments: realpow(x, y)\n");
                return 0;
            }
            next_token();
            if (!parse_expr(&arg2)) return 0;
            if (current_token.type != TOK_RPAREN) {
                printf("Error: expected ')'\n");
                return 0;
            }
            next_token();
            if (!apfc_is_real(&arg) || !apfc_is_real(&arg2)) {
                printf("Error: realpow requires real arguments\n");
                return 0;
            }
            /* Check if result would be complex (negative base with non-integer exponent) */
            if (apf_cmp_int(&arg.re, 0) < 0) {
                long y_int = apf_to_long(&arg2.re);
                apf y_check;
                apf_from_int(&y_check, y_int);
                if (!apf_eq(&arg2.re, &y_check)) {
                    printf("Error: realpow would produce complex result\n");
                    return 0;
                }
            }
            apfx_pow(&result->re, &arg.re, &arg2.re);
            apf_zero(&result->im);
            return 1;
        }
        /* complex(a, b) - create complex number from real and imaginary parts */
        if (str_eq(name, "complex")) {
            apfc arg2;
            if (current_token.type != TOK_COMMA) {
                printf("Error: complex requires two arguments: complex(re, im)\n");
                return 0;
            }
            next_token();
            if (!parse_expr(&arg2)) return 0;
            if (current_token.type != TOK_RPAREN) {
                printf("Error: expected ')'\n");
                return 0;
            }
            next_token();
            result->re = arg.re;
            result->im = arg2.re;
            return 1;
        }
        /* isequal(a, b) - test equality */
        if (str_eq(name, "isequal")) {
            apfc arg2;
            int eq;
            if (current_token.type != TOK_COMMA) {
                printf("Error: isequal requires two arguments\n");
                return 0;
            }
            next_token();
            if (!parse_expr(&arg2)) return 0;
            if (current_token.type != TOK_RPAREN) {
                printf("Error: expected ')'\n");
                return 0;
            }
            next_token();
            eq = apf_eq(&arg.re, &arg2.re) && apf_eq(&arg.im, &arg2.im);
            apf_from_int(&result->re, eq);
            apf_zero(&result->im);
            result_is_boolean = 1;
            return 1;
        }
#ifdef HAVE_BITWISE
        if (str_eq(name, "and") || str_eq(name, "band") || str_eq(name, "bitand")) {
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
        if (str_eq(name, "or") || str_eq(name, "bor") || str_eq(name, "bitor")) {
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
        if (str_eq(name, "xor") || str_eq(name, "bxor") || str_eq(name, "bitxor")) {
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
        if (str_eq(name, "lsl") || str_eq(name, "shl") || str_eq(name, "bitshift")) {
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
        /* Beta function B(a,b) = Gamma(a)*Gamma(b)/Gamma(a+b) */
        if (str_eq(name, "beta")) {
            apfc arg2;
            if (current_token.type != TOK_COMMA) {
                printf("Error: beta requires two arguments: beta(a, b)\n");
                return 0;
            }
            next_token();
            if (!parse_expr(&arg2)) return 0;
            if (current_token.type != TOK_RPAREN) {
                printf("Error: expected ')'\n");
                return 0;
            }
            next_token();
            apfx_beta(&result->re, &arg.re, &arg2.re);
            apf_zero(&result->im);
            return 1;
        }
        /* Stirling numbers of the second kind S(n,k) */
        if (str_eq(name, "stirling2") || str_eq(name, "stirling")) {
            apfc arg2;
            long n, k;
            if (current_token.type != TOK_COMMA) {
                printf("Error: stirling2 requires two arguments: stirling2(n, k)\n");
                return 0;
            }
            next_token();
            if (!parse_expr(&arg2)) return 0;
            if (current_token.type != TOK_RPAREN) {
                printf("Error: expected ')'\n");
                return 0;
            }
            next_token();
            n = apf_to_long(&arg.re);
            k = apf_to_long(&arg2.re);
            if (n < 0 || k < 0 || n > 100 || k > n) {
                printf("Error: stirling2 arguments out of range\n");
                return 0;
            }
            apf_stirling2(&result->re, n, k);
            apf_zero(&result->im);
            return 1;
        }
        /* Pochhammer symbol (rising factorial) (x)_n */
        if (str_eq(name, "pochhammer") || str_eq(name, "rising")) {
            apfc arg2;
            long n;
            if (current_token.type != TOK_COMMA) {
                printf("Error: pochhammer requires two arguments: pochhammer(x, n)\n");
                return 0;
            }
            next_token();
            if (!parse_expr(&arg2)) return 0;
            if (current_token.type != TOK_RPAREN) {
                printf("Error: expected ')'\n");
                return 0;
            }
            next_token();
            n = apf_to_long(&arg2.re);
            if (n < 0 || n > 100) {
                printf("Error: pochhammer n out of range (0-100)\n");
                return 0;
            }
            apfx_pochhammer(&result->re, &arg.re, n);
            apf_zero(&result->im);
            return 1;
        }
        /* Falling factorial x^(n) */
        if (str_eq(name, "falling") || str_eq(name, "ffactorial")) {
            apfc arg2;
            long n;
            if (current_token.type != TOK_COMMA) {
                printf("Error: falling requires two arguments: falling(x, n)\n");
                return 0;
            }
            next_token();
            if (!parse_expr(&arg2)) return 0;
            if (current_token.type != TOK_RPAREN) {
                printf("Error: expected ')'\n");
                return 0;
            }
            next_token();
            n = apf_to_long(&arg2.re);
            if (n < 0 || n > 100) {
                printf("Error: falling n out of range (0-100)\n");
                return 0;
            }
            apfx_falling(&result->re, &arg.re, n);
            apf_zero(&result->im);
            return 1;
        }
        /* Generalized harmonic number H_{n,m} */
        if (str_eq(name, "harmonic2") || str_eq(name, "genharmonic")) {
            apfc arg2;
            long n, m;
            if (current_token.type != TOK_COMMA) {
                printf("Error: harmonic2 requires two arguments: harmonic2(n, m)\n");
                return 0;
            }
            next_token();
            if (!parse_expr(&arg2)) return 0;
            if (current_token.type != TOK_RPAREN) {
                printf("Error: expected ')'\n");
                return 0;
            }
            next_token();
            n = apf_to_long(&arg.re);
            m = apf_to_long(&arg2.re);
            if (n < 0 || n > 10000) {
                printf("Error: harmonic2 n out of range\n");
                return 0;
            }
            apfx_harmonic_gen(&result->re, n, m);
            apf_zero(&result->im);
            return 1;
        }
        /* Sigma function with power parameter */
        if (str_eq(name, "sigmak") || str_eq(name, "divisorsum_k")) {
            apfc arg2;
            long n, k;
            if (current_token.type != TOK_COMMA) {
                printf("Error: sigmak requires two arguments: sigmak(n, k)\n");
                return 0;
            }
            next_token();
            if (!parse_expr(&arg2)) return 0;
            if (current_token.type != TOK_RPAREN) {
                printf("Error: expected ')'\n");
                return 0;
            }
            next_token();
            n = apf_to_long(&arg.re);
            k = (int)apf_to_long(&arg2.re);
            apf_from_int(&result->re, sigma_long(n, k));
            apf_zero(&result->im);
            return 1;
        }

        /* Single-argument functions - require ')' */
        if (current_token.type != TOK_RPAREN) {
            printf("Error: expected ')'\n");
            return 0;
        }
        next_token();

        /* Built-in functions */
        if (str_eq(name, "sin")) {
            apfc rad_arg;
            if (angle_mode != ANGLE_RAD && apfc_is_real(&arg)) {
                angle_to_rad(&rad_arg.re, &arg.re, angle_mode);
                apf_zero(&rad_arg.im);
                apfc_sin(result, &rad_arg);
            } else {
                apfc_sin(result, &arg);
            }
        } else if (str_eq(name, "cos")) {
            apfc rad_arg;
            if (angle_mode != ANGLE_RAD && apfc_is_real(&arg)) {
                angle_to_rad(&rad_arg.re, &arg.re, angle_mode);
                apf_zero(&rad_arg.im);
                apfc_cos(result, &rad_arg);
            } else {
                apfc_cos(result, &arg);
            }
        } else if (str_eq(name, "tan")) {
            apfc rad_arg;
            if (angle_mode != ANGLE_RAD && apfc_is_real(&arg)) {
                angle_to_rad(&rad_arg.re, &arg.re, angle_mode);
                apf_zero(&rad_arg.im);
                apfc_tan(result, &rad_arg);
            } else {
                apfc_tan(result, &arg);
            }
        } else if (str_eq(name, "sind")) {
            /* sind: sine of angle in degrees */
            apfc rad_arg;
            apf pi_val, d180;
            apfx_pi(&pi_val);
            apf_from_int(&d180, 180);
            apf_mul(&rad_arg.re, &arg.re, &pi_val);
            apf_div(&rad_arg.re, &rad_arg.re, &d180);
            apf_zero(&rad_arg.im);
            apfc_sin(result, &rad_arg);
        } else if (str_eq(name, "cosd")) {
            /* cosd: cosine of angle in degrees */
            apfc rad_arg;
            apf pi_val, d180;
            apfx_pi(&pi_val);
            apf_from_int(&d180, 180);
            apf_mul(&rad_arg.re, &arg.re, &pi_val);
            apf_div(&rad_arg.re, &rad_arg.re, &d180);
            apf_zero(&rad_arg.im);
            apfc_cos(result, &rad_arg);
        } else if (str_eq(name, "tand")) {
            /* tand: tangent of angle in degrees */
            apfc rad_arg;
            apf pi_val, d180;
            apfx_pi(&pi_val);
            apf_from_int(&d180, 180);
            apf_mul(&rad_arg.re, &arg.re, &pi_val);
            apf_div(&rad_arg.re, &rad_arg.re, &d180);
            apf_zero(&rad_arg.im);
            apfc_tan(result, &rad_arg);
        } else if (str_eq(name, "asind")) {
            /* asind: arcsine returning degrees */
            apf rad_result, pi_val, d180;
            if (!apfc_is_real(&arg)) {
                printf("Error: asind requires real argument\n");
                return 0;
            }
            apfx_asin(&rad_result, &arg.re);
            apfx_pi(&pi_val);
            apf_from_int(&d180, 180);
            apf_mul(&result->re, &rad_result, &d180);
            apf_div(&result->re, &result->re, &pi_val);
            apf_zero(&result->im);
        } else if (str_eq(name, "acosd")) {
            /* acosd: arccosine returning degrees */
            apf rad_result, pi_val, d180;
            if (!apfc_is_real(&arg)) {
                printf("Error: acosd requires real argument\n");
                return 0;
            }
            apfx_acos(&rad_result, &arg.re);
            apfx_pi(&pi_val);
            apf_from_int(&d180, 180);
            apf_mul(&result->re, &rad_result, &d180);
            apf_div(&result->re, &result->re, &pi_val);
            apf_zero(&result->im);
        } else if (str_eq(name, "atand")) {
            /* atand: arctangent returning degrees */
            apf rad_result, pi_val, d180;
            if (!apfc_is_real(&arg)) {
                printf("Error: atand requires real argument\n");
                return 0;
            }
            apfx_atan(&rad_result, &arg.re);
            apfx_pi(&pi_val);
            apf_from_int(&d180, 180);
            apf_mul(&result->re, &rad_result, &d180);
            apf_div(&result->re, &result->re, &pi_val);
            apf_zero(&result->im);
        } else if (str_eq(name, "secd")) {
            /* secd: secant of angle in degrees = 1/cos(x) */
            apfc rad_arg, cos_val;
            apf pi_val, d180, one;
            apfx_pi(&pi_val);
            apf_from_int(&d180, 180);
            apf_mul(&rad_arg.re, &arg.re, &pi_val);
            apf_div(&rad_arg.re, &rad_arg.re, &d180);
            apf_zero(&rad_arg.im);
            apfc_cos(&cos_val, &rad_arg);
            apf_from_int(&one, 1);
            apfc_from_real(result, &one);
            apfc_div(result, result, &cos_val);
        } else if (str_eq(name, "cscd")) {
            /* cscd: cosecant of angle in degrees = 1/sin(x) */
            apfc rad_arg, sin_val;
            apf pi_val, d180, one;
            apfx_pi(&pi_val);
            apf_from_int(&d180, 180);
            apf_mul(&rad_arg.re, &arg.re, &pi_val);
            apf_div(&rad_arg.re, &rad_arg.re, &d180);
            apf_zero(&rad_arg.im);
            apfc_sin(&sin_val, &rad_arg);
            apf_from_int(&one, 1);
            apfc_from_real(result, &one);
            apfc_div(result, result, &sin_val);
        } else if (str_eq(name, "cotd")) {
            /* cotd: cotangent of angle in degrees = cos(x)/sin(x) */
            apfc rad_arg, sin_val, cos_val;
            apf pi_val, d180;
            apfx_pi(&pi_val);
            apf_from_int(&d180, 180);
            apf_mul(&rad_arg.re, &arg.re, &pi_val);
            apf_div(&rad_arg.re, &rad_arg.re, &d180);
            apf_zero(&rad_arg.im);
            apfc_sin(&sin_val, &rad_arg);
            apfc_cos(&cos_val, &rad_arg);
            apfc_div(result, &cos_val, &sin_val);
        } else if (str_eq(name, "asecd")) {
            /* asecd: inverse secant returning degrees = acos(1/x) * 180/pi */
            apf rad_result, pi_val, d180, one;
            apfc inv_arg;
            if (!apfc_is_real(&arg)) {
                printf("Error: asecd requires real argument\n");
                return 0;
            }
            apf_from_int(&one, 1);
            apf_div(&inv_arg.re, &one, &arg.re);
            apfx_acos(&rad_result, &inv_arg.re);
            apfx_pi(&pi_val);
            apf_from_int(&d180, 180);
            apf_mul(&result->re, &rad_result, &d180);
            apf_div(&result->re, &result->re, &pi_val);
            apf_zero(&result->im);
        } else if (str_eq(name, "acscd")) {
            /* acscd: inverse cosecant returning degrees = asin(1/x) * 180/pi */
            apf rad_result, pi_val, d180, one;
            apfc inv_arg;
            if (!apfc_is_real(&arg)) {
                printf("Error: acscd requires real argument\n");
                return 0;
            }
            apf_from_int(&one, 1);
            apf_div(&inv_arg.re, &one, &arg.re);
            apfx_asin(&rad_result, &inv_arg.re);
            apfx_pi(&pi_val);
            apf_from_int(&d180, 180);
            apf_mul(&result->re, &rad_result, &d180);
            apf_div(&result->re, &result->re, &pi_val);
            apf_zero(&result->im);
        } else if (str_eq(name, "acotd")) {
            /* acotd: inverse cotangent returning degrees = atan(1/x) * 180/pi */
            apf rad_result, pi_val, d180, one;
            apfc inv_arg;
            if (!apfc_is_real(&arg)) {
                printf("Error: acotd requires real argument\n");
                return 0;
            }
            apf_from_int(&one, 1);
            apf_div(&inv_arg.re, &one, &arg.re);
            apfx_atan(&rad_result, &inv_arg.re);
            apfx_pi(&pi_val);
            apf_from_int(&d180, 180);
            apf_mul(&result->re, &rad_result, &d180);
            apf_div(&result->re, &result->re, &pi_val);
            apf_zero(&result->im);
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
        } else if (str_eq(name, "log2")) {
            /* log2(x) = ln(x) / ln(2) */
            apfc log_arg, ln2;
            apf two;
            apfc_log(&log_arg, &arg);
            apf_from_int(&two, 2);
            apfc_from_real(&ln2, &two);
            apfc_log(&ln2, &ln2);
            apfc_div(result, &log_arg, &ln2);
        } else if (str_eq(name, "exp2") || str_eq(name, "pow2")) {
            /* exp2(x) = 2^x */
            apf two;
            apf_from_int(&two, 2);
            if (apfc_is_real(&arg)) {
                apfx_pow(&result->re, &two, &arg.re);
                apf_zero(&result->im);
            } else {
                /* 2^(a+bi) = 2^a * (cos(b*ln2) + i*sin(b*ln2)) */
                apf mag, ln2, b_ln2;
                apfx_pow(&mag, &two, &arg.re);
                apf_from_int(&ln2, 2);
                apfx_log(&ln2, &ln2);
                apf_mul(&b_ln2, &arg.im, &ln2);
                apfx_cos(&result->re, &b_ln2);
                apfx_sin(&result->im, &b_ln2);
                apf_mul(&result->re, &result->re, &mag);
                apf_mul(&result->im, &result->im, &mag);
            }
        } else if (str_eq(name, "deg2rad")) {
            /* Convert degrees to radians: x * pi / 180 */
            apf pi_val, d180;
            if (!apfc_is_real(&arg)) {
                printf("Error: deg2rad requires real argument\n");
                return 0;
            }
            apfx_pi(&pi_val);
            apf_from_int(&d180, 180);
            apf_mul(&result->re, &arg.re, &pi_val);
            apf_div(&result->re, &result->re, &d180);
            apf_zero(&result->im);
        } else if (str_eq(name, "rad2deg")) {
            /* Convert radians to degrees: x * 180 / pi */
            apf pi_val, d180;
            if (!apfc_is_real(&arg)) {
                printf("Error: rad2deg requires real argument\n");
                return 0;
            }
            apfx_pi(&pi_val);
            apf_from_int(&d180, 180);
            apf_mul(&result->re, &arg.re, &d180);
            apf_div(&result->re, &result->re, &pi_val);
            apf_zero(&result->im);
        } else if (str_eq(name, "log1p")) {
            /* log1p(x) = log(1+x), accurate for small x */
            apf one, sum;
            if (!apfc_is_real(&arg)) {
                printf("Error: log1p requires real argument\n");
                return 0;
            }
            apf_from_int(&one, 1);
            apf_add(&sum, &one, &arg.re);
            apfx_log(&result->re, &sum);
            apf_zero(&result->im);
        } else if (str_eq(name, "expm1")) {
            /* expm1(x) = exp(x) - 1, accurate for small x */
            apf one;
            if (!apfc_is_real(&arg)) {
                printf("Error: expm1 requires real argument\n");
                return 0;
            }
            apfx_exp(&result->re, &arg.re);
            apf_from_int(&one, 1);
            apf_sub(&result->re, &result->re, &one);
            apf_zero(&result->im);
        } else if (str_eq(name, "wrapToPi")) {
            /* Wrap angle to [-pi, pi] */
            apf pi_val, two_pi, tmp;
            if (!apfc_is_real(&arg)) {
                printf("Error: wrap requires real argument\n");
                return 0;
            }
            apfx_pi(&pi_val);
            apf_add(&two_pi, &pi_val, &pi_val);
            /* result = x - 2*pi * floor((x + pi) / (2*pi)) */
            apf_add(&tmp, &arg.re, &pi_val);
            apf_div(&tmp, &tmp, &two_pi);
            apf_floor(&tmp, &tmp);
            apf_mul(&tmp, &tmp, &two_pi);
            apf_sub(&result->re, &arg.re, &tmp);
            apf_zero(&result->im);
        } else if (str_eq(name, "wrap360") || str_eq(name, "wrapTo360")) {
            /* Wrap angle to [0, 360) */
            apf d360, tmp;
            if (!apfc_is_real(&arg)) {
                printf("Error: wrap360 requires real argument\n");
                return 0;
            }
            apf_from_int(&d360, 360);
            apf_div(&tmp, &arg.re, &d360);
            apf_floor(&tmp, &tmp);
            apf_mul(&tmp, &tmp, &d360);
            apf_sub(&result->re, &arg.re, &tmp);
            /* Handle negative result */
            if (result->re.sign) {
                apf_add(&result->re, &result->re, &d360);
            }
            apf_zero(&result->im);
        } else if (str_eq(name, "wrap180") || str_eq(name, "wrapTo180")) {
            /* Wrap angle to [-180, 180] */
            apf d180, d360, tmp;
            if (!apfc_is_real(&arg)) {
                printf("Error: wrap180 requires real argument\n");
                return 0;
            }
            apf_from_int(&d180, 180);
            apf_from_int(&d360, 360);
            apf_add(&tmp, &arg.re, &d180);
            apf_div(&tmp, &tmp, &d360);
            apf_floor(&tmp, &tmp);
            apf_mul(&tmp, &tmp, &d360);
            apf_sub(&result->re, &arg.re, &tmp);
            apf_zero(&result->im);
        } else if (str_eq(name, "nextpow2")) {
            /* Next power of 2 >= |x| */
            long n, p;
            if (!apfc_is_real(&arg)) {
                printf("Error: nextpow2 requires real argument\n");
                return 0;
            }
            n = apf_to_long(&arg.re);
            if (n < 0) n = -n;
            if (n == 0) {
                apf_from_int(&result->re, 0);
            } else {
                p = 0;
                n--;
                while (n > 0) { n >>= 1; p++; }
                apf_from_int(&result->re, p);
            }
            apf_zero(&result->im);
        } else if (str_eq(name, "sqrt")) {
            apfc_sqrt(result, &arg);
        } else if (str_eq(name, "realsqrt")) {
            /* realsqrt: real square root, error if negative */
            if (!apfc_is_real(&arg)) {
                printf("Error: realsqrt requires real argument\n");
                return 0;
            }
            if (apf_cmp_int(&arg.re, 0) < 0) {
                printf("Error: realsqrt of negative number\n");
                return 0;
            }
            apf_sqrt(&result->re, &arg.re);
            apf_zero(&result->im);
        } else if (str_eq(name, "reallog")) {
            /* reallog: real logarithm, error if non-positive */
            if (!apfc_is_real(&arg)) {
                printf("Error: reallog requires real argument\n");
                return 0;
            }
            if (apf_cmp_int(&arg.re, 0) <= 0) {
                printf("Error: reallog of non-positive number\n");
                return 0;
            }
            apfx_log(&result->re, &arg.re);
            apf_zero(&result->im);
        } else if (str_eq(name, "cbrt")) {
            /* cbrt: cube root (preserves sign) */
            if (!apfc_is_real(&arg)) {
                printf("Error: cbrt requires real argument\n");
                return 0;
            }
            if (apf_cmp_int(&arg.re, 0) < 0) {
                apf pos, third;
                apf_neg(&pos, &arg.re);
                apf_from_str(&third, "0.333333333333333333333333333333333");
                apfx_pow(&result->re, &pos, &third);
                apf_neg(&result->re, &result->re);
            } else {
                apf third;
                apf_from_str(&third, "0.333333333333333333333333333333333");
                apfx_pow(&result->re, &arg.re, &third);
            }
            apf_zero(&result->im);
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
        } else if (str_eq(name, "arg") || str_eq(name, "angle") || str_eq(name, "phase")) {
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
        } else if (str_eq(name, "sech")) {
            /* sech(x) = 1/cosh(x) */
            apfc one, cosh_val;
            apf_from_int(&one.re, 1);
            apf_zero(&one.im);
            apfc_cosh(&cosh_val, &arg);
            apfc_div(result, &one, &cosh_val);
        } else if (str_eq(name, "csch")) {
            /* csch(x) = 1/sinh(x) */
            apfc one, sinh_val;
            apf_from_int(&one.re, 1);
            apf_zero(&one.im);
            apfc_sinh(&sinh_val, &arg);
            apfc_div(result, &one, &sinh_val);
        } else if (str_eq(name, "coth")) {
            /* coth(x) = cosh(x)/sinh(x) */
            apfc sinh_val, cosh_val;
            apfc_sinh(&sinh_val, &arg);
            apfc_cosh(&cosh_val, &arg);
            apfc_div(result, &cosh_val, &sinh_val);
        } else if (str_eq(name, "asech")) {
            /* asech(x) = acosh(1/x) */
            apf one, inv, tmp;
            apf_from_int(&one, 1);
            apf_div(&inv, &one, &arg.re);
            apfx_acosh(&tmp, &inv);
            result->re = tmp;
            apf_zero(&result->im);
        } else if (str_eq(name, "acsch")) {
            /* acsch(x) = asinh(1/x) */
            apf one, inv, tmp;
            apf_from_int(&one, 1);
            apf_div(&inv, &one, &arg.re);
            apfx_asinh(&tmp, &inv);
            result->re = tmp;
            apf_zero(&result->im);
        } else if (str_eq(name, "acoth")) {
            /* acoth(x) = atanh(1/x) */
            apf one, inv, tmp;
            apf_from_int(&one, 1);
            apf_div(&inv, &one, &arg.re);
            apfx_atanh(&tmp, &inv);
            result->re = tmp;
            apf_zero(&result->im);
        } else if (str_eq(name, "atan") || str_eq(name, "arctan")) {
            /* atan for real numbers only (complex atan not implemented)
             * But allow NaN through - should return NaN */
            if (!apfc_is_real(&arg) && arg.re.cls != APF_CLASS_NAN) {
                printf("Error: atan requires real argument\n");
                return 0;
            }
            apfx_atan(&result->re, &arg.re);
            if (angle_mode != ANGLE_RAD) {
                apf tmp;
                rad_to_angle(&tmp, &result->re, angle_mode);
                result->re = tmp;
            }
            apf_zero(&result->im);
        } else if (str_eq(name, "asin") || str_eq(name, "arcsin")) {
            /* asin for real numbers only, but allow NaN */
            if (!apfc_is_real(&arg) && arg.re.cls != APF_CLASS_NAN) {
                printf("Error: asin requires real argument\n");
                return 0;
            }
            apfx_asin(&result->re, &arg.re);
            if (angle_mode != ANGLE_RAD) {
                apf tmp;
                rad_to_angle(&tmp, &result->re, angle_mode);
                result->re = tmp;
            }
            apf_zero(&result->im);
        } else if (str_eq(name, "acos") || str_eq(name, "arccos")) {
            /* acos for real numbers only, but allow NaN */
            if (!apfc_is_real(&arg) && arg.re.cls != APF_CLASS_NAN) {
                printf("Error: acos requires real argument\n");
                return 0;
            }
            apfx_acos(&result->re, &arg.re);
            if (angle_mode != ANGLE_RAD) {
                apf tmp;
                rad_to_angle(&tmp, &result->re, angle_mode);
                result->re = tmp;
            }
            apf_zero(&result->im);
        } else if (str_eq(name, "tanh")) {
            /* tanh for real numbers only, but allow NaN */
            if (!apfc_is_real(&arg) && arg.re.cls != APF_CLASS_NAN) {
                printf("Error: tanh requires real argument\n");
                return 0;
            }
            apfx_tanh(&result->re, &arg.re);
            apf_zero(&result->im);
        } else if (str_eq(name, "asinh") || str_eq(name, "arcsinh")) {
            if (!apfc_is_real(&arg)) {
                printf("Error: asinh requires real argument\n");
                return 0;
            }
            apfx_asinh(&result->re, &arg.re);
            apf_zero(&result->im);
        } else if (str_eq(name, "acosh") || str_eq(name, "arccosh")) {
            if (!apfc_is_real(&arg)) {
                printf("Error: acosh requires real argument\n");
                return 0;
            }
            apfx_acosh(&result->re, &arg.re);
            apf_zero(&result->im);
        } else if (str_eq(name, "atanh") || str_eq(name, "arctanh")) {
            if (!apfc_is_real(&arg)) {
                printf("Error: atanh requires real argument\n");
                return 0;
            }
            apfx_atanh(&result->re, &arg.re);
            apf_zero(&result->im);
        } else if (str_eq(name, "sec")) {
            /* sec(x) = 1/cos(x) */
            if (!apfc_is_real(&arg)) {
                printf("Error: sec requires real argument\n");
                return 0;
            }
            if (angle_mode != ANGLE_RAD) {
                apf rad_angle;
                angle_to_rad(&rad_angle, &arg.re, angle_mode);
                apfx_sec(&result->re, &rad_angle);
            } else {
                apfx_sec(&result->re, &arg.re);
            }
            apf_zero(&result->im);
        } else if (str_eq(name, "csc") || str_eq(name, "cosec")) {
            /* csc(x) = 1/sin(x) */
            if (!apfc_is_real(&arg)) {
                printf("Error: csc requires real argument\n");
                return 0;
            }
            if (angle_mode != ANGLE_RAD) {
                apf rad_angle;
                angle_to_rad(&rad_angle, &arg.re, angle_mode);
                apfx_csc(&result->re, &rad_angle);
            } else {
                apfx_csc(&result->re, &arg.re);
            }
            apf_zero(&result->im);
        } else if (str_eq(name, "cot") || str_eq(name, "cotan")) {
            /* cot(x) = 1/tan(x) */
            if (!apfc_is_real(&arg)) {
                printf("Error: cot requires real argument\n");
                return 0;
            }
            if (angle_mode != ANGLE_RAD) {
                apf rad_angle;
                angle_to_rad(&rad_angle, &arg.re, angle_mode);
                apfx_cot(&result->re, &rad_angle);
            } else {
                apfx_cot(&result->re, &arg.re);
            }
            apf_zero(&result->im);
        } else if (str_eq(name, "asec") || str_eq(name, "arcsec")) {
            /* arcsec(x) = acos(1/x) */
            if (!apfc_is_real(&arg)) {
                printf("Error: asec requires real argument\n");
                return 0;
            }
            apfx_asec(&result->re, &arg.re);
            if (angle_mode != ANGLE_RAD) {
                apf tmp;
                rad_to_angle(&tmp, &result->re, angle_mode);
                result->re = tmp;
            }
            apf_zero(&result->im);
        } else if (str_eq(name, "acsc") || str_eq(name, "arccsc") || str_eq(name, "arccosec")) {
            /* arccsc(x) = asin(1/x) */
            if (!apfc_is_real(&arg)) {
                printf("Error: acsc requires real argument\n");
                return 0;
            }
            apfx_acsc(&result->re, &arg.re);
            if (angle_mode != ANGLE_RAD) {
                apf tmp;
                rad_to_angle(&tmp, &result->re, angle_mode);
                result->re = tmp;
            }
            apf_zero(&result->im);
        } else if (str_eq(name, "acot") || str_eq(name, "arccot") || str_eq(name, "arccotan")) {
            /* arccot(x) = atan(1/x) */
            if (!apfc_is_real(&arg)) {
                printf("Error: acot requires real argument\n");
                return 0;
            }
            apfx_acot(&result->re, &arg.re);
            if (angle_mode != ANGLE_RAD) {
                apf tmp;
                rad_to_angle(&tmp, &result->re, angle_mode);
                result->re = tmp;
            }
            apf_zero(&result->im);
        } else if (str_eq(name, "fact") || str_eq(name, "factorial")) {
            long n;
            if (!apfc_is_real(&arg)) {
                printf("Error: factorial requires real integer\n");
                return 0;
            }
            n = apf_to_long(&arg.re);
            if (n < 0 || n > 100000) {
                printf("Error: factorial out of range\n");
                return 0;
            }
            apfx_fact(&result->re, n);
            apf_zero(&result->im);
        } else if (str_eq(name, "factorial2") || str_eq(name, "fact2")) {
            /* Double factorial: n!! = n*(n-2)*(n-4)*... */
            long n, i;
            apf prod, tmp, step;
            if (!apfc_is_real(&arg)) {
                printf("Error: factorial2 requires real integer\n");
                return 0;
            }
            n = apf_to_long(&arg.re);
            if (n < 0 || n > 10000) {
                printf("Error: factorial2 out of range\n");
                return 0;
            }
            apf_from_int(&prod, 1);
            apf_from_int(&step, 2);
            for (i = n; i > 1; i -= 2) {
                apf_from_int(&tmp, i);
                apf_mul(&prod, &prod, &tmp);
            }
            result->re = prod;
            apf_zero(&result->im);
        } else if (str_eq(name, "fibonacci") || str_eq(name, "fib")) {
            /* Fibonacci number F(n) */
            long n, i;
            apf a, b, tmp;
            if (!apfc_is_real(&arg)) {
                printf("Error: fibonacci requires real integer\n");
                return 0;
            }
            n = apf_to_long(&arg.re);
            if (n < 0 || n > 10000) {
                printf("Error: fibonacci out of range\n");
                return 0;
            }
            if (n == 0) {
                apf_from_int(&result->re, 0);
            } else {
                apf_from_int(&a, 0);
                apf_from_int(&b, 1);
                for (i = 2; i <= n; i++) {
                    apf_add(&tmp, &a, &b);
                    a = b;
                    b = tmp;
                }
                result->re = b;
            }
            apf_zero(&result->im);
        } else if (str_eq(name, "lucas")) {
            /* Lucas number L(n): L(0)=2, L(1)=1, L(n)=L(n-1)+L(n-2) */
            long n, i;
            apf a, b, tmp;
            if (!apfc_is_real(&arg)) {
                printf("Error: lucas requires real integer\n");
                return 0;
            }
            n = apf_to_long(&arg.re);
            if (n < 0 || n > 10000) {
                printf("Error: lucas out of range\n");
                return 0;
            }
            if (n == 0) {
                apf_from_int(&result->re, 2);
            } else if (n == 1) {
                apf_from_int(&result->re, 1);
            } else {
                apf_from_int(&a, 2);
                apf_from_int(&b, 1);
                for (i = 2; i <= n; i++) {
                    apf_add(&tmp, &a, &b);
                    a = b;
                    b = tmp;
                }
                result->re = b;
            }
            apf_zero(&result->im);
        } else if (str_eq(name, "catalan")) {
            /* Catalan number C(n) = (2n)! / ((n+1)! * n!) */
            long n, i;
            apf num, tmp;
            if (!apfc_is_real(&arg)) {
                printf("Error: catalan requires real integer\n");
                return 0;
            }
            n = apf_to_long(&arg.re);
            if (n < 0 || n > 500) {
                printf("Error: catalan out of range\n");
                return 0;
            }
            /* C(n) = prod(k=2 to n) of (n+k)/k */
            apf_from_int(&num, 1);
            for (i = 2; i <= n; i++) {
                apf_from_int(&tmp, n + i);
                apf_mul(&num, &num, &tmp);
                apf_from_int(&tmp, i);
                apf_div(&num, &num, &tmp);
            }
            result->re = num;
            apf_zero(&result->im);
        } else if (str_eq(name, "gamma") || str_eq(name, "tgamma")) {
            /* gamma function: Gamma(x) */
            if (!apfc_is_real(&arg)) {
                printf("Error: gamma requires real argument\n");
                return 0;
            }
            apfx_tgamma(&result->re, &arg.re);
            apf_zero(&result->im);
        } else if (str_eq(name, "lgamma")) {
            /* log-gamma: ln(|Gamma(x)|) */
            if (!apfc_is_real(&arg)) {
                printf("Error: lgamma requires real argument\n");
                return 0;
            }
            apfx_lgamma(&result->re, &arg.re);
            apf_zero(&result->im);
        } else if (str_eq(name, "erf")) {
            /* error function */
            if (!apfc_is_real(&arg)) {
                printf("Error: erf requires real argument\n");
                return 0;
            }
            apfx_erf(&result->re, &arg.re);
            apf_zero(&result->im);
        } else if (str_eq(name, "erfc")) {
            /* complementary error function: 1-erf(x) */
            if (!apfc_is_real(&arg)) {
                printf("Error: erfc requires real argument\n");
                return 0;
            }
            apfx_erfc(&result->re, &arg.re);
            apf_zero(&result->im);
        } else if (str_eq(name, "digamma") || str_eq(name, "psi")) {
            /* digamma function: psi(x) = d/dx ln(Gamma(x)) */
            if (!apfc_is_real(&arg)) {
                printf("Error: digamma requires real argument\n");
                return 0;
            }
            apfx_digamma(&result->re, &arg.re);
            apf_zero(&result->im);
        } else if (str_eq(name, "zeta") || str_eq(name, "riemann_zeta")) {
            /* Riemann zeta function */
            if (!apfc_is_real(&arg)) {
                printf("Error: zeta requires real argument\n");
                return 0;
            }
            apfx_zeta(&result->re, &arg.re);
            apf_zero(&result->im);
        } else if (str_eq(name, "harmonic")) {
            /* Harmonic number H_n */
            long n = apf_to_long(&arg.re);
            if (n < 0 || n > 10000) {
                printf("Error: harmonic argument out of range\n");
                return 0;
            }
            apfx_harmonic(&result->re, n);
            apf_zero(&result->im);
        } else if (str_eq(name, "sinc")) {
            /* sinc(x) = sin(pi*x)/(pi*x), with sinc(0) = 1 */
            if (apf_is_zero(&arg.re)) {
                apf_from_int(&result->re, 1);
            } else {
                apf pi_x, sin_val;
                apfx_pi(&pi_x);
                apf_mul(&pi_x, &pi_x, &arg.re);
                apfx_sin(&sin_val, &pi_x);
                apf_div(&result->re, &sin_val, &pi_x);
            }
            apf_zero(&result->im);
        } else if (str_eq(name, "sigmoid") || str_eq(name, "logistic")) {
            /* sigmoid(x) = 1 / (1 + exp(-x)) */
            apf neg_x, exp_val, one, denom;
            apf_neg(&neg_x, &arg.re);
            apfx_exp(&exp_val, &neg_x);
            apf_from_int(&one, 1);
            apf_add(&denom, &one, &exp_val);
            apf_div(&result->re, &one, &denom);
            apf_zero(&result->im);
        } else if (str_eq(name, "softplus")) {
            /* softplus(x) = ln(1 + exp(x)) */
            apf exp_val, one, sum_val;
            apfx_exp(&exp_val, &arg.re);
            apf_from_int(&one, 1);
            apf_add(&sum_val, &one, &exp_val);
            apfx_log(&result->re, &sum_val);
            apf_zero(&result->im);
        } else if (str_eq(name, "step") || str_eq(name, "heaviside")) {
            /* Heaviside step: 0 if x<0, 0.5 if x=0, 1 if x>0 */
            if (apf_is_zero(&arg.re)) {
                apf_from_str(&result->re, "0.5");
            } else if (arg.re.sign) {
                apf_zero(&result->re);
            } else {
                apf_from_int(&result->re, 1);
            }
            apf_zero(&result->im);
        } else if (str_eq(name, "rect") || str_eq(name, "rectangle")) {
            /* rect(x) = 1 if |x| < 0.5, 0.5 if |x| = 0.5, 0 otherwise */
            apf abs_x, half;
            int cmp;
            apf_abs(&abs_x, &arg.re);
            apf_from_str(&half, "0.5");
            cmp = apf_cmp(&abs_x, &half);
            if (cmp < 0) {
                apf_from_int(&result->re, 1);
            } else if (cmp == 0) {
                apf_from_str(&result->re, "0.5");
            } else {
                apf_zero(&result->re);
            }
            apf_zero(&result->im);
        } else if (str_eq(name, "tri") || str_eq(name, "triangle")) {
            /* tri(x) = max(0, 1 - |x|) */
            apf abs_x, one, diff;
            apf_abs(&abs_x, &arg.re);
            apf_from_int(&one, 1);
            apf_sub(&diff, &one, &abs_x);
            if (diff.sign) {
                apf_zero(&result->re);
            } else {
                apf_copy(&result->re, &diff);
            }
            apf_zero(&result->im);
        } else if (str_eq(name, "floor")) {
            if (!apfc_is_real(&arg)) {
                printf("Error: floor requires real argument\n");
                return 0;
            }
            apf_floor(&result->re, &arg.re);
            apf_zero(&result->im);
        } else if (str_eq(name, "ceil")) {
            if (!apfc_is_real(&arg)) {
                printf("Error: ceil requires real argument\n");
                return 0;
            }
            apf_ceil(&result->re, &arg.re);
            apf_zero(&result->im);
        } else if (str_eq(name, "trunc") || str_eq(name, "int") || str_eq(name, "fix")) {
            if (!apfc_is_real(&arg)) {
                printf("Error: trunc requires real argument\n");
                return 0;
            }
            apf_trunc(&result->re, &arg.re);
            apf_zero(&result->im);
        } else if (str_eq(name, "sign") || str_eq(name, "signum")) {
            /* sign: -1 for negative, 0 for zero, 1 for positive */
            if (apf_is_zero(&arg.re)) {
                apf_zero(&result->re);
            } else if (arg.re.sign) {
                apf_from_int(&result->re, -1);
            } else {
                apf_from_int(&result->re, 1);
            }
            apf_zero(&result->im);
        } else if (str_eq(name, "round")) {
            /* round: round to nearest integer (ties away from zero) */
            apf half, tmp;
            if (!apfc_is_real(&arg)) {
                printf("Error: round requires real argument\n");
                return 0;
            }
            apf_from_str(&half, "0.5");
            if (arg.re.sign) {
                apf_sub(&tmp, &arg.re, &half);
                apf_ceil(&result->re, &tmp);
            } else {
                apf_add(&tmp, &arg.re, &half);
                apf_floor(&result->re, &tmp);
            }
            apf_zero(&result->im);
        } else if (str_eq(name, "isnan")) {
            /* isnan: check if value is NaN */
            int is_nan = (arg.re.cls == APF_CLASS_NAN);
            apf_from_int(&result->re, is_nan);
            apf_zero(&result->im);
            result_is_boolean = 1;
        } else if (str_eq(name, "isinf")) {
            /* isinf: check if value is infinite */
            int is_inf = (arg.re.cls == APF_CLASS_INF);
            apf_from_int(&result->re, is_inf);
            apf_zero(&result->im);
            result_is_boolean = 1;
        } else if (str_eq(name, "isfinite")) {
            /* isfinite: check if value is normal or zero (not NaN/Inf) */
            int is_finite = (arg.re.cls == APF_CLASS_NORMAL || arg.re.cls == APF_CLASS_ZERO);
            apf_from_int(&result->re, is_finite);
            apf_zero(&result->im);
            result_is_boolean = 1;
        } else if (str_eq(name, "isreal")) {
            /* isreal: check if imaginary part is zero */
            int is_real = apf_is_zero(&arg.im);
            apf_from_int(&result->re, is_real);
            apf_zero(&result->im);
            result_is_boolean = 1;
        } else if (str_eq(name, "isnumeric")) {
            /* isnumeric: always true (everything is numeric in calculator) */
            apf_from_int(&result->re, 1);
            apf_zero(&result->im);
            result_is_boolean = 1;
        } else if (str_eq(name, "isinteger")) {
            /* isinteger: check if value is an integer (real and no fractional part) */
            int is_int = 0;
            if (apfc_is_real(&arg) && arg.re.cls == APF_CLASS_NORMAL) {
                apf truncated, diff;
                apf_trunc(&truncated, &arg.re);
                apf_sub(&diff, &arg.re, &truncated);
                is_int = apf_is_zero(&diff);
            } else if (arg.re.cls == APF_CLASS_ZERO) {
                is_int = 1;
            }
            apf_from_int(&result->re, is_int);
            apf_zero(&result->im);
            result_is_boolean = 1;
        } else if (str_eq(name, "islogical")) {
            /* islogical: check if value is boolean (0 or 1) */
            int is_log = 0;
            if (apfc_is_real(&arg)) {
                if (arg.re.cls == APF_CLASS_ZERO || 
                    (arg.re.cls == APF_CLASS_NORMAL && apf_cmp_int(&arg.re, 1) == 0)) {
                    is_log = 1;
                }
            }
            apf_from_int(&result->re, is_log);
            apf_zero(&result->im);
            result_is_boolean = 1;
        } else if (str_eq(name, "ispositive")) {
            /* ispositive: check if value is strictly positive */
            int is_pos = apfc_is_real(&arg) && apf_cmp_int(&arg.re, 0) > 0;
            apf_from_int(&result->re, is_pos);
            apf_zero(&result->im);
            result_is_boolean = 1;
        } else if (str_eq(name, "isnegative")) {
            /* isnegative: check if value is strictly negative */
            int is_neg = apfc_is_real(&arg) && apf_cmp_int(&arg.re, 0) < 0;
            apf_from_int(&result->re, is_neg);
            apf_zero(&result->im);
            result_is_boolean = 1;
        } else if (str_eq(name, "frac")) {
            /* Fractional part: x - trunc(x) */
            apf truncated;
            if (!apfc_is_real(&arg)) {
                printf("Error: frac requires real argument\n");
                return 0;
            }
            apf_trunc(&truncated, &arg.re);
            apf_sub(&result->re, &arg.re, &truncated);
            apf_zero(&result->im);
        } else if (str_eq(name, "even")) {
            /* Returns 1 if integer is even, 0 otherwise */
            long n;
            if (!apfc_is_real(&arg)) {
                printf("Error: even requires real integer\n");
                return 0;
            }
            n = apf_to_long(&arg.re);
            apf_from_int(&result->re, (n % 2 == 0) ? 1 : 0);
            apf_zero(&result->im);
            result_is_boolean = 1;
        } else if (str_eq(name, "odd")) {
            /* Returns 1 if integer is odd, 0 otherwise */
            long n;
            if (!apfc_is_real(&arg)) {
                printf("Error: odd requires real integer\n");
                return 0;
            }
            n = apf_to_long(&arg.re);
            apf_from_int(&result->re, (n % 2 != 0) ? 1 : 0);
            apf_zero(&result->im);
            result_is_boolean = 1;
        } else if (str_eq(name, "isprime")) {
            /* Miller-Rabin primality test */
            long n, d, r, a, x, i, j;
            int witnesses[] = {2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37};
            int is_prime = 1;
            if (!apfc_is_real(&arg)) {
                printf("Error: isprime requires real integer\n");
                return 0;
            }
            n = apf_to_long(&arg.re);
            if (n < 2) {
                is_prime = 0;
            } else if (n == 2 || n == 3) {
                is_prime = 1;
            } else if (n % 2 == 0) {
                is_prime = 0;
            } else {
                /* Write n-1 as 2^r * d */
                d = n - 1;
                r = 0;
                while (d % 2 == 0) {
                    d /= 2;
                    r++;
                }
                /* Miller-Rabin witnesses */
                for (i = 0; i < 12 && is_prime; i++) {
                    a = witnesses[i];
                    if (a >= n) continue;
                    /* Compute x = a^d mod n */
                    x = 1;
                    {
                        long base = a % n;
                        long exp = d;
                        while (exp > 0) {
                            if (exp % 2 == 1) {
                                x = (x * base) % n;
                            }
                            base = (base * base) % n;
                            exp /= 2;
                        }
                    }
                    if (x == 1 || x == n - 1) continue;
                    for (j = 0; j < r - 1; j++) {
                        x = (x * x) % n;
                        if (x == n - 1) break;
                    }
                    if (x != n - 1) {
                        is_prime = 0;
                    }
                }
            }
            apf_from_int(&result->re, is_prime);
            apf_zero(&result->im);
            result_is_boolean = 1;
        } else if (str_eq(name, "mod")) {
            /* mod(a, b) - two argument function */
            printf("Error: mod requires two arguments: mod(a, b)\n");
            return 0;
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
        } else if (str_eq(name, "ieee") || str_eq(name, "explain") || str_eq(name, "breakdown")) {
            apfc_ieee_display(&arg);
            *result = arg;
#ifdef HAVE_GCD
        } else if (str_eq(name, "isprime")) {
            long n = apf_to_long(&arg.re);
            apf_from_int(&result->re, is_prime_long(n));
            apf_zero(&result->im);
            result_is_boolean = 1;
        } else if (str_eq(name, "nthprime")) {
            long n = apf_to_long(&arg.re);
            if (n < 1 || n > 100000) {
                printf("Error: nthprime argument out of range (1-100000)\n");
                return 0;
            }
            apf_from_int(&result->re, nthprime_long(n));
            apf_zero(&result->im);
        } else if (str_eq(name, "nextprime")) {
            long n = apf_to_long(&arg.re);
            apf_from_int(&result->re, nextprime_long(n));
            apf_zero(&result->im);
        } else if (str_eq(name, "prevprime")) {
            long n = apf_to_long(&arg.re);
            apf_from_int(&result->re, prevprime_long(n));
            apf_zero(&result->im);
        } else if (str_eq(name, "primepi")) {
            long n = apf_to_long(&arg.re);
            apf_from_int(&result->re, primepi_long(n));
            apf_zero(&result->im);
        } else if (str_eq(name, "radical")) {
            long n = apf_to_long(&arg.re);
            apf_from_int(&result->re, radical_long(n));
            apf_zero(&result->im);
        } else if (str_eq(name, "omega")) {
            long n = apf_to_long(&arg.re);
            apf_from_int(&result->re, omega_long(n));
            apf_zero(&result->im);
        } else if (str_eq(name, "bigomega") || str_eq(name, "Omega")) {
            long n = apf_to_long(&arg.re);
            apf_from_int(&result->re, bigomega_long(n));
            apf_zero(&result->im);
        } else if (str_eq(name, "issquarefree")) {
            long n = apf_to_long(&arg.re);
            apf_from_int(&result->re, issquarefree_long(n));
            apf_zero(&result->im);
            result_is_boolean = 1;
        } else if (str_eq(name, "moebius") || str_eq(name, "mobius") || str_eq(name, "mu")) {
            long n = apf_to_long(&arg.re);
            apf_from_int(&result->re, moebius_long(n));
            apf_zero(&result->im);
        } else if (str_eq(name, "sigma") || str_eq(name, "divisorsum")) {
            long n = apf_to_long(&arg.re);
            apf_from_int(&result->re, sigma_long(n, 1));
            apf_zero(&result->im);
        } else if (str_eq(name, "sigma0") || str_eq(name, "numdivisors") || str_eq(name, "tau")) {
            long n = apf_to_long(&arg.re);
            apf_from_int(&result->re, numdivisors_long(n));
            apf_zero(&result->im);
        } else if (str_eq(name, "digsum") || str_eq(name, "digitsum")) {
            long n = apf_to_long(&arg.re);
            apf_from_int(&result->re, digsum_long(n));
            apf_zero(&result->im);
        } else if (str_eq(name, "numdigits") || str_eq(name, "ndigits")) {
            long n = apf_to_long(&arg.re);
            apf_from_int(&result->re, numdigits_long(n));
            apf_zero(&result->im);
        } else if (str_eq(name, "digitalroot") || str_eq(name, "digroot")) {
            long n = apf_to_long(&arg.re);
            apf_from_int(&result->re, digitalroot_long(n));
            apf_zero(&result->im);
        } else if (str_eq(name, "reversedigits") || str_eq(name, "revdigits")) {
            long n = apf_to_long(&arg.re);
            apf_from_int(&result->re, reverse_long(n));
            apf_zero(&result->im);
        } else if (str_eq(name, "ispalindrome")) {
            long n = apf_to_long(&arg.re);
            apf_from_int(&result->re, ispalindrome_long(n));
            apf_zero(&result->im);
            result_is_boolean = 1;
        } else if (str_eq(name, "isperfect")) {
            long n = apf_to_long(&arg.re);
            apf_from_int(&result->re, isperfect_long(n));
            apf_zero(&result->im);
            result_is_boolean = 1;
        } else if (str_eq(name, "isabundant")) {
            long n = apf_to_long(&arg.re);
            apf_from_int(&result->re, isabundant_long(n));
            apf_zero(&result->im);
            result_is_boolean = 1;
        } else if (str_eq(name, "isdeficient")) {
            long n = apf_to_long(&arg.re);
            apf_from_int(&result->re, isdeficient_long(n));
            apf_zero(&result->im);
            result_is_boolean = 1;
        } else if (str_eq(name, "triangular") || str_eq(name, "tri_num")) {
            long n = apf_to_long(&arg.re);
            apf_from_int(&result->re, triangular_long(n));
            apf_zero(&result->im);
        } else if (str_eq(name, "istriangular")) {
            long n = apf_to_long(&arg.re);
            apf_from_int(&result->re, istriangular_long(n));
            apf_zero(&result->im);
            result_is_boolean = 1;
        } else if (str_eq(name, "isperfectsquare")) {
            long n = apf_to_long(&arg.re);
            apf_from_int(&result->re, issquare_long(n));
            apf_zero(&result->im);
            result_is_boolean = 1;
        } else if (str_eq(name, "pentagonal") || str_eq(name, "pent_num")) {
            long n = apf_to_long(&arg.re);
            apf_from_int(&result->re, pentagonal_long(n));
            apf_zero(&result->im);
        } else if (str_eq(name, "hexagonal") || str_eq(name, "hex_num")) {
            long n = apf_to_long(&arg.re);
            apf_from_int(&result->re, hexagonal_long(n));
            apf_zero(&result->im);
        } else if (str_eq(name, "ispower") || str_eq(name, "isperfectpower")) {
            long n = apf_to_long(&arg.re);
            apf_from_int(&result->re, ispower_long(n));
            apf_zero(&result->im);
            result_is_boolean = 1;
        } else if (str_eq(name, "subfactorial") || str_eq(name, "derangements")) {
            long n = apf_to_long(&arg.re);
            if (n < 0 || n > 100) {
                printf("Error: subfactorial out of range (0-100)\n");
                return 0;
            }
            apf_subfactorial(&result->re, n);
            apf_zero(&result->im);
        } else if (str_eq(name, "bell")) {
            long n = apf_to_long(&arg.re);
            if (n < 0 || n > 100) {
                printf("Error: bell out of range (0-100)\n");
                return 0;
            }
            apf_bell(&result->re, n);
            apf_zero(&result->im);
        } else if (str_eq(name, "partition") || str_eq(name, "partitions")) {
            long n = apf_to_long(&arg.re);
            if (n < 0 || n > 1000) {
                printf("Error: partition out of range (0-1000)\n");
                return 0;
            }
            apf_partition(&result->re, n);
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
        } else if (str_eq(name, "linspace")) {
            /* linspace(start, end, n) */
            int n, i;
            apf start, end_val, step, val;
            
            start = arg1.re;
            end_val = arg2.re;
            n = 10;  /* default */
            
            if (current_token.type == TOK_COMMA) {
                next_token();
                if (!parse_expr(&arg3)) return 0;
                n = apf_to_long(&arg3.re);
            }
            
            if (current_token.type != TOK_RPAREN) {
                printf("Error: expected ')'\n");
                return 0;
            }
            next_token();
            
            if (n < 1) n = 1;
            if (n > MAT_MAX_COLS) n = MAT_MAX_COLS;
            
            result->type = VAL_MATRIX;
            result->v.matrix.rows = 1;
            result->v.matrix.cols = n;
            
            if (n == 1) {
                apf_copy(&MAT_AT(&result->v.matrix, 0, 0).re, &end_val);
                apf_zero(&MAT_AT(&result->v.matrix, 0, 0).im);
            } else {
                apf n_minus_1;
                apf_from_int(&n_minus_1, n - 1);
                apf_sub(&step, &end_val, &start);
                apf_div(&step, &step, &n_minus_1);
                
                for (i = 0; i < n; i++) {
                    apf idx;
                    apf_from_int(&idx, i);
                    apf_mul(&val, &step, &idx);
                    apf_add(&MAT_AT(&result->v.matrix, 0, i).re, &start, &val);
                    apf_zero(&MAT_AT(&result->v.matrix, 0, i).im);
                }
            }
            return 1;
        } else if (str_eq(name, "logspace")) {
            /* logspace(a, b, n) - n logarithmically spaced points from 10^a to 10^b */
            int n, i;
            apf start_exp, end_exp, step, ten, exp_val;
            
            start_exp = arg1.re;
            end_exp = arg2.re;
            n = 50;  /* default */
            
            if (current_token.type == TOK_COMMA) {
                next_token();
                if (!parse_expr(&arg3)) return 0;
                n = apf_to_long(&arg3.re);
            }
            
            if (current_token.type != TOK_RPAREN) {
                printf("Error: expected ')'\n");
                return 0;
            }
            next_token();
            
            if (n < 1) n = 1;
            if (n > MAT_MAX_COLS) n = MAT_MAX_COLS;
            
            result->type = VAL_MATRIX;
            result->v.matrix.rows = 1;
            result->v.matrix.cols = n;
            
            apf_from_int(&ten, 10);
            
            if (n == 1) {
                apfx_pow(&MAT_AT(&result->v.matrix, 0, 0).re, &ten, &end_exp);
                apf_zero(&MAT_AT(&result->v.matrix, 0, 0).im);
            } else {
                apf n_minus_1;
                apf_from_int(&n_minus_1, n - 1);
                apf_sub(&step, &end_exp, &start_exp);
                apf_div(&step, &step, &n_minus_1);
                
                for (i = 0; i < n; i++) {
                    apf idx, cur_exp;
                    apf_from_int(&idx, i);
                    apf_mul(&cur_exp, &step, &idx);
                    apf_add(&exp_val, &start_exp, &cur_exp);
                    apfx_pow(&MAT_AT(&result->v.matrix, 0, i).re, &ten, &exp_val);
                    apf_zero(&MAT_AT(&result->v.matrix, 0, i).im);
                }
            }
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
    } else if (str_eq(name, "magic")) {
        /* magic(n) - magic square where all rows, cols, diagonals sum to same value */
        int n = rows;
        int i, j, num, r, c;
        if (n < 1) n = 1;
        if (n > MAT_MAX_ROWS) n = MAT_MAX_ROWS;
        mat_zero(&result->v.matrix, n, n);
        if (n == 1) {
            apf_from_int(&MAT_AT(&result->v.matrix, 0, 0).re, 1);
        } else if (n == 2) {
            /* 2x2 has no magic square, use simple version */
            apf_from_int(&MAT_AT(&result->v.matrix, 0, 0).re, 1);
            apf_from_int(&MAT_AT(&result->v.matrix, 0, 1).re, 2);
            apf_from_int(&MAT_AT(&result->v.matrix, 1, 0).re, 3);
            apf_from_int(&MAT_AT(&result->v.matrix, 1, 1).re, 4);
        } else if (n % 2 == 1) {
            /* Odd order: Siamese method */
            int newr, newc;
            r = 0; c = n / 2;
            for (num = 1; num <= n * n; num++) {
                apf_from_int(&MAT_AT(&result->v.matrix, r, c).re, num);
                newr = (r - 1 + n) % n;
                newc = (c + 1) % n;
                if (!apf_is_zero(&MAT_AT(&result->v.matrix, newr, newc).re)) {
                    r = (r + 1) % n;
                } else {
                    r = newr; c = newc;
                }
            }
        } else if (n % 4 == 0) {
            /* Doubly even (4k) */
            num = 1;
            for (i = 0; i < n; i++) {
                for (j = 0; j < n; j++) {
                    apf_from_int(&MAT_AT(&result->v.matrix, i, j).re, num++);
                }
            }
            for (i = 0; i < n; i++) {
                for (j = 0; j < n; j++) {
                    int ii = i % 4, jj = j % 4;
                    if ((ii == 0 || ii == 3) == (jj == 0 || jj == 3)) {
                        int val = n*n + 1 - apf_to_long(&MAT_AT(&result->v.matrix, i, j).re);
                        apf_from_int(&MAT_AT(&result->v.matrix, i, j).re, val);
                    }
                }
            }
        } else {
            /* Singly even (4k+2): LUX method simplified */
            int half = n / 2;
            int k = (n - 2) / 4;
            matrix_t sub;
            int newr, newc;
            /* Build magic square for half size recursively */
            mat_zero(&sub, half, half);
            r = 0; c = half / 2;
            for (num = 1; num <= half * half; num++) {
                apf_from_int(&MAT_AT(&sub, r, c).re, num);
                newr = (r - 1 + half) % half;
                newc = (c + 1) % half;
                if (!apf_is_zero(&MAT_AT(&sub, newr, newc).re)) {
                    r = (r + 1) % half;
                } else {
                    r = newr; c = newc;
                }
            }
            /* Fill quadrants */
            for (i = 0; i < half; i++) {
                for (j = 0; j < half; j++) {
                    int v = apf_to_long(&MAT_AT(&sub, i, j).re);
                    apf_from_int(&MAT_AT(&result->v.matrix, i, j).re, v);
                    apf_from_int(&MAT_AT(&result->v.matrix, i + half, j + half).re, v + half*half);
                    apf_from_int(&MAT_AT(&result->v.matrix, i, j + half).re, v + 2*half*half);
                    apf_from_int(&MAT_AT(&result->v.matrix, i + half, j).re, v + 3*half*half);
                }
            }
            /* Swap columns */
            for (i = 0; i < half; i++) {
                for (j = 0; j < k; j++) {
                    int jj = (i == half/2) ? (j + k) : j;
                    apfc tmp = MAT_AT(&result->v.matrix, i, jj);
                    MAT_AT(&result->v.matrix, i, jj) = MAT_AT(&result->v.matrix, i + half, jj);
                    MAT_AT(&result->v.matrix, i + half, jj) = tmp;
                }
                for (j = n - k + 1; j < n; j++) {
                    apfc tmp = MAT_AT(&result->v.matrix, i, j);
                    MAT_AT(&result->v.matrix, i, j) = MAT_AT(&result->v.matrix, i + half, j);
                    MAT_AT(&result->v.matrix, i + half, j) = tmp;
                }
            }
        }
    } else if (str_eq(name, "pascal")) {
        /* pascal(n) - Pascal's triangle as lower triangular matrix */
        int n = rows;
        int i, j;
        if (n < 1) n = 1;
        if (n > MAT_MAX_ROWS) n = MAT_MAX_ROWS;
        mat_zero(&result->v.matrix, n, n);
        for (i = 0; i < n; i++) {
            apf_from_int(&MAT_AT(&result->v.matrix, i, 0).re, 1);
            for (j = 1; j <= i; j++) {
                /* C(i,j) = C(i-1,j-1) + C(i-1,j) */
                apf_add(&MAT_AT(&result->v.matrix, i, j).re,
                       &MAT_AT(&result->v.matrix, i-1, j-1).re,
                       &MAT_AT(&result->v.matrix, i-1, j).re);
            }
        }
    } else if (str_eq(name, "hilb")) {
        /* hilb(n) - Hilbert matrix: H(i,j) = 1/(i+j-1) */
        int n = rows;
        int i, j;
        if (n < 1) n = 1;
        if (n > MAT_MAX_ROWS) n = MAT_MAX_ROWS;
        mat_zero(&result->v.matrix, n, n);
        for (i = 0; i < n; i++) {
            for (j = 0; j < n; j++) {
                apf denom, one;
                apf_from_int(&denom, i + j + 1);
                apf_from_int(&one, 1);
                apf_div(&MAT_AT(&result->v.matrix, i, j).re, &one, &denom);
            }
        }
    } else if (str_eq(name, "invhilb")) {
        /* invhilb(n) - inverse Hilbert matrix (exact integer entries) */
        int n = rows;
        int i, j;
        if (n < 1) n = 1;
        if (n > MAT_MAX_ROWS) n = MAT_MAX_ROWS;
        mat_zero(&result->v.matrix, n, n);
        for (i = 0; i < n; i++) {
            for (j = 0; j < n; j++) {
                /* Formula: H^-1(i,j) = (-1)^(i+j) * (i+j+1) * C(n+i,n-j-1) * C(n+j,n-i-1) * C(i+j,i)^2 */
                apf sign, term, c1, c2, c3, tmp;
                int ii = i + 1, jj = j + 1;  /* 1-based */
                apf_from_int(&sign, ((i + j) % 2 == 0) ? 1 : -1);
                apf_from_int(&term, ii + jj - 1);
                
                /* Compute binomial coefficients using nCr */
                /* C(n+i-1, n-j-1) = C(n+ii-1, n-jj) */
                {
                    int top = n + ii - 1, bot = n - jj;
                    apf numer, denom;
                    int k;
                    apf_from_int(&numer, 1);
                    apf_from_int(&denom, 1);
                    for (k = 0; k < bot && k < top; k++) {
                        apf t;
                        apf_from_int(&t, top - k);
                        apf_mul(&numer, &numer, &t);
                        apf_from_int(&t, k + 1);
                        apf_mul(&denom, &denom, &t);
                    }
                    apf_div(&c1, &numer, &denom);
                }
                {
                    int top = n + jj - 1, bot = n - ii;
                    apf numer, denom;
                    int k;
                    apf_from_int(&numer, 1);
                    apf_from_int(&denom, 1);
                    for (k = 0; k < bot && k < top; k++) {
                        apf t;
                        apf_from_int(&t, top - k);
                        apf_mul(&numer, &numer, &t);
                        apf_from_int(&t, k + 1);
                        apf_mul(&denom, &denom, &t);
                    }
                    apf_div(&c2, &numer, &denom);
                }
                {
                    int top = ii + jj - 2, bot = ii - 1;
                    apf numer, denom;
                    int k;
                    apf_from_int(&numer, 1);
                    apf_from_int(&denom, 1);
                    for (k = 0; k < bot; k++) {
                        apf t;
                        apf_from_int(&t, top - k);
                        apf_mul(&numer, &numer, &t);
                        apf_from_int(&t, k + 1);
                        apf_mul(&denom, &denom, &t);
                    }
                    apf_div(&c3, &numer, &denom);
                }
                
                apf_mul(&tmp, &sign, &term);
                apf_mul(&tmp, &tmp, &c1);
                apf_mul(&tmp, &tmp, &c2);
                apf_mul(&tmp, &tmp, &c3);
                apf_mul(&MAT_AT(&result->v.matrix, i, j).re, &tmp, &c3);
            }
        }
    }
    
    return 1;
}

/* Parse value factor (base element) */
static int parse_value_factor(value_t *result)
{
    /* Uses shared statics: pv_tmp_mat, pv_arg */
    
    if (current_token.type == TOK_LPAREN) {
        int leading_not = 0;
        next_token();
        
        /* Handle leading 'not' inside parens */
        while (current_token.type == TOK_NOT) {
            leading_not = !leading_not;
            next_token();
        }
        
        if (!parse_value_expr(result)) return 0;
        
        /* Handle comparisons and boolean operators inside parentheses */
        while (result->type == VAL_SCALAR && 
               (current_token.type == TOK_ASSIGN || current_token.type == TOK_EQUAL ||
                current_token.type == TOK_NE || current_token.type == TOK_LT ||
                current_token.type == TOK_LE || current_token.type == TOK_GT ||
                current_token.type == TOK_GE || current_token.type == TOK_APPROX ||
                current_token.type == TOK_AND || current_token.type == TOK_OR ||
                current_token.type == TOK_XOR)) {
            
            token_type_t op = current_token.type;
            
            if (op == TOK_AND || op == TOK_OR || op == TOK_XOR) {
                /* Boolean operator */
                int lhs_bool = !apf_is_zero(&result->v.scalar.re);
                int rhs_bool, bool_result;
                value_t rhs;
                int rhs_not = 0;
                
                next_token();
                while (current_token.type == TOK_NOT) {
                    rhs_not = !rhs_not;
                    next_token();
                }
                
                if (!parse_value_expr(&rhs)) return 0;
                if (rhs.type != VAL_SCALAR) {
                    printf("Error: boolean operators require scalars\n");
                    return 0;
                }
                
                /* Check for comparison on right side */
                if (current_token.type == TOK_ASSIGN || current_token.type == TOK_EQUAL ||
                    current_token.type == TOK_NE || current_token.type == TOK_LT ||
                    current_token.type == TOK_LE || current_token.type == TOK_GT ||
                    current_token.type == TOK_GE || current_token.type == TOK_APPROX) {
                    value_t rhs2;
                    int cmp_result = 0;
                    token_type_t cmp_op = current_token.type;
                    next_token();
                    if (!parse_value_expr(&rhs2)) return 0;
                    if (rhs2.type == VAL_SCALAR) {
                        int cmp = apf_cmp(&rhs.v.scalar.re, &rhs2.v.scalar.re);
                        int eq = apf_eq(&rhs.v.scalar.re, &rhs2.v.scalar.re) &&
                                 apf_eq(&rhs.v.scalar.im, &rhs2.v.scalar.im);
                        switch (cmp_op) {
                            case TOK_ASSIGN: case TOK_EQUAL: cmp_result = eq; break;
                            case TOK_NE: cmp_result = !eq; break;
                            case TOK_LT: cmp_result = (cmp < 0); break;
                            case TOK_LE: cmp_result = (cmp <= 0); break;
                            case TOK_GT: cmp_result = (cmp > 0); break;
                            case TOK_GE: cmp_result = (cmp >= 0); break;
                            case TOK_APPROX: {
                                apf abs_l, abs_r, up, lo, fh, fl;
                                apf_abs(&abs_l, &rhs.v.scalar.re);
                                apf_abs(&abs_r, &rhs2.v.scalar.re);
                                apf_from_str(&fh, "1.0500001");
                                apf_from_str(&fl, "0.9499999");
                                apf_mul(&up, &abs_l, &fh);
                                apf_mul(&lo, &abs_l, &fl);
                                if (apf_is_zero(&rhs.v.scalar.re)) cmp_result = apf_is_zero(&rhs2.v.scalar.re);
                                else if (rhs.v.scalar.re.sign != rhs2.v.scalar.re.sign) cmp_result = 0;
                                else cmp_result = (apf_cmp(&abs_r, &lo) >= 0 && apf_cmp(&abs_r, &up) <= 0);
                                break;
                            }
                            default: cmp_result = 0; break;
                        }
                    }
                    apf_from_int(&rhs.v.scalar.re, cmp_result);
                    apf_zero(&rhs.v.scalar.im);
                }
                
                rhs_bool = !apf_is_zero(&rhs.v.scalar.re);
                if (rhs_not) rhs_bool = !rhs_bool;
                
                switch (op) {
                    case TOK_AND: bool_result = lhs_bool && rhs_bool; break;
                    case TOK_OR:  bool_result = lhs_bool || rhs_bool; break;
                    case TOK_XOR: bool_result = lhs_bool != rhs_bool; break;
                    default: bool_result = 0; break;
                }
                apf_from_int(&result->v.scalar.re, bool_result);
                apf_zero(&result->v.scalar.im);
            } else {
                /* Comparison operator */
                value_t rhs;
                int cmp_result = 0;
                next_token();
                if (!parse_value_expr(&rhs)) return 0;
                if (rhs.type == VAL_SCALAR) {
                    int cmp = apf_cmp(&result->v.scalar.re, &rhs.v.scalar.re);
                    int eq = apf_eq(&result->v.scalar.re, &rhs.v.scalar.re) &&
                             apf_eq(&result->v.scalar.im, &rhs.v.scalar.im);
                    switch (op) {
                        case TOK_ASSIGN: case TOK_EQUAL: cmp_result = eq; break;
                        case TOK_NE: cmp_result = !eq; break;
                        case TOK_LT: cmp_result = (cmp < 0); break;
                        case TOK_LE: cmp_result = (cmp <= 0); break;
                        case TOK_GT: cmp_result = (cmp > 0); break;
                        case TOK_GE: cmp_result = (cmp >= 0); break;
                        case TOK_APPROX: {
                            apf abs_l, abs_r, up, lo, fh, fl;
                            apf_abs(&abs_l, &result->v.scalar.re);
                            apf_abs(&abs_r, &rhs.v.scalar.re);
                            apf_from_str(&fh, "1.0500001");
                            apf_from_str(&fl, "0.9499999");
                            apf_mul(&up, &abs_l, &fh);
                            apf_mul(&lo, &abs_l, &fl);
                            if (apf_is_zero(&result->v.scalar.re)) cmp_result = apf_is_zero(&rhs.v.scalar.re);
                            else if (result->v.scalar.re.sign != rhs.v.scalar.re.sign) cmp_result = 0;
                            else cmp_result = (apf_cmp(&abs_r, &lo) >= 0 && apf_cmp(&abs_r, &up) <= 0);
                            break;
                        }
                        default: cmp_result = 0; break;
                    }
                }
                apf_from_int(&result->v.scalar.re, cmp_result);
                apf_zero(&result->v.scalar.im);
            }
        }
        
        /* Apply leading 'not' */
        if (leading_not && result->type == VAL_SCALAR) {
            int val = apf_is_zero(&result->v.scalar.re) ? 1 : 0;
            apf_from_int(&result->v.scalar.re, val);
            apf_zero(&result->v.scalar.im);
        }
        
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
        char name[32];
        int var_idx;
        const char *saved_pos = input_ptr;
        token_t saved_tok = current_token;
        
        strncpy(name, current_token.func_name, sizeof(name)-1);
        name[sizeof(name)-1] = '\0';
        var_idx = get_var_index(name);
        
        /* Check for matrix variable with possible indexing */
        if (var_idx >= 0 && is_var_matrix(var_idx)) {
            next_token();
            if (current_token.type == TOK_LPAREN) {
                /* Matrix indexing: M(i,j) or v(i) */
                value_t mat_val;
                apfc idx1, idx2;
                int has_idx2 = 0;
                int idx1_is_colon = 0, idx2_is_colon = 0;
                int i, j;
                
                get_var_value(var_idx, &mat_val);
                next_token();
                
                /* Check for colon (all elements) */
                if (current_token.type == TOK_COLON) {
                    idx1_is_colon = 1;
                    next_token();
                } else {
                    if (!parse_expr(&idx1)) return 0;
                }
                
                if (current_token.type == TOK_COMMA) {
                    next_token();
                    has_idx2 = 1;
                    if (current_token.type == TOK_COLON) {
                        idx2_is_colon = 1;
                        next_token();
                    } else {
                        if (!parse_expr(&idx2)) return 0;
                    }
                }
                
                if (current_token.type != TOK_RPAREN) {
                    printf("Error: expected ')'\n");
                    return 0;
                }
                next_token();
                
                if (mat_val.type != VAL_MATRIX) {
                    printf("Error: indexing requires matrix\n");
                    return 0;
                }
                
                /* Handle different indexing cases */
                if (has_idx2) {
                    /* M(i,j) or M(:,j) or M(i,:) or M(:,:) */
                    if (idx1_is_colon && idx2_is_colon) {
                        /* M(:,:) - entire matrix */
                        *result = mat_val;
                    } else if (idx1_is_colon) {
                        /* M(:,j) - column j */
                        j = apf_to_long(&idx2.re) - 1;  /* 1-indexed */
                        if (j < 0 || j >= mat_val.v.matrix.cols) {
                            printf("Error: column index out of range\n");
                            return 0;
                        }
                        result->type = VAL_MATRIX;
                        result->v.matrix.rows = mat_val.v.matrix.rows;
                        result->v.matrix.cols = 1;
                        for (i = 0; i < mat_val.v.matrix.rows; i++) {
                            MAT_AT(&result->v.matrix, i, 0) = MAT_AT(&mat_val.v.matrix, i, j);
                        }
                    } else if (idx2_is_colon) {
                        /* M(i,:) - row i */
                        i = apf_to_long(&idx1.re) - 1;  /* 1-indexed */
                        if (i < 0 || i >= mat_val.v.matrix.rows) {
                            printf("Error: row index out of range\n");
                            return 0;
                        }
                        result->type = VAL_MATRIX;
                        result->v.matrix.rows = 1;
                        result->v.matrix.cols = mat_val.v.matrix.cols;
                        for (j = 0; j < mat_val.v.matrix.cols; j++) {
                            MAT_AT(&result->v.matrix, 0, j) = MAT_AT(&mat_val.v.matrix, i, j);
                        }
                    } else {
                        /* M(i,j) - single element */
                        i = apf_to_long(&idx1.re) - 1;
                        j = apf_to_long(&idx2.re) - 1;
                        if (i < 0 || i >= mat_val.v.matrix.rows ||
                            j < 0 || j >= mat_val.v.matrix.cols) {
                            printf("Error: index out of range\n");
                            return 0;
                        }
                        result->type = VAL_SCALAR;
                        result->v.scalar = MAT_AT(&mat_val.v.matrix, i, j);
                    }
                } else {
                    /* v(i) - linear indexing */
                    int idx;
                    if (idx1_is_colon) {
                        /* v(:) - linearize */
                        *result = mat_val;
                    } else {
                        idx = apf_to_long(&idx1.re) - 1;  /* 1-indexed */
                        if (mat_val.v.matrix.rows == 1) {
                            /* Row vector */
                            if (idx < 0 || idx >= mat_val.v.matrix.cols) {
                                printf("Error: index out of range\n");
                                return 0;
                            }
                            result->type = VAL_SCALAR;
                            result->v.scalar = MAT_AT(&mat_val.v.matrix, 0, idx);
                        } else if (mat_val.v.matrix.cols == 1) {
                            /* Column vector */
                            if (idx < 0 || idx >= mat_val.v.matrix.rows) {
                                printf("Error: index out of range\n");
                                return 0;
                            }
                            result->type = VAL_SCALAR;
                            result->v.scalar = MAT_AT(&mat_val.v.matrix, idx, 0);
                        } else {
                            /* Matrix linear indexing (column-major) */
                            int total = mat_val.v.matrix.rows * mat_val.v.matrix.cols;
                            int r, c;
                            if (idx < 0 || idx >= total) {
                                printf("Error: index out of range\n");
                                return 0;
                            }
                            c = idx / mat_val.v.matrix.rows;
                            r = idx % mat_val.v.matrix.rows;
                            result->type = VAL_SCALAR;
                            result->v.scalar = MAT_AT(&mat_val.v.matrix, r, c);
                        }
                    }
                }
                return 1;
            }
            /* No indexing - return whole matrix */
            get_var_value(var_idx, result);
            return 1;
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
        
#ifdef HAVE_NAMED_VARS
        /* Check for named variable (multi-character) with possible indexing */
        {
            value_t named_val;
            if (get_named_var(name, &named_val)) {
                next_token();
                if (current_token.type == TOK_LPAREN && named_val.type == VAL_MATRIX) {
                    /* Matrix indexing for named variable */
                    apfc idx1, idx2;
                    int has_idx2 = 0;
                    int idx1_is_colon = 0, idx2_is_colon = 0;
                    int i, j;
                    
                    next_token();
                    
                    if (current_token.type == TOK_COLON) {
                        idx1_is_colon = 1;
                        next_token();
                    } else {
                        if (!parse_expr(&idx1)) return 0;
                    }
                    
                    if (current_token.type == TOK_COMMA) {
                        next_token();
                        has_idx2 = 1;
                        if (current_token.type == TOK_COLON) {
                            idx2_is_colon = 1;
                            next_token();
                        } else {
                            if (!parse_expr(&idx2)) return 0;
                        }
                    }
                    
                    if (current_token.type != TOK_RPAREN) {
                        printf("Error: expected ')'\n");
                        return 0;
                    }
                    next_token();
                    
                    /* Handle different indexing cases */
                    if (has_idx2) {
                        if (idx1_is_colon && idx2_is_colon) {
                            *result = named_val;
                        } else if (idx1_is_colon) {
                            j = apf_to_long(&idx2.re) - 1;
                            if (j < 0 || j >= named_val.v.matrix.cols) {
                                printf("Error: column index out of range\n");
                                return 0;
                            }
                            result->type = VAL_MATRIX;
                            result->v.matrix.rows = named_val.v.matrix.rows;
                            result->v.matrix.cols = 1;
                            for (i = 0; i < named_val.v.matrix.rows; i++) {
                                MAT_AT(&result->v.matrix, i, 0) = MAT_AT(&named_val.v.matrix, i, j);
                            }
                        } else if (idx2_is_colon) {
                            i = apf_to_long(&idx1.re) - 1;
                            if (i < 0 || i >= named_val.v.matrix.rows) {
                                printf("Error: row index out of range\n");
                                return 0;
                            }
                            result->type = VAL_MATRIX;
                            result->v.matrix.rows = 1;
                            result->v.matrix.cols = named_val.v.matrix.cols;
                            for (j = 0; j < named_val.v.matrix.cols; j++) {
                                MAT_AT(&result->v.matrix, 0, j) = MAT_AT(&named_val.v.matrix, i, j);
                            }
                        } else {
                            i = apf_to_long(&idx1.re) - 1;
                            j = apf_to_long(&idx2.re) - 1;
                            if (i < 0 || i >= named_val.v.matrix.rows ||
                                j < 0 || j >= named_val.v.matrix.cols) {
                                printf("Error: index out of range\n");
                                return 0;
                            }
                            result->type = VAL_SCALAR;
                            result->v.scalar = MAT_AT(&named_val.v.matrix, i, j);
                        }
                    } else {
                        int idx;
                        if (idx1_is_colon) {
                            *result = named_val;
                        } else {
                            idx = apf_to_long(&idx1.re) - 1;
                            if (named_val.v.matrix.rows == 1) {
                                if (idx < 0 || idx >= named_val.v.matrix.cols) {
                                    printf("Error: index out of range\n");
                                    return 0;
                                }
                                result->type = VAL_SCALAR;
                                result->v.scalar = MAT_AT(&named_val.v.matrix, 0, idx);
                            } else if (named_val.v.matrix.cols == 1) {
                                if (idx < 0 || idx >= named_val.v.matrix.rows) {
                                    printf("Error: index out of range\n");
                                    return 0;
                                }
                                result->type = VAL_SCALAR;
                                result->v.scalar = MAT_AT(&named_val.v.matrix, idx, 0);
                            } else {
                                int total = named_val.v.matrix.rows * named_val.v.matrix.cols;
                                int r, c;
                                if (idx < 0 || idx >= total) {
                                    printf("Error: index out of range\n");
                                    return 0;
                                }
                                c = idx / named_val.v.matrix.rows;
                                r = idx % named_val.v.matrix.rows;
                                result->type = VAL_SCALAR;
                                result->v.scalar = MAT_AT(&named_val.v.matrix, r, c);
                            }
                        }
                    }
                    return 1;
                }
                if (current_token.type != TOK_LPAREN) {
                    *result = named_val;
                    return 1;
                }
                input_ptr = saved_pos;
                current_token = saved_tok;
            }
        }
#endif
        
        next_token();
        if (current_token.type == TOK_LPAREN) {
            if (str_eq(name, "zeros") || str_eq(name, "ones") || 
                str_eq(name, "eye") || str_eq(name, "identity") ||
                str_eq(name, "rand") || str_eq(name, "randn") || str_eq(name, "randi") ||
                str_eq(name, "linspace") || str_eq(name, "logspace")) {
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
                    apfc saved_s = pv_arg.v.scalar;
                    mat_zero(&pv_arg.v.matrix, 1, 1);
                    MAT_AT(&pv_arg.v.matrix, 0, 0) = saved_s;
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
            
            /* Cholesky decomposition: chol(A) returns L where A = L*L' */
            if (str_eq(name, "chol") || str_eq(name, "cholesky")) {
                next_token();
                if (!parse_value(&pv_arg)) return 0;
                if (current_token.type != TOK_RPAREN) {
                    printf("Error: expected ')'\n");
                    return 0;
                }
                next_token();
                
                if (pv_arg.type != VAL_MATRIX) {
                    apfc saved_s = pv_arg.v.scalar;
                    mat_zero(&pv_arg.v.matrix, 1, 1);
                    MAT_AT(&pv_arg.v.matrix, 0, 0) = saved_s;
                    pv_arg.type = VAL_MATRIX;
                }
                
                result->type = VAL_MATRIX;
                if (!mat_chol(&result->v.matrix, &pv_arg.v.matrix)) {
                    return 0;
                }
                return 1;
            }
            
            /* Null space: null(A) returns basis for null space */
            if (str_eq(name, "null") || str_eq(name, "nullspace")) {
                next_token();
                if (!parse_value(&pv_arg)) return 0;
                if (current_token.type != TOK_RPAREN) {
                    printf("Error: expected ')'\n");
                    return 0;
                }
                next_token();
                
                if (pv_arg.type != VAL_MATRIX) {
                    apfc saved_s = pv_arg.v.scalar;
                    mat_zero(&pv_arg.v.matrix, 1, 1);
                    MAT_AT(&pv_arg.v.matrix, 0, 0) = saved_s;
                    pv_arg.type = VAL_MATRIX;
                }
                
                result->type = VAL_MATRIX;
                mat_null(&result->v.matrix, &pv_arg.v.matrix);
                return 1;
            }
            
            /* SVD: svd(A) returns singular values as vector */
            if (str_eq(name, "svd")) {
                matrix_t U, S, V;
                int i, min_dim;
                
                next_token();
                if (!parse_value(&pv_arg)) return 0;
                if (current_token.type != TOK_RPAREN) {
                    printf("Error: expected ')'\n");
                    return 0;
                }
                next_token();
                
                if (pv_arg.type != VAL_MATRIX) {
                    apfc saved_s = pv_arg.v.scalar;
                    mat_zero(&pv_arg.v.matrix, 1, 1);
                    MAT_AT(&pv_arg.v.matrix, 0, 0) = saved_s;
                    pv_arg.type = VAL_MATRIX;
                }
                
                if (!mat_svd(&U, &S, &V, &pv_arg.v.matrix)) {
                    printf("Error: SVD failed\n");
                    return 0;
                }
                
                /* Extract diagonal of S as column vector */
                min_dim = (S.rows < S.cols) ? S.rows : S.cols;
                mat_zero(&result->v.matrix, min_dim, 1);
                for (i = 0; i < min_dim; i++) {
                    MAT_AT(&result->v.matrix, i, 0) = MAT_AT(&S, i, i);
                }
                result->type = VAL_MATRIX;
                return 1;
            }
            
            /* QR decomposition: qr(A) returns R */
            if (str_eq(name, "qr")) {
                matrix_t Q, R;
                
                next_token();
                if (!parse_value(&pv_arg)) return 0;
                if (current_token.type != TOK_RPAREN) {
                    printf("Error: expected ')'\n");
                    return 0;
                }
                next_token();
                
                if (pv_arg.type != VAL_MATRIX) {
                    apfc saved_s = pv_arg.v.scalar;
                    mat_zero(&pv_arg.v.matrix, 1, 1);
                    MAT_AT(&pv_arg.v.matrix, 0, 0) = saved_s;
                    pv_arg.type = VAL_MATRIX;
                }
                
                if (!mat_qr(&Q, &R, &pv_arg.v.matrix)) {
                    printf("Error: QR decomposition failed\n");
                    return 0;
                }
                
                mat_copy(&result->v.matrix, &R);
                result->type = VAL_MATRIX;
                return 1;
            }
            
            /* QR decomposition: qrq(A) returns Q */
            if (str_eq(name, "qrq")) {
                matrix_t Q, R;
                
                next_token();
                if (!parse_value(&pv_arg)) return 0;
                if (current_token.type != TOK_RPAREN) {
                    printf("Error: expected ')'\n");
                    return 0;
                }
                next_token();
                
                if (pv_arg.type != VAL_MATRIX) {
                    apfc saved_s = pv_arg.v.scalar;
                    mat_zero(&pv_arg.v.matrix, 1, 1);
                    MAT_AT(&pv_arg.v.matrix, 0, 0) = saved_s;
                    pv_arg.type = VAL_MATRIX;
                }
                
                if (!mat_qr(&Q, &R, &pv_arg.v.matrix)) {
                    printf("Error: QR decomposition failed\n");
                    return 0;
                }
                
                mat_copy(&result->v.matrix, &Q);
                result->type = VAL_MATRIX;
                return 1;
            }
            
            /* Schur decomposition: schur(A) returns T */
            if (str_eq(name, "schur")) {
                matrix_t Q, T;
                
                next_token();
                if (!parse_value(&pv_arg)) return 0;
                if (current_token.type != TOK_RPAREN) {
                    printf("Error: expected ')'\n");
                    return 0;
                }
                next_token();
                
                if (pv_arg.type != VAL_MATRIX) {
                    apfc saved_s = pv_arg.v.scalar;
                    mat_zero(&pv_arg.v.matrix, 1, 1);
                    MAT_AT(&pv_arg.v.matrix, 0, 0) = saved_s;
                    pv_arg.type = VAL_MATRIX;
                }
                
                if (!mat_schur(&Q, &T, &pv_arg.v.matrix)) {
                    printf("Error: Schur decomposition failed\n");
                    return 0;
                }
                
                mat_copy(&result->v.matrix, &T);
                result->type = VAL_MATRIX;
                return 1;
            }
            
            /* Schur decomposition: schurq(A) returns Q */
            if (str_eq(name, "schurq")) {
                matrix_t Q, T;
                
                next_token();
                if (!parse_value(&pv_arg)) return 0;
                if (current_token.type != TOK_RPAREN) {
                    printf("Error: expected ')'\n");
                    return 0;
                }
                next_token();
                
                if (pv_arg.type != VAL_MATRIX) {
                    apfc saved_s = pv_arg.v.scalar;
                    mat_zero(&pv_arg.v.matrix, 1, 1);
                    MAT_AT(&pv_arg.v.matrix, 0, 0) = saved_s;
                    pv_arg.type = VAL_MATRIX;
                }
                
                if (!mat_schur(&Q, &T, &pv_arg.v.matrix)) {
                    printf("Error: Schur decomposition failed\n");
                    return 0;
                }
                
                mat_copy(&result->v.matrix, &Q);
                result->type = VAL_MATRIX;
                return 1;
            }
            
            /* Array manipulation functions that return matrices */
            if (str_eq(name, "fliplr") || str_eq(name, "flipud") || str_eq(name, "flip") ||
                str_eq(name, "sort") || str_eq(name, "cumsum") || str_eq(name, "cumprod") ||
                str_eq(name, "cummax") || str_eq(name, "cummin") ||
                str_eq(name, "diff") || str_eq(name, "find") || str_eq(name, "unique") ||
                str_eq(name, "rot90") || str_eq(name, "normalize") ||
                str_eq(name, "triu") || str_eq(name, "tril") || str_eq(name, "squeeze") ||
                str_eq(name, "zscore") || str_eq(name, "rescale")) {
                int r, c;
                
                next_token();
                if (!parse_value(&pv_arg)) return 0;
                if (current_token.type != TOK_RPAREN) {
                    printf("Error: expected ')'\n");
                    return 0;
                }
                next_token();
                
                if (pv_arg.type != VAL_MATRIX) {
                    apfc saved_s = pv_arg.v.scalar;
                    mat_zero(&pv_arg.v.matrix, 1, 1);
                    MAT_AT(&pv_arg.v.matrix, 0, 0) = saved_s;
                    pv_arg.type = VAL_MATRIX;
                }
                
                result->type = VAL_MATRIX;
                
                if (str_eq(name, "fliplr")) {
                    result->v.matrix.rows = pv_arg.v.matrix.rows;
                    result->v.matrix.cols = pv_arg.v.matrix.cols;
                    for (r = 0; r < pv_arg.v.matrix.rows; r++) {
                        for (c = 0; c < pv_arg.v.matrix.cols; c++) {
                            MAT_AT(&result->v.matrix, r, pv_arg.v.matrix.cols - 1 - c) = 
                                MAT_AT(&pv_arg.v.matrix, r, c);
                        }
                    }
                } else if (str_eq(name, "flipud")) {
                    result->v.matrix.rows = pv_arg.v.matrix.rows;
                    result->v.matrix.cols = pv_arg.v.matrix.cols;
                    for (r = 0; r < pv_arg.v.matrix.rows; r++) {
                        for (c = 0; c < pv_arg.v.matrix.cols; c++) {
                            MAT_AT(&result->v.matrix, pv_arg.v.matrix.rows - 1 - r, c) = 
                                MAT_AT(&pv_arg.v.matrix, r, c);
                        }
                    }
                } else if (str_eq(name, "flip")) {
                    /* flip: flip vector (like fliplr for row vec, flipud for col vec) */
                    result->v.matrix.rows = pv_arg.v.matrix.rows;
                    result->v.matrix.cols = pv_arg.v.matrix.cols;
                    if (pv_arg.v.matrix.rows == 1) {
                        /* Row vector - flip horizontally */
                        for (c = 0; c < pv_arg.v.matrix.cols; c++) {
                            MAT_AT(&result->v.matrix, 0, pv_arg.v.matrix.cols - 1 - c) = 
                                MAT_AT(&pv_arg.v.matrix, 0, c);
                        }
                    } else {
                        /* Column vector or matrix - flip vertically */
                        for (r = 0; r < pv_arg.v.matrix.rows; r++) {
                            for (c = 0; c < pv_arg.v.matrix.cols; c++) {
                                MAT_AT(&result->v.matrix, pv_arg.v.matrix.rows - 1 - r, c) = 
                                    MAT_AT(&pv_arg.v.matrix, r, c);
                            }
                        }
                    }
                } else if (str_eq(name, "rot90")) {
                    /* rot90: rotate matrix 90 degrees counterclockwise */
                    result->v.matrix.rows = pv_arg.v.matrix.cols;
                    result->v.matrix.cols = pv_arg.v.matrix.rows;
                    for (r = 0; r < pv_arg.v.matrix.rows; r++) {
                        for (c = 0; c < pv_arg.v.matrix.cols; c++) {
                            MAT_AT(&result->v.matrix, pv_arg.v.matrix.cols - 1 - c, r) = 
                                MAT_AT(&pv_arg.v.matrix, r, c);
                        }
                    }
                } else if (str_eq(name, "normalize")) {
                    /* normalize: normalize vector to unit length */
                    apf sum_sq, norm_val, sq;
                    result->v.matrix.rows = pv_arg.v.matrix.rows;
                    result->v.matrix.cols = pv_arg.v.matrix.cols;
                    /* Compute norm */
                    apf_zero(&sum_sq);
                    for (r = 0; r < pv_arg.v.matrix.rows; r++) {
                        for (c = 0; c < pv_arg.v.matrix.cols; c++) {
                            apf_mul(&sq, &MAT_AT(&pv_arg.v.matrix, r, c).re, &MAT_AT(&pv_arg.v.matrix, r, c).re);
                            apf_add(&sum_sq, &sum_sq, &sq);
                        }
                    }
                    apf_sqrt(&norm_val, &sum_sq);
                    /* Divide each element by norm */
                    for (r = 0; r < pv_arg.v.matrix.rows; r++) {
                        for (c = 0; c < pv_arg.v.matrix.cols; c++) {
                            apf_div(&MAT_AT(&result->v.matrix, r, c).re, 
                                   &MAT_AT(&pv_arg.v.matrix, r, c).re, &norm_val);
                            apf_div(&MAT_AT(&result->v.matrix, r, c).im,
                                   &MAT_AT(&pv_arg.v.matrix, r, c).im, &norm_val);
                        }
                    }
                } else if (str_eq(name, "unique")) {
                    /* unique: return sorted unique elements */
                    int n, i, j, count;
                    n = pv_arg.v.matrix.rows * pv_arg.v.matrix.cols;
                    /* First sort */
                    mat_copy(&pv_tmp_mat, &pv_arg.v.matrix);
                    for (i = 0; i < n - 1; i++) {
                        for (j = 0; j < n - 1 - i; j++) {
                            int r1 = j % pv_tmp_mat.rows;
                            int c1 = j / pv_tmp_mat.rows;
                            int r2 = (j+1) % pv_tmp_mat.rows;
                            int c2 = (j+1) / pv_tmp_mat.rows;
                            if (apf_cmp(&MAT_AT(&pv_tmp_mat, r1, c1).re,
                                       &MAT_AT(&pv_tmp_mat, r2, c2).re) > 0) {
                                apfc tmp = MAT_AT(&pv_tmp_mat, r1, c1);
                                MAT_AT(&pv_tmp_mat, r1, c1) = MAT_AT(&pv_tmp_mat, r2, c2);
                                MAT_AT(&pv_tmp_mat, r2, c2) = tmp;
                            }
                        }
                    }
                    /* Then remove duplicates */
                    result->v.matrix.rows = 1;
                    count = 0;
                    for (i = 0; i < n; i++) {
                        int ri = i % pv_tmp_mat.rows;
                        int ci = i / pv_tmp_mat.rows;
                        int is_dup = 0;
                        for (j = 0; j < count; j++) {
                            if (apf_eq(&MAT_AT(&pv_tmp_mat, ri, ci).re, &MAT_AT(&result->v.matrix, 0, j).re) &&
                                apf_eq(&MAT_AT(&pv_tmp_mat, ri, ci).im, &MAT_AT(&result->v.matrix, 0, j).im)) {
                                is_dup = 1;
                                break;
                            }
                        }
                        if (!is_dup && count < MAT_MAX_COLS) {
                            MAT_AT(&result->v.matrix, 0, count) = MAT_AT(&pv_tmp_mat, ri, ci);
                            count++;
                        }
                    }
                    result->v.matrix.cols = count;
                } else if (str_eq(name, "sort")) {
                    int n, i, j;
                    mat_copy(&result->v.matrix, &pv_arg.v.matrix);
                    n = result->v.matrix.rows * result->v.matrix.cols;
                    for (i = 0; i < n - 1; i++) {
                        for (j = 0; j < n - 1 - i; j++) {
                            int r1 = j % result->v.matrix.rows;
                            int c1 = j / result->v.matrix.rows;
                            int r2 = (j+1) % result->v.matrix.rows;
                            int c2 = (j+1) / result->v.matrix.rows;
                            if (apf_cmp(&MAT_AT(&result->v.matrix, r1, c1).re,
                                       &MAT_AT(&result->v.matrix, r2, c2).re) > 0) {
                                apfc tmp = MAT_AT(&result->v.matrix, r1, c1);
                                MAT_AT(&result->v.matrix, r1, c1) = MAT_AT(&result->v.matrix, r2, c2);
                                MAT_AT(&result->v.matrix, r2, c2) = tmp;
                            }
                        }
                    }
                } else if (str_eq(name, "cumsum")) {
                    int n, i;
                    apfc running;
                    mat_copy(&result->v.matrix, &pv_arg.v.matrix);
                    n = result->v.matrix.rows * result->v.matrix.cols;
                    apf_zero(&running.re);
                    apf_zero(&running.im);
                    for (i = 0; i < n; i++) {
                        int ri = i % result->v.matrix.rows;
                        int ci = i / result->v.matrix.rows;
                        apfc_add(&running, &running, &MAT_AT(&result->v.matrix, ri, ci));
                        MAT_AT(&result->v.matrix, ri, ci) = running;
                    }
                } else if (str_eq(name, "cumprod")) {
                    int n, i;
                    apfc running;
                    mat_copy(&result->v.matrix, &pv_arg.v.matrix);
                    n = result->v.matrix.rows * result->v.matrix.cols;
                    apf_from_int(&running.re, 1);
                    apf_zero(&running.im);
                    for (i = 0; i < n; i++) {
                        int ri = i % result->v.matrix.rows;
                        int ci = i / result->v.matrix.rows;
                        apfc_mul(&running, &running, &MAT_AT(&result->v.matrix, ri, ci));
                        MAT_AT(&result->v.matrix, ri, ci) = running;
                    }
                } else if (str_eq(name, "diff")) {
                    int n, i;
                    n = pv_arg.v.matrix.rows * pv_arg.v.matrix.cols;
                    if (n < 2) {
                        result->v.matrix.rows = 0;
                        result->v.matrix.cols = 0;
                    } else {
                        if (pv_arg.v.matrix.rows == 1) {
                            result->v.matrix.rows = 1;
                            result->v.matrix.cols = n - 1;
                        } else {
                            result->v.matrix.rows = n - 1;
                            result->v.matrix.cols = 1;
                        }
                        for (i = 0; i < n - 1; i++) {
                            int r1 = i % pv_arg.v.matrix.rows;
                            int c1 = i / pv_arg.v.matrix.rows;
                            int r2 = (i+1) % pv_arg.v.matrix.rows;
                            int c2 = (i+1) / pv_arg.v.matrix.rows;
                            int dr = i % result->v.matrix.rows;
                            int dc = i / result->v.matrix.rows;
                            apfc_sub(&MAT_AT(&result->v.matrix, dr, dc),
                                    &MAT_AT(&pv_arg.v.matrix, r2, c2),
                                    &MAT_AT(&pv_arg.v.matrix, r1, c1));
                        }
                    }
                } else if (str_eq(name, "find")) {
                    int n, i, count = 0;
                    n = pv_arg.v.matrix.rows * pv_arg.v.matrix.cols;
                    result->v.matrix.rows = 1;
                    result->v.matrix.cols = 0;
                    for (i = 0; i < n && count < MAT_MAX_COLS; i++) {
                        int ri = i % pv_arg.v.matrix.rows;
                        int ci = i / pv_arg.v.matrix.rows;
                        if (!apf_is_zero(&MAT_AT(&pv_arg.v.matrix, ri, ci).re) ||
                            !apf_is_zero(&MAT_AT(&pv_arg.v.matrix, ri, ci).im)) {
                            apf_from_int(&MAT_AT(&result->v.matrix, 0, count).re, i + 1);
                            apf_zero(&MAT_AT(&result->v.matrix, 0, count).im);
                            count++;
                        }
                    }
                    result->v.matrix.cols = count;
                } else if (str_eq(name, "triu")) {
                    /* triu: upper triangular part */
                    result->v.matrix.rows = pv_arg.v.matrix.rows;
                    result->v.matrix.cols = pv_arg.v.matrix.cols;
                    for (r = 0; r < pv_arg.v.matrix.rows; r++) {
                        for (c = 0; c < pv_arg.v.matrix.cols; c++) {
                            if (c >= r) {
                                MAT_AT(&result->v.matrix, r, c) = MAT_AT(&pv_arg.v.matrix, r, c);
                            } else {
                                apf_zero(&MAT_AT(&result->v.matrix, r, c).re);
                                apf_zero(&MAT_AT(&result->v.matrix, r, c).im);
                            }
                        }
                    }
                } else if (str_eq(name, "tril")) {
                    /* tril: lower triangular part */
                    result->v.matrix.rows = pv_arg.v.matrix.rows;
                    result->v.matrix.cols = pv_arg.v.matrix.cols;
                    for (r = 0; r < pv_arg.v.matrix.rows; r++) {
                        for (c = 0; c < pv_arg.v.matrix.cols; c++) {
                            if (c <= r) {
                                MAT_AT(&result->v.matrix, r, c) = MAT_AT(&pv_arg.v.matrix, r, c);
                            } else {
                                apf_zero(&MAT_AT(&result->v.matrix, r, c).re);
                                apf_zero(&MAT_AT(&result->v.matrix, r, c).im);
                            }
                        }
                    }
                } else if (str_eq(name, "squeeze")) {
                    /* squeeze: remove singleton dimensions (trivial for 2D) */
                    mat_copy(&result->v.matrix, &pv_arg.v.matrix);
                } else if (str_eq(name, "zscore")) {
                    /* zscore: standardize to z-scores (x - mean) / std */
                    apf mean_val, std_val, sum_val, sum_sq, tmp, n_val, var_val;
                    int n = pv_arg.v.matrix.rows * pv_arg.v.matrix.cols;
                    
                    result->v.matrix.rows = pv_arg.v.matrix.rows;
                    result->v.matrix.cols = pv_arg.v.matrix.cols;
                    
                    /* Compute mean */
                    apf_zero(&sum_val);
                    for (r = 0; r < pv_arg.v.matrix.rows; r++) {
                        for (c = 0; c < pv_arg.v.matrix.cols; c++) {
                            apf_add(&sum_val, &sum_val, &MAT_AT(&pv_arg.v.matrix, r, c).re);
                        }
                    }
                    apf_from_int(&n_val, n);
                    apf_div(&mean_val, &sum_val, &n_val);
                    
                    /* Compute std */
                    apf_zero(&sum_sq);
                    for (r = 0; r < pv_arg.v.matrix.rows; r++) {
                        for (c = 0; c < pv_arg.v.matrix.cols; c++) {
                            apf_sub(&tmp, &MAT_AT(&pv_arg.v.matrix, r, c).re, &mean_val);
                            apf_mul(&tmp, &tmp, &tmp);
                            apf_add(&sum_sq, &sum_sq, &tmp);
                        }
                    }
                    if (n > 1) {
                        apf_from_int(&n_val, n - 1);
                        apf_div(&var_val, &sum_sq, &n_val);
                        apf_sqrt(&std_val, &var_val);
                    } else {
                        apf_from_int(&std_val, 1);
                    }
                    
                    /* Apply z-score */
                    for (r = 0; r < pv_arg.v.matrix.rows; r++) {
                        for (c = 0; c < pv_arg.v.matrix.cols; c++) {
                            apf_sub(&MAT_AT(&result->v.matrix, r, c).re, 
                                   &MAT_AT(&pv_arg.v.matrix, r, c).re, &mean_val);
                            apf_div(&MAT_AT(&result->v.matrix, r, c).re,
                                   &MAT_AT(&result->v.matrix, r, c).re, &std_val);
                            apf_zero(&MAT_AT(&result->v.matrix, r, c).im);
                        }
                    }
                } else if (str_eq(name, "cummax")) {
                    /* cummax: cumulative maximum */
                    apf current_max;
                    int first = 1;
                    
                    result->v.matrix.rows = pv_arg.v.matrix.rows;
                    result->v.matrix.cols = pv_arg.v.matrix.cols;
                    
                    for (r = 0; r < pv_arg.v.matrix.rows; r++) {
                        for (c = 0; c < pv_arg.v.matrix.cols; c++) {
                            if (first) {
                                apf_copy(&current_max, &MAT_AT(&pv_arg.v.matrix, r, c).re);
                                first = 0;
                            } else if (apf_cmp(&MAT_AT(&pv_arg.v.matrix, r, c).re, &current_max) > 0) {
                                apf_copy(&current_max, &MAT_AT(&pv_arg.v.matrix, r, c).re);
                            }
                            apf_copy(&MAT_AT(&result->v.matrix, r, c).re, &current_max);
                            apf_zero(&MAT_AT(&result->v.matrix, r, c).im);
                        }
                    }
                } else if (str_eq(name, "cummin")) {
                    /* cummin: cumulative minimum */
                    apf current_min;
                    int first = 1;
                    
                    result->v.matrix.rows = pv_arg.v.matrix.rows;
                    result->v.matrix.cols = pv_arg.v.matrix.cols;
                    
                    for (r = 0; r < pv_arg.v.matrix.rows; r++) {
                        for (c = 0; c < pv_arg.v.matrix.cols; c++) {
                            if (first) {
                                apf_copy(&current_min, &MAT_AT(&pv_arg.v.matrix, r, c).re);
                                first = 0;
                            } else if (apf_cmp(&MAT_AT(&pv_arg.v.matrix, r, c).re, &current_min) < 0) {
                                apf_copy(&current_min, &MAT_AT(&pv_arg.v.matrix, r, c).re);
                            }
                            apf_copy(&MAT_AT(&result->v.matrix, r, c).re, &current_min);
                            apf_zero(&MAT_AT(&result->v.matrix, r, c).im);
                        }
                    }
                } else if (str_eq(name, "rescale")) {
                    /* rescale: scale to [0, 1] */
                    apf min_val, max_val, range_val, tmp;
                    int first = 1;
                    
                    result->v.matrix.rows = pv_arg.v.matrix.rows;
                    result->v.matrix.cols = pv_arg.v.matrix.cols;
                    
                    /* Find min and max */
                    for (r = 0; r < pv_arg.v.matrix.rows; r++) {
                        for (c = 0; c < pv_arg.v.matrix.cols; c++) {
                            if (first) {
                                apf_copy(&min_val, &MAT_AT(&pv_arg.v.matrix, r, c).re);
                                apf_copy(&max_val, &MAT_AT(&pv_arg.v.matrix, r, c).re);
                                first = 0;
                            } else {
                                if (apf_cmp(&MAT_AT(&pv_arg.v.matrix, r, c).re, &min_val) < 0) {
                                    apf_copy(&min_val, &MAT_AT(&pv_arg.v.matrix, r, c).re);
                                }
                                if (apf_cmp(&MAT_AT(&pv_arg.v.matrix, r, c).re, &max_val) > 0) {
                                    apf_copy(&max_val, &MAT_AT(&pv_arg.v.matrix, r, c).re);
                                }
                            }
                        }
                    }
                    
                    /* Compute range */
                    apf_sub(&range_val, &max_val, &min_val);
                    
                    /* Scale to [0, 1] */
                    for (r = 0; r < pv_arg.v.matrix.rows; r++) {
                        for (c = 0; c < pv_arg.v.matrix.cols; c++) {
                            if (apf_is_zero(&range_val)) {
                                apf_from_int(&MAT_AT(&result->v.matrix, r, c).re, 0);
                            } else {
                                apf_sub(&tmp, &MAT_AT(&pv_arg.v.matrix, r, c).re, &min_val);
                                apf_div(&MAT_AT(&result->v.matrix, r, c).re, &tmp, &range_val);
                            }
                            apf_zero(&MAT_AT(&result->v.matrix, r, c).im);
                        }
                    }
                }
                return 1;
            }
            
            /* circshift(v, n) - circular shift */
            if (str_eq(name, "circshift")) {
                int len, i, shift;
                apfc shift_arg;
                
                next_token();
                if (!parse_value(&pv_arg)) return 0;
                if (pv_arg.type != VAL_MATRIX) {
                    apfc saved_s = pv_arg.v.scalar;
                    mat_zero(&pv_arg.v.matrix, 1, 1);
                    MAT_AT(&pv_arg.v.matrix, 0, 0) = saved_s;
                    pv_arg.type = VAL_MATRIX;
                }
                if (current_token.type != TOK_COMMA) {
                    printf("Error: circshift requires (v, n)\n");
                    return 0;
                }
                next_token();
                if (!parse_expr(&shift_arg)) return 0;
                shift = apf_to_long(&shift_arg.re);
                if (current_token.type != TOK_RPAREN) {
                    printf("Error: expected ')'\n");
                    return 0;
                }
                next_token();
                
                result->type = VAL_MATRIX;
                len = pv_arg.v.matrix.rows * pv_arg.v.matrix.cols;
                result->v.matrix.rows = pv_arg.v.matrix.rows;
                result->v.matrix.cols = pv_arg.v.matrix.cols;
                
                if (len > 0) {
                    /* Normalize shift to [0, len) */
                    shift = shift % len;
                    if (shift < 0) shift += len;
                    
                    for (i = 0; i < len; i++) {
                        int src_i = (i - shift + len) % len;
                        int dst_r = i % result->v.matrix.rows;
                        int dst_c = i / result->v.matrix.rows;
                        int src_r = src_i % pv_arg.v.matrix.rows;
                        int src_c = src_i / pv_arg.v.matrix.rows;
                        MAT_AT(&result->v.matrix, dst_r, dst_c) = MAT_AT(&pv_arg.v.matrix, src_r, src_c);
                    }
                }
                return 1;
            }
            
            /* union(A, B) - set union of unique sorted elements */
            if (str_eq(name, "union")) {
                value_t arg2;
                int na, nb, i, j, count;
                
                next_token();
                if (!parse_value(&pv_arg)) return 0;
                if (pv_arg.type != VAL_MATRIX) {
                    apfc saved_s = pv_arg.v.scalar;
                    mat_zero(&pv_arg.v.matrix, 1, 1);
                    MAT_AT(&pv_arg.v.matrix, 0, 0) = saved_s;
                    pv_arg.type = VAL_MATRIX;
                }
                if (current_token.type != TOK_COMMA) {
                    printf("Error: union requires (A, B)\n");
                    return 0;
                }
                next_token();
                if (!parse_value(&arg2)) return 0;
                if (arg2.type != VAL_MATRIX) {
                    apfc saved_s2 = arg2.v.scalar;
                    mat_zero(&arg2.v.matrix, 1, 1);
                    MAT_AT(&arg2.v.matrix, 0, 0) = saved_s2;
                    arg2.type = VAL_MATRIX;
                }
                if (current_token.type != TOK_RPAREN) {
                    printf("Error: expected ')'\n");
                    return 0;
                }
                next_token();
                
                result->type = VAL_MATRIX;
                result->v.matrix.rows = 1;
                na = pv_arg.v.matrix.rows * pv_arg.v.matrix.cols;
                nb = arg2.v.matrix.rows * arg2.v.matrix.cols;
                
                /* Collect all unique elements */
                count = 0;
                for (i = 0; i < na + nb && count < MAT_MAX_COLS; i++) {
                    apfc val;
                    int is_dup = 0;
                    if (i < na) {
                        int ri = i % pv_arg.v.matrix.rows;
                        int ci = i / pv_arg.v.matrix.rows;
                        val = MAT_AT(&pv_arg.v.matrix, ri, ci);
                    } else {
                        int idx = i - na;
                        int ri = idx % arg2.v.matrix.rows;
                        int ci = idx / arg2.v.matrix.rows;
                        val = MAT_AT(&arg2.v.matrix, ri, ci);
                    }
                    for (j = 0; j < count; j++) {
                        if (apf_eq(&val.re, &MAT_AT(&result->v.matrix, 0, j).re)) {
                            is_dup = 1;
                            break;
                        }
                    }
                    if (!is_dup) {
                        MAT_AT(&result->v.matrix, 0, count) = val;
                        count++;
                    }
                }
                result->v.matrix.cols = count;
                /* Sort result */
                for (i = 0; i < count - 1; i++) {
                    for (j = 0; j < count - 1 - i; j++) {
                        if (apf_cmp(&MAT_AT(&result->v.matrix, 0, j).re,
                                   &MAT_AT(&result->v.matrix, 0, j+1).re) > 0) {
                            apfc tmp = MAT_AT(&result->v.matrix, 0, j);
                            MAT_AT(&result->v.matrix, 0, j) = MAT_AT(&result->v.matrix, 0, j+1);
                            MAT_AT(&result->v.matrix, 0, j+1) = tmp;
                        }
                    }
                }
                return 1;
            }
            
            /* intersect(A, B) - set intersection */
            if (str_eq(name, "intersect")) {
                value_t arg2;
                int na, nb, i, j, count;
                
                next_token();
                if (!parse_value(&pv_arg)) return 0;
                if (pv_arg.type != VAL_MATRIX) {
                    apfc saved_s = pv_arg.v.scalar;
                    mat_zero(&pv_arg.v.matrix, 1, 1);
                    MAT_AT(&pv_arg.v.matrix, 0, 0) = saved_s;
                    pv_arg.type = VAL_MATRIX;
                }
                if (current_token.type != TOK_COMMA) {
                    printf("Error: intersect requires (A, B)\n");
                    return 0;
                }
                next_token();
                if (!parse_value(&arg2)) return 0;
                if (arg2.type != VAL_MATRIX) {
                    apfc saved_s2 = arg2.v.scalar;
                    mat_zero(&arg2.v.matrix, 1, 1);
                    MAT_AT(&arg2.v.matrix, 0, 0) = saved_s2;
                    arg2.type = VAL_MATRIX;
                }
                if (current_token.type != TOK_RPAREN) {
                    printf("Error: expected ')'\n");
                    return 0;
                }
                next_token();
                
                result->type = VAL_MATRIX;
                result->v.matrix.rows = 1;
                na = pv_arg.v.matrix.rows * pv_arg.v.matrix.cols;
                nb = arg2.v.matrix.rows * arg2.v.matrix.cols;
                
                count = 0;
                for (i = 0; i < na && count < MAT_MAX_COLS; i++) {
                    int ri = i % pv_arg.v.matrix.rows;
                    int ci = i / pv_arg.v.matrix.rows;
                    apfc val = MAT_AT(&pv_arg.v.matrix, ri, ci);
                    int in_b = 0, already = 0;
                    
                    /* Check if in B */
                    for (j = 0; j < nb; j++) {
                        int rj = j % arg2.v.matrix.rows;
                        int cj = j / arg2.v.matrix.rows;
                        if (apf_eq(&val.re, &MAT_AT(&arg2.v.matrix, rj, cj).re)) {
                            in_b = 1;
                            break;
                        }
                    }
                    /* Check if already in result */
                    for (j = 0; j < count; j++) {
                        if (apf_eq(&val.re, &MAT_AT(&result->v.matrix, 0, j).re)) {
                            already = 1;
                            break;
                        }
                    }
                    if (in_b && !already) {
                        MAT_AT(&result->v.matrix, 0, count) = val;
                        count++;
                    }
                }
                result->v.matrix.cols = count;
                /* Sort result */
                for (i = 0; i < count - 1; i++) {
                    for (j = 0; j < count - 1 - i; j++) {
                        if (apf_cmp(&MAT_AT(&result->v.matrix, 0, j).re,
                                   &MAT_AT(&result->v.matrix, 0, j+1).re) > 0) {
                            apfc tmp = MAT_AT(&result->v.matrix, 0, j);
                            MAT_AT(&result->v.matrix, 0, j) = MAT_AT(&result->v.matrix, 0, j+1);
                            MAT_AT(&result->v.matrix, 0, j+1) = tmp;
                        }
                    }
                }
                return 1;
            }
            
            /* setdiff(A, B) - set difference (elements in A but not in B) */
            if (str_eq(name, "setdiff")) {
                value_t arg2;
                int na, nb, i, j, count;
                
                next_token();
                if (!parse_value(&pv_arg)) return 0;
                if (pv_arg.type != VAL_MATRIX) {
                    apfc saved_s = pv_arg.v.scalar;
                    mat_zero(&pv_arg.v.matrix, 1, 1);
                    MAT_AT(&pv_arg.v.matrix, 0, 0) = saved_s;
                    pv_arg.type = VAL_MATRIX;
                }
                if (current_token.type != TOK_COMMA) {
                    printf("Error: setdiff requires (A, B)\n");
                    return 0;
                }
                next_token();
                if (!parse_value(&arg2)) return 0;
                if (arg2.type != VAL_MATRIX) {
                    apfc saved_s2 = arg2.v.scalar;
                    mat_zero(&arg2.v.matrix, 1, 1);
                    MAT_AT(&arg2.v.matrix, 0, 0) = saved_s2;
                    arg2.type = VAL_MATRIX;
                }
                if (current_token.type != TOK_RPAREN) {
                    printf("Error: expected ')'\n");
                    return 0;
                }
                next_token();
                
                result->type = VAL_MATRIX;
                result->v.matrix.rows = 1;
                na = pv_arg.v.matrix.rows * pv_arg.v.matrix.cols;
                nb = arg2.v.matrix.rows * arg2.v.matrix.cols;
                
                count = 0;
                for (i = 0; i < na && count < MAT_MAX_COLS; i++) {
                    int ri = i % pv_arg.v.matrix.rows;
                    int ci = i / pv_arg.v.matrix.rows;
                    apfc val = MAT_AT(&pv_arg.v.matrix, ri, ci);
                    int in_b = 0, already = 0;
                    
                    /* Check if in B */
                    for (j = 0; j < nb; j++) {
                        int rj = j % arg2.v.matrix.rows;
                        int cj = j / arg2.v.matrix.rows;
                        if (apf_eq(&val.re, &MAT_AT(&arg2.v.matrix, rj, cj).re)) {
                            in_b = 1;
                            break;
                        }
                    }
                    /* Check if already in result */
                    for (j = 0; j < count; j++) {
                        if (apf_eq(&val.re, &MAT_AT(&result->v.matrix, 0, j).re)) {
                            already = 1;
                            break;
                        }
                    }
                    if (!in_b && !already) {
                        MAT_AT(&result->v.matrix, 0, count) = val;
                        count++;
                    }
                }
                result->v.matrix.cols = count;
                /* Sort result */
                for (i = 0; i < count - 1; i++) {
                    for (j = 0; j < count - 1 - i; j++) {
                        if (apf_cmp(&MAT_AT(&result->v.matrix, 0, j).re,
                                   &MAT_AT(&result->v.matrix, 0, j+1).re) > 0) {
                            apfc tmp = MAT_AT(&result->v.matrix, 0, j);
                            MAT_AT(&result->v.matrix, 0, j) = MAT_AT(&result->v.matrix, 0, j+1);
                            MAT_AT(&result->v.matrix, 0, j+1) = tmp;
                        }
                    }
                }
                return 1;
            }
            
            /* setxor(A, B) - set symmetric difference */
            if (str_eq(name, "setxor")) {
                value_t arg2;
                int na, nb, i, j, count;
                
                next_token();
                if (!parse_value(&pv_arg)) return 0;
                if (pv_arg.type != VAL_MATRIX) {
                    apfc saved_s = pv_arg.v.scalar;
                    mat_zero(&pv_arg.v.matrix, 1, 1);
                    MAT_AT(&pv_arg.v.matrix, 0, 0) = saved_s;
                    pv_arg.type = VAL_MATRIX;
                }
                if (current_token.type != TOK_COMMA) {
                    printf("Error: setxor requires (A, B)\n");
                    return 0;
                }
                next_token();
                if (!parse_value(&arg2)) return 0;
                if (arg2.type != VAL_MATRIX) {
                    apfc saved_s2 = arg2.v.scalar;
                    mat_zero(&arg2.v.matrix, 1, 1);
                    MAT_AT(&arg2.v.matrix, 0, 0) = saved_s2;
                    arg2.type = VAL_MATRIX;
                }
                if (current_token.type != TOK_RPAREN) {
                    printf("Error: expected ')'\n");
                    return 0;
                }
                next_token();
                
                result->type = VAL_MATRIX;
                result->v.matrix.rows = 1;
                na = pv_arg.v.matrix.rows * pv_arg.v.matrix.cols;
                nb = arg2.v.matrix.rows * arg2.v.matrix.cols;
                
                count = 0;
                /* Add elements from A not in B */
                for (i = 0; i < na && count < MAT_MAX_COLS; i++) {
                    int ri = i % pv_arg.v.matrix.rows;
                    int ci = i / pv_arg.v.matrix.rows;
                    apfc val = MAT_AT(&pv_arg.v.matrix, ri, ci);
                    int in_b = 0, already = 0;
                    
                    for (j = 0; j < nb; j++) {
                        int rj = j % arg2.v.matrix.rows;
                        int cj = j / arg2.v.matrix.rows;
                        if (apf_eq(&val.re, &MAT_AT(&arg2.v.matrix, rj, cj).re)) {
                            in_b = 1;
                            break;
                        }
                    }
                    for (j = 0; j < count; j++) {
                        if (apf_eq(&val.re, &MAT_AT(&result->v.matrix, 0, j).re)) {
                            already = 1;
                            break;
                        }
                    }
                    if (!in_b && !already) {
                        MAT_AT(&result->v.matrix, 0, count) = val;
                        count++;
                    }
                }
                /* Add elements from B not in A */
                for (i = 0; i < nb && count < MAT_MAX_COLS; i++) {
                    int ri = i % arg2.v.matrix.rows;
                    int ci = i / arg2.v.matrix.rows;
                    apfc val = MAT_AT(&arg2.v.matrix, ri, ci);
                    int in_a = 0, already = 0;
                    
                    for (j = 0; j < na; j++) {
                        int rj = j % pv_arg.v.matrix.rows;
                        int cj = j / pv_arg.v.matrix.rows;
                        if (apf_eq(&val.re, &MAT_AT(&pv_arg.v.matrix, rj, cj).re)) {
                            in_a = 1;
                            break;
                        }
                    }
                    for (j = 0; j < count; j++) {
                        if (apf_eq(&val.re, &MAT_AT(&result->v.matrix, 0, j).re)) {
                            already = 1;
                            break;
                        }
                    }
                    if (!in_a && !already) {
                        MAT_AT(&result->v.matrix, 0, count) = val;
                        count++;
                    }
                }
                result->v.matrix.cols = count;
                /* Sort result */
                for (i = 0; i < count - 1; i++) {
                    for (j = 0; j < count - 1 - i; j++) {
                        if (apf_cmp(&MAT_AT(&result->v.matrix, 0, j).re,
                                   &MAT_AT(&result->v.matrix, 0, j+1).re) > 0) {
                            apfc tmp = MAT_AT(&result->v.matrix, 0, j);
                            MAT_AT(&result->v.matrix, 0, j) = MAT_AT(&result->v.matrix, 0, j+1);
                            MAT_AT(&result->v.matrix, 0, j+1) = tmp;
                        }
                    }
                }
                return 1;
            }
            
            /* maxk(v, k) - k largest elements */
            if (str_eq(name, "maxk")) {
                int k, n, i, j;
                apfc k_arg;
                
                next_token();
                if (!parse_value(&pv_arg)) return 0;
                if (pv_arg.type != VAL_MATRIX) {
                    apfc saved_s = pv_arg.v.scalar;
                    mat_zero(&pv_arg.v.matrix, 1, 1);
                    MAT_AT(&pv_arg.v.matrix, 0, 0) = saved_s;
                    pv_arg.type = VAL_MATRIX;
                }
                if (current_token.type != TOK_COMMA) {
                    printf("Error: maxk requires (v, k)\n");
                    return 0;
                }
                next_token();
                if (!parse_expr(&k_arg)) return 0;
                k = apf_to_long(&k_arg.re);
                if (current_token.type != TOK_RPAREN) {
                    printf("Error: expected ')'\n");
                    return 0;
                }
                next_token();
                
                n = pv_arg.v.matrix.rows * pv_arg.v.matrix.cols;
                if (k > n) k = n;
                if (k < 1) k = 1;
                
                /* Sort a copy descending */
                mat_copy(&pv_tmp_mat, &pv_arg.v.matrix);
                for (i = 0; i < n - 1; i++) {
                    for (j = 0; j < n - 1 - i; j++) {
                        int r1 = j % pv_tmp_mat.rows;
                        int c1 = j / pv_tmp_mat.rows;
                        int r2 = (j+1) % pv_tmp_mat.rows;
                        int c2 = (j+1) / pv_tmp_mat.rows;
                        if (apf_cmp(&MAT_AT(&pv_tmp_mat, r1, c1).re,
                                   &MAT_AT(&pv_tmp_mat, r2, c2).re) < 0) {
                            apfc tmp = MAT_AT(&pv_tmp_mat, r1, c1);
                            MAT_AT(&pv_tmp_mat, r1, c1) = MAT_AT(&pv_tmp_mat, r2, c2);
                            MAT_AT(&pv_tmp_mat, r2, c2) = tmp;
                        }
                    }
                }
                
                result->type = VAL_MATRIX;
                result->v.matrix.rows = 1;
                result->v.matrix.cols = k;
                for (i = 0; i < k; i++) {
                    int ri = i % pv_tmp_mat.rows;
                    int ci = i / pv_tmp_mat.rows;
                    MAT_AT(&result->v.matrix, 0, i) = MAT_AT(&pv_tmp_mat, ri, ci);
                }
                return 1;
            }
            
            /* mink(v, k) - k smallest elements */
            if (str_eq(name, "mink")) {
                int k, n, i, j;
                apfc k_arg;
                
                next_token();
                if (!parse_value(&pv_arg)) return 0;
                if (pv_arg.type != VAL_MATRIX) {
                    apfc saved_s = pv_arg.v.scalar;
                    mat_zero(&pv_arg.v.matrix, 1, 1);
                    MAT_AT(&pv_arg.v.matrix, 0, 0) = saved_s;
                    pv_arg.type = VAL_MATRIX;
                }
                if (current_token.type != TOK_COMMA) {
                    printf("Error: mink requires (v, k)\n");
                    return 0;
                }
                next_token();
                if (!parse_expr(&k_arg)) return 0;
                k = apf_to_long(&k_arg.re);
                if (current_token.type != TOK_RPAREN) {
                    printf("Error: expected ')'\n");
                    return 0;
                }
                next_token();
                
                n = pv_arg.v.matrix.rows * pv_arg.v.matrix.cols;
                if (k > n) k = n;
                if (k < 1) k = 1;
                
                /* Sort a copy ascending */
                mat_copy(&pv_tmp_mat, &pv_arg.v.matrix);
                for (i = 0; i < n - 1; i++) {
                    for (j = 0; j < n - 1 - i; j++) {
                        int r1 = j % pv_tmp_mat.rows;
                        int c1 = j / pv_tmp_mat.rows;
                        int r2 = (j+1) % pv_tmp_mat.rows;
                        int c2 = (j+1) / pv_tmp_mat.rows;
                        if (apf_cmp(&MAT_AT(&pv_tmp_mat, r1, c1).re,
                                   &MAT_AT(&pv_tmp_mat, r2, c2).re) > 0) {
                            apfc tmp = MAT_AT(&pv_tmp_mat, r1, c1);
                            MAT_AT(&pv_tmp_mat, r1, c1) = MAT_AT(&pv_tmp_mat, r2, c2);
                            MAT_AT(&pv_tmp_mat, r2, c2) = tmp;
                        }
                    }
                }
                
                result->type = VAL_MATRIX;
                result->v.matrix.rows = 1;
                result->v.matrix.cols = k;
                for (i = 0; i < k; i++) {
                    int ri = i % pv_tmp_mat.rows;
                    int ci = i / pv_tmp_mat.rows;
                    MAT_AT(&result->v.matrix, 0, i) = MAT_AT(&pv_tmp_mat, ri, ci);
                }
                return 1;
            }
            
            /* prctile(v, p) - percentile (0-100) */
            if (str_eq(name, "prctile")) {
                int n, i, j, idx;
                apfc p_arg;
                apf p_val, pos_val, n_minus_1, hundred, frac, one_minus_frac, idx_apf;
                
                next_token();
                if (!parse_value(&pv_arg)) return 0;
                if (pv_arg.type != VAL_MATRIX) {
                    apfc saved_s = pv_arg.v.scalar;
                    mat_zero(&pv_arg.v.matrix, 1, 1);
                    MAT_AT(&pv_arg.v.matrix, 0, 0) = saved_s;
                    pv_arg.type = VAL_MATRIX;
                }
                if (current_token.type != TOK_COMMA) {
                    printf("Error: prctile requires (v, p)\n");
                    return 0;
                }
                next_token();
                if (!parse_expr(&p_arg)) return 0;
                apf_copy(&p_val, &p_arg.re);
                if (current_token.type != TOK_RPAREN) {
                    printf("Error: expected ')'\n");
                    return 0;
                }
                next_token();
                
                /* Clamp p to [0, 100] */
                {
                    apf zero_apf, hundred_apf;
                    apf_zero(&zero_apf);
                    apf_from_int(&hundred_apf, 100);
                    if (apf_cmp(&p_val, &zero_apf) < 0) apf_zero(&p_val);
                    if (apf_cmp(&p_val, &hundred_apf) > 0) apf_from_int(&p_val, 100);
                }
                
                n = pv_arg.v.matrix.rows * pv_arg.v.matrix.cols;
                
                /* Sort a copy */
                mat_copy(&pv_tmp_mat, &pv_arg.v.matrix);
                for (i = 0; i < n - 1; i++) {
                    for (j = 0; j < n - 1 - i; j++) {
                        int r1 = j % pv_tmp_mat.rows;
                        int c1 = j / pv_tmp_mat.rows;
                        int r2 = (j+1) % pv_tmp_mat.rows;
                        int c2 = (j+1) / pv_tmp_mat.rows;
                        if (apf_cmp(&MAT_AT(&pv_tmp_mat, r1, c1).re,
                                   &MAT_AT(&pv_tmp_mat, r2, c2).re) > 0) {
                            apfc tmp = MAT_AT(&pv_tmp_mat, r1, c1);
                            MAT_AT(&pv_tmp_mat, r1, c1) = MAT_AT(&pv_tmp_mat, r2, c2);
                            MAT_AT(&pv_tmp_mat, r2, c2) = tmp;
                        }
                    }
                }
                
                /* Linear interpolation: pos = (p / 100) * (n - 1) */
                apf_from_int(&hundred, 100);
                apf_from_int(&n_minus_1, n - 1);
                apf_div(&pos_val, &p_val, &hundred);
                apf_mul(&pos_val, &pos_val, &n_minus_1);
                
                idx = apf_to_long(&pos_val);
                if (idx < 0) idx = 0;
                if (idx >= n - 1) {
                    int ri = (n-1) % pv_tmp_mat.rows;
                    int ci = (n-1) / pv_tmp_mat.rows;
                    result->type = VAL_SCALAR;
                    result->v.scalar = MAT_AT(&pv_tmp_mat, ri, ci);
                } else {
                    int r1 = idx % pv_tmp_mat.rows;
                    int c1 = idx / pv_tmp_mat.rows;
                    int r2 = (idx+1) % pv_tmp_mat.rows;
                    int c2 = (idx+1) / pv_tmp_mat.rows;
                    
                    /* frac = pos - idx */
                    apf_from_int(&idx_apf, idx);
                    apf_sub(&frac, &pos_val, &idx_apf);
                    
                    apf_from_int(&one_minus_frac, 1);
                    apf_sub(&one_minus_frac, &one_minus_frac, &frac);
                    result->type = VAL_SCALAR;
                    apf_mul(&result->v.scalar.re, &MAT_AT(&pv_tmp_mat, r1, c1).re, &one_minus_frac);
                    {
                        apf tmp;
                        apf_mul(&tmp, &MAT_AT(&pv_tmp_mat, r2, c2).re, &frac);
                        apf_add(&result->v.scalar.re, &result->v.scalar.re, &tmp);
                    }
                    apf_zero(&result->v.scalar.im);
                }
                return 1;
            }
            
            /* iqr(v) - interquartile range (Q3 - Q1) */
            if (str_eq(name, "iqr")) {
                int n, i, j;
                apf q1, q3;
                apf pos_val, n_minus_1, quarter, three_quarters, idx_apf;
                int idx;
                
                next_token();
                if (!parse_value(&pv_arg)) return 0;
                if (current_token.type != TOK_RPAREN) {
                    printf("Error: expected ')'\n");
                    return 0;
                }
                next_token();
                
                if (pv_arg.type != VAL_MATRIX) {
                    apfc saved_s = pv_arg.v.scalar;
                    mat_zero(&pv_arg.v.matrix, 1, 1);
                    MAT_AT(&pv_arg.v.matrix, 0, 0) = saved_s;
                    pv_arg.type = VAL_MATRIX;
                }
                
                n = pv_arg.v.matrix.rows * pv_arg.v.matrix.cols;
                
                /* Sort a copy */
                mat_copy(&pv_tmp_mat, &pv_arg.v.matrix);
                for (i = 0; i < n - 1; i++) {
                    for (j = 0; j < n - 1 - i; j++) {
                        int r1 = j % pv_tmp_mat.rows;
                        int c1 = j / pv_tmp_mat.rows;
                        int r2 = (j+1) % pv_tmp_mat.rows;
                        int c2 = (j+1) / pv_tmp_mat.rows;
                        if (apf_cmp(&MAT_AT(&pv_tmp_mat, r1, c1).re,
                                   &MAT_AT(&pv_tmp_mat, r2, c2).re) > 0) {
                            apfc tmp = MAT_AT(&pv_tmp_mat, r1, c1);
                            MAT_AT(&pv_tmp_mat, r1, c1) = MAT_AT(&pv_tmp_mat, r2, c2);
                            MAT_AT(&pv_tmp_mat, r2, c2) = tmp;
                        }
                    }
                }
                
                apf_from_int(&n_minus_1, n - 1);
                
                /* Q1 (25th percentile): pos = 0.25 * (n - 1) */
                apf_from_str(&quarter, "0.25");
                apf_mul(&pos_val, &quarter, &n_minus_1);
                idx = apf_to_long(&pos_val);
                if (idx < 0) idx = 0;
                if (idx >= n - 1) {
                    int ri = (n-1) % pv_tmp_mat.rows;
                    int ci = (n-1) / pv_tmp_mat.rows;
                    q1 = MAT_AT(&pv_tmp_mat, ri, ci).re;
                } else {
                    int r1 = idx % pv_tmp_mat.rows;
                    int c1 = idx / pv_tmp_mat.rows;
                    int r2 = (idx+1) % pv_tmp_mat.rows;
                    int c2 = (idx+1) / pv_tmp_mat.rows;
                    apf frac, one_minus, tmp;
                    apf_from_int(&idx_apf, idx);
                    apf_sub(&frac, &pos_val, &idx_apf);
                    apf_from_int(&one_minus, 1);
                    apf_sub(&one_minus, &one_minus, &frac);
                    apf_mul(&q1, &MAT_AT(&pv_tmp_mat, r1, c1).re, &one_minus);
                    apf_mul(&tmp, &MAT_AT(&pv_tmp_mat, r2, c2).re, &frac);
                    apf_add(&q1, &q1, &tmp);
                }
                
                /* Q3 (75th percentile): pos = 0.75 * (n - 1) */
                apf_from_str(&three_quarters, "0.75");
                apf_mul(&pos_val, &three_quarters, &n_minus_1);
                idx = apf_to_long(&pos_val);
                if (idx < 0) idx = 0;
                if (idx >= n - 1) {
                    int ri = (n-1) % pv_tmp_mat.rows;
                    int ci = (n-1) / pv_tmp_mat.rows;
                    q3 = MAT_AT(&pv_tmp_mat, ri, ci).re;
                } else {
                    int r1 = idx % pv_tmp_mat.rows;
                    int c1 = idx / pv_tmp_mat.rows;
                    int r2 = (idx+1) % pv_tmp_mat.rows;
                    int c2 = (idx+1) / pv_tmp_mat.rows;
                    apf frac, one_minus, tmp;
                    apf_from_int(&idx_apf, idx);
                    apf_sub(&frac, &pos_val, &idx_apf);
                    apf_from_int(&one_minus, 1);
                    apf_sub(&one_minus, &one_minus, &frac);
                    apf_mul(&q3, &MAT_AT(&pv_tmp_mat, r1, c1).re, &one_minus);
                    apf_mul(&tmp, &MAT_AT(&pv_tmp_mat, r2, c2).re, &frac);
                    apf_add(&q3, &q3, &tmp);
                }
                
                result->type = VAL_SCALAR;
                apf_sub(&result->v.scalar.re, &q3, &q1);
                apf_zero(&result->v.scalar.im);
                return 1;
            }
            
            /* mode(v) - most frequent value */
            if (str_eq(name, "mode")) {
                int n, i, j, max_count, curr_count;
                apfc mode_val;
                
                next_token();
                if (!parse_value(&pv_arg)) return 0;
                if (current_token.type != TOK_RPAREN) {
                    printf("Error: expected ')'\n");
                    return 0;
                }
                next_token();
                
                if (pv_arg.type != VAL_MATRIX) {
                    result->type = VAL_SCALAR;
                    result->v.scalar = pv_arg.v.scalar;
                    return 1;
                }
                
                n = pv_arg.v.matrix.rows * pv_arg.v.matrix.cols;
                if (n == 0) {
                    result->type = VAL_SCALAR;
                    apf_zero(&result->v.scalar.re);
                    apf_zero(&result->v.scalar.im);
                    return 1;
                }
                
                max_count = 0;
                apf_zero(&mode_val.re);
                apf_zero(&mode_val.im);
                
                for (i = 0; i < n; i++) {
                    int ri = i % pv_arg.v.matrix.rows;
                    int ci = i / pv_arg.v.matrix.rows;
                    apfc val = MAT_AT(&pv_arg.v.matrix, ri, ci);
                    curr_count = 0;
                    
                    for (j = 0; j < n; j++) {
                        int rj = j % pv_arg.v.matrix.rows;
                        int cj = j / pv_arg.v.matrix.rows;
                        if (apf_eq(&val.re, &MAT_AT(&pv_arg.v.matrix, rj, cj).re)) {
                            curr_count++;
                        }
                    }
                    
                    if (curr_count > max_count || 
                        (curr_count == max_count && apf_cmp(&val.re, &mode_val.re) < 0)) {
                        max_count = curr_count;
                        mode_val = val;
                    }
                }
                
                result->type = VAL_SCALAR;
                result->v.scalar = mode_val;
                return 1;
            }
            
            /* sortrows(M) - sort matrix rows by first column */
            if (str_eq(name, "sortrows")) {
                int i, j, c;
                
                next_token();
                if (!parse_value(&pv_arg)) return 0;
                if (current_token.type != TOK_RPAREN) {
                    printf("Error: expected ')'\n");
                    return 0;
                }
                next_token();
                
                if (pv_arg.type != VAL_MATRIX) {
                    apfc saved_s = pv_arg.v.scalar;
                    mat_zero(&pv_arg.v.matrix, 1, 1);
                    MAT_AT(&pv_arg.v.matrix, 0, 0) = saved_s;
                    pv_arg.type = VAL_MATRIX;
                }
                
                result->type = VAL_MATRIX;
                mat_copy(&result->v.matrix, &pv_arg.v.matrix);
                
                /* Bubble sort rows by first column */
                for (i = 0; i < result->v.matrix.rows - 1; i++) {
                    for (j = 0; j < result->v.matrix.rows - 1 - i; j++) {
                        if (apf_cmp(&MAT_AT(&result->v.matrix, j, 0).re,
                                   &MAT_AT(&result->v.matrix, j+1, 0).re) > 0) {
                            /* Swap rows j and j+1 */
                            for (c = 0; c < result->v.matrix.cols; c++) {
                                apfc tmp = MAT_AT(&result->v.matrix, j, c);
                                MAT_AT(&result->v.matrix, j, c) = MAT_AT(&result->v.matrix, j+1, c);
                                MAT_AT(&result->v.matrix, j+1, c) = tmp;
                            }
                        }
                    }
                }
                return 1;
            }
            
            /* horzcat(A, B) - horizontal concatenation */
            if (str_eq(name, "horzcat")) {
                value_t arg2;
                int r, c;
                
                next_token();
                if (!parse_value(&pv_arg)) return 0;
                if (pv_arg.type != VAL_MATRIX) {
                    apfc saved_s = pv_arg.v.scalar;
                    mat_zero(&pv_arg.v.matrix, 1, 1);
                    MAT_AT(&pv_arg.v.matrix, 0, 0) = saved_s;
                    pv_arg.type = VAL_MATRIX;
                }
                if (current_token.type != TOK_COMMA) {
                    printf("Error: horzcat requires (A, B)\n");
                    return 0;
                }
                next_token();
                if (!parse_value(&arg2)) return 0;
                if (arg2.type != VAL_MATRIX) {
                    apfc saved_s2 = arg2.v.scalar;
                    mat_zero(&arg2.v.matrix, 1, 1);
                    MAT_AT(&arg2.v.matrix, 0, 0) = saved_s2;
                    arg2.type = VAL_MATRIX;
                }
                if (current_token.type != TOK_RPAREN) {
                    printf("Error: expected ')'\n");
                    return 0;
                }
                next_token();
                
                if (pv_arg.v.matrix.rows != arg2.v.matrix.rows) {
                    printf("Error: horzcat requires same number of rows\n");
                    return 0;
                }
                
                result->type = VAL_MATRIX;
                result->v.matrix.rows = pv_arg.v.matrix.rows;
                result->v.matrix.cols = pv_arg.v.matrix.cols + arg2.v.matrix.cols;
                
                for (r = 0; r < pv_arg.v.matrix.rows; r++) {
                    for (c = 0; c < pv_arg.v.matrix.cols; c++) {
                        MAT_AT(&result->v.matrix, r, c) = MAT_AT(&pv_arg.v.matrix, r, c);
                    }
                    for (c = 0; c < arg2.v.matrix.cols; c++) {
                        MAT_AT(&result->v.matrix, r, pv_arg.v.matrix.cols + c) = MAT_AT(&arg2.v.matrix, r, c);
                    }
                }
                return 1;
            }
            
            /* vertcat(A, B) - vertical concatenation */
            if (str_eq(name, "vertcat")) {
                value_t arg2;
                int r, c;
                
                next_token();
                if (!parse_value(&pv_arg)) return 0;
                if (pv_arg.type != VAL_MATRIX) {
                    apfc saved_s = pv_arg.v.scalar;
                    mat_zero(&pv_arg.v.matrix, 1, 1);
                    MAT_AT(&pv_arg.v.matrix, 0, 0) = saved_s;
                    pv_arg.type = VAL_MATRIX;
                }
                if (current_token.type != TOK_COMMA) {
                    printf("Error: vertcat requires (A, B)\n");
                    return 0;
                }
                next_token();
                if (!parse_value(&arg2)) return 0;
                if (arg2.type != VAL_MATRIX) {
                    apfc saved_s2 = arg2.v.scalar;
                    mat_zero(&arg2.v.matrix, 1, 1);
                    MAT_AT(&arg2.v.matrix, 0, 0) = saved_s2;
                    arg2.type = VAL_MATRIX;
                }
                if (current_token.type != TOK_RPAREN) {
                    printf("Error: expected ')'\n");
                    return 0;
                }
                next_token();
                
                if (pv_arg.v.matrix.cols != arg2.v.matrix.cols) {
                    printf("Error: vertcat requires same number of columns\n");
                    return 0;
                }
                
                result->type = VAL_MATRIX;
                result->v.matrix.rows = pv_arg.v.matrix.rows + arg2.v.matrix.rows;
                result->v.matrix.cols = pv_arg.v.matrix.cols;
                
                for (r = 0; r < pv_arg.v.matrix.rows; r++) {
                    for (c = 0; c < pv_arg.v.matrix.cols; c++) {
                        MAT_AT(&result->v.matrix, r, c) = MAT_AT(&pv_arg.v.matrix, r, c);
                    }
                }
                for (r = 0; r < arg2.v.matrix.rows; r++) {
                    for (c = 0; c < arg2.v.matrix.cols; c++) {
                        MAT_AT(&result->v.matrix, pv_arg.v.matrix.rows + r, c) = MAT_AT(&arg2.v.matrix, r, c);
                    }
                }
                return 1;
            }
            
            /* Element-wise comparison functions: eq, ne, lt, le, gt, ge */
            if (str_eq(name, "eq") || str_eq(name, "ne") || str_eq(name, "lt") ||
                str_eq(name, "le") || str_eq(name, "gt") || str_eq(name, "ge")) {
                value_t arg2;
                int r, c;
                
                next_token();
                if (!parse_value(&pv_arg)) return 0;
                if (current_token.type != TOK_COMMA) {
                    printf("Error: %s requires two arguments\n", name);
                    return 0;
                }
                next_token();
                if (!parse_value(&arg2)) return 0;
                if (current_token.type != TOK_RPAREN) {
                    printf("Error: expected ')'\n");
                    return 0;
                }
                next_token();
                
                /* Convert scalars to 1x1 matrices */
                if (pv_arg.type == VAL_SCALAR) {
                    apfc saved = pv_arg.v.scalar;
                    mat_zero(&pv_arg.v.matrix, 1, 1);
                    MAT_AT(&pv_arg.v.matrix, 0, 0) = saved;
                    pv_arg.type = VAL_MATRIX;
                }
                if (arg2.type == VAL_SCALAR) {
                    apfc saved = arg2.v.scalar;
                    mat_zero(&arg2.v.matrix, 1, 1);
                    MAT_AT(&arg2.v.matrix, 0, 0) = saved;
                    arg2.type = VAL_MATRIX;
                }
                
                /* Check dimensions match */
                if (pv_arg.v.matrix.rows != arg2.v.matrix.rows ||
                    pv_arg.v.matrix.cols != arg2.v.matrix.cols) {
                    printf("Error: dimensions must match\n");
                    return 0;
                }
                
                result->type = VAL_MATRIX;
                result->v.matrix.rows = pv_arg.v.matrix.rows;
                result->v.matrix.cols = pv_arg.v.matrix.cols;
                
                for (r = 0; r < pv_arg.v.matrix.rows; r++) {
                    for (c = 0; c < pv_arg.v.matrix.cols; c++) {
                        int cmp = apf_cmp(&MAT_AT(&pv_arg.v.matrix, r, c).re,
                                         &MAT_AT(&arg2.v.matrix, r, c).re);
                        int res = 0;
                        if (str_eq(name, "eq")) res = (cmp == 0);
                        else if (str_eq(name, "ne")) res = (cmp != 0);
                        else if (str_eq(name, "lt")) res = (cmp < 0);
                        else if (str_eq(name, "le")) res = (cmp <= 0);
                        else if (str_eq(name, "gt")) res = (cmp > 0);
                        else if (str_eq(name, "ge")) res = (cmp >= 0);
                        apf_from_int(&MAT_AT(&result->v.matrix, r, c).re, res);
                        apf_zero(&MAT_AT(&result->v.matrix, r, c).im);
                    }
                }
                return 1;
            }
            
            /* Element-wise operations: times, rdivide, ldivide, plus, minus */
            if (str_eq(name, "times") || str_eq(name, "rdivide") || str_eq(name, "ldivide") ||
                str_eq(name, "plus") || str_eq(name, "minus")) {
                value_t arg2;
                int r, c;
                
                next_token();
                if (!parse_value(&pv_arg)) return 0;
                if (current_token.type != TOK_COMMA) {
                    printf("Error: %s requires two arguments\n", name);
                    return 0;
                }
                next_token();
                if (!parse_value(&arg2)) return 0;
                if (current_token.type != TOK_RPAREN) {
                    printf("Error: expected ')'\n");
                    return 0;
                }
                next_token();
                
                /* Convert scalars to 1x1 matrices */
                if (pv_arg.type == VAL_SCALAR) {
                    apfc saved = pv_arg.v.scalar;
                    mat_zero(&pv_arg.v.matrix, 1, 1);
                    MAT_AT(&pv_arg.v.matrix, 0, 0) = saved;
                    pv_arg.type = VAL_MATRIX;
                }
                if (arg2.type == VAL_SCALAR) {
                    apfc saved = arg2.v.scalar;
                    mat_zero(&arg2.v.matrix, 1, 1);
                    MAT_AT(&arg2.v.matrix, 0, 0) = saved;
                    arg2.type = VAL_MATRIX;
                }
                
                /* Check dimensions match */
                if (pv_arg.v.matrix.rows != arg2.v.matrix.rows ||
                    pv_arg.v.matrix.cols != arg2.v.matrix.cols) {
                    printf("Error: dimensions must match\n");
                    return 0;
                }
                
                result->type = VAL_MATRIX;
                result->v.matrix.rows = pv_arg.v.matrix.rows;
                result->v.matrix.cols = pv_arg.v.matrix.cols;
                
                for (r = 0; r < pv_arg.v.matrix.rows; r++) {
                    for (c = 0; c < pv_arg.v.matrix.cols; c++) {
                        if (str_eq(name, "times")) {
                            apfc_mul(&MAT_AT(&result->v.matrix, r, c),
                                    &MAT_AT(&pv_arg.v.matrix, r, c),
                                    &MAT_AT(&arg2.v.matrix, r, c));
                        } else if (str_eq(name, "rdivide")) {
                            apfc_div(&MAT_AT(&result->v.matrix, r, c),
                                    &MAT_AT(&pv_arg.v.matrix, r, c),
                                    &MAT_AT(&arg2.v.matrix, r, c));
                        } else if (str_eq(name, "ldivide")) {
                            apfc_div(&MAT_AT(&result->v.matrix, r, c),
                                    &MAT_AT(&arg2.v.matrix, r, c),
                                    &MAT_AT(&pv_arg.v.matrix, r, c));
                        } else if (str_eq(name, "plus")) {
                            apfc_add(&MAT_AT(&result->v.matrix, r, c),
                                    &MAT_AT(&pv_arg.v.matrix, r, c),
                                    &MAT_AT(&arg2.v.matrix, r, c));
                        } else if (str_eq(name, "minus")) {
                            apfc_sub(&MAT_AT(&result->v.matrix, r, c),
                                    &MAT_AT(&pv_arg.v.matrix, r, c),
                                    &MAT_AT(&arg2.v.matrix, r, c));
                        }
                    }
                }
                return 1;
            }
            
            /* conv(a, b) - convolution of two vectors */
            if (str_eq(name, "conv")) {
                value_t arg2;
                int na, nb, nc, i, j;
                
                next_token();
                if (!parse_value(&pv_arg)) return 0;
                if (current_token.type != TOK_COMMA) {
                    printf("Error: conv requires two arguments\n");
                    return 0;
                }
                next_token();
                if (!parse_value(&arg2)) return 0;
                if (current_token.type != TOK_RPAREN) {
                    printf("Error: expected ')'\n");
                    return 0;
                }
                next_token();
                
                /* Convert scalars to 1x1 matrices */
                if (pv_arg.type == VAL_SCALAR) {
                    apfc saved = pv_arg.v.scalar;
                    mat_zero(&pv_arg.v.matrix, 1, 1);
                    MAT_AT(&pv_arg.v.matrix, 0, 0) = saved;
                    pv_arg.type = VAL_MATRIX;
                }
                if (arg2.type == VAL_SCALAR) {
                    apfc saved = arg2.v.scalar;
                    mat_zero(&arg2.v.matrix, 1, 1);
                    MAT_AT(&arg2.v.matrix, 0, 0) = saved;
                    arg2.type = VAL_MATRIX;
                }
                
                na = pv_arg.v.matrix.rows * pv_arg.v.matrix.cols;
                nb = arg2.v.matrix.rows * arg2.v.matrix.cols;
                nc = na + nb - 1;
                
                if (nc > MAT_MAX_ROWS * MAT_MAX_COLS) {
                    printf("Error: result too large\n");
                    return 0;
                }
                
                result->type = VAL_MATRIX;
                mat_zero(&result->v.matrix, 1, nc);
                
                /* Convolution: c[k] = sum_j a[j] * b[k-j] */
                for (i = 0; i < nc; i++) {
                    apf_zero(&MAT_AT(&result->v.matrix, 0, i).re);
                    apf_zero(&MAT_AT(&result->v.matrix, 0, i).im);
                    for (j = 0; j < na; j++) {
                        int bj = i - j;
                        if (bj >= 0 && bj < nb) {
                            apfc prod;
                            int ra = j % pv_arg.v.matrix.rows;
                            int ca = j / pv_arg.v.matrix.rows;
                            int rb = bj % arg2.v.matrix.rows;
                            int cb = bj / arg2.v.matrix.rows;
                            apfc_mul(&prod, &MAT_AT(&pv_arg.v.matrix, ra, ca),
                                    &MAT_AT(&arg2.v.matrix, rb, cb));
                            apfc_add(&MAT_AT(&result->v.matrix, 0, i),
                                    &MAT_AT(&result->v.matrix, 0, i), &prod);
                        }
                    }
                }
                return 1;
            }
            
            /* deconv(u, v) - deconvolution / polynomial division */
            /* Returns [q, r] where u = conv(v, q) + r */
            if (str_eq(name, "deconv")) {
                value_t arg2;
                int nu, nv, nq, i, j;
                
                next_token();
                if (!parse_value(&pv_arg)) return 0;
                if (current_token.type != TOK_COMMA) {
                    printf("Error: deconv requires two arguments\n");
                    return 0;
                }
                next_token();
                if (!parse_value(&arg2)) return 0;
                if (current_token.type != TOK_RPAREN) {
                    printf("Error: expected ')'\n");
                    return 0;
                }
                next_token();
                
                if (pv_arg.type == VAL_SCALAR) {
                    apfc saved = pv_arg.v.scalar;
                    mat_zero(&pv_arg.v.matrix, 1, 1);
                    MAT_AT(&pv_arg.v.matrix, 0, 0) = saved;
                    pv_arg.type = VAL_MATRIX;
                }
                if (arg2.type == VAL_SCALAR) {
                    apfc saved = arg2.v.scalar;
                    mat_zero(&arg2.v.matrix, 1, 1);
                    MAT_AT(&arg2.v.matrix, 0, 0) = saved;
                    arg2.type = VAL_MATRIX;
                }
                
                nu = pv_arg.v.matrix.rows * pv_arg.v.matrix.cols;
                nv = arg2.v.matrix.rows * arg2.v.matrix.cols;
                
                if (nv > nu) {
                    /* Quotient is 0, remainder is u */
                    result->type = VAL_MATRIX;
                    mat_copy(&result->v.matrix, &pv_arg.v.matrix);
                    return 1;
                }
                
                nq = nu - nv + 1;
                
                /* Copy u to working array */
                mat_copy(&pv_tmp_mat, &pv_arg.v.matrix);
                
                result->type = VAL_MATRIX;
                mat_zero(&result->v.matrix, 1, nq);
                
                /* Polynomial long division */
                for (i = 0; i < nq; i++) {
                    int ru = i % pv_tmp_mat.rows;
                    int cu = i / pv_tmp_mat.rows;
                    int rv = 0, cv = 0;  /* First element of v */
                    apf_div(&MAT_AT(&result->v.matrix, 0, i).re,
                           &MAT_AT(&pv_tmp_mat, ru, cu).re,
                           &MAT_AT(&arg2.v.matrix, rv, cv).re);
                    apf_zero(&MAT_AT(&result->v.matrix, 0, i).im);
                    
                    /* Subtract q[i] * v from u */
                    for (j = 0; j < nv; j++) {
                        int ruj = (i + j) % pv_tmp_mat.rows;
                        int cuj = (i + j) / pv_tmp_mat.rows;
                        int rvj = j % arg2.v.matrix.rows;
                        int cvj = j / arg2.v.matrix.rows;
                        apf prod;
                        apf_mul(&prod, &MAT_AT(&result->v.matrix, 0, i).re,
                               &MAT_AT(&arg2.v.matrix, rvj, cvj).re);
                        apf_sub(&MAT_AT(&pv_tmp_mat, ruj, cuj).re,
                               &MAT_AT(&pv_tmp_mat, ruj, cuj).re, &prod);
                    }
                }
                return 1;
            }
            
            /* corrcoef(X, Y) - correlation coefficient */
            if (str_eq(name, "corrcoef")) {
                value_t arg2;
                apf mean_x, mean_y, sum_x, sum_y, sum_xy, sum_xx, sum_yy;
                apf diff_x, diff_y, prod, n_val;
                int i, n;
                
                next_token();
                if (!parse_value(&pv_arg)) return 0;
                if (current_token.type != TOK_COMMA) {
                    printf("Error: corrcoef requires two arguments\n");
                    return 0;
                }
                next_token();
                if (!parse_value(&arg2)) return 0;
                if (current_token.type != TOK_RPAREN) {
                    printf("Error: expected ')'\n");
                    return 0;
                }
                next_token();
                
                /* Convert scalars to 1x1 matrices */
                if (pv_arg.type == VAL_SCALAR) {
                    apfc saved = pv_arg.v.scalar;
                    mat_zero(&pv_arg.v.matrix, 1, 1);
                    MAT_AT(&pv_arg.v.matrix, 0, 0) = saved;
                    pv_arg.type = VAL_MATRIX;
                }
                if (arg2.type == VAL_SCALAR) {
                    apfc saved = arg2.v.scalar;
                    mat_zero(&arg2.v.matrix, 1, 1);
                    MAT_AT(&arg2.v.matrix, 0, 0) = saved;
                    arg2.type = VAL_MATRIX;
                }
                
                n = pv_arg.v.matrix.rows * pv_arg.v.matrix.cols;
                if (n != arg2.v.matrix.rows * arg2.v.matrix.cols) {
                    printf("Error: vectors must be same length\n");
                    return 0;
                }
                
                /* Compute means */
                apf_zero(&sum_x);
                apf_zero(&sum_y);
                for (i = 0; i < n; i++) {
                    int ri = i % pv_arg.v.matrix.rows;
                    int ci = i / pv_arg.v.matrix.rows;
                    int ri2 = i % arg2.v.matrix.rows;
                    int ci2 = i / arg2.v.matrix.rows;
                    apf_add(&sum_x, &sum_x, &MAT_AT(&pv_arg.v.matrix, ri, ci).re);
                    apf_add(&sum_y, &sum_y, &MAT_AT(&arg2.v.matrix, ri2, ci2).re);
                }
                apf_from_int(&n_val, n);
                apf_div(&mean_x, &sum_x, &n_val);
                apf_div(&mean_y, &sum_y, &n_val);
                
                /* Compute covariance and variances */
                apf_zero(&sum_xy);
                apf_zero(&sum_xx);
                apf_zero(&sum_yy);
                for (i = 0; i < n; i++) {
                    int ri = i % pv_arg.v.matrix.rows;
                    int ci = i / pv_arg.v.matrix.rows;
                    int ri2 = i % arg2.v.matrix.rows;
                    int ci2 = i / arg2.v.matrix.rows;
                    apf_sub(&diff_x, &MAT_AT(&pv_arg.v.matrix, ri, ci).re, &mean_x);
                    apf_sub(&diff_y, &MAT_AT(&arg2.v.matrix, ri2, ci2).re, &mean_y);
                    apf_mul(&prod, &diff_x, &diff_y);
                    apf_add(&sum_xy, &sum_xy, &prod);
                    apf_mul(&prod, &diff_x, &diff_x);
                    apf_add(&sum_xx, &sum_xx, &prod);
                    apf_mul(&prod, &diff_y, &diff_y);
                    apf_add(&sum_yy, &sum_yy, &prod);
                }
                
                /* r = sum_xy / sqrt(sum_xx * sum_yy) */
                apf_mul(&prod, &sum_xx, &sum_yy);
                apf_sqrt(&prod, &prod);
                if (apf_is_zero(&prod)) {
                    apf_from_int(&result->v.scalar.re, 0);
                } else {
                    apf_div(&result->v.scalar.re, &sum_xy, &prod);
                }
                apf_zero(&result->v.scalar.im);
                result->type = VAL_SCALAR;
                return 1;
            }
            
            /* movmean(v, k) - moving average with window k */
            if (str_eq(name, "movmean") || str_eq(name, "movsum")) {
                apfc k_val;
                int n, k, i, j, half;
                int is_mean = str_eq(name, "movmean");
                
                next_token();
                if (!parse_value(&pv_arg)) return 0;
                if (current_token.type != TOK_COMMA) {
                    printf("Error: %s requires (v, k)\n", name);
                    return 0;
                }
                next_token();
                if (!parse_expr(&k_val)) return 0;
                if (current_token.type != TOK_RPAREN) {
                    printf("Error: expected ')'\n");
                    return 0;
                }
                next_token();
                
                if (pv_arg.type == VAL_SCALAR) {
                    apfc saved = pv_arg.v.scalar;
                    mat_zero(&pv_arg.v.matrix, 1, 1);
                    MAT_AT(&pv_arg.v.matrix, 0, 0) = saved;
                    pv_arg.type = VAL_MATRIX;
                }
                
                n = pv_arg.v.matrix.rows * pv_arg.v.matrix.cols;
                k = apf_to_long(&k_val.re);
                if (k < 1) k = 1;
                if (k > n) k = n;
                half = k / 2;
                
                result->type = VAL_MATRIX;
                mat_zero(&result->v.matrix, pv_arg.v.matrix.rows, pv_arg.v.matrix.cols);
                
                for (i = 0; i < n; i++) {
                    apf sum, count;
                    int ri = i % pv_arg.v.matrix.rows;
                    int ci = i / pv_arg.v.matrix.rows;
                    int start = i - half;
                    int end = start + k;
                    int cnt = 0;
                    
                    if (start < 0) start = 0;
                    if (end > n) end = n;
                    
                    apf_zero(&sum);
                    for (j = start; j < end; j++) {
                        int rj = j % pv_arg.v.matrix.rows;
                        int cj = j / pv_arg.v.matrix.rows;
                        apf_add(&sum, &sum, &MAT_AT(&pv_arg.v.matrix, rj, cj).re);
                        cnt++;
                    }
                    
                    if (is_mean && cnt > 0) {
                        apf_from_int(&count, cnt);
                        apf_div(&MAT_AT(&result->v.matrix, ri, ci).re, &sum, &count);
                    } else {
                        apf_copy(&MAT_AT(&result->v.matrix, ri, ci).re, &sum);
                    }
                    apf_zero(&MAT_AT(&result->v.matrix, ri, ci).im);
                }
                return 1;
            }
            
            /* kron(A, B) - Kronecker product */
            if (str_eq(name, "kron")) {
                value_t arg2;
                int ra, ca, rb, cb, i, j, m, n;
                
                next_token();
                if (!parse_value(&pv_arg)) return 0;
                if (current_token.type != TOK_COMMA) {
                    printf("Error: kron requires two arguments\n");
                    return 0;
                }
                next_token();
                if (!parse_value(&arg2)) return 0;
                if (current_token.type != TOK_RPAREN) {
                    printf("Error: expected ')'\n");
                    return 0;
                }
                next_token();
                
                if (pv_arg.type == VAL_SCALAR) {
                    apfc saved = pv_arg.v.scalar;
                    mat_zero(&pv_arg.v.matrix, 1, 1);
                    MAT_AT(&pv_arg.v.matrix, 0, 0) = saved;
                    pv_arg.type = VAL_MATRIX;
                }
                if (arg2.type == VAL_SCALAR) {
                    apfc saved = arg2.v.scalar;
                    mat_zero(&arg2.v.matrix, 1, 1);
                    MAT_AT(&arg2.v.matrix, 0, 0) = saved;
                    arg2.type = VAL_MATRIX;
                }
                
                ra = pv_arg.v.matrix.rows;
                ca = pv_arg.v.matrix.cols;
                rb = arg2.v.matrix.rows;
                cb = arg2.v.matrix.cols;
                
                if (ra * rb > MAT_MAX_ROWS || ca * cb > MAT_MAX_COLS) {
                    printf("Error: result too large\n");
                    return 0;
                }
                
                result->type = VAL_MATRIX;
                mat_zero(&result->v.matrix, ra * rb, ca * cb);
                
                for (i = 0; i < ra; i++) {
                    for (j = 0; j < ca; j++) {
                        for (m = 0; m < rb; m++) {
                            for (n = 0; n < cb; n++) {
                                apfc_mul(&MAT_AT(&result->v.matrix, i*rb + m, j*cb + n),
                                        &MAT_AT(&pv_arg.v.matrix, i, j),
                                        &MAT_AT(&arg2.v.matrix, m, n));
                            }
                        }
                    }
                }
                return 1;
            }
            
            /* interp1(x, y, xi) - linear interpolation */
            if (str_eq(name, "interp1")) {
                value_t arg2, arg3;
                int n, i, j;
                
                next_token();
                if (!parse_value(&pv_arg)) return 0;
                if (current_token.type != TOK_COMMA) {
                    printf("Error: interp1 requires (x, y, xi)\n");
                    return 0;
                }
                next_token();
                if (!parse_value(&arg2)) return 0;
                if (current_token.type != TOK_COMMA) {
                    printf("Error: interp1 requires (x, y, xi)\n");
                    return 0;
                }
                next_token();
                if (!parse_value(&arg3)) return 0;
                if (current_token.type != TOK_RPAREN) {
                    printf("Error: expected ')'\n");
                    return 0;
                }
                next_token();
                
                /* Convert scalars to matrices */
                if (pv_arg.type == VAL_SCALAR) {
                    apfc saved = pv_arg.v.scalar;
                    mat_zero(&pv_arg.v.matrix, 1, 1);
                    MAT_AT(&pv_arg.v.matrix, 0, 0) = saved;
                    pv_arg.type = VAL_MATRIX;
                }
                if (arg2.type == VAL_SCALAR) {
                    apfc saved = arg2.v.scalar;
                    mat_zero(&arg2.v.matrix, 1, 1);
                    MAT_AT(&arg2.v.matrix, 0, 0) = saved;
                    arg2.type = VAL_MATRIX;
                }
                
                n = pv_arg.v.matrix.rows * pv_arg.v.matrix.cols;
                
                if (arg3.type == VAL_SCALAR) {
                    /* Single query point */
                    apf xi = arg3.v.scalar.re;
                    result->type = VAL_SCALAR;
                    
                    /* Find interval containing xi */
                    for (i = 0; i < n - 1; i++) {
                        int r1 = i % pv_arg.v.matrix.rows;
                        int c1 = i / pv_arg.v.matrix.rows;
                        int r2 = (i+1) % pv_arg.v.matrix.rows;
                        int c2 = (i+1) / pv_arg.v.matrix.rows;
                        
                        if (apf_cmp(&xi, &MAT_AT(&pv_arg.v.matrix, r1, c1).re) >= 0 &&
                            apf_cmp(&xi, &MAT_AT(&pv_arg.v.matrix, r2, c2).re) <= 0) {
                            /* Interpolate: y = y1 + (y2-y1) * (xi-x1) / (x2-x1) */
                            apf x1 = MAT_AT(&pv_arg.v.matrix, r1, c1).re;
                            apf x2 = MAT_AT(&pv_arg.v.matrix, r2, c2).re;
                            apf y1 = MAT_AT(&arg2.v.matrix, r1, c1).re;
                            apf y2 = MAT_AT(&arg2.v.matrix, r2, c2).re;
                            apf dx, dy, t;
                            apf_sub(&dx, &x2, &x1);
                            apf_sub(&dy, &y2, &y1);
                            apf_sub(&t, &xi, &x1);
                            apf_div(&t, &t, &dx);
                            apf_mul(&t, &t, &dy);
                            apf_add(&result->v.scalar.re, &y1, &t);
                            apf_zero(&result->v.scalar.im);
                            return 1;
                        }
                    }
                    /* Extrapolate using first/last interval */
                    if (n >= 2) {
                        int r1 = 0, c1 = 0;
                        int r2 = 1 % pv_arg.v.matrix.rows;
                        int c2 = 1 / pv_arg.v.matrix.rows;
                        apf x1 = MAT_AT(&pv_arg.v.matrix, r1, c1).re;
                        apf x2 = MAT_AT(&pv_arg.v.matrix, r2, c2).re;
                        apf y1 = MAT_AT(&arg2.v.matrix, r1, c1).re;
                        apf y2 = MAT_AT(&arg2.v.matrix, r2, c2).re;
                        apf dx, dy, t;
                        apf_sub(&dx, &x2, &x1);
                        apf_sub(&dy, &y2, &y1);
                        apf_sub(&t, &xi, &x1);
                        apf_div(&t, &t, &dx);
                        apf_mul(&t, &t, &dy);
                        apf_add(&result->v.scalar.re, &y1, &t);
                    } else {
                        apf_copy(&result->v.scalar.re, &MAT_AT(&arg2.v.matrix, 0, 0).re);
                    }
                    apf_zero(&result->v.scalar.im);
                } else {
                    /* Vector of query points */
                    int nq = arg3.v.matrix.rows * arg3.v.matrix.cols;
                    result->type = VAL_MATRIX;
                    mat_zero(&result->v.matrix, arg3.v.matrix.rows, arg3.v.matrix.cols);
                    
                    for (j = 0; j < nq; j++) {
                        int rq = j % arg3.v.matrix.rows;
                        int cq = j / arg3.v.matrix.rows;
                        apf xi = MAT_AT(&arg3.v.matrix, rq, cq).re;
                        int found = 0;
                        
                        for (i = 0; i < n - 1 && !found; i++) {
                            int r1 = i % pv_arg.v.matrix.rows;
                            int c1 = i / pv_arg.v.matrix.rows;
                            int r2 = (i+1) % pv_arg.v.matrix.rows;
                            int c2 = (i+1) / pv_arg.v.matrix.rows;
                            
                            if (apf_cmp(&xi, &MAT_AT(&pv_arg.v.matrix, r1, c1).re) >= 0 &&
                                apf_cmp(&xi, &MAT_AT(&pv_arg.v.matrix, r2, c2).re) <= 0) {
                                apf x1 = MAT_AT(&pv_arg.v.matrix, r1, c1).re;
                                apf x2 = MAT_AT(&pv_arg.v.matrix, r2, c2).re;
                                apf y1 = MAT_AT(&arg2.v.matrix, r1, c1).re;
                                apf y2 = MAT_AT(&arg2.v.matrix, r2, c2).re;
                                apf dx, dy, t;
                                apf_sub(&dx, &x2, &x1);
                                apf_sub(&dy, &y2, &y1);
                                apf_sub(&t, &xi, &x1);
                                apf_div(&t, &t, &dx);
                                apf_mul(&t, &t, &dy);
                                apf_add(&MAT_AT(&result->v.matrix, rq, cq).re, &y1, &t);
                                found = 1;
                            }
                        }
                        if (!found && n >= 2) {
                            /* Extrapolate */
                            int r1 = 0, c1 = 0;
                            int r2 = 1 % pv_arg.v.matrix.rows;
                            int c2 = 1 / pv_arg.v.matrix.rows;
                            apf x1 = MAT_AT(&pv_arg.v.matrix, r1, c1).re;
                            apf x2 = MAT_AT(&pv_arg.v.matrix, r2, c2).re;
                            apf y1 = MAT_AT(&arg2.v.matrix, r1, c1).re;
                            apf y2 = MAT_AT(&arg2.v.matrix, r2, c2).re;
                            apf dx, dy, t;
                            apf_sub(&dx, &x2, &x1);
                            apf_sub(&dy, &y2, &y1);
                            apf_sub(&t, &xi, &x1);
                            apf_div(&t, &t, &dx);
                            apf_mul(&t, &t, &dy);
                            apf_add(&MAT_AT(&result->v.matrix, rq, cq).re, &y1, &t);
                        }
                    }
                }
                return 1;
            }
            
            /* trapz(y) or trapz(x, y) - trapezoidal integration */
            if (str_eq(name, "trapz")) {
                int n, i;
                apf sum, h, yi, yi1, avg;
                int has_x = 0;
                value_t arg2;
                
                next_token();
                if (!parse_value(&pv_arg)) return 0;
                
                if (current_token.type == TOK_COMMA) {
                    has_x = 1;
                    next_token();
                    if (!parse_value(&arg2)) return 0;
                }
                
                if (current_token.type != TOK_RPAREN) {
                    printf("Error: expected ')'\n");
                    return 0;
                }
                next_token();
                
                if (pv_arg.type == VAL_SCALAR) {
                    apfc saved = pv_arg.v.scalar;
                    mat_zero(&pv_arg.v.matrix, 1, 1);
                    MAT_AT(&pv_arg.v.matrix, 0, 0) = saved;
                    pv_arg.type = VAL_MATRIX;
                }
                
                if (has_x) {
                    /* trapz(x, y): use x values for spacing */
                    if (arg2.type == VAL_SCALAR) {
                        apfc saved = arg2.v.scalar;
                        mat_zero(&arg2.v.matrix, 1, 1);
                        MAT_AT(&arg2.v.matrix, 0, 0) = saved;
                        arg2.type = VAL_MATRIX;
                    }
                    n = arg2.v.matrix.rows * arg2.v.matrix.cols;
                } else {
                    n = pv_arg.v.matrix.rows * pv_arg.v.matrix.cols;
                }
                
                apf_zero(&sum);
                for (i = 0; i < n - 1; i++) {
                    int ri = i % pv_arg.v.matrix.rows;
                    int ci = i / pv_arg.v.matrix.rows;
                    int ri1 = (i+1) % pv_arg.v.matrix.rows;
                    int ci1 = (i+1) / pv_arg.v.matrix.rows;
                    
                    if (has_x) {
                        int xi = i % arg2.v.matrix.rows;
                        int xci = i / arg2.v.matrix.rows;
                        int xi1 = (i+1) % arg2.v.matrix.rows;
                        int xci1 = (i+1) / arg2.v.matrix.rows;
                        apf_sub(&h, &MAT_AT(&pv_arg.v.matrix, xi1, xci1).re,
                               &MAT_AT(&pv_arg.v.matrix, xi, xci).re);
                        apf_copy(&yi, &MAT_AT(&arg2.v.matrix, ri, ci).re);
                        apf_copy(&yi1, &MAT_AT(&arg2.v.matrix, ri1, ci1).re);
                    } else {
                        apf_from_int(&h, 1);
                        apf_copy(&yi, &MAT_AT(&pv_arg.v.matrix, ri, ci).re);
                        apf_copy(&yi1, &MAT_AT(&pv_arg.v.matrix, ri1, ci1).re);
                    }
                    
                    apf_add(&avg, &yi, &yi1);
                    apf_from_int(&yi, 2);
                    apf_div(&avg, &avg, &yi);
                    apf_mul(&avg, &avg, &h);
                    apf_add(&sum, &sum, &avg);
                }
                
                result->type = VAL_SCALAR;
                apf_copy(&result->v.scalar.re, &sum);
                apf_zero(&result->v.scalar.im);
                return 1;
            }
            
            /* gradient(y) - numerical gradient */
            if (str_eq(name, "gradient")) {
                int n, i;
                
                next_token();
                if (!parse_value(&pv_arg)) return 0;
                if (current_token.type != TOK_RPAREN) {
                    printf("Error: expected ')'\n");
                    return 0;
                }
                next_token();
                
                if (pv_arg.type == VAL_SCALAR) {
                    result->type = VAL_SCALAR;
                    apf_zero(&result->v.scalar.re);
                    apf_zero(&result->v.scalar.im);
                    return 1;
                }
                
                n = pv_arg.v.matrix.rows * pv_arg.v.matrix.cols;
                result->type = VAL_MATRIX;
                mat_zero(&result->v.matrix, pv_arg.v.matrix.rows, pv_arg.v.matrix.cols);
                
                for (i = 0; i < n; i++) {
                    int ri = i % pv_arg.v.matrix.rows;
                    int ci = i / pv_arg.v.matrix.rows;
                    apf grad;
                    
                    if (i == 0) {
                        /* Forward difference */
                        int ri1 = 1 % pv_arg.v.matrix.rows;
                        int ci1 = 1 / pv_arg.v.matrix.rows;
                        apf_sub(&grad, &MAT_AT(&pv_arg.v.matrix, ri1, ci1).re,
                               &MAT_AT(&pv_arg.v.matrix, ri, ci).re);
                    } else if (i == n - 1) {
                        /* Backward difference */
                        int rim1 = (n-2) % pv_arg.v.matrix.rows;
                        int cim1 = (n-2) / pv_arg.v.matrix.rows;
                        apf_sub(&grad, &MAT_AT(&pv_arg.v.matrix, ri, ci).re,
                               &MAT_AT(&pv_arg.v.matrix, rim1, cim1).re);
                    } else {
                        /* Central difference */
                        int ri1 = (i+1) % pv_arg.v.matrix.rows;
                        int ci1 = (i+1) / pv_arg.v.matrix.rows;
                        int rim1 = (i-1) % pv_arg.v.matrix.rows;
                        int cim1 = (i-1) / pv_arg.v.matrix.rows;
                        apf two;
                        apf_sub(&grad, &MAT_AT(&pv_arg.v.matrix, ri1, ci1).re,
                               &MAT_AT(&pv_arg.v.matrix, rim1, cim1).re);
                        apf_from_int(&two, 2);
                        apf_div(&grad, &grad, &two);
                    }
                    apf_copy(&MAT_AT(&result->v.matrix, ri, ci).re, &grad);
                    apf_zero(&MAT_AT(&result->v.matrix, ri, ci).im);
                }
                return 1;
            }
            
            /* linsolve(A, b) or mldivide(A, b) - solve linear system Ax = b */
            if (str_eq(name, "linsolve") || str_eq(name, "mldivide")) {
                value_t arg2;
                
                next_token();
                if (!parse_value(&pv_arg)) return 0;
                if (current_token.type != TOK_COMMA) {
                    printf("Error: linsolve requires (A, b)\n");
                    return 0;
                }
                next_token();
                if (!parse_value(&arg2)) return 0;
                if (current_token.type != TOK_RPAREN) {
                    printf("Error: expected ')'\n");
                    return 0;
                }
                next_token();
                
                if (pv_arg.type == VAL_SCALAR) {
                    apfc saved = pv_arg.v.scalar;
                    mat_zero(&pv_arg.v.matrix, 1, 1);
                    MAT_AT(&pv_arg.v.matrix, 0, 0) = saved;
                    pv_arg.type = VAL_MATRIX;
                }
                if (arg2.type == VAL_SCALAR) {
                    apfc saved = arg2.v.scalar;
                    mat_zero(&arg2.v.matrix, 1, 1);
                    MAT_AT(&arg2.v.matrix, 0, 0) = saved;
                    arg2.type = VAL_MATRIX;
                }
                
                result->type = VAL_MATRIX;
                if (!mat_solve(&result->v.matrix, &pv_arg.v.matrix, &arg2.v.matrix)) {
                    printf("Error: singular matrix\n");
                    return 0;
                }
                return 1;
            }
            
            /* vander(v, n) - Vandermonde matrix */
            if (str_eq(name, "vander")) {
                int n, cols, i, j;
                apfc n_arg;
                
                next_token();
                if (!parse_value(&pv_arg)) return 0;
                
                if (pv_arg.type == VAL_SCALAR) {
                    apfc saved = pv_arg.v.scalar;
                    mat_zero(&pv_arg.v.matrix, 1, 1);
                    MAT_AT(&pv_arg.v.matrix, 0, 0) = saved;
                    pv_arg.type = VAL_MATRIX;
                }
                
                n = pv_arg.v.matrix.rows * pv_arg.v.matrix.cols;
                cols = n;
                
                if (current_token.type == TOK_COMMA) {
                    next_token();
                    if (!parse_expr(&n_arg)) return 0;
                    cols = apf_to_long(&n_arg.re);
                    if (cols < 1) cols = 1;
                    if (cols > MAT_MAX_COLS) cols = MAT_MAX_COLS;
                }
                
                if (current_token.type != TOK_RPAREN) {
                    printf("Error: expected ')'\n");
                    return 0;
                }
                next_token();
                
                result->type = VAL_MATRIX;
                mat_zero(&result->v.matrix, n, cols);
                
                /* V(i,j) = v(i)^(cols-1-j) */
                for (i = 0; i < n; i++) {
                    int ri = i % pv_arg.v.matrix.rows;
                    int ci = i / pv_arg.v.matrix.rows;
                    apf base, power;
                    apf_copy(&base, &MAT_AT(&pv_arg.v.matrix, ri, ci).re);
                    
                    for (j = 0; j < cols; j++) {
                        int exp = cols - 1 - j;
                        if (exp == 0) {
                            apf_from_int(&MAT_AT(&result->v.matrix, i, j).re, 1);
                        } else {
                            apf_from_int(&power, 1);
                            apf_copy(&power, &base);
                            {
                                int k;
                                for (k = 1; k < exp; k++) {
                                    apf_mul(&power, &power, &base);
                                }
                            }
                            apf_copy(&MAT_AT(&result->v.matrix, i, j).re, &power);
                        }
                        apf_zero(&MAT_AT(&result->v.matrix, i, j).im);
                    }
                }
                return 1;
            }
            
            /* toeplitz(c) or toeplitz(c, r) - Toeplitz matrix */
            if (str_eq(name, "toeplitz")) {
                int nc, nr, i, j;
                value_t arg2;
                int has_r = 0;
                
                next_token();
                if (!parse_value(&pv_arg)) return 0;
                
                if (pv_arg.type == VAL_SCALAR) {
                    apfc saved = pv_arg.v.scalar;
                    mat_zero(&pv_arg.v.matrix, 1, 1);
                    MAT_AT(&pv_arg.v.matrix, 0, 0) = saved;
                    pv_arg.type = VAL_MATRIX;
                }
                
                if (current_token.type == TOK_COMMA) {
                    has_r = 1;
                    next_token();
                    if (!parse_value(&arg2)) return 0;
                    if (arg2.type == VAL_SCALAR) {
                        apfc saved = arg2.v.scalar;
                        mat_zero(&arg2.v.matrix, 1, 1);
                        MAT_AT(&arg2.v.matrix, 0, 0) = saved;
                        arg2.type = VAL_MATRIX;
                    }
                }
                
                if (current_token.type != TOK_RPAREN) {
                    printf("Error: expected ')'\n");
                    return 0;
                }
                next_token();
                
                nc = pv_arg.v.matrix.rows * pv_arg.v.matrix.cols;
                nr = has_r ? (arg2.v.matrix.rows * arg2.v.matrix.cols) : nc;
                
                if (nc > MAT_MAX_ROWS || nr > MAT_MAX_COLS) {
                    printf("Error: result too large\n");
                    return 0;
                }
                
                result->type = VAL_MATRIX;
                mat_zero(&result->v.matrix, nc, nr);
                
                for (i = 0; i < nc; i++) {
                    for (j = 0; j < nr; j++) {
                        int idx = i - j;
                        if (idx >= 0) {
                            /* From first column (c) */
                            int ri = idx % pv_arg.v.matrix.rows;
                            int ci = idx / pv_arg.v.matrix.rows;
                            MAT_AT(&result->v.matrix, i, j) = MAT_AT(&pv_arg.v.matrix, ri, ci);
                        } else {
                            /* From first row (r or c conjugate) */
                            idx = -idx;
                            if (has_r) {
                                int ri = idx % arg2.v.matrix.rows;
                                int ci = idx / arg2.v.matrix.rows;
                                MAT_AT(&result->v.matrix, i, j) = MAT_AT(&arg2.v.matrix, ri, ci);
                            } else {
                                int ri = idx % pv_arg.v.matrix.rows;
                                int ci = idx / pv_arg.v.matrix.rows;
                                MAT_AT(&result->v.matrix, i, j) = MAT_AT(&pv_arg.v.matrix, ri, ci);
                            }
                        }
                    }
                }
                return 1;
            }
            
            /* hankel(c) or hankel(c, r) - Hankel matrix */
            if (str_eq(name, "hankel")) {
                int nc, nr, i, j;
                value_t arg2;
                int has_r = 0;
                
                next_token();
                if (!parse_value(&pv_arg)) return 0;
                
                if (pv_arg.type == VAL_SCALAR) {
                    apfc saved = pv_arg.v.scalar;
                    mat_zero(&pv_arg.v.matrix, 1, 1);
                    MAT_AT(&pv_arg.v.matrix, 0, 0) = saved;
                    pv_arg.type = VAL_MATRIX;
                }
                
                if (current_token.type == TOK_COMMA) {
                    has_r = 1;
                    next_token();
                    if (!parse_value(&arg2)) return 0;
                    if (arg2.type == VAL_SCALAR) {
                        apfc saved = arg2.v.scalar;
                        mat_zero(&arg2.v.matrix, 1, 1);
                        MAT_AT(&arg2.v.matrix, 0, 0) = saved;
                        arg2.type = VAL_MATRIX;
                    }
                }
                
                if (current_token.type != TOK_RPAREN) {
                    printf("Error: expected ')'\n");
                    return 0;
                }
                next_token();
                
                nc = pv_arg.v.matrix.rows * pv_arg.v.matrix.cols;
                nr = has_r ? (arg2.v.matrix.rows * arg2.v.matrix.cols) : nc;
                
                if (nc > MAT_MAX_ROWS || nr > MAT_MAX_COLS) {
                    printf("Error: result too large\n");
                    return 0;
                }
                
                result->type = VAL_MATRIX;
                mat_zero(&result->v.matrix, nc, nr);
                
                for (i = 0; i < nc; i++) {
                    for (j = 0; j < nr; j++) {
                        int idx = i + j;
                        if (idx < nc) {
                            /* From c */
                            int ri = idx % pv_arg.v.matrix.rows;
                            int ci = idx / pv_arg.v.matrix.rows;
                            MAT_AT(&result->v.matrix, i, j) = MAT_AT(&pv_arg.v.matrix, ri, ci);
                        } else if (has_r) {
                            /* From r */
                            idx = idx - nc + 1;
                            if (idx < nr) {
                                int ri = idx % arg2.v.matrix.rows;
                                int ci = idx / arg2.v.matrix.rows;
                                MAT_AT(&result->v.matrix, i, j) = MAT_AT(&arg2.v.matrix, ri, ci);
                            }
                        }
                    }
                }
                return 1;
            }
            
            /* blkdiag(A, B) - block diagonal matrix */
            if (str_eq(name, "blkdiag")) {
                value_t arg2;
                int ra, ca, rb, cb, i, j;
                
                next_token();
                if (!parse_value(&pv_arg)) return 0;
                if (current_token.type != TOK_COMMA) {
                    printf("Error: blkdiag requires two arguments\n");
                    return 0;
                }
                next_token();
                if (!parse_value(&arg2)) return 0;
                if (current_token.type != TOK_RPAREN) {
                    printf("Error: expected ')'\n");
                    return 0;
                }
                next_token();
                
                if (pv_arg.type == VAL_SCALAR) {
                    apfc saved = pv_arg.v.scalar;
                    mat_zero(&pv_arg.v.matrix, 1, 1);
                    MAT_AT(&pv_arg.v.matrix, 0, 0) = saved;
                    pv_arg.type = VAL_MATRIX;
                }
                if (arg2.type == VAL_SCALAR) {
                    apfc saved = arg2.v.scalar;
                    mat_zero(&arg2.v.matrix, 1, 1);
                    MAT_AT(&arg2.v.matrix, 0, 0) = saved;
                    arg2.type = VAL_MATRIX;
                }
                
                ra = pv_arg.v.matrix.rows;
                ca = pv_arg.v.matrix.cols;
                rb = arg2.v.matrix.rows;
                cb = arg2.v.matrix.cols;
                
                if (ra + rb > MAT_MAX_ROWS || ca + cb > MAT_MAX_COLS) {
                    printf("Error: result too large\n");
                    return 0;
                }
                
                result->type = VAL_MATRIX;
                mat_zero(&result->v.matrix, ra + rb, ca + cb);
                
                /* Copy A to top-left */
                for (i = 0; i < ra; i++) {
                    for (j = 0; j < ca; j++) {
                        MAT_AT(&result->v.matrix, i, j) = MAT_AT(&pv_arg.v.matrix, i, j);
                    }
                }
                /* Copy B to bottom-right */
                for (i = 0; i < rb; i++) {
                    for (j = 0; j < cb; j++) {
                        MAT_AT(&result->v.matrix, ra + i, ca + j) = MAT_AT(&arg2.v.matrix, i, j);
                    }
                }
                return 1;
            }
            
            /* sub2ind([rows, cols], r, c) - convert subscripts to linear index */
            if (str_eq(name, "sub2ind")) {
                apfc row_arg, col_arg;
                int rows, r_idx, c_idx;
                
                next_token();
                if (!parse_value(&pv_arg)) return 0;
                if (current_token.type != TOK_COMMA) {
                    printf("Error: sub2ind requires (dims, row, col)\n");
                    return 0;
                }
                next_token();
                if (!parse_expr(&row_arg)) return 0;
                if (current_token.type != TOK_COMMA) {
                    printf("Error: sub2ind requires (dims, row, col)\n");
                    return 0;
                }
                next_token();
                if (!parse_expr(&col_arg)) return 0;
                if (current_token.type != TOK_RPAREN) {
                    printf("Error: expected ')'\n");
                    return 0;
                }
                next_token();
                
                if (pv_arg.type == VAL_MATRIX) {
                    rows = apf_to_long(&MAT_AT(&pv_arg.v.matrix, 0, 0).re);
                } else {
                    rows = apf_to_long(&pv_arg.v.scalar.re);
                }
                
                r_idx = apf_to_long(&row_arg.re);
                c_idx = apf_to_long(&col_arg.re);
                
                /* MATLAB uses 1-based indexing, linear index is column-major */
                result->type = VAL_SCALAR;
                apf_from_int(&result->v.scalar.re, (c_idx - 1) * rows + r_idx);
                apf_zero(&result->v.scalar.im);
                return 1;
            }
            
            /* ind2sub([rows, cols], idx) - convert linear index to subscripts */
            if (str_eq(name, "ind2sub")) {
                apfc idx_arg;
                int rows, idx, r_idx, c_idx;
                
                next_token();
                if (!parse_value(&pv_arg)) return 0;
                if (current_token.type != TOK_COMMA) {
                    printf("Error: ind2sub requires (dims, idx)\n");
                    return 0;
                }
                next_token();
                if (!parse_expr(&idx_arg)) return 0;
                if (current_token.type != TOK_RPAREN) {
                    printf("Error: expected ')'\n");
                    return 0;
                }
                next_token();
                
                if (pv_arg.type == VAL_MATRIX) {
                    rows = apf_to_long(&MAT_AT(&pv_arg.v.matrix, 0, 0).re);
                } else {
                    rows = apf_to_long(&pv_arg.v.scalar.re);
                }
                
                idx = apf_to_long(&idx_arg.re);
                
                /* Convert 1-based linear index to row, col (column-major) */
                c_idx = (idx - 1) / rows + 1;
                r_idx = (idx - 1) % rows + 1;
                
                result->type = VAL_MATRIX;
                mat_zero(&result->v.matrix, 1, 2);
                apf_from_int(&MAT_AT(&result->v.matrix, 0, 0).re, r_idx);
                apf_from_int(&MAT_AT(&result->v.matrix, 0, 1).re, c_idx);
                return 1;
            }
            
            /* clamp(x, lo, hi) - limit value to range [lo, hi] */
            if (str_eq(name, "clamp")) {
                apfc lo, hi;
                int i, j;
                
                next_token();
                if (!parse_value(&pv_arg)) return 0;
                if (current_token.type != TOK_COMMA) {
                    printf("Error: clamp requires (x, lo, hi)\n");
                    return 0;
                }
                next_token();
                if (!parse_expr(&lo)) return 0;
                if (current_token.type != TOK_COMMA) {
                    printf("Error: clamp requires (x, lo, hi)\n");
                    return 0;
                }
                next_token();
                if (!parse_expr(&hi)) return 0;
                if (current_token.type != TOK_RPAREN) {
                    printf("Error: expected ')'\n");
                    return 0;
                }
                next_token();
                
                if (pv_arg.type == VAL_SCALAR) {
                    result->type = VAL_SCALAR;
                    if (apf_cmp(&pv_arg.v.scalar.re, &lo.re) < 0) {
                        apf_copy(&result->v.scalar.re, &lo.re);
                    } else if (apf_cmp(&pv_arg.v.scalar.re, &hi.re) > 0) {
                        apf_copy(&result->v.scalar.re, &hi.re);
                    } else {
                        apf_copy(&result->v.scalar.re, &pv_arg.v.scalar.re);
                    }
                    apf_zero(&result->v.scalar.im);
                } else {
                    result->type = VAL_MATRIX;
                    mat_zero(&result->v.matrix, pv_arg.v.matrix.rows, pv_arg.v.matrix.cols);
                    for (i = 0; i < pv_arg.v.matrix.rows; i++) {
                        for (j = 0; j < pv_arg.v.matrix.cols; j++) {
                            if (apf_cmp(&MAT_AT(&pv_arg.v.matrix, i, j).re, &lo.re) < 0) {
                                apf_copy(&MAT_AT(&result->v.matrix, i, j).re, &lo.re);
                            } else if (apf_cmp(&MAT_AT(&pv_arg.v.matrix, i, j).re, &hi.re) > 0) {
                                apf_copy(&MAT_AT(&result->v.matrix, i, j).re, &hi.re);
                            } else {
                                apf_copy(&MAT_AT(&result->v.matrix, i, j).re, &MAT_AT(&pv_arg.v.matrix, i, j).re);
                            }
                            apf_zero(&MAT_AT(&result->v.matrix, i, j).im);
                        }
                    }
                }
                return 1;
            }
            
            /* isapprox(a, b) or isapprox(a, b, tol) - approximate equality */
            if (str_eq(name, "isapprox") || str_eq(name, "approxeq")) {
                value_t arg2;
                apfc tol;
                apf diff, abs_diff;
                int equal = 1;
                
                next_token();
                if (!parse_value(&pv_arg)) return 0;
                if (current_token.type != TOK_COMMA) {
                    printf("Error: %s requires at least two arguments\n", name);
                    return 0;
                }
                next_token();
                if (!parse_value(&arg2)) return 0;
                
                /* Default tolerance */
                apf_from_str(&tol.re, "1e-10");
                apf_zero(&tol.im);
                
                if (current_token.type == TOK_COMMA) {
                    next_token();
                    if (!parse_expr(&tol)) return 0;
                }
                
                if (current_token.type != TOK_RPAREN) {
                    printf("Error: expected ')'\n");
                    return 0;
                }
                next_token();
                
                /* Compare values */
                if (pv_arg.type == VAL_SCALAR && arg2.type == VAL_SCALAR) {
                    apf_sub(&diff, &pv_arg.v.scalar.re, &arg2.v.scalar.re);
                    apf_abs(&abs_diff, &diff);
                    equal = (apf_cmp(&abs_diff, &tol.re) <= 0);
                } else {
                    /* Matrix comparison */
                    int i, j;
                    if (pv_arg.type == VAL_SCALAR) {
                        apfc saved = pv_arg.v.scalar;
                        mat_zero(&pv_arg.v.matrix, 1, 1);
                        MAT_AT(&pv_arg.v.matrix, 0, 0) = saved;
                        pv_arg.type = VAL_MATRIX;
                    }
                    if (arg2.type == VAL_SCALAR) {
                        apfc saved = arg2.v.scalar;
                        mat_zero(&arg2.v.matrix, 1, 1);
                        MAT_AT(&arg2.v.matrix, 0, 0) = saved;
                        arg2.type = VAL_MATRIX;
                    }
                    
                    if (pv_arg.v.matrix.rows != arg2.v.matrix.rows ||
                        pv_arg.v.matrix.cols != arg2.v.matrix.cols) {
                        equal = 0;
                    } else {
                        for (i = 0; i < pv_arg.v.matrix.rows && equal; i++) {
                            for (j = 0; j < pv_arg.v.matrix.cols && equal; j++) {
                                apf_sub(&diff, &MAT_AT(&pv_arg.v.matrix, i, j).re,
                                       &MAT_AT(&arg2.v.matrix, i, j).re);
                                apf_abs(&abs_diff, &diff);
                                if (apf_cmp(&abs_diff, &tol.re) > 0) {
                                    equal = 0;
                                }
                            }
                        }
                    }
                }
                
                result->type = VAL_SCALAR;
                apf_from_int(&result->v.scalar.re, equal);
                apf_zero(&result->v.scalar.im);
                return 1;
            }
            
            /* colon(a, b) or colon(a, d, b) - explicit colon operator */
            if (str_eq(name, "colon")) {
                apfc a, b, d;
                int n, i;
                
                next_token();
                if (!parse_expr(&a)) return 0;
                if (current_token.type != TOK_COMMA) {
                    printf("Error: colon requires (a, b) or (a, d, b)\n");
                    return 0;
                }
                next_token();
                if (!parse_expr(&b)) return 0;
                
                /* Default step of 1 */
                apf_from_int(&d.re, 1);
                apf_zero(&d.im);
                
                if (current_token.type == TOK_COMMA) {
                    /* colon(a, d, b) */
                    d = b;
                    next_token();
                    if (!parse_expr(&b)) return 0;
                }
                
                if (current_token.type != TOK_RPAREN) {
                    printf("Error: expected ')'\n");
                    return 0;
                }
                next_token();
                
                /* Compute number of elements */
                {
                    apf range, n_apf;
                    apf_sub(&range, &b.re, &a.re);
                    apf_div(&n_apf, &range, &d.re);
                    n = apf_to_long(&n_apf) + 1;
                    if (n < 1) n = 0;
                    if (n > MAT_MAX_COLS) n = MAT_MAX_COLS;
                }
                
                result->type = VAL_MATRIX;
                mat_zero(&result->v.matrix, 1, n);
                
                for (i = 0; i < n; i++) {
                    apf step;
                    apf_from_int(&step, i);
                    apf_mul(&step, &step, &d.re);
                    apf_add(&MAT_AT(&result->v.matrix, 0, i).re, &a.re, &step);
                    apf_zero(&MAT_AT(&result->v.matrix, 0, i).im);
                }
                return 1;
            }
            
            /* cat(dim, A, B) - concatenate arrays */
            if (str_eq(name, "cat")) {
                value_t arg2;
                apfc dim_arg;
                int dim, i, j;
                
                next_token();
                if (!parse_expr(&dim_arg)) return 0;
                if (current_token.type != TOK_COMMA) {
                    printf("Error: cat requires (dim, A, B)\n");
                    return 0;
                }
                next_token();
                if (!parse_value(&pv_arg)) return 0;
                if (current_token.type != TOK_COMMA) {
                    printf("Error: cat requires (dim, A, B)\n");
                    return 0;
                }
                next_token();
                if (!parse_value(&arg2)) return 0;
                if (current_token.type != TOK_RPAREN) {
                    printf("Error: expected ')'\n");
                    return 0;
                }
                next_token();
                
                dim = apf_to_long(&dim_arg.re);
                
                if (pv_arg.type == VAL_SCALAR) {
                    apfc saved = pv_arg.v.scalar;
                    mat_zero(&pv_arg.v.matrix, 1, 1);
                    MAT_AT(&pv_arg.v.matrix, 0, 0) = saved;
                    pv_arg.type = VAL_MATRIX;
                }
                if (arg2.type == VAL_SCALAR) {
                    apfc saved = arg2.v.scalar;
                    mat_zero(&arg2.v.matrix, 1, 1);
                    MAT_AT(&arg2.v.matrix, 0, 0) = saved;
                    arg2.type = VAL_MATRIX;
                }
                
                result->type = VAL_MATRIX;
                
                if (dim == 1) {
                    /* Concatenate vertically (rows) */
                    if (pv_arg.v.matrix.cols != arg2.v.matrix.cols) {
                        printf("Error: column dimensions must match\n");
                        return 0;
                    }
                    mat_zero(&result->v.matrix, 
                            pv_arg.v.matrix.rows + arg2.v.matrix.rows,
                            pv_arg.v.matrix.cols);
                    for (i = 0; i < pv_arg.v.matrix.rows; i++) {
                        for (j = 0; j < pv_arg.v.matrix.cols; j++) {
                            MAT_AT(&result->v.matrix, i, j) = MAT_AT(&pv_arg.v.matrix, i, j);
                        }
                    }
                    for (i = 0; i < arg2.v.matrix.rows; i++) {
                        for (j = 0; j < arg2.v.matrix.cols; j++) {
                            MAT_AT(&result->v.matrix, pv_arg.v.matrix.rows + i, j) = MAT_AT(&arg2.v.matrix, i, j);
                        }
                    }
                } else {
                    /* Concatenate horizontally (columns) */
                    if (pv_arg.v.matrix.rows != arg2.v.matrix.rows) {
                        printf("Error: row dimensions must match\n");
                        return 0;
                    }
                    mat_zero(&result->v.matrix,
                            pv_arg.v.matrix.rows,
                            pv_arg.v.matrix.cols + arg2.v.matrix.cols);
                    for (i = 0; i < pv_arg.v.matrix.rows; i++) {
                        for (j = 0; j < pv_arg.v.matrix.cols; j++) {
                            MAT_AT(&result->v.matrix, i, j) = MAT_AT(&pv_arg.v.matrix, i, j);
                        }
                    }
                    for (i = 0; i < arg2.v.matrix.rows; i++) {
                        for (j = 0; j < arg2.v.matrix.cols; j++) {
                            MAT_AT(&result->v.matrix, i, pv_arg.v.matrix.cols + j) = MAT_AT(&arg2.v.matrix, i, j);
                        }
                    }
                }
                return 1;
            }
            
            /* linspace2(a, b, n) - linearly spaced vector (explicit 3-arg) */
            if (str_eq(name, "linspace2")) {
                apfc a, b, n_arg;
                int n, i;
                apf step;
                
                next_token();
                if (!parse_expr(&a)) return 0;
                if (current_token.type != TOK_COMMA) {
                    printf("Error: linspace2 requires (a, b, n)\n");
                    return 0;
                }
                next_token();
                if (!parse_expr(&b)) return 0;
                if (current_token.type != TOK_COMMA) {
                    printf("Error: linspace2 requires (a, b, n)\n");
                    return 0;
                }
                next_token();
                if (!parse_expr(&n_arg)) return 0;
                if (current_token.type != TOK_RPAREN) {
                    printf("Error: expected ')'\n");
                    return 0;
                }
                next_token();
                
                n = apf_to_long(&n_arg.re);
                if (n < 1) n = 1;
                if (n > MAT_MAX_COLS) n = MAT_MAX_COLS;
                
                result->type = VAL_MATRIX;
                mat_zero(&result->v.matrix, 1, n);
                
                if (n == 1) {
                    MAT_AT(&result->v.matrix, 0, 0).re = b.re;
                } else {
                    apf n_minus_1;
                    apf_from_int(&n_minus_1, n - 1);
                    apf_sub(&step, &b.re, &a.re);
                    apf_div(&step, &step, &n_minus_1);
                    
                    for (i = 0; i < n; i++) {
                        apf idx;
                        apf_from_int(&idx, i);
                        apf_mul(&idx, &idx, &step);
                        apf_add(&MAT_AT(&result->v.matrix, 0, i).re, &a.re, &idx);
                        apf_zero(&MAT_AT(&result->v.matrix, 0, i).im);
                    }
                }
                return 1;
            }
            
            /* magic(n) - n-by-n magic square */
            if (str_eq(name, "magic")) {
                int n, i, j, r, c, num;
                apfc n_arg;
                
                next_token();
                if (!parse_expr(&n_arg)) return 0;
                if (current_token.type != TOK_RPAREN) {
                    printf("Error: expected ')'\n");
                    return 0;
                }
                next_token();
                
                n = apf_to_long(&n_arg.re);
                if (n < 1 || n > MAT_MAX_ROWS) {
                    printf("Error: magic square size must be 1-%d\n", MAT_MAX_ROWS);
                    return 0;
                }
                
                result->type = VAL_MATRIX;
                mat_zero(&result->v.matrix, n, n);
                
                if (n == 1) {
                    apf_from_int(&MAT_AT(&result->v.matrix, 0, 0).re, 1);
                } else if (n == 2) {
                    /* No magic square for n=2, return [1,3;4,2] */
                    apf_from_int(&MAT_AT(&result->v.matrix, 0, 0).re, 1);
                    apf_from_int(&MAT_AT(&result->v.matrix, 0, 1).re, 3);
                    apf_from_int(&MAT_AT(&result->v.matrix, 1, 0).re, 4);
                    apf_from_int(&MAT_AT(&result->v.matrix, 1, 1).re, 2);
                } else if (n % 2 == 1) {
                    /* Siamese method for odd n */
                    int nr, nc;
                    r = 0;
                    c = n / 2;
                    for (num = 1; num <= n * n; num++) {
                        apf_from_int(&MAT_AT(&result->v.matrix, r, c).re, num);
                        nr = (r - 1 + n) % n;
                        nc = (c + 1) % n;
                        if (apf_to_long(&MAT_AT(&result->v.matrix, nr, nc).re) != 0) {
                            r = (r + 1) % n;
                        } else {
                            r = nr;
                            c = nc;
                        }
                    }
                } else if (n % 4 == 0) {
                    /* Doubly even (divisible by 4) */
                    int ii, jj;
                    long v;
                    num = 1;
                    for (i = 0; i < n; i++) {
                        for (j = 0; j < n; j++) {
                            apf_from_int(&MAT_AT(&result->v.matrix, i, j).re, num);
                            num++;
                        }
                    }
                    /* Swap on diagonals */
                    for (i = 0; i < n; i++) {
                        for (j = 0; j < n; j++) {
                            ii = i % 4;
                            jj = j % 4;
                            if ((ii == jj) || (ii + jj == 3)) {
                                v = apf_to_long(&MAT_AT(&result->v.matrix, i, j).re);
                                apf_from_int(&MAT_AT(&result->v.matrix, i, j).re, n*n + 1 - v);
                            }
                        }
                    }
                } else {
                    /* Singly even (n=6, 10, 14...) - LUX method simplified */
                    int m = n / 2;
                    int k = (n - 2) / 4;
                    /* Build quadrants */
                    for (i = 0; i < m; i++) {
                        for (j = 0; j < m; j++) {
                            /* Siamese on m-by-m */
                            /* Simplified: just fill with pattern */
                            apf_from_int(&MAT_AT(&result->v.matrix, i, j).re, 
                                        (i * m + j + 1));
                            apf_from_int(&MAT_AT(&result->v.matrix, i, j + m).re, 
                                        (i * m + j + 1) + 2*m*m);
                            apf_from_int(&MAT_AT(&result->v.matrix, i + m, j).re, 
                                        (i * m + j + 1) + 3*m*m);
                            apf_from_int(&MAT_AT(&result->v.matrix, i + m, j + m).re, 
                                        (i * m + j + 1) + m*m);
                        }
                    }
                    /* Swap columns for magic property */
                    for (i = 0; i < m; i++) {
                        for (j = 0; j < k; j++) {
                            apfc tmp = MAT_AT(&result->v.matrix, i, j);
                            MAT_AT(&result->v.matrix, i, j) = MAT_AT(&result->v.matrix, i + m, j);
                            MAT_AT(&result->v.matrix, i + m, j) = tmp;
                        }
                    }
                }
                return 1;
            }
            
            /* pascal(n) - n-by-n Pascal matrix (binomial coefficients) */
            if (str_eq(name, "pascal")) {
                int n, i, j;
                apfc n_arg;
                
                next_token();
                if (!parse_expr(&n_arg)) return 0;
                if (current_token.type != TOK_RPAREN) {
                    printf("Error: expected ')'\n");
                    return 0;
                }
                next_token();
                
                n = apf_to_long(&n_arg.re);
                if (n < 1 || n > MAT_MAX_ROWS) {
                    printf("Error: pascal size must be 1-%d\n", MAT_MAX_ROWS);
                    return 0;
                }
                
                result->type = VAL_MATRIX;
                mat_zero(&result->v.matrix, n, n);
                
                /* P(i,j) = C(i+j-2, i-1) for 1-based indexing */
                /* P(i,j) = C(i+j, i) for 0-based */
                for (i = 0; i < n; i++) {
                    apf_from_int(&MAT_AT(&result->v.matrix, i, 0).re, 1);
                    apf_from_int(&MAT_AT(&result->v.matrix, 0, i).re, 1);
                }
                for (i = 1; i < n; i++) {
                    for (j = 1; j < n; j++) {
                        apf sum;
                        apf_add(&sum, 
                               &MAT_AT(&result->v.matrix, i-1, j).re,
                               &MAT_AT(&result->v.matrix, i, j-1).re);
                        MAT_AT(&result->v.matrix, i, j).re = sum;
                    }
                }
                return 1;
            }
            
            /* hilb(n) - n-by-n Hilbert matrix */
            if (str_eq(name, "hilb")) {
                int n, i, j;
                apfc n_arg;
                
                next_token();
                if (!parse_expr(&n_arg)) return 0;
                if (current_token.type != TOK_RPAREN) {
                    printf("Error: expected ')'\n");
                    return 0;
                }
                next_token();
                
                n = apf_to_long(&n_arg.re);
                if (n < 1 || n > MAT_MAX_ROWS) {
                    printf("Error: hilb size must be 1-%d\n", MAT_MAX_ROWS);
                    return 0;
                }
                
                result->type = VAL_MATRIX;
                mat_zero(&result->v.matrix, n, n);
                
                /* H(i,j) = 1 / (i + j - 1) for 1-based = 1/(i+j+1) for 0-based */
                for (i = 0; i < n; i++) {
                    for (j = 0; j < n; j++) {
                        apf num, denom;
                        apf_from_int(&num, 1);
                        apf_from_int(&denom, i + j + 1);
                        apf_div(&MAT_AT(&result->v.matrix, i, j).re, &num, &denom);
                    }
                }
                return 1;
            }
            
            /* compan(p) - companion matrix for polynomial p */
            if (str_eq(name, "compan")) {
                int n, i;
                apf leading;
                
                next_token();
                if (!parse_value(&pv_arg)) return 0;
                if (current_token.type != TOK_RPAREN) {
                    printf("Error: expected ')'\n");
                    return 0;
                }
                next_token();
                
                if (pv_arg.type == VAL_SCALAR) {
                    result->type = VAL_MATRIX;
                    mat_zero(&result->v.matrix, 0, 0);
                    return 1;
                }
                
                n = pv_arg.v.matrix.rows * pv_arg.v.matrix.cols - 1;
                if (n < 1) {
                    result->type = VAL_MATRIX;
                    mat_zero(&result->v.matrix, 0, 0);
                    return 1;
                }
                
                result->type = VAL_MATRIX;
                mat_zero(&result->v.matrix, n, n);
                
                /* Get leading coefficient */
                apf_copy(&leading, &MAT_AT(&pv_arg.v.matrix, 0, 0).re);
                
                /* First row: -p(2:end) / p(1) */
                for (i = 0; i < n; i++) {
                    int ri = (i + 1) % pv_arg.v.matrix.rows;
                    int ci = (i + 1) / pv_arg.v.matrix.rows;
                    apf_div(&MAT_AT(&result->v.matrix, 0, i).re,
                           &MAT_AT(&pv_arg.v.matrix, ri, ci).re, &leading);
                    apf_neg(&MAT_AT(&result->v.matrix, 0, i).re, &MAT_AT(&result->v.matrix, 0, i).re);
                }
                
                /* Sub-diagonal of ones */
                for (i = 0; i < n - 1; i++) {
                    apf_from_int(&MAT_AT(&result->v.matrix, i + 1, i).re, 1);
                }
                
                return 1;
            }
            
            /* ===== SESSION 3: 40 NEW FUNCTIONS ===== */
            
            /* === Coordinate Conversions === */
            
            /* cart2pol(x, y) - Cartesian to polar [theta, r] */
            if (str_eq(name, "cart2pol")) {
                apfc x_arg, y_arg;
                apf theta, r, x2, y2, sum;
                
                next_token();
                if (!parse_expr(&x_arg)) return 0;
                if (current_token.type != TOK_COMMA) {
                    printf("Error: cart2pol requires (x, y)\n");
                    return 0;
                }
                next_token();
                if (!parse_expr(&y_arg)) return 0;
                if (current_token.type != TOK_RPAREN) {
                    printf("Error: expected ')'\n");
                    return 0;
                }
                next_token();
                
                apfx_atan2(&theta, &y_arg.re, &x_arg.re);
                apf_mul(&x2, &x_arg.re, &x_arg.re);
                apf_mul(&y2, &y_arg.re, &y_arg.re);
                apf_add(&sum, &x2, &y2);
                apf_sqrt(&r, &sum);
                
                result->type = VAL_MATRIX;
                mat_zero(&result->v.matrix, 1, 2);
                MAT_AT(&result->v.matrix, 0, 0).re = theta;
                MAT_AT(&result->v.matrix, 0, 1).re = r;
                return 1;
            }
            
            /* pol2cart(theta, r) - Polar to Cartesian [x, y] */
            if (str_eq(name, "pol2cart")) {
                apfc theta_arg, r_arg;
                apf x, y, sin_t, cos_t;
                
                next_token();
                if (!parse_expr(&theta_arg)) return 0;
                if (current_token.type != TOK_COMMA) {
                    printf("Error: pol2cart requires (theta, r)\n");
                    return 0;
                }
                next_token();
                if (!parse_expr(&r_arg)) return 0;
                if (current_token.type != TOK_RPAREN) {
                    printf("Error: expected ')'\n");
                    return 0;
                }
                next_token();
                
                apfx_sincos(&sin_t, &cos_t, &theta_arg.re);
                apf_mul(&x, &r_arg.re, &cos_t);
                apf_mul(&y, &r_arg.re, &sin_t);
                
                result->type = VAL_MATRIX;
                mat_zero(&result->v.matrix, 1, 2);
                MAT_AT(&result->v.matrix, 0, 0).re = x;
                MAT_AT(&result->v.matrix, 0, 1).re = y;
                return 1;
            }
            
            /* cart2sph(x, y, z) - Cartesian to spherical [azimuth, elevation, r] */
            if (str_eq(name, "cart2sph")) {
                apfc x_arg, y_arg, z_arg;
                apf azimuth, elevation, r, x2, y2, z2, hypotxy, sum;
                
                next_token();
                if (!parse_expr(&x_arg)) return 0;
                if (current_token.type != TOK_COMMA) {
                    printf("Error: cart2sph requires (x, y, z)\n");
                    return 0;
                }
                next_token();
                if (!parse_expr(&y_arg)) return 0;
                if (current_token.type != TOK_COMMA) {
                    printf("Error: cart2sph requires (x, y, z)\n");
                    return 0;
                }
                next_token();
                if (!parse_expr(&z_arg)) return 0;
                if (current_token.type != TOK_RPAREN) {
                    printf("Error: expected ')'\n");
                    return 0;
                }
                next_token();
                
                /* azimuth = atan2(y, x) */
                apfx_atan2(&azimuth, &y_arg.re, &x_arg.re);
                
                /* hypotxy = sqrt(x^2 + y^2) */
                apf_mul(&x2, &x_arg.re, &x_arg.re);
                apf_mul(&y2, &y_arg.re, &y_arg.re);
                apf_add(&hypotxy, &x2, &y2);
                apf_sqrt(&hypotxy, &hypotxy);
                
                /* elevation = atan2(z, hypotxy) */
                apfx_atan2(&elevation, &z_arg.re, &hypotxy);
                
                /* r = sqrt(x^2 + y^2 + z^2) */
                apf_mul(&z2, &z_arg.re, &z_arg.re);
                apf_add(&sum, &x2, &y2);
                apf_add(&sum, &sum, &z2);
                apf_sqrt(&r, &sum);
                
                result->type = VAL_MATRIX;
                mat_zero(&result->v.matrix, 1, 3);
                MAT_AT(&result->v.matrix, 0, 0).re = azimuth;
                MAT_AT(&result->v.matrix, 0, 1).re = elevation;
                MAT_AT(&result->v.matrix, 0, 2).re = r;
                return 1;
            }
            
            /* sph2cart(azimuth, elevation, r) - Spherical to Cartesian [x, y, z] */
            if (str_eq(name, "sph2cart")) {
                apfc az_arg, el_arg, r_arg;
                apf x, y, z, sin_az, cos_az, sin_el, cos_el, rcosel;
                
                next_token();
                if (!parse_expr(&az_arg)) return 0;
                if (current_token.type != TOK_COMMA) {
                    printf("Error: sph2cart requires (azimuth, elevation, r)\n");
                    return 0;
                }
                next_token();
                if (!parse_expr(&el_arg)) return 0;
                if (current_token.type != TOK_COMMA) {
                    printf("Error: sph2cart requires (azimuth, elevation, r)\n");
                    return 0;
                }
                next_token();
                if (!parse_expr(&r_arg)) return 0;
                if (current_token.type != TOK_RPAREN) {
                    printf("Error: expected ')'\n");
                    return 0;
                }
                next_token();
                
                apfx_sincos(&sin_az, &cos_az, &az_arg.re);
                apfx_sincos(&sin_el, &cos_el, &el_arg.re);
                
                /* x = r * cos(el) * cos(az) */
                apf_mul(&rcosel, &r_arg.re, &cos_el);
                apf_mul(&x, &rcosel, &cos_az);
                
                /* y = r * cos(el) * sin(az) */
                apf_mul(&y, &rcosel, &sin_az);
                
                /* z = r * sin(el) */
                apf_mul(&z, &r_arg.re, &sin_el);
                
                result->type = VAL_MATRIX;
                mat_zero(&result->v.matrix, 1, 3);
                MAT_AT(&result->v.matrix, 0, 0).re = x;
                MAT_AT(&result->v.matrix, 0, 1).re = y;
                MAT_AT(&result->v.matrix, 0, 2).re = z;
                return 1;
            }
            
            /* === Array Manipulation === */
            
            /* fliplr(M) - Flip left-right */
            if (str_eq(name, "fliplr")) {
                int i, j;
                
                next_token();
                if (!parse_value(&pv_arg)) return 0;
                if (current_token.type != TOK_RPAREN) {
                    printf("Error: expected ')'\n");
                    return 0;
                }
                next_token();
                
                if (pv_arg.type == VAL_SCALAR) {
                    result->type = VAL_SCALAR;
                    result->v.scalar = pv_arg.v.scalar;
                    return 1;
                }
                
                result->type = VAL_MATRIX;
                mat_zero(&result->v.matrix, pv_arg.v.matrix.rows, pv_arg.v.matrix.cols);
                for (i = 0; i < pv_arg.v.matrix.rows; i++) {
                    for (j = 0; j < pv_arg.v.matrix.cols; j++) {
                        MAT_AT(&result->v.matrix, i, pv_arg.v.matrix.cols - 1 - j) =
                            MAT_AT(&pv_arg.v.matrix, i, j);
                    }
                }
                return 1;
            }
            
            /* flipud(M) - Flip up-down */
            if (str_eq(name, "flipud")) {
                int i, j;
                
                next_token();
                if (!parse_value(&pv_arg)) return 0;
                if (current_token.type != TOK_RPAREN) {
                    printf("Error: expected ')'\n");
                    return 0;
                }
                next_token();
                
                if (pv_arg.type == VAL_SCALAR) {
                    result->type = VAL_SCALAR;
                    result->v.scalar = pv_arg.v.scalar;
                    return 1;
                }
                
                result->type = VAL_MATRIX;
                mat_zero(&result->v.matrix, pv_arg.v.matrix.rows, pv_arg.v.matrix.cols);
                for (i = 0; i < pv_arg.v.matrix.rows; i++) {
                    for (j = 0; j < pv_arg.v.matrix.cols; j++) {
                        MAT_AT(&result->v.matrix, pv_arg.v.matrix.rows - 1 - i, j) =
                            MAT_AT(&pv_arg.v.matrix, i, j);
                    }
                }
                return 1;
            }
            
            /* circshift(v, k) - Circular shift */
            if (str_eq(name, "circshift")) {
                apfc k_arg;
                int k, n, i;
                
                next_token();
                if (!parse_value(&pv_arg)) return 0;
                if (current_token.type != TOK_COMMA) {
                    printf("Error: circshift requires (v, k)\n");
                    return 0;
                }
                next_token();
                if (!parse_expr(&k_arg)) return 0;
                if (current_token.type != TOK_RPAREN) {
                    printf("Error: expected ')'\n");
                    return 0;
                }
                next_token();
                
                if (pv_arg.type == VAL_SCALAR) {
                    result->type = VAL_SCALAR;
                    result->v.scalar = pv_arg.v.scalar;
                    return 1;
                }
                
                n = pv_arg.v.matrix.rows * pv_arg.v.matrix.cols;
                k = apf_to_long(&k_arg.re);
                k = ((k % n) + n) % n;
                
                result->type = VAL_MATRIX;
                mat_zero(&result->v.matrix, pv_arg.v.matrix.rows, pv_arg.v.matrix.cols);
                for (i = 0; i < n; i++) {
                    int src_r = i % pv_arg.v.matrix.rows;
                    int src_c = i / pv_arg.v.matrix.rows;
                    int dst_i = (i + k) % n;
                    int dst_r = dst_i % result->v.matrix.rows;
                    int dst_c = dst_i / result->v.matrix.rows;
                    MAT_AT(&result->v.matrix, dst_r, dst_c) = MAT_AT(&pv_arg.v.matrix, src_r, src_c);
                }
                return 1;
            }
            
            /* diff(v) - Differences (discrete derivative) */
            if (str_eq(name, "diff")) {
                int i, n;
                
                next_token();
                if (!parse_value(&pv_arg)) return 0;
                if (current_token.type != TOK_RPAREN) {
                    printf("Error: expected ')'\n");
                    return 0;
                }
                next_token();
                
                if (pv_arg.type == VAL_SCALAR) {
                    result->type = VAL_MATRIX;
                    mat_zero(&result->v.matrix, 0, 0);
                    return 1;
                }
                
                n = pv_arg.v.matrix.rows * pv_arg.v.matrix.cols;
                if (n < 2) {
                    result->type = VAL_MATRIX;
                    mat_zero(&result->v.matrix, 0, 0);
                    return 1;
                }
                
                result->type = VAL_MATRIX;
                mat_zero(&result->v.matrix, 1, n - 1);
                for (i = 0; i < n - 1; i++) {
                    int r1 = i % pv_arg.v.matrix.rows;
                    int c1 = i / pv_arg.v.matrix.rows;
                    int r2 = (i + 1) % pv_arg.v.matrix.rows;
                    int c2 = (i + 1) / pv_arg.v.matrix.rows;
                    apf_sub(&MAT_AT(&result->v.matrix, 0, i).re,
                           &MAT_AT(&pv_arg.v.matrix, r2, c2).re,
                           &MAT_AT(&pv_arg.v.matrix, r1, c1).re);
                }
                return 1;
            }
            
            /* === HP/TI Calculator Functions === */
            
            /* frac(x) - Fractional part */
            if (str_eq(name, "frac")) {
                apf trunc_val;
                int i, j;
                
                next_token();
                if (!parse_value(&pv_arg)) return 0;
                if (current_token.type != TOK_RPAREN) {
                    printf("Error: expected ')'\n");
                    return 0;
                }
                next_token();
                
                if (pv_arg.type == VAL_SCALAR) {
                    result->type = VAL_SCALAR;
                    apf_trunc(&trunc_val, &pv_arg.v.scalar.re);
                    apf_sub(&result->v.scalar.re, &pv_arg.v.scalar.re, &trunc_val);
                    apf_zero(&result->v.scalar.im);
                } else {
                    result->type = VAL_MATRIX;
                    mat_zero(&result->v.matrix, pv_arg.v.matrix.rows, pv_arg.v.matrix.cols);
                    for (i = 0; i < pv_arg.v.matrix.rows; i++) {
                        for (j = 0; j < pv_arg.v.matrix.cols; j++) {
                            apf_trunc(&trunc_val, &MAT_AT(&pv_arg.v.matrix, i, j).re);
                            apf_sub(&MAT_AT(&result->v.matrix, i, j).re,
                                   &MAT_AT(&pv_arg.v.matrix, i, j).re, &trunc_val);
                        }
                    }
                }
                return 1;
            }
            
            /* dms(d, m, s) - DMS to decimal degrees */
            if (str_eq(name, "dms") || str_eq(name, "dms2deg")) {
                apfc d_arg, m_arg, s_arg;
                apf sixty, temp;
                
                next_token();
                if (!parse_expr(&d_arg)) return 0;
                if (current_token.type != TOK_COMMA) {
                    printf("Error: dms requires (d, m, s)\n");
                    return 0;
                }
                next_token();
                if (!parse_expr(&m_arg)) return 0;
                if (current_token.type != TOK_COMMA) {
                    printf("Error: dms requires (d, m, s)\n");
                    return 0;
                }
                next_token();
                if (!parse_expr(&s_arg)) return 0;
                if (current_token.type != TOK_RPAREN) {
                    printf("Error: expected ')'\n");
                    return 0;
                }
                next_token();
                
                apf_from_int(&sixty, 60);
                apf_div(&temp, &s_arg.re, &sixty);
                apf_add(&temp, &temp, &m_arg.re);
                apf_div(&temp, &temp, &sixty);
                apf_add(&result->v.scalar.re, &temp, &d_arg.re);
                apf_zero(&result->v.scalar.im);
                result->type = VAL_SCALAR;
                return 1;
            }
            
            /* hms(h, m, s) - HMS to decimal hours */
            if (str_eq(name, "hms") || str_eq(name, "hms2hr")) {
                apfc h_arg, m_arg, s_arg;
                apf sixty, temp;
                
                next_token();
                if (!parse_expr(&h_arg)) return 0;
                if (current_token.type != TOK_COMMA) {
                    printf("Error: hms requires (h, m, s)\n");
                    return 0;
                }
                next_token();
                if (!parse_expr(&m_arg)) return 0;
                if (current_token.type != TOK_COMMA) {
                    printf("Error: hms requires (h, m, s)\n");
                    return 0;
                }
                next_token();
                if (!parse_expr(&s_arg)) return 0;
                if (current_token.type != TOK_RPAREN) {
                    printf("Error: expected ')'\n");
                    return 0;
                }
                next_token();
                
                apf_from_int(&sixty, 60);
                apf_div(&temp, &s_arg.re, &sixty);
                apf_add(&temp, &temp, &m_arg.re);
                apf_div(&temp, &temp, &sixty);
                apf_add(&result->v.scalar.re, &temp, &h_arg.re);
                apf_zero(&result->v.scalar.im);
                result->type = VAL_SCALAR;
                return 1;
            }
            
            /* deg2dms(d) - Decimal degrees to [d, m, s] */
            if (str_eq(name, "deg2dms")) {
                apf deg, min_val, sec, sixty, abs_deg, frac_val;
                int sign;
                
                next_token();
                if (!parse_expr(&pv_arg.v.scalar)) return 0;
                if (current_token.type != TOK_RPAREN) {
                    printf("Error: expected ')'\n");
                    return 0;
                }
                next_token();
                
                sign = pv_arg.v.scalar.re.sign;
                apf_abs(&abs_deg, &pv_arg.v.scalar.re);
                apf_from_int(&sixty, 60);
                
                apf_trunc(&deg, &abs_deg);
                apf_sub(&frac_val, &abs_deg, &deg);
                apf_mul(&frac_val, &frac_val, &sixty);
                apf_trunc(&min_val, &frac_val);
                apf_sub(&frac_val, &frac_val, &min_val);
                apf_mul(&sec, &frac_val, &sixty);
                
                if (sign) apf_neg(&deg, &deg);
                
                result->type = VAL_MATRIX;
                mat_zero(&result->v.matrix, 1, 3);
                MAT_AT(&result->v.matrix, 0, 0).re = deg;
                MAT_AT(&result->v.matrix, 0, 1).re = min_val;
                MAT_AT(&result->v.matrix, 0, 2).re = sec;
                return 1;
            }
            
            /* hr2hms(h) - Decimal hours to [h, m, s] */
            if (str_eq(name, "hr2hms")) {
                apf hr, min_val, sec, sixty, abs_hr, frac_val;
                int sign;
                
                next_token();
                if (!parse_expr(&pv_arg.v.scalar)) return 0;
                if (current_token.type != TOK_RPAREN) {
                    printf("Error: expected ')'\n");
                    return 0;
                }
                next_token();
                
                sign = pv_arg.v.scalar.re.sign;
                apf_abs(&abs_hr, &pv_arg.v.scalar.re);
                apf_from_int(&sixty, 60);
                
                apf_trunc(&hr, &abs_hr);
                apf_sub(&frac_val, &abs_hr, &hr);
                apf_mul(&frac_val, &frac_val, &sixty);
                apf_trunc(&min_val, &frac_val);
                apf_sub(&frac_val, &frac_val, &min_val);
                apf_mul(&sec, &frac_val, &sixty);
                
                if (sign) apf_neg(&hr, &hr);
                
                result->type = VAL_MATRIX;
                mat_zero(&result->v.matrix, 1, 3);
                MAT_AT(&result->v.matrix, 0, 0).re = hr;
                MAT_AT(&result->v.matrix, 0, 1).re = min_val;
                MAT_AT(&result->v.matrix, 0, 2).re = sec;
                return 1;
            }
            
            /* percent(x, p) - p% of x */
            if (str_eq(name, "percent") || str_eq(name, "pct")) {
                apfc x_arg, p_arg;
                apf hundred;
                
                next_token();
                if (!parse_expr(&x_arg)) return 0;
                if (current_token.type != TOK_COMMA) {
                    printf("Error: percent requires (x, p)\n");
                    return 0;
                }
                next_token();
                if (!parse_expr(&p_arg)) return 0;
                if (current_token.type != TOK_RPAREN) {
                    printf("Error: expected ')'\n");
                    return 0;
                }
                next_token();
                
                apf_from_int(&hundred, 100);
                apf_mul(&result->v.scalar.re, &x_arg.re, &p_arg.re);
                apf_div(&result->v.scalar.re, &result->v.scalar.re, &hundred);
                apf_zero(&result->v.scalar.im);
                result->type = VAL_SCALAR;
                return 1;
            }
            
            /* pctchange(old, new) - Percent change */
            if (str_eq(name, "pctchange") || str_eq(name, "percentchange")) {
                apfc old_arg, new_arg;
                apf diff, hundred;
                
                next_token();
                if (!parse_expr(&old_arg)) return 0;
                if (current_token.type != TOK_COMMA) {
                    printf("Error: pctchange requires (old, new)\n");
                    return 0;
                }
                next_token();
                if (!parse_expr(&new_arg)) return 0;
                if (current_token.type != TOK_RPAREN) {
                    printf("Error: expected ')'\n");
                    return 0;
                }
                next_token();
                
                apf_sub(&diff, &new_arg.re, &old_arg.re);
                apf_div(&diff, &diff, &old_arg.re);
                apf_from_int(&hundred, 100);
                apf_mul(&result->v.scalar.re, &diff, &hundred);
                apf_zero(&result->v.scalar.im);
                result->type = VAL_SCALAR;
                return 1;
            }
            
            /* === Utility Functions === */
            
            /* lerp(a, b, t) - Linear interpolation */
            if (str_eq(name, "lerp") || str_eq(name, "mix")) {
                apfc a_arg, b_arg, t_arg;
                apf diff, scaled;
                
                next_token();
                if (!parse_expr(&a_arg)) return 0;
                if (current_token.type != TOK_COMMA) {
                    printf("Error: lerp requires (a, b, t)\n");
                    return 0;
                }
                next_token();
                if (!parse_expr(&b_arg)) return 0;
                if (current_token.type != TOK_COMMA) {
                    printf("Error: lerp requires (a, b, t)\n");
                    return 0;
                }
                next_token();
                if (!parse_expr(&t_arg)) return 0;
                if (current_token.type != TOK_RPAREN) {
                    printf("Error: expected ')'\n");
                    return 0;
                }
                next_token();
                
                apf_sub(&diff, &b_arg.re, &a_arg.re);
                apf_mul(&scaled, &diff, &t_arg.re);
                apf_add(&result->v.scalar.re, &a_arg.re, &scaled);
                apf_zero(&result->v.scalar.im);
                result->type = VAL_SCALAR;
                return 1;
            }
            
            /* wrap(x, lo, hi) - Wrap to range [lo, hi) */
            if (str_eq(name, "wrap")) {
                apfc lo, hi;
                apf range, offset, n_periods, floor_n;
                
                next_token();
                if (!parse_expr(&pv_arg.v.scalar)) return 0;
                if (current_token.type != TOK_COMMA) {
                    printf("Error: wrap requires (x, lo, hi)\n");
                    return 0;
                }
                next_token();
                if (!parse_expr(&lo)) return 0;
                if (current_token.type != TOK_COMMA) {
                    printf("Error: wrap requires (x, lo, hi)\n");
                    return 0;
                }
                next_token();
                if (!parse_expr(&hi)) return 0;
                if (current_token.type != TOK_RPAREN) {
                    printf("Error: expected ')'\n");
                    return 0;
                }
                next_token();
                
                apf_sub(&range, &hi.re, &lo.re);
                apf_sub(&offset, &pv_arg.v.scalar.re, &lo.re);
                apf_div(&n_periods, &offset, &range);
                apf_floor(&floor_n, &n_periods);
                apf_mul(&floor_n, &floor_n, &range);
                apf_sub(&result->v.scalar.re, &offset, &floor_n);
                apf_add(&result->v.scalar.re, &result->v.scalar.re, &lo.re);
                apf_zero(&result->v.scalar.im);
                result->type = VAL_SCALAR;
                return 1;
            }
            
            /* === Number Theory === */
            
            /* totient(n) / eulerphi(n) - Euler's totient */
            if (str_eq(name, "totient") || str_eq(name, "eulerphi")) {
                long n, res, temp_n, p;
                
                next_token();
                if (!parse_expr(&pv_arg.v.scalar)) return 0;
                if (current_token.type != TOK_RPAREN) {
                    printf("Error: expected ')'\n");
                    return 0;
                }
                next_token();
                
                n = apf_to_long(&pv_arg.v.scalar.re);
                if (n < 1) n = 1;
                
                res = n;
                temp_n = n;
                
                for (p = 2; p * p <= temp_n; p++) {
                    if (temp_n % p == 0) {
                        while (temp_n % p == 0) temp_n /= p;
                        res -= res / p;
                    }
                }
                if (temp_n > 1) res -= res / temp_n;
                
                apf_from_int(&result->v.scalar.re, res);
                apf_zero(&result->v.scalar.im);
                result->type = VAL_SCALAR;
                return 1;
            }
            
            /* divisors(n) - All divisors as row vector */
            if (str_eq(name, "divisors")) {
                long n, i, count;
                long divs[100];
                int j;
                
                next_token();
                if (!parse_expr(&pv_arg.v.scalar)) return 0;
                if (current_token.type != TOK_RPAREN) {
                    printf("Error: expected ')'\n");
                    return 0;
                }
                next_token();
                
                n = apf_to_long(&pv_arg.v.scalar.re);
                if (n < 1) n = 1;
                
                count = 0;
                for (i = 1; i * i <= n && count < 100; i++) {
                    if (n % i == 0) {
                        divs[count++] = i;
                        if (i != n / i && count < 100) divs[count++] = n / i;
                    }
                }
                
                /* Sort */
                for (i = 0; i < count - 1; i++) {
                    for (j = (int)i + 1; j < count; j++) {
                        if (divs[j] < divs[i]) {
                            long tmp = divs[i]; divs[i] = divs[j]; divs[j] = tmp;
                        }
                    }
                }
                
                result->type = VAL_MATRIX;
                if (count > MAT_MAX_COLS) count = MAT_MAX_COLS;
                mat_zero(&result->v.matrix, 1, (int)count);
                for (i = 0; i < count; i++) {
                    apf_from_int(&MAT_AT(&result->v.matrix, 0, (int)i).re, divs[i]);
                }
                return 1;
            }
            
            /* mobius(n) - Möbius function */
            if (str_eq(name, "mobius") || str_eq(name, "moebius")) {
                long n, temp_n, p, num_factors, res;
                
                next_token();
                if (!parse_expr(&pv_arg.v.scalar)) return 0;
                if (current_token.type != TOK_RPAREN) {
                    printf("Error: expected ')'\n");
                    return 0;
                }
                next_token();
                
                n = apf_to_long(&pv_arg.v.scalar.re);
                if (n < 1) {
                    apf_zero(&result->v.scalar.re);
                    apf_zero(&result->v.scalar.im);
                    result->type = VAL_SCALAR;
                    return 1;
                }
                if (n == 1) {
                    apf_from_int(&result->v.scalar.re, 1);
                    apf_zero(&result->v.scalar.im);
                    result->type = VAL_SCALAR;
                    return 1;
                }
                
                temp_n = n;
                num_factors = 0;
                
                for (p = 2; p * p <= temp_n; p++) {
                    if (temp_n % p == 0) {
                        temp_n /= p;
                        num_factors++;
                        if (temp_n % p == 0) {
                            apf_zero(&result->v.scalar.re);
                            apf_zero(&result->v.scalar.im);
                            result->type = VAL_SCALAR;
                            return 1;
                        }
                    }
                }
                if (temp_n > 1) num_factors++;
                
                res = (num_factors % 2 == 0) ? 1 : -1;
                apf_from_int(&result->v.scalar.re, res);
                apf_zero(&result->v.scalar.im);
                result->type = VAL_SCALAR;
                return 1;
            }
            
            /* === Statistical Distributions === */
            
            /* normpdf(x) or normpdf(x, mu, sigma) */
            if (str_eq(name, "normpdf") || str_eq(name, "normalpdf")) {
                apfc x_arg;
                apf mu, sigma, z, z2, exp_term, coef, sqrt_2pi, half;
                
                next_token();
                if (!parse_expr(&x_arg)) return 0;
                
                apf_zero(&mu);
                apf_from_int(&sigma, 1);
                
                if (current_token.type == TOK_COMMA) {
                    next_token();
                    if (!parse_expr(&pv_arg.v.scalar)) return 0;
                    mu = pv_arg.v.scalar.re;
                    if (current_token.type == TOK_COMMA) {
                        next_token();
                        if (!parse_expr(&pv_arg.v.scalar)) return 0;
                        sigma = pv_arg.v.scalar.re;
                    }
                }
                
                if (current_token.type != TOK_RPAREN) {
                    printf("Error: expected ')'\n");
                    return 0;
                }
                next_token();
                
                apf_sub(&z, &x_arg.re, &mu);
                apf_div(&z, &z, &sigma);
                apf_mul(&z2, &z, &z);
                apf_from_str(&half, "0.5");
                apf_mul(&z2, &z2, &half);
                apf_neg(&z2, &z2);
                apfx_exp(&exp_term, &z2);
                
                apf_from_str(&sqrt_2pi, "2.506628274631");
                apf_mul(&coef, &sigma, &sqrt_2pi);
                apf_div(&result->v.scalar.re, &exp_term, &coef);
                apf_zero(&result->v.scalar.im);
                result->type = VAL_SCALAR;
                return 1;
            }
            
            /* normcdf(x) or normcdf(x, mu, sigma) */
            if (str_eq(name, "normcdf") || str_eq(name, "normalcdf")) {
                apfc x_arg;
                apf mu, sigma, z, erf_val, one, sqrt2, half;
                
                next_token();
                if (!parse_expr(&x_arg)) return 0;
                
                apf_zero(&mu);
                apf_from_int(&sigma, 1);
                
                if (current_token.type == TOK_COMMA) {
                    next_token();
                    if (!parse_expr(&pv_arg.v.scalar)) return 0;
                    mu = pv_arg.v.scalar.re;
                    if (current_token.type == TOK_COMMA) {
                        next_token();
                        if (!parse_expr(&pv_arg.v.scalar)) return 0;
                        sigma = pv_arg.v.scalar.re;
                    }
                }
                
                if (current_token.type != TOK_RPAREN) {
                    printf("Error: expected ')'\n");
                    return 0;
                }
                next_token();
                
                apf_sub(&z, &x_arg.re, &mu);
                apf_from_str(&sqrt2, "1.41421356237");
                apf_mul(&sqrt2, &sqrt2, &sigma);
                apf_div(&z, &z, &sqrt2);
                apfx_erf(&erf_val, &z);
                
                apf_from_int(&one, 1);
                apf_add(&erf_val, &one, &erf_val);
                apf_from_str(&half, "0.5");
                apf_mul(&result->v.scalar.re, &erf_val, &half);
                apf_zero(&result->v.scalar.im);
                result->type = VAL_SCALAR;
                return 1;
            }
            
            /* norminv(p) - inverse standard normal CDF (quantile function) */
            /* Uses Beasley-Springer-Moro approximation */
            if (str_eq(name, "norminv") || str_eq(name, "normalinv") || str_eq(name, "probit")) {
                apf p, t, c0, c1, c2, d1, d2, d3, num, den, half, one, zero_apf, one_apf;
                int p_leq_half;
                
                next_token();
                if (!parse_expr(&pv_arg.v.scalar)) return 0;
                if (current_token.type != TOK_RPAREN) {
                    printf("Error: expected ')'\n");
                    return 0;
                }
                next_token();
                
                apf_copy(&p, &pv_arg.v.scalar.re);
                apf_zero(&zero_apf);
                apf_from_int(&one_apf, 1);
                
                /* Check bounds: p <= 0 or p >= 1 */
                if (apf_cmp(&p, &zero_apf) <= 0) {
                    apf_set_inf(&result->v.scalar.re, 1);  /* -Inf */
                    apf_zero(&result->v.scalar.im);
                    result->type = VAL_SCALAR;
                    return 1;
                }
                if (apf_cmp(&p, &one_apf) >= 0) {
                    apf_set_inf(&result->v.scalar.re, 0);  /* +Inf */
                    apf_zero(&result->v.scalar.im);
                    result->type = VAL_SCALAR;
                    return 1;
                }
                
                /* Rational approximation for central region */
                apf_from_str(&half, "0.5");
                apf_from_int(&one, 1);
                
                p_leq_half = (apf_cmp(&p, &half) <= 0);
                if (!p_leq_half) {
                    apf_sub(&p, &one, &p);
                }
                
                /* t = sqrt(-2*ln(p)) for p < 0.5 */
                {
                    apf ln_p, neg2, two;
                    apf_from_int(&two, 2);
                    apfx_log(&ln_p, &p);
                    apf_neg(&neg2, &two);
                    apf_mul(&t, &neg2, &ln_p);
                    apf_sqrt(&t, &t);
                }
                
                /* Approximation coefficients (Hastings) */
                apf_from_str(&c0, "2.515517");
                apf_from_str(&c1, "0.802853");
                apf_from_str(&c2, "0.010328");
                apf_from_str(&d1, "1.432788");
                apf_from_str(&d2, "0.189269");
                apf_from_str(&d3, "0.001308");
                
                /* num = c0 + c1*t + c2*t^2 */
                {
                    apf t2, tmp1, tmp2;
                    apf_mul(&t2, &t, &t);
                    apf_mul(&tmp1, &c1, &t);
                    apf_mul(&tmp2, &c2, &t2);
                    apf_add(&num, &c0, &tmp1);
                    apf_add(&num, &num, &tmp2);
                }
                
                /* den = 1 + d1*t + d2*t^2 + d3*t^3 */
                {
                    apf t2, t3, tmp1, tmp2, tmp3;
                    apf_mul(&t2, &t, &t);
                    apf_mul(&t3, &t2, &t);
                    apf_mul(&tmp1, &d1, &t);
                    apf_mul(&tmp2, &d2, &t2);
                    apf_mul(&tmp3, &d3, &t3);
                    apf_add(&den, &one, &tmp1);
                    apf_add(&den, &den, &tmp2);
                    apf_add(&den, &den, &tmp3);
                }
                
                /* result = t - num/den */
                apf_div(&num, &num, &den);
                apf_sub(&result->v.scalar.re, &t, &num);
                
                /* Negate if original p <= 0.5 */
                if (p_leq_half) {
                    apf_neg(&result->v.scalar.re, &result->v.scalar.re);
                }
                
                apf_zero(&result->v.scalar.im);
                result->type = VAL_SCALAR;
                return 1;
            }
            
            /* === Trig Extensions === */
            
            /* sinpi(x) - sin(pi*x), exact at integers */
            if (str_eq(name, "sinpi")) {
                apf pi_x, trunc_x, frac_x;
                
                next_token();
                if (!parse_expr(&pv_arg.v.scalar)) return 0;
                if (current_token.type != TOK_RPAREN) {
                    printf("Error: expected ')'\n");
                    return 0;
                }
                next_token();
                
                apf_trunc(&trunc_x, &pv_arg.v.scalar.re);
                apf_sub(&frac_x, &pv_arg.v.scalar.re, &trunc_x);
                
                if (apf_is_zero(&frac_x)) {
                    apf_zero(&result->v.scalar.re);
                } else {
                    apfx_pi(&pi_x);
                    apf_mul(&pi_x, &pi_x, &pv_arg.v.scalar.re);
                    apfx_sin(&result->v.scalar.re, &pi_x);
                }
                apf_zero(&result->v.scalar.im);
                result->type = VAL_SCALAR;
                return 1;
            }
            
            /* cospi(x) - cos(pi*x), exact at half-integers */
            if (str_eq(name, "cospi")) {
                apf pi_x, frac_x, half, temp;
                
                next_token();
                if (!parse_expr(&pv_arg.v.scalar)) return 0;
                if (current_token.type != TOK_RPAREN) {
                    printf("Error: expected ')'\n");
                    return 0;
                }
                next_token();
                
                apf_from_str(&half, "0.5");
                apf_sub(&frac_x, &pv_arg.v.scalar.re, &half);
                apf_trunc(&temp, &frac_x);
                apf_sub(&frac_x, &frac_x, &temp);
                
                if (apf_is_zero(&frac_x)) {
                    apf_zero(&result->v.scalar.re);
                } else {
                    apfx_pi(&pi_x);
                    apf_mul(&pi_x, &pi_x, &pv_arg.v.scalar.re);
                    apfx_cos(&result->v.scalar.re, &pi_x);
                }
                apf_zero(&result->v.scalar.im);
                result->type = VAL_SCALAR;
                return 1;
            }
            
            /* === Logical Operations === */
            
            /* nand(a, b) - NOT AND */
            if (str_eq(name, "nand")) {
                apfc a_arg, b_arg;
                int a_val, b_val;
                
                next_token();
                if (!parse_expr(&a_arg)) return 0;
                if (current_token.type != TOK_COMMA) {
                    printf("Error: nand requires (a, b)\n");
                    return 0;
                }
                next_token();
                if (!parse_expr(&b_arg)) return 0;
                if (current_token.type != TOK_RPAREN) {
                    printf("Error: expected ')'\n");
                    return 0;
                }
                next_token();
                
                a_val = !apf_is_zero(&a_arg.re);
                b_val = !apf_is_zero(&b_arg.re);
                apf_from_int(&result->v.scalar.re, !(a_val && b_val));
                apf_zero(&result->v.scalar.im);
                result->type = VAL_SCALAR;
                return 1;
            }
            
            /* nor(a, b) - NOT OR */
            if (str_eq(name, "nor")) {
                apfc a_arg, b_arg;
                int a_val, b_val;
                
                next_token();
                if (!parse_expr(&a_arg)) return 0;
                if (current_token.type != TOK_COMMA) {
                    printf("Error: nor requires (a, b)\n");
                    return 0;
                }
                next_token();
                if (!parse_expr(&b_arg)) return 0;
                if (current_token.type != TOK_RPAREN) {
                    printf("Error: expected ')'\n");
                    return 0;
                }
                next_token();
                
                a_val = !apf_is_zero(&a_arg.re);
                b_val = !apf_is_zero(&b_arg.re);
                apf_from_int(&result->v.scalar.re, !(a_val || b_val));
                apf_zero(&result->v.scalar.im);
                result->type = VAL_SCALAR;
                return 1;
            }
            
            /* implies(a, b) - Logical implication */
            if (str_eq(name, "implies") || str_eq(name, "imply")) {
                apfc a_arg, b_arg;
                int a_val, b_val;
                
                next_token();
                if (!parse_expr(&a_arg)) return 0;
                if (current_token.type != TOK_COMMA) {
                    printf("Error: implies requires (a, b)\n");
                    return 0;
                }
                next_token();
                if (!parse_expr(&b_arg)) return 0;
                if (current_token.type != TOK_RPAREN) {
                    printf("Error: expected ')'\n");
                    return 0;
                }
                next_token();
                
                a_val = !apf_is_zero(&a_arg.re);
                b_val = !apf_is_zero(&b_arg.re);
                apf_from_int(&result->v.scalar.re, !a_val || b_val);
                apf_zero(&result->v.scalar.im);
                result->type = VAL_SCALAR;
                return 1;
            }
            
            /* === Equation Solving === */
            
            /* quadratic(a, b, c) - Solve ax^2+bx+c=0 */
            if (str_eq(name, "quadratic") || str_eq(name, "quadroots")) {
                apfc a_arg, b_arg, c_arg;
                apf disc, sqrt_disc, neg_b, two_a, four_ac, x1, x2;
                
                next_token();
                if (!parse_expr(&a_arg)) return 0;
                if (current_token.type != TOK_COMMA) {
                    printf("Error: quadratic requires (a, b, c)\n");
                    return 0;
                }
                next_token();
                if (!parse_expr(&b_arg)) return 0;
                if (current_token.type != TOK_COMMA) {
                    printf("Error: quadratic requires (a, b, c)\n");
                    return 0;
                }
                next_token();
                if (!parse_expr(&c_arg)) return 0;
                if (current_token.type != TOK_RPAREN) {
                    printf("Error: expected ')'\n");
                    return 0;
                }
                next_token();
                
                /* disc = b^2 - 4ac */
                apf_mul(&disc, &b_arg.re, &b_arg.re);
                apf_mul(&four_ac, &a_arg.re, &c_arg.re);
                apf_from_int(&two_a, 4);
                apf_mul(&four_ac, &four_ac, &two_a);
                apf_sub(&disc, &disc, &four_ac);
                
                apf_neg(&neg_b, &b_arg.re);
                apf_from_int(&two_a, 2);
                apf_mul(&two_a, &two_a, &a_arg.re);
                
                result->type = VAL_MATRIX;
                mat_zero(&result->v.matrix, 1, 2);
                
                if (disc.sign) {
                    /* Complex roots */
                    apf real_part, imag_part, abs_disc;
                    apf_div(&real_part, &neg_b, &two_a);
                    apf_abs(&abs_disc, &disc);
                    apf_sqrt(&imag_part, &abs_disc);
                    apf_div(&imag_part, &imag_part, &two_a);
                    
                    MAT_AT(&result->v.matrix, 0, 0).re = real_part;
                    MAT_AT(&result->v.matrix, 0, 0).im = imag_part;
                    MAT_AT(&result->v.matrix, 0, 1).re = real_part;
                    apf_neg(&MAT_AT(&result->v.matrix, 0, 1).im, &imag_part);
                } else {
                    apf_sqrt(&sqrt_disc, &disc);
                    apf_add(&x1, &neg_b, &sqrt_disc);
                    apf_div(&x1, &x1, &two_a);
                    apf_sub(&x2, &neg_b, &sqrt_disc);
                    apf_div(&x2, &x2, &two_a);
                    
                    MAT_AT(&result->v.matrix, 0, 0).re = x1;
                    MAT_AT(&result->v.matrix, 0, 1).re = x2;
                }
                return 1;
            }
            
            /* === Additional Utility Functions === */
            
            /* clip(x, lo, hi) - Alias for clamp */
            if (str_eq(name, "clip")) {
                apfc lo, hi;
                int i, j;
                
                next_token();
                if (!parse_value(&pv_arg)) return 0;
                if (current_token.type != TOK_COMMA) {
                    printf("Error: clip requires (x, lo, hi)\n");
                    return 0;
                }
                next_token();
                if (!parse_expr(&lo)) return 0;
                if (current_token.type != TOK_COMMA) {
                    printf("Error: clip requires (x, lo, hi)\n");
                    return 0;
                }
                next_token();
                if (!parse_expr(&hi)) return 0;
                if (current_token.type != TOK_RPAREN) {
                    printf("Error: expected ')'\n");
                    return 0;
                }
                next_token();
                
                if (pv_arg.type == VAL_SCALAR) {
                    result->type = VAL_SCALAR;
                    if (apf_cmp(&pv_arg.v.scalar.re, &lo.re) < 0) {
                        apf_copy(&result->v.scalar.re, &lo.re);
                    } else if (apf_cmp(&pv_arg.v.scalar.re, &hi.re) > 0) {
                        apf_copy(&result->v.scalar.re, &hi.re);
                    } else {
                        apf_copy(&result->v.scalar.re, &pv_arg.v.scalar.re);
                    }
                    apf_zero(&result->v.scalar.im);
                } else {
                    result->type = VAL_MATRIX;
                    mat_zero(&result->v.matrix, pv_arg.v.matrix.rows, pv_arg.v.matrix.cols);
                    for (i = 0; i < pv_arg.v.matrix.rows; i++) {
                        for (j = 0; j < pv_arg.v.matrix.cols; j++) {
                            if (apf_cmp(&MAT_AT(&pv_arg.v.matrix, i, j).re, &lo.re) < 0) {
                                apf_copy(&MAT_AT(&result->v.matrix, i, j).re, &lo.re);
                            } else if (apf_cmp(&MAT_AT(&pv_arg.v.matrix, i, j).re, &hi.re) > 0) {
                                apf_copy(&MAT_AT(&result->v.matrix, i, j).re, &hi.re);
                            } else {
                                apf_copy(&MAT_AT(&result->v.matrix, i, j).re, 
                                        &MAT_AT(&pv_arg.v.matrix, i, j).re);
                            }
                        }
                    }
                }
                return 1;
            }
            
            /* remap(x, a1, b1, a2, b2) - Map from [a1,b1] to [a2,b2] */
            if (str_eq(name, "remap") || str_eq(name, "map")) {
                apfc a1, b1, a2, b2;
                apf range1, range2, t, mapped;
                
                next_token();
                if (!parse_expr(&pv_arg.v.scalar)) return 0;
                if (current_token.type != TOK_COMMA) {
                    printf("Error: remap requires (x, a1, b1, a2, b2)\n");
                    return 0;
                }
                next_token();
                if (!parse_expr(&a1)) return 0;
                if (current_token.type != TOK_COMMA) return 0;
                next_token();
                if (!parse_expr(&b1)) return 0;
                if (current_token.type != TOK_COMMA) return 0;
                next_token();
                if (!parse_expr(&a2)) return 0;
                if (current_token.type != TOK_COMMA) return 0;
                next_token();
                if (!parse_expr(&b2)) return 0;
                if (current_token.type != TOK_RPAREN) {
                    printf("Error: expected ')'\n");
                    return 0;
                }
                next_token();
                
                /* t = (x - a1) / (b1 - a1) */
                apf_sub(&range1, &b1.re, &a1.re);
                apf_sub(&t, &pv_arg.v.scalar.re, &a1.re);
                apf_div(&t, &t, &range1);
                
                /* result = a2 + t * (b2 - a2) */
                apf_sub(&range2, &b2.re, &a2.re);
                apf_mul(&mapped, &t, &range2);
                apf_add(&result->v.scalar.re, &a2.re, &mapped);
                apf_zero(&result->v.scalar.im);
                result->type = VAL_SCALAR;
                return 1;
            }
            
            /* dist(x1, y1, x2, y2) - Euclidean distance */
            if (str_eq(name, "dist") || str_eq(name, "distance")) {
                apfc x1, y1, x2, y2;
                apf dx, dy, dx2, dy2, sum;
                
                next_token();
                if (!parse_expr(&x1)) return 0;
                if (current_token.type != TOK_COMMA) return 0;
                next_token();
                if (!parse_expr(&y1)) return 0;
                if (current_token.type != TOK_COMMA) return 0;
                next_token();
                if (!parse_expr(&x2)) return 0;
                if (current_token.type != TOK_COMMA) return 0;
                next_token();
                if (!parse_expr(&y2)) return 0;
                if (current_token.type != TOK_RPAREN) {
                    printf("Error: expected ')'\n");
                    return 0;
                }
                next_token();
                
                apf_sub(&dx, &x2.re, &x1.re);
                apf_sub(&dy, &y2.re, &y1.re);
                apf_mul(&dx2, &dx, &dx);
                apf_mul(&dy2, &dy, &dy);
                apf_add(&sum, &dx2, &dy2);
                apf_sqrt(&result->v.scalar.re, &sum);
                apf_zero(&result->v.scalar.im);
                result->type = VAL_SCALAR;
                return 1;
            }
            
            /* manhattan(x1, y1, x2, y2) - Manhattan distance */
            if (str_eq(name, "manhattan")) {
                apfc x1, y1, x2, y2;
                apf dx, dy, abs_dx, abs_dy;
                
                next_token();
                if (!parse_expr(&x1)) return 0;
                if (current_token.type != TOK_COMMA) return 0;
                next_token();
                if (!parse_expr(&y1)) return 0;
                if (current_token.type != TOK_COMMA) return 0;
                next_token();
                if (!parse_expr(&x2)) return 0;
                if (current_token.type != TOK_COMMA) return 0;
                next_token();
                if (!parse_expr(&y2)) return 0;
                if (current_token.type != TOK_RPAREN) {
                    printf("Error: expected ')'\n");
                    return 0;
                }
                next_token();
                
                apf_sub(&dx, &x2.re, &x1.re);
                apf_sub(&dy, &y2.re, &y1.re);
                apf_abs(&abs_dx, &dx);
                apf_abs(&abs_dy, &dy);
                apf_add(&result->v.scalar.re, &abs_dx, &abs_dy);
                apf_zero(&result->v.scalar.im);
                result->type = VAL_SCALAR;
                return 1;
            }
            
            /* midpoint(x1, y1, x2, y2) - Midpoint [x, y] */
            if (str_eq(name, "midpoint")) {
                apfc x1, y1, x2, y2;
                apf two, mx, my;
                
                next_token();
                if (!parse_expr(&x1)) return 0;
                if (current_token.type != TOK_COMMA) return 0;
                next_token();
                if (!parse_expr(&y1)) return 0;
                if (current_token.type != TOK_COMMA) return 0;
                next_token();
                if (!parse_expr(&x2)) return 0;
                if (current_token.type != TOK_COMMA) return 0;
                next_token();
                if (!parse_expr(&y2)) return 0;
                if (current_token.type != TOK_RPAREN) {
                    printf("Error: expected ')'\n");
                    return 0;
                }
                next_token();
                
                apf_from_int(&two, 2);
                apf_add(&mx, &x1.re, &x2.re);
                apf_div(&mx, &mx, &two);
                apf_add(&my, &y1.re, &y2.re);
                apf_div(&my, &my, &two);
                
                result->type = VAL_MATRIX;
                mat_zero(&result->v.matrix, 1, 2);
                MAT_AT(&result->v.matrix, 0, 0).re = mx;
                MAT_AT(&result->v.matrix, 0, 1).re = my;
                return 1;
            }
            
            /* golden() - Golden ratio phi = (1+sqrt(5))/2 */
            if (str_eq(name, "golden") || str_eq(name, "phi")) {
                apf one, five, sqrt5, two;
                
                next_token();
                if (current_token.type != TOK_RPAREN) {
                    printf("Error: expected ')'\n");
                    return 0;
                }
                next_token();
                
                apf_from_int(&one, 1);
                apf_from_int(&five, 5);
                apf_sqrt(&sqrt5, &five);
                apf_add(&result->v.scalar.re, &one, &sqrt5);
                apf_from_int(&two, 2);
                apf_div(&result->v.scalar.re, &result->v.scalar.re, &two);
                apf_zero(&result->v.scalar.im);
                result->type = VAL_SCALAR;
                return 1;
            }
            
            /* silver() - Silver ratio delta_S = 1+sqrt(2) */
            if (str_eq(name, "silver")) {
                apf one, two, sqrt2;
                
                next_token();
                if (current_token.type != TOK_RPAREN) {
                    printf("Error: expected ')'\n");
                    return 0;
                }
                next_token();
                
                apf_from_int(&one, 1);
                apf_from_int(&two, 2);
                apf_sqrt(&sqrt2, &two);
                apf_add(&result->v.scalar.re, &one, &sqrt2);
                apf_zero(&result->v.scalar.im);
                result->type = VAL_SCALAR;
                return 1;
            }
            
            /* cov(X) - covariance (variance for vector) */
            if (str_eq(name, "cov")) {
                next_token();
                if (!parse_value(&pv_arg)) return 0;
                if (current_token.type != TOK_RPAREN) {
                    printf("Error: expected ')'\n");
                    return 0;
                }
                next_token();
                
                if (pv_arg.type == VAL_SCALAR) {
                    apf_zero(&result->v.scalar.re);
                    apf_zero(&result->v.scalar.im);
                    result->type = VAL_SCALAR;
                } else {
                    /* For a vector, cov = var */
                    mat_var(&result->v.scalar, &pv_arg.v.matrix);
                    result->type = VAL_SCALAR;
                }
                return 1;
            }
            
            /* reshape(M, rows, cols) */
            if (str_eq(name, "reshape")) {
                int new_rows, new_cols, total, i;
                apfc dim1, dim2;
                
                next_token();
                if (!parse_value(&pv_arg)) return 0;
                if (pv_arg.type != VAL_MATRIX) {
                    apfc saved_s = pv_arg.v.scalar;
                    mat_zero(&pv_arg.v.matrix, 1, 1);
                    MAT_AT(&pv_arg.v.matrix, 0, 0) = saved_s;
                    pv_arg.type = VAL_MATRIX;
                }
                if (current_token.type != TOK_COMMA) {
                    printf("Error: reshape requires (M, rows, cols)\n");
                    return 0;
                }
                next_token();
                if (!parse_expr(&dim1)) return 0;
                new_rows = apf_to_long(&dim1.re);
                if (current_token.type != TOK_COMMA) {
                    printf("Error: reshape requires (M, rows, cols)\n");
                    return 0;
                }
                next_token();
                if (!parse_expr(&dim2)) return 0;
                new_cols = apf_to_long(&dim2.re);
                if (current_token.type != TOK_RPAREN) {
                    printf("Error: expected ')'\n");
                    return 0;
                }
                next_token();
                
                total = pv_arg.v.matrix.rows * pv_arg.v.matrix.cols;
                if (new_rows * new_cols != total) {
                    printf("Error: reshape size mismatch\n");
                    return 0;
                }
                
                result->type = VAL_MATRIX;
                result->v.matrix.rows = new_rows;
                result->v.matrix.cols = new_cols;
                for (i = 0; i < total; i++) {
                    int old_c = i / pv_arg.v.matrix.rows;
                    int old_r = i % pv_arg.v.matrix.rows;
                    int new_c = i / new_rows;
                    int new_r = i % new_rows;
                    MAT_AT(&result->v.matrix, new_r, new_c) = 
                        MAT_AT(&pv_arg.v.matrix, old_r, old_c);
                }
                return 1;
            }
            
            /* repmat(M, rows, cols) - replicate matrix */
            if (str_eq(name, "repmat")) {
                int rep_rows, rep_cols, r, c, rr, rc;
                apfc dim1, dim2;
                
                next_token();
                if (!parse_value(&pv_arg)) return 0;
                if (pv_arg.type != VAL_MATRIX) {
                    apfc saved_s = pv_arg.v.scalar;
                    mat_zero(&pv_arg.v.matrix, 1, 1);
                    MAT_AT(&pv_arg.v.matrix, 0, 0) = saved_s;
                    pv_arg.type = VAL_MATRIX;
                }
                if (current_token.type != TOK_COMMA) {
                    printf("Error: repmat requires (M, rows, cols)\n");
                    return 0;
                }
                next_token();
                if (!parse_expr(&dim1)) return 0;
                rep_rows = apf_to_long(&dim1.re);
                if (current_token.type != TOK_COMMA) {
                    printf("Error: repmat requires (M, rows, cols)\n");
                    return 0;
                }
                next_token();
                if (!parse_expr(&dim2)) return 0;
                rep_cols = apf_to_long(&dim2.re);
                if (current_token.type != TOK_RPAREN) {
                    printf("Error: expected ')'\n");
                    return 0;
                }
                next_token();
                
                if (pv_arg.v.matrix.rows * rep_rows > MAT_MAX_ROWS ||
                    pv_arg.v.matrix.cols * rep_cols > MAT_MAX_COLS) {
                    printf("Error: repmat result too large\n");
                    return 0;
                }
                
                result->type = VAL_MATRIX;
                result->v.matrix.rows = pv_arg.v.matrix.rows * rep_rows;
                result->v.matrix.cols = pv_arg.v.matrix.cols * rep_cols;
                
                for (rr = 0; rr < rep_rows; rr++) {
                    for (rc = 0; rc < rep_cols; rc++) {
                        for (r = 0; r < pv_arg.v.matrix.rows; r++) {
                            for (c = 0; c < pv_arg.v.matrix.cols; c++) {
                                MAT_AT(&result->v.matrix, 
                                      rr * pv_arg.v.matrix.rows + r,
                                      rc * pv_arg.v.matrix.cols + c) = 
                                    MAT_AT(&pv_arg.v.matrix, r, c);
                            }
                        }
                    }
                }
                return 1;
            }
            
            /* polyval(p, x) - evaluate polynomial at x */
            if (str_eq(name, "polyval")) {
                apfc x_val;
                int n, i;
                
                next_token();
                if (!parse_value(&pv_arg)) return 0;
                if (pv_arg.type != VAL_MATRIX) {
                    apfc saved_s = pv_arg.v.scalar;
                    mat_zero(&pv_arg.v.matrix, 1, 1);
                    MAT_AT(&pv_arg.v.matrix, 0, 0) = saved_s;
                    pv_arg.type = VAL_MATRIX;
                }
                if (current_token.type != TOK_COMMA) {
                    printf("Error: polyval requires (p, x)\n");
                    return 0;
                }
                next_token();
                if (!parse_expr(&x_val)) return 0;
                if (current_token.type != TOK_RPAREN) {
                    printf("Error: expected ')'\n");
                    return 0;
                }
                next_token();
                
                /* Horner's method: p = [a_n ... a_1 a_0] */
                n = pv_arg.v.matrix.rows * pv_arg.v.matrix.cols;
                result->type = VAL_SCALAR;
                apf_zero(&result->v.scalar.re);
                apf_zero(&result->v.scalar.im);
                for (i = 0; i < n; i++) {
                    int ri = i % pv_arg.v.matrix.rows;
                    int ci = i / pv_arg.v.matrix.rows;
                    apfc_mul(&result->v.scalar, &result->v.scalar, &x_val);
                    apfc_add(&result->v.scalar, &result->v.scalar, &MAT_AT(&pv_arg.v.matrix, ri, ci));
                }
                return 1;
            }
            
            /* polyder(p) - polynomial derivative coefficients */
            if (str_eq(name, "polyder")) {
                int n, i;
                
                next_token();
                if (!parse_value(&pv_arg)) return 0;
                if (current_token.type != TOK_RPAREN) {
                    printf("Error: expected ')'\n");
                    return 0;
                }
                next_token();
                
                if (pv_arg.type != VAL_MATRIX) {
                    apfc saved_s = pv_arg.v.scalar;
                    mat_zero(&pv_arg.v.matrix, 1, 1);
                    MAT_AT(&pv_arg.v.matrix, 0, 0) = saved_s;
                    pv_arg.type = VAL_MATRIX;
                }
                
                n = pv_arg.v.matrix.rows * pv_arg.v.matrix.cols;
                result->type = VAL_MATRIX;
                if (n <= 1) {
                    result->v.matrix.rows = 1;
                    result->v.matrix.cols = 1;
                    apf_zero(&MAT_AT(&result->v.matrix, 0, 0).re);
                    apf_zero(&MAT_AT(&result->v.matrix, 0, 0).im);
                } else {
                    result->v.matrix.rows = 1;
                    result->v.matrix.cols = n - 1;
                    for (i = 0; i < n - 1; i++) {
                        int ri = i % pv_arg.v.matrix.rows;
                        int ci = i / pv_arg.v.matrix.rows;
                        apf coef;
                        apf_from_int(&coef, n - 1 - i);
                        apf_mul(&MAT_AT(&result->v.matrix, 0, i).re, &MAT_AT(&pv_arg.v.matrix, ri, ci).re, &coef);
                        apf_mul(&MAT_AT(&result->v.matrix, 0, i).im, &MAT_AT(&pv_arg.v.matrix, ri, ci).im, &coef);
                    }
                }
                return 1;
            }
            
            /* polyint(p) - polynomial integral coefficients (constant = 0) */
            if (str_eq(name, "polyint")) {
                int n, i;
                
                next_token();
                if (!parse_value(&pv_arg)) return 0;
                if (current_token.type != TOK_RPAREN) {
                    printf("Error: expected ')'\n");
                    return 0;
                }
                next_token();
                
                if (pv_arg.type != VAL_MATRIX) {
                    apfc saved_s = pv_arg.v.scalar;
                    mat_zero(&pv_arg.v.matrix, 1, 1);
                    MAT_AT(&pv_arg.v.matrix, 0, 0) = saved_s;
                    pv_arg.type = VAL_MATRIX;
                }
                
                n = pv_arg.v.matrix.rows * pv_arg.v.matrix.cols;
                result->type = VAL_MATRIX;
                result->v.matrix.rows = 1;
                result->v.matrix.cols = n + 1;
                for (i = 0; i < n; i++) {
                    int ri = i % pv_arg.v.matrix.rows;
                    int ci = i / pv_arg.v.matrix.rows;
                    apf coef;
                    apf_from_int(&coef, n - i);
                    apf_div(&MAT_AT(&result->v.matrix, 0, i).re, &MAT_AT(&pv_arg.v.matrix, ri, ci).re, &coef);
                    apf_div(&MAT_AT(&result->v.matrix, 0, i).im, &MAT_AT(&pv_arg.v.matrix, ri, ci).im, &coef);
                }
                /* Constant term = 0 */
                apf_zero(&MAT_AT(&result->v.matrix, 0, n).re);
                apf_zero(&MAT_AT(&result->v.matrix, 0, n).im);
                return 1;
            }
            
            /* nnz(v) - count nonzero elements */
            if (str_eq(name, "nnz")) {
                int n, i, count;
                
                next_token();
                if (!parse_value(&pv_arg)) return 0;
                if (current_token.type != TOK_RPAREN) {
                    printf("Error: expected ')'\n");
                    return 0;
                }
                next_token();
                
                if (pv_arg.type != VAL_MATRIX) {
                    count = (!apf_is_zero(&pv_arg.v.scalar.re) || !apf_is_zero(&pv_arg.v.scalar.im)) ? 1 : 0;
                } else {
                    n = pv_arg.v.matrix.rows * pv_arg.v.matrix.cols;
                    count = 0;
                    for (i = 0; i < n; i++) {
                        int ri = i % pv_arg.v.matrix.rows;
                        int ci = i / pv_arg.v.matrix.rows;
                        if (!apf_is_zero(&MAT_AT(&pv_arg.v.matrix, ri, ci).re) ||
                            !apf_is_zero(&MAT_AT(&pv_arg.v.matrix, ri, ci).im)) {
                            count++;
                        }
                    }
                }
                result->type = VAL_SCALAR;
                apf_from_int(&result->v.scalar.re, count);
                apf_zero(&result->v.scalar.im);
                return 1;
            }
            
            /* randperm(n) - random permutation of 1:n */
            if (str_eq(name, "randperm")) {
                int n, i, j;
                apfc tmp;
                
                next_token();
                if (!parse_value(&pv_arg)) return 0;
                if (current_token.type != TOK_RPAREN) {
                    printf("Error: expected ')'\n");
                    return 0;
                }
                next_token();
                
                if (pv_arg.type == VAL_SCALAR) {
                    n = apf_to_long(&pv_arg.v.scalar.re);
                } else {
                    n = apf_to_long(&MAT_AT(&pv_arg.v.matrix, 0, 0).re);
                }
                
                if (n < 1) n = 1;
                if (n > MAT_MAX_COLS) n = MAT_MAX_COLS;
                
                result->type = VAL_MATRIX;
                result->v.matrix.rows = 1;
                result->v.matrix.cols = n;
                
                /* Initialize with 1:n */
                for (i = 0; i < n; i++) {
                    apf_from_int(&MAT_AT(&result->v.matrix, 0, i).re, i + 1);
                    apf_zero(&MAT_AT(&result->v.matrix, 0, i).im);
                }
                
                /* Fisher-Yates shuffle */
                for (i = n - 1; i > 0; i--) {
                    j = rand() % (i + 1);
                    tmp = MAT_AT(&result->v.matrix, 0, i);
                    MAT_AT(&result->v.matrix, 0, i) = MAT_AT(&result->v.matrix, 0, j);
                    MAT_AT(&result->v.matrix, 0, j) = tmp;
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
                if (n < 0 || n > 100000) {
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

/* Parse value term (*, /, %, \, .*, ./) */
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

/* Parse a value (scalar or matrix) with optional colon range */
static int parse_value_range(value_t *result)
{
    if (!parse_value_expr(result)) return 0;
    
    /* Check for colon notation: start:end or start:step:end */
    if (current_token.type == TOK_COLON && result->type == VAL_SCALAR) {
        apf start, end_val, step;
        int n, i;
        
        start = result->v.scalar.re;
        
        next_token();
        if (!parse_value_expr(result)) return 0;
        
        if (result->type != VAL_SCALAR) {
            printf("Error: range requires scalar values\n");
            return 0;
        }
        
        if (current_token.type == TOK_COLON) {
            /* start:step:end format */
            step = result->v.scalar.re;
            next_token();
            if (!parse_value_expr(result)) return 0;
            if (result->type != VAL_SCALAR) {
                printf("Error: range requires scalar values\n");
                return 0;
            }
            end_val = result->v.scalar.re;
        } else {
            /* start:end format (step = 1) */
            end_val = result->v.scalar.re;
            apf_from_int(&step, 1);
        }
        
        /* Calculate number of elements */
        {
            apf diff, count;
            apf_sub(&diff, &end_val, &start);
            apf_div(&count, &diff, &step);
            n = apf_to_long(&count) + 1;
        }
        
        if (n < 1) n = 0;
        if (n > MAT_MAX_COLS) n = MAT_MAX_COLS;
        
        /* Create row vector */
        result->type = VAL_MATRIX;
        result->v.matrix.rows = 1;
        result->v.matrix.cols = n;
        
        for (i = 0; i < n; i++) {
            apf idx, val;
            apf_from_int(&idx, i);
            apf_mul(&val, &step, &idx);
            apf_add(&MAT_AT(&result->v.matrix, 0, i).re, &start, &val);
            apf_zero(&MAT_AT(&result->v.matrix, 0, i).im);
        }
    }
    
    return 1;
}

/* Parse a value (scalar or matrix) - main entry point */
int parse_value(value_t *result)
{
    return parse_value_range(result);
}
