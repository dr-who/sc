/* parser.c - Expression parser with built-in functions
 * Recursive descent parser with precedence climbing
 * C89 compliant for Watcom C / DOS
 * 16-bit clean: int is 16-bit, long is 32-bit
 */
#define _XOPEN_SOURCE 600
#include <math.h>
#ifndef INFINITY
#define INFINITY (1.0/0.0)
#endif

/* C89-compatible lgamma approximation */
static double p_lgamma(double x)
{
    static double c[7] = {
        1.000000000190015,
        76.18009172947146,
        -86.50532032941677,
        24.01409824083091,
        -1.231739572450155,
        0.1208650973866179e-2,
        -0.5395239384953e-5
    };
    double g = 5.0;
    double sum = c[0];
    int i;
    double t;
    
    if (x <= 0 && x == floor(x)) return 1e308;
    
    for (i = 1; i < 7; i++) sum += c[i] / (x + i - 1);
    t = x + g - 0.5;
    return 0.5 * log(2 * 3.14159265358979323846) + (x - 0.5) * log(t) - t + log(sum);
}
#define lgamma p_lgamma

#include "sc.h"
#include "scalar_funcs.h"
#include "matrix_funcs.h"
#include "decomp_funcs.h"
#include "table.h"
#include "apf_native.h"

#ifdef HAVE_MATRIX
#include "ml.h"
#endif

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

/* Pre-allocated data for static matrix temporaries */
#define PARSER_TMP_SIZE 256
static apfc pv_tmp_data[PARSER_TMP_SIZE];
static apfc pv_tmp_data2[PARSER_TMP_SIZE];
static apfc pv_arg_data[PARSER_TMP_SIZE];

/* Result buffer pool for parse_value results */
#define RESULT_BUF_COUNT 8
#define RESULT_BUF_SIZE 256
static apfc result_bufs[RESULT_BUF_COUNT][RESULT_BUF_SIZE];
static int result_buf_idx = 0;

/* Get next result buffer (cycles through pool) */
static apfc *get_result_buf(void)
{
    apfc *buf = result_bufs[result_buf_idx];
    result_buf_idx = (result_buf_idx + 1) % RESULT_BUF_COUNT;
    return buf;
}

/* Ensure result matrix has allocated data - uses pool for small, arena for large */
int ensure_result_matrix(matrix_t *m, int rows, int cols)
{
    int n = rows * cols;
    m->rows = rows;
    m->cols = cols;
    if (n <= RESULT_BUF_SIZE) {
        m->data = get_result_buf();
        return 1;
    } else {
        m->data = (apfc *)mat_arena_alloc((size_t)n * sizeof(apfc));
        return m->data != NULL;
    }
}

/* Initialize parser temporaries - call once at startup */
void parser_init_temps(void)
{
    pv_tmp_mat.data = pv_tmp_data;
    pv_tmp_mat.rows = 0;
    pv_tmp_mat.cols = 0;
    
    pv_tmp_mat2.data = pv_tmp_data2;
    pv_tmp_mat2.rows = 0;
    pv_tmp_mat2.cols = 0;
    
    pv_arg.v.matrix.data = pv_arg_data;
    pv_arg.v.matrix.rows = 0;
    pv_arg.v.matrix.cols = 0;
    
    result_buf_idx = 0;
}

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
        char name[32];
        int func_idx, var_idx;
        strncpy(name, current_token.func_name, sizeof(name)-1);
        name[sizeof(name)-1] = '\0';


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
        
        /* Named scalar variable lookup (multi-character names like "span") */
        if (var_idx < 0 && current_token.type != TOK_LPAREN) {
            extern int get_named_var(const char *name, value_t *val);
            value_t named_val;
            if (get_named_var(name, &named_val)) {
                if (named_val.type == VAL_SCALAR) {
                    *result = named_val.v.scalar;
                    return 1;
                } else if (named_val.type == VAL_MATRIX && 
                           named_val.v.matrix.rows == 1 && named_val.v.matrix.cols == 1) {
                    *result = MAT_AT(&named_val.v.matrix, 0, 0);
                    return 1;
                }
                /* It's a named matrix but not 1x1 - fall through to error */
            }
        }
        
        /* Single-letter matrix variable with indexing: t(1), a(2,3) */
        if (var_idx >= 0 && current_token.type == TOK_LPAREN && is_var_matrix(var_idx)) {
            int mi;
            for (mi = 0; mi < MAX_MATRIX_VARS; mi++) {
                if (matrix_vars[mi].defined && matrix_vars[mi].name == 'a' + var_idx) {
                    matrix_t *mat = &matrix_vars[mi].val;
                    apfc idx_val;
                    long idx;
                    int n;
                    
                    next_token();  /* skip ( */
                    if (!parse_expr(&idx_val)) return 0;
                    
                    /* Check for second index (2D indexing) */
                    if (current_token.type == TOK_COMMA) {
                        long row, col;
                        apfc idx2_val;
                        next_token();
                        if (!parse_expr(&idx2_val)) return 0;
                        if (current_token.type != TOK_RPAREN) {
                            printf("Error: expected ')'\n");
                            return 0;
                        }
                        next_token();
                        row = apf_to_long(&idx_val.re) - 1;
                        col = apf_to_long(&idx2_val.re) - 1;
                        if (row < 0 || row >= mat->rows || col < 0 || col >= mat->cols) {
                            printf("Error: index out of range\n");
                            return 0;
                        }
                        *result = MAT_AT(mat, row, col);
                        return 1;
                    }
                    
                    if (current_token.type != TOK_RPAREN) {
                        printf("Error: expected ')' after index\n");
                        return 0;
                    }
                    next_token();
                    
                    idx = apf_to_long(&idx_val.re) - 1;  /* 1-based to 0-based */
                    n = mat->rows * mat->cols;
                    if (idx < 0 || idx >= n) {
                        printf("Error: index %ld out of range (1-%d)\n", idx+1, n);
                        return 0;
                    }
                    /* Linear indexing */
                    {
                        int ri = idx % mat->rows;
                        int ci = idx / mat->rows;
                        *result = MAT_AT(mat, ri, ci);
                    }
                    return 1;
                }
            }
        }
        
        /* Named variable lookup with indexing: varname(i) */
        if (current_token.type == TOK_LPAREN) {
            extern int get_named_var(const char *name, value_t *val);
            value_t named_val;
            if (get_named_var(name, &named_val)) {
                /* It's a named variable - handle indexing */
                apfc idx_val;
                long idx;
                next_token();  /* skip ( */
                if (!parse_expr(&idx_val)) return 0;
                if (current_token.type != TOK_RPAREN) {
                    printf("Error: expected ')' after index\n");
                    return 0;
                }
                next_token();
                
                idx = apf_to_long(&idx_val.re) - 1;  /* 1-based to 0-based */
                
                if (named_val.type == VAL_MATRIX) {
                    int n = named_val.v.matrix.rows * named_val.v.matrix.cols;
                    if (idx < 0 || idx >= n) {
                        printf("Error: index %ld out of range (1-%d)\n", idx+1, n);
                        return 0;
                    }
                    /* Linear indexing */
                    {
                        int ri = idx % named_val.v.matrix.rows;
                        int ci = idx / named_val.v.matrix.rows;
                        *result = MAT_AT(&named_val.v.matrix, ri, ci);
                    }
                } else {
                    /* Scalar - index 1 returns the value */
                    if (idx != 0) {
                        printf("Error: scalar indexed with %ld\n", idx+1);
                        return 0;
                    }
                    *result = named_val.v.scalar;
                }
                return 1;
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
        if (str_eq(name, "i") || str_eq(name, "j")) {
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
        if (str_eq(name, "epsp") || str_eq(name, "epsn")) {
            /* Native epsilon for current AP precision: 2^(-AP_BITS+1) */
            apf two;
            apf_from_int(&two, 2);
            apf_from_int(&result->re, -(AP_BITS - 1));
            apfx_pow(&result->re, &two, &result->re);
            apf_zero(&result->im);
            return 1;
        }
        if (str_eq(name, "precision") || str_eq(name, "apbits")) {
            /* Report current precision in bits */
            apf_from_int(&result->re, AP_BITS);
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
        if (str_eq(name, "phi") || str_eq(name, "golden")) {
            /* Golden ratio: (1+sqrt(5))/2 */
            apf five, sqrtfive, one, two;
            apf_from_int(&five, 5);
            apf_sqrt(&sqrtfive, &five);
            apf_from_int(&one, 1);
            apf_add(&result->re, &one, &sqrtfive);
            apf_from_int(&two, 2);
            apf_div(&result->re, &result->re, &two);
            apf_zero(&result->im);
            return 1;
        }

        /* Random functions */
        if (str_eq(name, "rand") || str_eq(name, "randn") || str_eq(name, "randi")) {
            return parse_random_func(result, name);
        }

        /* planet("name", t_days) - Get planet position */
        if (str_eq(name, "planet")) {
            extern void fn_planet(const char *name, double t_days, double *x, double *y);
            char planet_name[32];
            double t_days = 0, px, py;
            int ni = 0;
            
            if (current_token.type != TOK_LPAREN) {
                printf("Error: expected '(' after 'planet'\n");
                return 0;
            }
            
            /* Check for string literal */
            while (*input_ptr == ' ') input_ptr++;
            if (*input_ptr == '\'' || *input_ptr == '"') {
                char delim = *input_ptr++;
                while (*input_ptr && *input_ptr != delim && ni < 31) {
                    planet_name[ni++] = *input_ptr++;
                }
                planet_name[ni] = '\0';
                if (*input_ptr == delim) input_ptr++;
                while (*input_ptr == ' ') input_ptr++;
                
                /* Check for comma and second arg */
                if (*input_ptr == ',') {
                    input_ptr++;
                    next_token();
                    if (!parse_expr(&arg)) return 0;
                    t_days = apf_to_double(&arg.re);
                } else {
                    next_token();
                }
                
                if (current_token.type != TOK_RPAREN) {
                    printf("Error: expected ')'\n");
                    return 0;
                }
                next_token();
                
                fn_planet(planet_name, t_days, &px, &py);
                apf_from_double(&result->re, px);
                apf_from_double(&result->im, py);
                return 1;
            } else {
                printf("Error: planet() requires quoted string: planet(\"earth\", 0)\n");
                return 0;
            }
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
        
        /* fprintf - formatted print with %d, %f, %s, \n, etc. */
        if (str_eq(name, "fprintf") || str_eq(name, "printf")) {
            char fmt[256];
            int fi = 0;
            const char *p;
            
            if (current_token.type != TOK_LPAREN) {
                printf("Error: expected '(' after '%s'\n", name);
                return 0;
            }
            while (*input_ptr == ' ') input_ptr++;
            if (*input_ptr != '\'' && *input_ptr != '"') {
                printf("Error: fprintf requires format string\n");
                return 0;
            }
            
            {
                char delim = *input_ptr++;
                
                /* Read format string */
                while (*input_ptr && *input_ptr != delim && fi < 255) {
                    if (*input_ptr == '\\' && *(input_ptr+1)) {
                        input_ptr++;
                        if (*input_ptr == 'n') { fmt[fi++] = '\n'; input_ptr++; }
                        else if (*input_ptr == 't') { fmt[fi++] = '\t'; input_ptr++; }
                        else if (*input_ptr == '\\') { fmt[fi++] = '\\'; input_ptr++; }
                        else { fmt[fi++] = *input_ptr++; }
                    } else {
                        fmt[fi++] = *input_ptr++;
                    }
                }
                fmt[fi] = '\0';
                if (*input_ptr == delim) input_ptr++;
            }
            
            /* Skip comma before first arg if present */
            while (*input_ptr == ' ') input_ptr++;
            if (*input_ptr == ',') input_ptr++;
            next_token();  /* Prepare for argument parsing */
            
            /* Process format string, parsing args as needed */
            p = fmt;
            while (*p) {
                if (*p == '%' && *(p+1)) {
                    char spec = *(p+1);  /* Save the format specifier */
                    p++;  /* Skip the '%' */
                    if (spec == '%') {
                        putchar('%');
                        p++;  /* Skip the second '%' */
                    } else if (spec == 'd' || spec == 'i') {
                        /* Integer arg */
                        if (!parse_expr(&arg)) return 0;
                        printf("%ld", apf_to_long(&arg.re));
                        p++;
                        /* Skip comma before next arg */
                        if (current_token.type == TOK_COMMA) next_token();
                    } else if (spec == 'f' || spec == 'g' || spec == 'e') {
                        /* Float arg */
                        if (!parse_expr(&arg)) return 0;
                        printf("%g", apf_to_double(&arg.re));
                        p++;
                        if (current_token.type == TOK_COMMA) next_token();
                    } else if (spec == '.') {
                        /* Format like %.2f */
                        int prec = 0;
                        p++;  /* Skip the '.' */
                        while (*p >= '0' && *p <= '9') {
                            prec = prec * 10 + (*p - '0');
                            p++;
                        }
                        if (*p == 'f' || *p == 'g' || *p == 'e') {
                            if (!parse_expr(&arg)) return 0;
                            printf("%.*f", prec, apf_to_double(&arg.re));
                            p++;
                            if (current_token.type == TOK_COMMA) next_token();
                        } else {
                            putchar('%');
                            putchar('.');
                        }
                    } else if (spec == 's') {
                        /* String arg - skip for now */
                        if (!parse_expr(&arg)) return 0;
                        printf("<str>");
                        p++;
                        if (current_token.type == TOK_COMMA) next_token();
                    } else {
                        putchar('%');
                        putchar(spec);
                        p++;
                    }
                } else {
                    putchar(*p++);
                }
            }
            
            if (current_token.type != TOK_RPAREN) {
                printf("Error: expected ')' in fprintf\n");
                return 0;
            }
            next_token();
            apf_zero(&result->re);
            apf_zero(&result->im);
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

        /* size(A, dim) - returns rows (dim=1) or cols (dim=2) */
        if (str_eq(name, "size")) {
            value_t local_arg, dim_arg;
            if (current_token.type != TOK_LPAREN) {
                printf("Error: expected '(' after function 'size'\n");
                return 0;
            }
            next_token();
            if (!parse_value(&local_arg)) return 0;
            
            if (current_token.type == TOK_COMMA) {
                /* Two-argument form: size(A, dim) */
                int dim;
                next_token();
                if (!parse_value(&dim_arg)) return 0;
                if (current_token.type != TOK_RPAREN) {
                    printf("Error: expected ')'\n");
                    return 0;
                }
                next_token();
                
                dim = (dim_arg.type == VAL_SCALAR) ? (int)apf_to_long(&dim_arg.v.scalar.re) : 1;
                
                if (local_arg.type == VAL_SCALAR) {
                    apf_from_int(&result->re, 1);  /* scalar has size 1 in all dims */
                } else if (dim == 1) {
                    apf_from_int(&result->re, local_arg.v.matrix.rows);
                } else if (dim == 2) {
                    apf_from_int(&result->re, local_arg.v.matrix.cols);
                } else {
                    apf_from_int(&result->re, 1);  /* higher dims are 1 */
                }
                apf_zero(&result->im);
                return 1;
            }
            
            /* Single-argument form: size(A) returns total elements */
            if (current_token.type != TOK_RPAREN) {
                printf("Error: expected ')'\n");
                return 0;
            }
            next_token();
            
            if (local_arg.type == VAL_SCALAR) {
                apf_from_int(&result->re, 1);
            } else {
                apf_from_int(&result->re, local_arg.v.matrix.rows * local_arg.v.matrix.cols);
            }
            apf_zero(&result->im);
            return 1;
        }
        
        /* Matrix functions that return scalar */
        if (str_eq(name, "det") || str_eq(name, "trace") || str_eq(name, "tr") ||
            str_eq(name, "norm") || str_eq(name, "length") ||
            str_eq(name, "rows") || str_eq(name, "cols") || str_eq(name, "numel") ||
            str_eq(name, "height") || str_eq(name, "width") ||
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
            
            /* Handle two-argument min/max: min(a,b), max(a,b) */
            if ((str_eq(name, "min") || str_eq(name, "max")) && current_token.type == TOK_COMMA) {
                value_t arg2;
                int cmp;
                next_token();
                if (!parse_value(&arg2)) return 0;
                if (current_token.type != TOK_RPAREN) {
                    printf("Error: expected ')'\n");
                    return 0;
                }
                next_token();
                /* Both must be scalars for two-arg version */
                if (local_arg.type == VAL_SCALAR && arg2.type == VAL_SCALAR) {
                    cmp = apf_cmp(&local_arg.v.scalar.re, &arg2.v.scalar.re);
                    if (str_eq(name, "max")) {
                        *result = (cmp >= 0) ? local_arg.v.scalar : arg2.v.scalar;
                    } else {
                        *result = (cmp <= 0) ? local_arg.v.scalar : arg2.v.scalar;
                    }
                    return 1;
                }
                printf("Error: two-argument min/max requires scalars\n");
                return 0;
            }
            
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
                } else if (str_eq(name, "numel") || str_eq(name, "length")) {
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
                } else if (str_eq(name, "numel")) {
                    apf_from_int(&result->re, m->rows * m->cols);
                    apf_zero(&result->im);
                } else if (str_eq(name, "length")) {
                    apf_from_int(&result->re, mat_length(m));
                    apf_zero(&result->im);
                } else if (str_eq(name, "rows") || str_eq(name, "height")) {
                    apf_from_int(&result->re, m->rows);
                    apf_zero(&result->im);
                } else if (str_eq(name, "cols") || str_eq(name, "width")) {
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
                matrix_t L = {0,0,NULL}, U = {0,0,NULL};
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
                    matrix_t inv_m = {0, 0, NULL};
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
        /* round(x, n) - round to n decimal places */
        if (str_eq(name, "round") && current_token.type == TOK_COMMA) {
            apfc arg2;
            long n;
            apf scale, half, tmp, scaled;
            
            next_token();
            if (!parse_expr(&arg2)) return 0;
            if (current_token.type != TOK_RPAREN) {
                printf("Error: expected ')'\n");
                return 0;
            }
            next_token();
            
            if (!apfc_is_real(&arg)) {
                printf("Error: round requires real argument\n");
                return 0;
            }
            
            n = apf_to_long(&arg2.re);
            
            /* scale = 10^n */
            apf_from_int(&scale, 1);
            {
                apf ten;
                long i;
                apf_from_int(&ten, 10);
                for (i = 0; i < n; i++) {
                    apf_mul(&tmp, &scale, &ten);
                    apf_copy(&scale, &tmp);
                }
            }
            
            /* scaled = x * scale */
            apf_mul(&scaled, &arg.re, &scale);
            
            /* Round scaled value */
            apf_from_str(&half, "0.5");
            if (scaled.sign) {
                apf_sub(&tmp, &scaled, &half);
                apf_ceil(&scaled, &tmp);
            } else {
                apf_add(&tmp, &scaled, &half);
                apf_floor(&scaled, &tmp);
            }
            
            /* result = scaled / scale */
            apf_div(&result->re, &scaled, &scale);
            apf_zero(&result->im);
            return 1;
        }
        /* dateshift(t, 'start', 'month') - shift datetime to period boundary */
        if (str_eq(name, "dateshift")) {
            /* For compatibility, we accept dateshift(t, 'start', 'month'|'year'|'day') */
            /* Skip the 'start' argument and check the unit */
            int is_month = 0, is_year = 0, is_day = 0;
            long ts, days_since_epoch, year, month, i;
            long new_days;
            static const int days_in_month[] = {31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31};
            apf secs_per_day, tmp;
            
            /* Skip comma and 'start' string */
            if (current_token.type == TOK_COMMA) {
                next_token();
                if (current_token.type == TOK_STRING) {
                    /* Skip 'start' */
                    next_token();
                }
            }
            /* Get unit string */
            if (current_token.type == TOK_COMMA) {
                next_token();
                if (current_token.type == TOK_STRING) {
                    if (str_eq(current_token.str_value, "month")) is_month = 1;
                    else if (str_eq(current_token.str_value, "year")) is_year = 1;
                    else if (str_eq(current_token.str_value, "day")) is_day = 1;
                    next_token();
                }
            }
            
            if (current_token.type != TOK_RPAREN) {
                printf("Error: expected ')'\n");
                return 0;
            }
            next_token();
            
            ts = apf_to_long(&arg.re);
            days_since_epoch = ts / 86400;
            
            /* Find year and month */
            year = 1970;
            while (1) {
                int year_days = 365;
                if ((year % 4 == 0 && year % 100 != 0) || (year % 400 == 0)) year_days = 366;
                if (days_since_epoch < year_days) break;
                days_since_epoch -= year_days;
                year++;
            }
            month = 1;
            for (i = 0; i < 12; i++) {
                int dm = days_in_month[i];
                if (i == 1 && ((year % 4 == 0 && year % 100 != 0) || (year % 400 == 0))) dm++;
                if (days_since_epoch < dm) break;
                days_since_epoch -= dm;
                month++;
            }
            
            /* Calculate new_days based on unit */
            new_days = 0;
            if (is_year) {
                for (i = 1970; i < year; i++) {
                    new_days += 365;
                    if ((i % 4 == 0 && i % 100 != 0) || (i % 400 == 0)) new_days++;
                }
            } else if (is_month || (!is_year && !is_day)) {
                /* Default to month */
                for (i = 1970; i < year; i++) {
                    new_days += 365;
                    if ((i % 4 == 0 && i % 100 != 0) || (i % 400 == 0)) new_days++;
                }
                for (i = 1; i < month; i++) {
                    new_days += days_in_month[i-1];
                    if (i == 2 && ((year % 4 == 0 && year % 100 != 0) || (year % 400 == 0))) new_days++;
                }
            } else if (is_day) {
                new_days = ts / 86400;
            }
            
            apf_from_int(&result->re, new_days);
            apf_from_int(&secs_per_day, 86400);
            apf_mul(&tmp, &result->re, &secs_per_day);
            apf_copy(&result->re, &tmp);
            apf_zero(&result->im);
            return 1;
        }
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
        
        /* Try 2-arg table dispatch for functions like atan2, gcd, lcm, mod */
        if (current_token.type == TOK_COMMA) {
            ScalarFuncEntry *sf2_entry = scalar_func_lookup(name);
            if (sf2_entry && (sf2_entry->type == SF_TRIG_INV_2 || 
                              sf2_entry->type == SF_TRIG_INV_2_DEG ||
                              sf2_entry->type == SF_REAL_2 ||
                              sf2_entry->type == SF_COMPLEX_2 ||
                              sf2_entry->type == SF_INT_2 ||
                              sf2_entry->type == SF_APF_INT ||
                              sf2_entry->type == SF_INT_APF ||
                              sf2_entry->type == SF_MAKE_COMPLEX)) {
                apfc arg2;
                next_token();
                if (!parse_expr(&arg2)) return 0;
                if (current_token.type != TOK_RPAREN) {
                    printf("Error: expected ')'\n");
                    return 0;
                }
                next_token();
                return scalar_func_eval2(sf2_entry, result, &arg, &arg2);
            }
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
        
        /* mandelbrot_iter(x, y, max_iter) - Mandelbrot iteration count */
        if (str_eq(name, "mandelbrot_iter") || str_eq(name, "mandel")) {
            extern void fn_mandelbrot_iter(apfc *result, const apfc *cx, const apfc *cy, const apfc *max_iter);
            apfc cy, max_iter_val;
            
            if (current_token.type != TOK_COMMA) {
                printf("Error: mandelbrot_iter requires (x, y, max_iter)\n");
                return 0;
            }
            next_token();
            if (!parse_expr(&cy)) return 0;
            
            if (current_token.type != TOK_COMMA) {
                /* Default max_iter = 100 */
                apf_from_int(&max_iter_val.re, 100);
                apf_zero(&max_iter_val.im);
            } else {
                next_token();
                if (!parse_expr(&max_iter_val)) return 0;
            }
            
            if (current_token.type != TOK_RPAREN) {
                printf("Error: expected ')'\n");
                return 0;
            }
            next_token();
            
            fn_mandelbrot_iter(result, &arg, &cy, &max_iter_val);
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
        /* datenum(y,m,d) or datenum(y,m,d,h,mi,s) - convert to Unix timestamp */
        if (str_eq(name, "datenum") || str_eq(name, "datetime")) {
            apfc arg2, arg3, arg4, arg5, arg6;
            int year, month, day, hour = 0, minute = 0, sec = 0;
            long days, i;
            static const int days_in_month[] = {31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31};
            apf tmp, tmp2;
            
            /* Parse month */
            if (current_token.type != TOK_COMMA) {
                printf("Error: datenum requires (year,month,day) or (year,month,day,hour,min,sec)\n");
                return 0;
            }
            next_token();
            if (!parse_expr(&arg2)) return 0;
            
            /* Parse day */
            if (current_token.type != TOK_COMMA) {
                printf("Error: datenum requires at least 3 arguments\n");
                return 0;
            }
            next_token();
            if (!parse_expr(&arg3)) return 0;
            
            /* Optional: hour, minute, second */
            if (current_token.type == TOK_COMMA) {
                next_token();
                if (!parse_expr(&arg4)) return 0;
                hour = (int)apf_to_long(&arg4.re);
                
                if (current_token.type == TOK_COMMA) {
                    next_token();
                    if (!parse_expr(&arg5)) return 0;
                    minute = (int)apf_to_long(&arg5.re);
                    
                    if (current_token.type == TOK_COMMA) {
                        next_token();
                        if (!parse_expr(&arg6)) return 0;
                        sec = (int)apf_to_long(&arg6.re);
                    }
                }
            }
            
            if (current_token.type != TOK_RPAREN) {
                printf("Error: expected ')'\n");
                return 0;
            }
            next_token();
            
            year = (int)apf_to_long(&arg.re);
            month = (int)apf_to_long(&arg2.re);
            day = (int)apf_to_long(&arg3.re);
            
            /* Calculate days since 1970-01-01 (Unix epoch) */
            days = 0;
            for (i = 1970; i < year; i++) {
                days += 365;
                if ((i % 4 == 0 && i % 100 != 0) || (i % 400 == 0)) days += 1;
            }
            for (i = year; i < 1970; i++) {
                days -= 365;
                if ((i % 4 == 0 && i % 100 != 0) || (i % 400 == 0)) days -= 1;
            }
            for (i = 1; i < month; i++) {
                days += days_in_month[i-1];
                if (i == 2 && ((year % 4 == 0 && year % 100 != 0) || (year % 400 == 0))) days += 1;
            }
            days += (day - 1);
            
            /* Convert to seconds: days * 86400 + hour * 3600 + min * 60 + sec */
            apf_from_int(&result->re, days);
            apf_from_int(&tmp, 86400);
            apf_mul(&tmp2, &result->re, &tmp);
            
            apf_from_int(&tmp, hour * 3600 + minute * 60 + sec);
            apf_add(&result->re, &tmp2, &tmp);
            apf_zero(&result->im);
            return 1;
        }
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

        /* Try table-driven dispatch first (handles ~50 functions) */
        {
            ScalarFuncEntry *sf_entry = scalar_func_lookup(name);
            if (sf_entry) {
                return scalar_func_eval(sf_entry, result, &arg);
            }
        }

        /* Built-in functions (fallback for functions not in table) */
        if (str_eq(name, "roots2")) {
            static char buf1[256], buf2[256];
            apfc r1, r2;
            apfc_sqrt(&r1, &arg);
            apfc_neg(&r2, &r1);
            apfc_to_str(buf1, sizeof(buf1), &r1, 0);
            apfc_to_str(buf2, sizeof(buf2), &r2, 0);
            printf("roots: %s, %s\n", buf1, buf2);
            *result = r1;
        } else if (str_eq(name, "startofmonth") || str_eq(name, "som")) {
            /* startofmonth: get Unix timestamp for start of month */
            long ts, days_since_epoch, year, month, i;
            long new_days;
            static const int days_in_month[] = {31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31};
            apf secs_per_day, tmp;
            
            ts = apf_to_long(&arg.re);
            days_since_epoch = ts / 86400;
            
            /* Convert to year/month/day */
            year = 1970;
            while (1) {
                int year_days = 365;
                if ((year % 4 == 0 && year % 100 != 0) || (year % 400 == 0)) year_days = 366;
                if (days_since_epoch < year_days) break;
                days_since_epoch -= year_days;
                year++;
            }
            month = 1;
            for (i = 0; i < 12; i++) {
                int dm = days_in_month[i];
                if (i == 1 && ((year % 4 == 0 && year % 100 != 0) || (year % 400 == 0))) dm++;
                if (days_since_epoch < dm) break;
                days_since_epoch -= dm;
                month++;
            }
            
            /* Calculate days from 1970 to start of this month */
            new_days = 0;
            for (i = 1970; i < year; i++) {
                new_days += 365;
                if ((i % 4 == 0 && i % 100 != 0) || (i % 400 == 0)) new_days++;
            }
            for (i = 1; i < month; i++) {
                new_days += days_in_month[i-1];
                if (i == 2 && ((year % 4 == 0 && year % 100 != 0) || (year % 400 == 0))) new_days++;
            }
            
            apf_from_int(&result->re, new_days);
            apf_from_int(&secs_per_day, 86400);
            apf_mul(&tmp, &result->re, &secs_per_day);
            apf_copy(&result->re, &tmp);
            apf_zero(&result->im);
        } else if (str_eq(name, "startofyear") || str_eq(name, "soy")) {
            /* startofyear: get Unix timestamp for start of year */
            long ts, days_since_epoch, year, i;
            long new_days;
            apf secs_per_day, tmp;
            
            ts = apf_to_long(&arg.re);
            days_since_epoch = ts / 86400;
            
            /* Find year */
            year = 1970;
            while (1) {
                int year_days = 365;
                if ((year % 4 == 0 && year % 100 != 0) || (year % 400 == 0)) year_days = 366;
                if (days_since_epoch < year_days) break;
                days_since_epoch -= year_days;
                year++;
            }
            
            /* Calculate days from 1970 to start of this year */
            new_days = 0;
            for (i = 1970; i < year; i++) {
                new_days += 365;
                if ((i % 4 == 0 && i % 100 != 0) || (i % 400 == 0)) new_days++;
            }
            
            apf_from_int(&result->re, new_days);
            apf_from_int(&secs_per_day, 86400);
            apf_mul(&tmp, &result->re, &secs_per_day);
            apf_copy(&result->re, &tmp);
            apf_zero(&result->im);
        } else if (str_eq(name, "startofday") || str_eq(name, "sod")) {
            /* startofday: get Unix timestamp for start of day (midnight) */
            long ts, days;
            apf secs_per_day, tmp;
            
            ts = apf_to_long(&arg.re);
            days = ts / 86400;
            
            apf_from_int(&result->re, days);
            apf_from_int(&secs_per_day, 86400);
            apf_mul(&tmp, &result->re, &secs_per_day);
            apf_copy(&result->re, &tmp);
            apf_zero(&result->im);
        } else if (str_eq(name, "year")) {
            /* year: extract year from Unix timestamp */
            long ts, days_since_epoch, year;
            
            ts = apf_to_long(&arg.re);
            days_since_epoch = ts / 86400;
            
            year = 1970;
            while (1) {
                int year_days = 365;
                if ((year % 4 == 0 && year % 100 != 0) || (year % 400 == 0)) year_days = 366;
                if (days_since_epoch < year_days) break;
                days_since_epoch -= year_days;
                year++;
            }
            
            apf_from_int(&result->re, year);
            apf_zero(&result->im);
        } else if (str_eq(name, "month")) {
            /* month: extract month (1-12) from Unix timestamp */
            long ts, days_since_epoch, year, month, i;
            static const int days_in_month[] = {31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31};
            
            ts = apf_to_long(&arg.re);
            days_since_epoch = ts / 86400;
            
            year = 1970;
            while (1) {
                int year_days = 365;
                if ((year % 4 == 0 && year % 100 != 0) || (year % 400 == 0)) year_days = 366;
                if (days_since_epoch < year_days) break;
                days_since_epoch -= year_days;
                year++;
            }
            month = 1;
            for (i = 0; i < 12; i++) {
                int dm = days_in_month[i];
                if (i == 1 && ((year % 4 == 0 && year % 100 != 0) || (year % 400 == 0))) dm++;
                if (days_since_epoch < dm) break;
                days_since_epoch -= dm;
                month++;
            }
            
            apf_from_int(&result->re, month);
            apf_zero(&result->im);
        } else if (str_eq(name, "day")) {
            /* day: extract day of month (1-31) from Unix timestamp */
            long ts, days_since_epoch, year, i;
            static const int days_in_month[] = {31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31};
            
            ts = apf_to_long(&arg.re);
            days_since_epoch = ts / 86400;
            
            year = 1970;
            while (1) {
                int year_days = 365;
                if ((year % 4 == 0 && year % 100 != 0) || (year % 400 == 0)) year_days = 366;
                if (days_since_epoch < year_days) break;
                days_since_epoch -= year_days;
                year++;
            }
            for (i = 0; i < 12; i++) {
                int dm = days_in_month[i];
                if (i == 1 && ((year % 4 == 0 && year % 100 != 0) || (year % 400 == 0))) dm++;
                if (days_since_epoch < dm) break;
                days_since_epoch -= dm;
            }
            
            apf_from_int(&result->re, days_since_epoch + 1);
            apf_zero(&result->im);
        } else if (str_eq(name, "hour")) {
            /* hour: extract hour (0-23) from Unix timestamp */
            long ts = apf_to_long(&arg.re);
            apf_from_int(&result->re, (ts % 86400) / 3600);
            apf_zero(&result->im);
        } else if (str_eq(name, "minute")) {
            /* minute: extract minute (0-59) from Unix timestamp */
            long ts = apf_to_long(&arg.re);
            apf_from_int(&result->re, (ts % 3600) / 60);
            apf_zero(&result->im);
        } else if (str_eq(name, "second")) {
            /* second: extract second (0-59) from Unix timestamp */
            long ts = apf_to_long(&arg.re);
            apf_from_int(&result->re, ts % 60);
            apf_zero(&result->im);
        } else if (str_eq(name, "rng")) {
            /* rng: set random number generator seed */
            extern void mat_rand_seed(unsigned long seed);
            long seed = apf_to_long(&arg.re);
            mat_rand_seed((unsigned long)seed);
            apf_from_int(&result->re, seed);
            apf_zero(&result->im);
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
        } else if (str_eq(name, "nthprime")) {
            long n = apf_to_long(&arg.re);
            if (n < 1 || n > 100000) {
                printf("Error: nthprime argument out of range (1-100000)\n");
                return 0;
            }
            apf_from_int(&result->re, nthprime_long(n));
            apf_zero(&result->im);
        } else if (str_eq(name, "primepi")) {
            long n = apf_to_long(&arg.re);
            apf_from_int(&result->re, primepi_long(n));
            apf_zero(&result->im);
        } else if (str_eq(name, "radical")) {
            long n = apf_to_long(&arg.re);
            apf_from_int(&result->re, radical_long(n));
            apf_zero(&result->im);
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
        } else if (str_eq(name, "ispower") || str_eq(name, "isperfectpower")) {
            long n = apf_to_long(&arg.re);
            apf_from_int(&result->re, ispower_long(n));
            apf_zero(&result->im);
            result_is_boolean = 1;
#endif
#ifdef HAVE_BITWISE
#endif
        /* Duration functions - convert to seconds */
        } else if (str_eq(name, "seconds") || str_eq(name, "sec")) {
            apf_copy(&result->re, &arg.re);
            apf_zero(&result->im);
        } else if (str_eq(name, "minutes") || str_eq(name, "min")) {
            apf sixty;
            apf_from_int(&sixty, 60);
            apf_mul(&result->re, &arg.re, &sixty);
            apf_zero(&result->im);
        } else if (str_eq(name, "hours") || str_eq(name, "hr")) {
            apf factor;
            apf_from_int(&factor, 3600);
            apf_mul(&result->re, &arg.re, &factor);
            apf_zero(&result->im);
        } else if (str_eq(name, "days")) {
            apf factor;
            apf_from_int(&factor, 86400);
            apf_mul(&result->re, &arg.re, &factor);
            apf_zero(&result->im);
        } else if (str_eq(name, "milliseconds") || str_eq(name, "ms")) {
            apf factor;
            apf_from_int(&factor, 1000);
            apf_div(&result->re, &arg.re, &factor);
            apf_zero(&result->im);
        } else {
            printf("Error: unknown function '%s'\n", name);
            return 0;
        }
        return 1;
    }

    /* Matrix literal */
    if (current_token.type == TOK_LBRACKET) {
        matrix_t mat = {0, 0, NULL};
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
    
    /* Allocate matrix from arena */
    mat_zero(result, row, cols_per_row > 0 ? cols_per_row : 1);
    if (!result->data) {
        printf("Error: cannot allocate matrix\n");
        return 0;
    }
    
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
    
    /* Handle randi([lo,hi], m, n) syntax */
    if (str_eq(name, "randi") && current_token.type == TOK_LBRACKET) {
        apfc lo_val, hi_val, rows_val, cols_val;
        next_token();  /* skip [ */
        if (!parse_expr(&lo_val)) return 0;
        if (current_token.type != TOK_COMMA && current_token.type != TOK_SEMI) {
            printf("Error: randi([lo,hi],...) requires two values in range\n");
            return 0;
        }
        next_token();
        if (!parse_expr(&hi_val)) return 0;
        if (current_token.type != TOK_RBRACKET) {
            printf("Error: expected ']'\n");
            return 0;
        }
        next_token();
        
        imin = apf_to_long(&lo_val.re);
        imax = apf_to_long(&hi_val.re);
        
        /* Get dimensions */
        if (current_token.type != TOK_COMMA) {
            printf("Error: randi([lo,hi], m, n) requires dimensions\n");
            return 0;
        }
        next_token();
        if (!parse_expr(&rows_val)) return 0;
        rows = apf_to_long(&rows_val.re);
        
        if (current_token.type == TOK_COMMA) {
            next_token();
            if (!parse_expr(&cols_val)) return 0;
            cols = apf_to_long(&cols_val.re);
        } else {
            cols = 1;  /* Column vector by default */
        }
        
        if (current_token.type != TOK_RPAREN) {
            printf("Error: expected ')'\n");
            return 0;
        }
        next_token();
        
        result->type = VAL_MATRIX;
        mat_randi_range(&result->v.matrix, rows, cols, imin, imax);
        return 1;
    }
    
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
            mat_zero(&result->v.matrix, 1, n);
            if (!result->v.matrix.data) return 0;
            
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
            mat_zero(&result->v.matrix, 1, n);
            if (!result->v.matrix.data) return 0;
            
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
            matrix_t sub = {0, 0, NULL};
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
                /* Matrix indexing: M(i,j), M(1:5,:), M(:,2), etc. */
                value_t mat_val;
                int has_idx2 = 0;
                int idx1_is_colon = 0, idx2_is_colon = 0;
                int idx1_is_range = 0, idx2_is_range = 0;
                int idx1_start = 0, idx1_end = 0, idx1_step = 1;
                int idx2_start = 0, idx2_end = 0, idx2_step = 1;
                int i, j, ri, rj;
                apfc tmp;
                
                get_var_value(var_idx, &mat_val);
                next_token();
                
                /* Parse first index: :, scalar, or range (start:end or start:step:end) */
                if (current_token.type == TOK_COLON) {
                    idx1_is_colon = 1;
                    next_token();
                } else {
                    if (!parse_expr(&tmp)) return 0;
                    idx1_start = apf_to_long(&tmp.re) - 1;  /* Convert to 0-indexed */
                    idx1_end = idx1_start;
                    
                    /* Check for range: expr:end or expr:step:end */
                    if (current_token.type == TOK_COLON) {
                        next_token();
                        if (!parse_expr(&tmp)) return 0;
                        idx1_end = apf_to_long(&tmp.re) - 1;
                        idx1_is_range = 1;
                        
                        /* Check for step: expr:step:end */
                        if (current_token.type == TOK_COLON) {
                            idx1_step = idx1_end - idx1_start + 1;  /* Previous "end" was actually step */
                            next_token();
                            if (!parse_expr(&tmp)) return 0;
                            idx1_end = apf_to_long(&tmp.re) - 1;
                        }
                    }
                }
                
                /* Check for second index */
                if (current_token.type == TOK_COMMA) {
                    next_token();
                    has_idx2 = 1;
                    
                    if (current_token.type == TOK_COLON) {
                        idx2_is_colon = 1;
                        next_token();
                    } else {
                        if (!parse_expr(&tmp)) return 0;
                        idx2_start = apf_to_long(&tmp.re) - 1;
                        idx2_end = idx2_start;
                        
                        if (current_token.type == TOK_COLON) {
                            next_token();
                            if (!parse_expr(&tmp)) return 0;
                            idx2_end = apf_to_long(&tmp.re) - 1;
                            idx2_is_range = 1;
                            
                            if (current_token.type == TOK_COLON) {
                                idx2_step = idx2_end - idx2_start + 1;
                                next_token();
                                if (!parse_expr(&tmp)) return 0;
                                idx2_end = apf_to_long(&tmp.re) - 1;
                            }
                        }
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
                
                /* Set defaults for colon indices */
                if (idx1_is_colon) {
                    idx1_start = 0;
                    idx1_end = mat_val.v.matrix.rows - 1;
                    idx1_is_range = 1;
                }
                if (idx2_is_colon) {
                    idx2_start = 0;
                    idx2_end = mat_val.v.matrix.cols - 1;
                    idx2_is_range = 1;
                }
                
                /* Handle different indexing cases */
                if (has_idx2) {
                    int nrows, ncols;
                    
                    /* Validate indices */
                    if (idx1_start < 0 || idx1_end >= mat_val.v.matrix.rows ||
                        idx2_start < 0 || idx2_end >= mat_val.v.matrix.cols) {
                        printf("Error: index out of range\n");
                        return 0;
                    }
                    
                    /* Calculate result dimensions */
                    nrows = (idx1_end - idx1_start) / idx1_step + 1;
                    ncols = (idx2_end - idx2_start) / idx2_step + 1;
                    
                    if (nrows == 1 && ncols == 1 && !idx1_is_range && !idx2_is_range && !idx1_is_colon && !idx2_is_colon) {
                        /* Single element: return scalar */
                        result->type = VAL_SCALAR;
                        result->v.scalar = MAT_AT(&mat_val.v.matrix, idx1_start, idx2_start);
                    } else {
                        /* Submatrix */
                        result->type = VAL_MATRIX;
                        mat_zero(&result->v.matrix, nrows, ncols);
                        if (!result->v.matrix.data) return 0;
                        
                        ri = 0;
                        for (i = idx1_start; i <= idx1_end; i += idx1_step) {
                            rj = 0;
                            for (j = idx2_start; j <= idx2_end; j += idx2_step) {
                                MAT_AT(&result->v.matrix, ri, rj) = MAT_AT(&mat_val.v.matrix, i, j);
                                rj++;
                            }
                            ri++;
                        }
                    }
                } else {
                    /* v(i) or v(1:5) - linear indexing for vectors */
                    if (idx1_is_colon || idx1_is_range) {
                        int n = (idx1_end - idx1_start) / idx1_step + 1;
                        result->type = VAL_MATRIX;
                        
                        if (mat_val.v.matrix.rows == 1) {
                            /* Row vector */
                            mat_zero(&result->v.matrix, 1, n);
                            if (!result->v.matrix.data) return 0;
                            ri = 0;
                            for (i = idx1_start; i <= idx1_end; i += idx1_step) {
                                if (i < 0 || i >= mat_val.v.matrix.cols) {
                                    printf("Error: index out of range\n");
                                    return 0;
                                }
                                MAT_AT(&result->v.matrix, 0, ri) = MAT_AT(&mat_val.v.matrix, 0, i);
                                ri++;
                            }
                        } else {
                            /* Column vector or general matrix linear indexing */
                            mat_zero(&result->v.matrix, n, 1);
                            if (!result->v.matrix.data) return 0;
                            ri = 0;
                            for (i = idx1_start; i <= idx1_end; i += idx1_step) {
                                if (i < 0 || i >= mat_val.v.matrix.rows * mat_val.v.matrix.cols) {
                                    printf("Error: index out of range\n");
                                    return 0;
                                }
                                /* Linear indexing: column-major order */
                                MAT_AT(&result->v.matrix, ri, 0) = mat_val.v.matrix.data[i];
                                ri++;
                            }
                        }
                    } else {
                        /* Single element */
                        int idx = idx1_start;
                        if (mat_val.v.matrix.rows == 1) {
                            if (idx < 0 || idx >= mat_val.v.matrix.cols) {
                                printf("Error: index out of range\n");
                                return 0;
                            }
                            result->type = VAL_SCALAR;
                            result->v.scalar = MAT_AT(&mat_val.v.matrix, 0, idx);
                        } else if (mat_val.v.matrix.cols == 1) {
                            if (idx < 0 || idx >= mat_val.v.matrix.rows) {
                                printf("Error: index out of range\n");
                                return 0;
                            }
                            result->type = VAL_SCALAR;
                            result->v.scalar = MAT_AT(&mat_val.v.matrix, idx, 0);
                        } else {
                            /* Linear indexing into 2D matrix */
                            if (idx < 0 || idx >= mat_val.v.matrix.rows * mat_val.v.matrix.cols) {
                                printf("Error: index out of range\n");
                                return 0;
                            }
                            result->type = VAL_SCALAR;
                            result->v.scalar = mat_val.v.matrix.data[idx];
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
        /* Check for table.column access syntax or table display */
        {
            table_t *tbl;
            
            /* Check if this is a table */
            
            tbl = table_get(name);
            
            if (tbl) {
                if (*input_ptr == '.') {
                    char col_name[32];
                    int ni = 0;
                    matrix_t *col;
                    
                    input_ptr++;  /* skip . */
                    
                    /* Read column name */
                    while ((*input_ptr >= 'a' && *input_ptr <= 'z') ||
                           (*input_ptr >= 'A' && *input_ptr <= 'Z') ||
                           (*input_ptr >= '0' && *input_ptr <= '9') ||
                           *input_ptr == '_') {
                        if (ni < 31) col_name[ni++] = *input_ptr;
                        input_ptr++;
                    }
                    col_name[ni] = '\0';
                    next_token();
                    
                    col = table_get_column(tbl, col_name);
                    if (col) {
                        result->type = VAL_MATRIX;
                        mat_copy(&result->v.matrix, col);
                        return 1;
                    } else {
                        printf("Error: unknown column '%s' in table '%s'\n", col_name, name);
                        return 0;
                    }
                } else {
                    /* Table without dot - display it */
                    extern void table_print(const table_t *t);
                    next_token();
                    table_print(tbl);
                    /* Return empty result */
                    result->type = VAL_SCALAR;
                    apf_zero(&result->v.scalar.re);
                    apf_zero(&result->v.scalar.im);
                    return 1;
                }
            }
        }
        
        /* Check for named variable (multi-character) with possible indexing */
        {
            value_t named_val;
            if (get_named_var(name, &named_val)) {
                next_token();
                if (current_token.type == TOK_LPAREN && named_val.type == VAL_MATRIX) {
                    /* Matrix indexing for named variable: M(i,j), M(1:5,:), M(:,2), etc. */
                    int has_idx2 = 0;
                    int idx1_is_colon = 0, idx2_is_colon = 0;
                    int idx1_is_range = 0, idx2_is_range = 0;
                    int idx1_start = 0, idx1_end = 0, idx1_step = 1;
                    int idx2_start = 0, idx2_end = 0, idx2_step = 1;
                    int i, j, ri, rj;
                    apfc tmp;
                    
                    next_token();
                    
                    /* Parse first index */
                    if (current_token.type == TOK_COLON) {
                        idx1_is_colon = 1;
                        next_token();
                    } else {
                        if (!parse_expr(&tmp)) return 0;
                        idx1_start = apf_to_long(&tmp.re) - 1;
                        idx1_end = idx1_start;
                        
                        if (current_token.type == TOK_COLON) {
                            next_token();
                            if (!parse_expr(&tmp)) return 0;
                            idx1_end = apf_to_long(&tmp.re) - 1;
                            idx1_is_range = 1;
                            
                            if (current_token.type == TOK_COLON) {
                                idx1_step = idx1_end - idx1_start + 1;
                                next_token();
                                if (!parse_expr(&tmp)) return 0;
                                idx1_end = apf_to_long(&tmp.re) - 1;
                            }
                        }
                    }
                    
                    /* Parse second index */
                    if (current_token.type == TOK_COMMA) {
                        next_token();
                        has_idx2 = 1;
                        
                        if (current_token.type == TOK_COLON) {
                            idx2_is_colon = 1;
                            next_token();
                        } else {
                            if (!parse_expr(&tmp)) return 0;
                            idx2_start = apf_to_long(&tmp.re) - 1;
                            idx2_end = idx2_start;
                            
                            if (current_token.type == TOK_COLON) {
                                next_token();
                                if (!parse_expr(&tmp)) return 0;
                                idx2_end = apf_to_long(&tmp.re) - 1;
                                idx2_is_range = 1;
                                
                                if (current_token.type == TOK_COLON) {
                                    idx2_step = idx2_end - idx2_start + 1;
                                    next_token();
                                    if (!parse_expr(&tmp)) return 0;
                                    idx2_end = apf_to_long(&tmp.re) - 1;
                                }
                            }
                        }
                    }
                    
                    if (current_token.type != TOK_RPAREN) {
                        printf("Error: expected ')'\n");
                        return 0;
                    }
                    next_token();
                    
                    /* Set defaults for colon indices */
                    if (idx1_is_colon) {
                        idx1_start = 0;
                        idx1_end = named_val.v.matrix.rows - 1;
                        idx1_is_range = 1;
                    }
                    if (idx2_is_colon) {
                        idx2_start = 0;
                        idx2_end = named_val.v.matrix.cols - 1;
                        idx2_is_range = 1;
                    }
                    
                    /* Handle different indexing cases */
                    if (has_idx2) {
                        int nrows, ncols;
                        
                        if (idx1_start < 0 || idx1_end >= named_val.v.matrix.rows ||
                            idx2_start < 0 || idx2_end >= named_val.v.matrix.cols) {
                            printf("Error: index out of range\n");
                            return 0;
                        }
                        
                        nrows = (idx1_end - idx1_start) / idx1_step + 1;
                        ncols = (idx2_end - idx2_start) / idx2_step + 1;
                        
                        if (nrows == 1 && ncols == 1 && !idx1_is_range && !idx2_is_range && !idx1_is_colon && !idx2_is_colon) {
                            result->type = VAL_SCALAR;
                            result->v.scalar = MAT_AT(&named_val.v.matrix, idx1_start, idx2_start);
                        } else {
                            result->type = VAL_MATRIX;
                            mat_zero(&result->v.matrix, nrows, ncols);
                            if (!result->v.matrix.data) return 0;
                            
                            ri = 0;
                            for (i = idx1_start; i <= idx1_end; i += idx1_step) {
                                rj = 0;
                                for (j = idx2_start; j <= idx2_end; j += idx2_step) {
                                    MAT_AT(&result->v.matrix, ri, rj) = MAT_AT(&named_val.v.matrix, i, j);
                                    rj++;
                                }
                                ri++;
                            }
                        }
                    } else {
                        /* Single index or range */
                        if (idx1_is_colon || idx1_is_range) {
                            int n = (idx1_end - idx1_start) / idx1_step + 1;
                            result->type = VAL_MATRIX;
                            
                            if (named_val.v.matrix.rows == 1) {
                                mat_zero(&result->v.matrix, 1, n);
                                if (!result->v.matrix.data) return 0;
                                ri = 0;
                                for (i = idx1_start; i <= idx1_end; i += idx1_step) {
                                    if (i < 0 || i >= named_val.v.matrix.cols) {
                                        printf("Error: index out of range\n");
                                        return 0;
                                    }
                                    MAT_AT(&result->v.matrix, 0, ri) = MAT_AT(&named_val.v.matrix, 0, i);
                                    ri++;
                                }
                            } else {
                                mat_zero(&result->v.matrix, n, 1);
                                if (!result->v.matrix.data) return 0;
                                ri = 0;
                                for (i = idx1_start; i <= idx1_end; i += idx1_step) {
                                    if (i < 0 || i >= named_val.v.matrix.rows * named_val.v.matrix.cols) {
                                        printf("Error: index out of range\n");
                                        return 0;
                                    }
                                    MAT_AT(&result->v.matrix, ri, 0) = named_val.v.matrix.data[i];
                                    ri++;
                                }
                            }
                        } else {
                            int idx = idx1_start;
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
            /* datetime({'string';'string';...}) - cell array of datetime strings */
            if (str_eq(name, "datetime")) {
                extern int datetime_parse(apf *result, const char *str);
                const char *check_ptr = input_ptr;
                
                /* Skip whitespace and check for { */
                while (*check_ptr == ' ') check_ptr++;
                if (*check_ptr == '{') {
                    /* This is a cell array of strings */
                    int count = 0;
                    int max_times = 100;
                    apf *time_vals;
                    int i;
                    
                    next_token();  /* skip ( */
                    input_ptr++;   /* skip { */
                    next_token();
                    
                    time_vals = (apf *)mat_arena_alloc(max_times * sizeof(apf));
                    if (!time_vals) {
                        printf("Error: memory allocation failed\n");
                        return 0;
                    }
                    
                    /* Parse strings */
                    while (count < max_times) {
                        if (current_token.type == TOK_STRING) {
                            if (!datetime_parse(&time_vals[count], current_token.str_value)) {
                                printf("Error: cannot parse datetime '%s'\n", current_token.str_value);
                                return 0;
                            }
                            count++;
                            next_token();
                            
                            /* Check for separator */
                            if (current_token.type == TOK_SEMI || current_token.type == TOK_COMMA) {
                                next_token();
                            } else {
                                break;
                            }
                        } else {
                            break;
                        }
                    }
                    
                    /* Skip } */
                    while (*input_ptr == ' ') input_ptr++;
                    if (*input_ptr == '}') {
                        input_ptr++;
                        next_token();
                    }
                    
                    /* Skip ) */
                    if (current_token.type == TOK_RPAREN) {
                        next_token();
                    }
                    
                    if (count == 0) {
                        printf("Error: no datetime strings found\n");
                        return 0;
                    }
                    
                    /* Return as column vector */
                    result->type = VAL_MATRIX;
                    mat_zero(&result->v.matrix, count, 1);
                    for (i = 0; i < count; i++) {
                        apf_copy(&MAT_AT(&result->v.matrix, i, 0).re, &time_vals[i]);
                        apf_zero(&MAT_AT(&result->v.matrix, i, 0).im);
                    }
                    return 1;
                }
            }
            
            /* categorical(["str"; "str"; ...]) - creates categorical array */
            if (str_eq(name, "categorical")) {
                const char *check_ptr = input_ptr;
                
                /* Skip whitespace and check for [ */
                while (*check_ptr == ' ') check_ptr++;
                if (*check_ptr == '[') {
                    /* This is a string array */
                    static char cat_labels[256][32];
                    int count = 0;
                    int i;
                    
                    next_token();  /* skip ( */
                    input_ptr++;   /* skip [ */
                    next_token();
                    
                    /* Parse strings */
                    while (count < 256) {
                        if (current_token.type == TOK_STRING) {
                            size_t slen = strlen(current_token.str_value);
                            if (slen > 31) slen = 31;
                            memcpy(cat_labels[count], current_token.str_value, slen);
                            cat_labels[count][slen] = '\0';
                            count++;
                            next_token();
                            
                            /* Check for separator */
                            if (current_token.type == TOK_SEMI || current_token.type == TOK_COMMA) {
                                next_token();
                            } else {
                                break;
                            }
                        } else {
                            break;
                        }
                    }
                    
                    /* Skip ] */
                    while (*input_ptr == ' ') input_ptr++;
                    if (*input_ptr == ']') {
                        input_ptr++;
                        next_token();
                    }
                    
                    /* Skip ) */
                    if (current_token.type == TOK_RPAREN) {
                        next_token();
                    }
                    
                    if (count == 0) {
                        printf("Error: no strings found in categorical\n");
                        return 0;
                    }
                    
                    /* Store as special categorical result - use indices as values */
                    /* Return indices as column vector for now */
                    result->type = VAL_MATRIX;
                    mat_zero(&result->v.matrix, count, 1);
                    
                    /* Build category lookup and store indices */
                    {
                        static char unique_cats[64][32];
                        int n_unique = 0;
                        
                        for (i = 0; i < count; i++) {
                            int idx = -1, j;
                            /* Find or add category */
                            for (j = 0; j < n_unique; j++) {
                                if (strcmp(unique_cats[j], cat_labels[i]) == 0) {
                                    idx = j;
                                    break;
                                }
                            }
                            if (idx < 0) {
                                size_t clen = strlen(cat_labels[i]);
                                if (clen > 31) clen = 31;
                                idx = n_unique;
                                memcpy(unique_cats[n_unique], cat_labels[i], clen);
                                unique_cats[n_unique][clen] = '\0';
                                n_unique++;
                            }
                            apf_from_int(&MAT_AT(&result->v.matrix, i, 0).re, idx + 1);
                            apf_zero(&MAT_AT(&result->v.matrix, i, 0).im);
                        }
                        
                        /* Print categories */
                        printf("Categories: ");
                        for (i = 0; i < n_unique; i++) {
                            printf("%s", unique_cats[i]);
                            if (i < n_unique - 1) printf(", ");
                        }
                        printf("\n");
                    }
                    return 1;
                }
            }
            
            /* table(var1, var2, ...) - create table from variables */
            if (str_eq(name, "table")) {
                extern int get_var_value(int idx, value_t *val);
                extern int get_named_var(const char *name, value_t *val);
                
                static char col_names[16][32];
                static value_t col_vals[16];
                static categorical_t *col_cats[16];  /* Pointers to categorical data */
                static int col_is_cat[16];  /* Is this column categorical? */
                int n_cols = 0;
                int n_rows = 0;
                int i, j;
                
                next_token();  /* skip ( */
                
                /* Parse variable names */
                while (n_cols < 16) {
                    if (current_token.type == TOK_FUNC) {
                        char vname[32];
                        size_t vlen = strlen(current_token.func_name);
                        void *cat;
                        
                        if (vlen > 31) vlen = 31;
                        memcpy(vname, current_token.func_name, vlen);
                        vname[vlen] = '\0';
                        memcpy(col_names[n_cols], vname, vlen);
                        col_names[n_cols][vlen] = '\0';
                        col_is_cat[n_cols] = 0;
                        col_cats[n_cols] = NULL;
                        
                        /* First check if it's a categorical */
                        cat = cat_get(vname);
                        if (cat) {
                            col_is_cat[n_cols] = 1;
                            col_cats[n_cols] = cat;
                            col_vals[n_cols].type = VAL_SCALAR;  /* Mark as handled */
                            if (n_rows == 0) n_rows = cat_num_elements(cat);
                            n_cols++;
                        }
                        /* Try to get numeric variable value */
                        else if (strlen(vname) == 1 && vname[0] >= 'a' && vname[0] <= 'z') {
                            if (get_var_value(vname[0] - 'a', &col_vals[n_cols])) {
                                if (col_vals[n_cols].type == VAL_MATRIX) {
                                    if (n_rows == 0) n_rows = col_vals[n_cols].v.matrix.rows;
                                }
                                n_cols++;
                            } else {
                                printf("Error: variable '%s' not found\n", vname);
                                return 0;
                            }
                        } else if (get_named_var(vname, &col_vals[n_cols])) {
                            if (col_vals[n_cols].type == VAL_MATRIX) {
                                if (n_rows == 0) n_rows = col_vals[n_cols].v.matrix.rows;
                            }
                            n_cols++;
                        } else {
                            printf("Error: variable '%s' not found\n", vname);
                            return 0;
                        }
                        
                        next_token();
                        if (current_token.type == TOK_COMMA) {
                            next_token();
                        } else {
                            break;
                        }
                    } else {
                        break;
                    }
                }
                
                if (current_token.type == TOK_RPAREN) {
                    next_token();
                }
                
                if (n_cols == 0) {
                    printf("Error: table requires at least one variable\n");
                    return 0;
                }
                
                /* Print table header */
                printf("\nT=%dx%d table\n", n_rows, n_cols);
                
                /* Calculate column widths */
                {
                    int widths[16];
                    char buf[64];
                    
                    for (j = 0; j < n_cols; j++) {
                        widths[j] = (int)strlen(col_names[j]);
                        if (col_is_cat[j] && col_cats[j]) {
                            int cw = cat_max_width(col_cats[j]);
                            if (cw > widths[j]) widths[j] = cw;
                        } else if (col_vals[j].type == VAL_MATRIX) {
                            for (i = 0; i < col_vals[j].v.matrix.rows; i++) {
                                int len;
                                apf_to_str(buf, sizeof(buf), &MAT_AT(&col_vals[j].v.matrix, i, 0).re, 6);
                                len = (int)strlen(buf);
                                if (len > widths[j]) widths[j] = len;
                            }
                        }
                        if (widths[j] < 6) widths[j] = 6;
                        widths[j] += 2;
                    }
                    
                    /* Print header */
                    printf("    ");
                    for (j = 0; j < n_cols; j++) {
                        printf("%*s", widths[j], col_names[j]);
                    }
                    printf("\n    ");
                    for (j = 0; j < n_cols; j++) {
                        for (i = 0; i < widths[j]; i++) printf("_");
                    }
                    printf("\n\n");
                    
                    /* Print rows */
                    for (i = 0; i < n_rows; i++) {
                        printf("    ");
                        for (j = 0; j < n_cols; j++) {
                            if (col_is_cat[j] && col_cats[j]) {
                                const char *lbl = cat_get_label(col_cats[j], i);
                                printf("%*s", widths[j], lbl ? lbl : "");
                            } else if (col_vals[j].type == VAL_MATRIX && i < col_vals[j].v.matrix.rows) {
                                apf_to_str(buf, sizeof(buf), &MAT_AT(&col_vals[j].v.matrix, i, 0).re, 6);
                                printf("%*s", widths[j], buf);
                            } else if (col_vals[j].type == VAL_SCALAR) {
                                apf_to_str(buf, sizeof(buf), &col_vals[j].v.scalar.re, 6);
                                printf("%*s", widths[j], buf);
                            }
                        }
                        printf("\n");
                    }
                }
                
                /* Return empty for now - table is a display function */
                result->type = VAL_SCALAR;
                apf_from_int(&result->v.scalar.re, n_rows);
                apf_zero(&result->v.scalar.im);
                return 1;
            }
            
            if (str_eq(name, "zeros") || str_eq(name, "ones") || 
                str_eq(name, "eye") || str_eq(name, "identity") ||
                str_eq(name, "rand") || str_eq(name, "randn") || str_eq(name, "randi") ||
                str_eq(name, "linspace") || str_eq(name, "logspace")) {
                input_ptr = saved_pos;
                current_token = saved_tok;
                next_token();
                return parse_matrix_func(result, name);
            }
            
            /* MATLAB-style: mean/sum/std/var on matrix returns column-wise results */
            if (str_eq(name, "mean") || str_eq(name, "sum") || 
                str_eq(name, "std") || str_eq(name, "sd") || str_eq(name, "var")) {
                next_token();
                if (!parse_value(&pv_arg)) return 0;
                if (current_token.type != TOK_RPAREN) {
                    printf("Error: expected ')'\n");
                    return 0;
                }
                next_token();
                
                if (pv_arg.type == VAL_MATRIX && pv_arg.v.matrix.cols > 1 && pv_arg.v.matrix.rows > 1) {
                    /* Multi-column matrix: return column-wise results */
                    extern void mat_mean_cols(matrix_t *r, const matrix_t *m);
                    extern void mat_std_cols(matrix_t *r, const matrix_t *m);
                    extern void mat_var_cols(matrix_t *r, const matrix_t *m);
                    extern void mat_col_sums(matrix_t *r, const matrix_t *m);
                    
                    result->type = VAL_MATRIX;
                    if (str_eq(name, "mean")) {
                        mat_mean_cols(&result->v.matrix, &pv_arg.v.matrix);
                    } else if (str_eq(name, "sum")) {
                        mat_col_sums(&result->v.matrix, &pv_arg.v.matrix);
                    } else if (str_eq(name, "std") || str_eq(name, "sd")) {
                        mat_std_cols(&result->v.matrix, &pv_arg.v.matrix);
                    } else if (str_eq(name, "var")) {
                        mat_var_cols(&result->v.matrix, &pv_arg.v.matrix);
                    }
                    return 1;
                } else {
                    /* Vector or scalar: return overall result */
                    extern void mat_mean(apfc *r, const matrix_t *m);
                    extern void mat_sum(apfc *r, const matrix_t *m);
                    extern void mat_std(apfc *r, const matrix_t *m);
                    extern void mat_var(apfc *r, const matrix_t *m);
                    
                    result->type = VAL_SCALAR;
                    if (pv_arg.type == VAL_SCALAR) {
                        if (str_eq(name, "mean") || str_eq(name, "sum")) {
                            result->v.scalar = pv_arg.v.scalar;
                        } else {
                            /* std/var of scalar is 0 */
                            apf_zero(&result->v.scalar.re);
                            apf_zero(&result->v.scalar.im);
                        }
                    } else {
                        if (str_eq(name, "mean")) {
                            mat_mean(&result->v.scalar, &pv_arg.v.matrix);
                        } else if (str_eq(name, "sum")) {
                            mat_sum(&result->v.scalar, &pv_arg.v.matrix);
                        } else if (str_eq(name, "std") || str_eq(name, "sd")) {
                            mat_std(&result->v.scalar, &pv_arg.v.matrix);
                        } else if (str_eq(name, "var")) {
                            mat_var(&result->v.scalar, &pv_arg.v.matrix);
                        }
                    }
                    return 1;
                }
            }
            
            /* MATLAB-style: mean/sum/std/var on multi-column matrix returns column-wise results */
            if (str_eq(name, "mean") || str_eq(name, "sum") || str_eq(name, "std") || 
                str_eq(name, "sd") || str_eq(name, "var")) {
                extern void mat_mean_cols(matrix_t *r, const matrix_t *m);
                extern void mat_col_sums(matrix_t *r, const matrix_t *m);
                extern void mat_std_cols(matrix_t *r, const matrix_t *m);
                extern void mat_var_cols(matrix_t *r, const matrix_t *m);
                
                next_token();
                if (!parse_value(&pv_arg)) return 0;
                if (current_token.type != TOK_RPAREN) {
                    printf("Error: expected ')'\n");
                    return 0;
                }
                next_token();
                
                if (pv_arg.type == VAL_MATRIX && pv_arg.v.matrix.cols > 1 && pv_arg.v.matrix.rows > 1) {
                    /* Multi-column matrix: return column-wise results */
                    result->type = VAL_MATRIX;
                    if (str_eq(name, "mean")) {
                        mat_mean_cols(&result->v.matrix, &pv_arg.v.matrix);
                    } else if (str_eq(name, "sum")) {
                        mat_col_sums(&result->v.matrix, &pv_arg.v.matrix);
                    } else if (str_eq(name, "std") || str_eq(name, "sd")) {
                        mat_std_cols(&result->v.matrix, &pv_arg.v.matrix);
                    } else if (str_eq(name, "var")) {
                        mat_var_cols(&result->v.matrix, &pv_arg.v.matrix);
                    }
                    return 1;
                } else if (pv_arg.type == VAL_MATRIX) {
                    /* Vector or single-column: return scalar */
                    result->type = VAL_SCALAR;
                    if (str_eq(name, "mean")) {
                        extern void mat_mean(apfc *r, const matrix_t *m);
                        mat_mean(&result->v.scalar, &pv_arg.v.matrix);
                    } else if (str_eq(name, "sum")) {
                        extern void mat_sum(apfc *r, const matrix_t *m);
                        mat_sum(&result->v.scalar, &pv_arg.v.matrix);
                    } else if (str_eq(name, "std") || str_eq(name, "sd")) {
                        extern void mat_std(apfc *r, const matrix_t *m);
                        mat_std(&result->v.scalar, &pv_arg.v.matrix);
                    } else if (str_eq(name, "var")) {
                        extern void mat_var(apfc *r, const matrix_t *m);
                        mat_var(&result->v.scalar, &pv_arg.v.matrix);
                    }
                    return 1;
                } else {
                    /* Scalar input */
                    result->type = VAL_SCALAR;
                    result->v.scalar = pv_arg.v.scalar;
                    if (str_eq(name, "std") || str_eq(name, "sd") || str_eq(name, "var")) {
                        apf_zero(&result->v.scalar.re);
                        apf_zero(&result->v.scalar.im);
                    }
                    return 1;
                }
            }
            
            /* Try table-driven matrix dispatch first */
            {
                MatFuncEntry *mf_entry = matrix_func_lookup(name);
                if (mf_entry) {
                    mf_value_t mf_result, mf_arg;
                    next_token();
                    if (!parse_value(&pv_arg)) return 0;
                    if (current_token.type != TOK_RPAREN) {
                        printf("Error: expected ')'\n");
                        return 0;
                    }
                    next_token();
                    
                    /* Convert to mf_value_t */
                    mf_arg.type = pv_arg.type;
                    if (pv_arg.type == VAL_SCALAR) {
                        mf_arg.v.scalar = pv_arg.v.scalar;
                    } else {
                        mat_copy(&mf_arg.v.matrix, &pv_arg.v.matrix);
                    }
                    
                    if (!matrix_func_eval1(mf_entry, &mf_result, &mf_arg)) {
                        return 0;
                    }
                    
                    /* Convert result back to value_t */
                    result->type = mf_result.type;
                    if (mf_result.type == MF_VAL_SCALAR) {
                        result->v.scalar = mf_result.v.scalar;
                    } else {
                        mat_copy(&result->v.matrix, &mf_result.v.matrix);
                    }
                    return 1;
                }
            }
            
            /* diag has dual behavior: creates matrix from vector or extracts diagonal */
            if (str_eq(name, "diag")) {
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
                if (mat_is_vector(&pv_arg.v.matrix)) {
                    mat_diag_create(&result->v.matrix, &pv_arg.v.matrix);
                } else {
                    mat_diag_extract(&result->v.matrix, &pv_arg.v.matrix);
                }
                return 1;
            }
            
            /* ============================================================
             * DECOMPOSITION FUNCTIONS - use table dispatch
             * Single-return calls like svd(A), qr(A), lu(A), chol(A), eig(A)
             * Multi-return [U,S,V]=svd(A) handled in main.c
             * ============================================================ */
            {
                DecompFuncEntry *dfe = decomp_func_lookup(name);
                if (dfe) {
                    DecompResult dresult;
                    matrix_t *input2 = NULL;
                    value_t arg2;
                    
                    next_token();
                    if (!parse_value(&pv_arg)) return 0;
                    
                    /* Ensure input is matrix */
                    if (pv_arg.type == VAL_SCALAR) {
                        mat_zero(&pv_arg.v.matrix, 1, 1);
                        MAT_AT(&pv_arg.v.matrix, 0, 0) = pv_arg.v.scalar;
                    }
                    
                    /* Check for second argument (e.g., kmeans(X, k)) */
                    if (current_token.type == TOK_COMMA && dfe->max_inputs > 1) {
                        next_token();
                        if (!parse_value(&arg2)) return 0;
                        if (arg2.type == VAL_SCALAR) {
                            mat_zero(&arg2.v.matrix, 1, 1);
                            MAT_AT(&arg2.v.matrix, 0, 0) = arg2.v.scalar;
                        }
                        input2 = &arg2.v.matrix;
                    }
                    
                    if (current_token.type != TOK_RPAREN) {
                        printf("Error: expected ')'\n");
                        return 0;
                    }
                    next_token();
                    
                    /* Execute with single output */
                    if (!decomp_func_exec(name, 1, &dresult, &pv_arg.v.matrix, input2)) {
                        printf("Error: %s failed\n", name);
                        return 0;
                    }
                    
                    result->type = VAL_MATRIX;
                    mat_copy(&result->v.matrix, &dresult.outputs[0]);
                    return 1;
                }
            }
            
            /* PCA reduce: pcareduce(X, k) */
            if (str_eq(name, "pcareduce")) {
                value_t k_arg;
                int k;
                
                next_token();
                if (!parse_value(&pv_arg)) return 0;
                if (current_token.type != TOK_COMMA) {
                    printf("Error: expected ',' in pcareduce(X, k)\n");
                    return 0;
                }
                next_token();
                if (!parse_value(&k_arg)) return 0;
                if (current_token.type != TOK_RPAREN) {
                    printf("Error: expected ')'\n");
                    return 0;
                }
                next_token();
                
                if (pv_arg.type != VAL_MATRIX) {
                    printf("Error: pcareduce requires matrix argument\n");
                    return 0;
                }
                
                k = (k_arg.type == VAL_SCALAR) ? apf_to_long(&k_arg.v.scalar.re) : 2;
                
                if (!mat_pca_reduce(&result->v.matrix, &pv_arg.v.matrix, k)) {
                    printf("Error: PCA reduce failed\n");
                    return 0;
                }
                result->type = VAL_MATRIX;
                return 1;
            }
            
            /* K-means: kmeans(X, k) returns cluster indices */
            if (str_eq(name, "kmeans")) {
                value_t k_arg;
                matrix_t idx = {0,0,NULL}, centroids = {0,0,NULL};
                int k;
                
                next_token();
                if (!parse_value(&pv_arg)) return 0;
                if (current_token.type != TOK_COMMA) {
                    printf("Error: expected ',' in kmeans(X, k)\n");
                    return 0;
                }
                next_token();
                if (!parse_value(&k_arg)) return 0;
                if (current_token.type != TOK_RPAREN) {
                    printf("Error: expected ')'\n");
                    return 0;
                }
                next_token();
                
                if (pv_arg.type != VAL_MATRIX) {
                    printf("Error: kmeans requires matrix argument\n");
                    return 0;
                }
                
                k = (k_arg.type == VAL_SCALAR) ? apf_to_long(&k_arg.v.scalar.re) : 2;
                
                if (!mat_kmeans(&idx, &centroids, &pv_arg.v.matrix, k)) {
                    printf("Error: kmeans failed\n");
                    return 0;
                }
                
                mat_copy(&result->v.matrix, &idx);
                result->type = VAL_MATRIX;
                return 1;
            }
            
            /* Pairwise distance: pdist(X) */
            if (str_eq(name, "pdist")) {
                next_token();
                if (!parse_value(&pv_arg)) return 0;
                if (current_token.type != TOK_RPAREN) {
                    printf("Error: expected ')'\n");
                    return 0;
                }
                next_token();
                
                if (pv_arg.type != VAL_MATRIX) {
                    printf("Error: pdist requires matrix argument\n");
                    return 0;
                }
                
                if (!mat_pdist(&result->v.matrix, &pv_arg.v.matrix)) {
                    printf("Error: pdist failed\n");
                    return 0;
                }
                result->type = VAL_MATRIX;
                return 1;
            }
            
            /* Silhouette score: silhouette(X, idx) */
            if (str_eq(name, "silhouette")) {
                value_t idx_arg;
                
                next_token();
                if (!parse_value(&pv_arg)) return 0;
                if (current_token.type != TOK_COMMA) {
                    printf("Error: expected ',' in silhouette(X, idx)\n");
                    return 0;
                }
                next_token();
                if (!parse_value(&idx_arg)) return 0;
                if (current_token.type != TOK_RPAREN) {
                    printf("Error: expected ')'\n");
                    return 0;
                }
                next_token();
                
                if (pv_arg.type != VAL_MATRIX || idx_arg.type != VAL_MATRIX) {
                    printf("Error: silhouette requires matrix arguments\n");
                    return 0;
                }
                
                mat_silhouette(&result->v.scalar, &pv_arg.v.matrix, &idx_arg.v.matrix);
                result->type = VAL_SCALAR;
                return 1;
            }
            
            /* Covariance matrix: cov(X) */
            if (str_eq(name, "cov")) {
                next_token();
                if (!parse_value(&pv_arg)) return 0;
                if (current_token.type != TOK_RPAREN) {
                    printf("Error: expected ')'\n");
                    return 0;
                }
                next_token();
                
                if (pv_arg.type != VAL_MATRIX) {
                    printf("Error: cov requires matrix argument\n");
                    return 0;
                }
                mat_cov(&result->v.matrix, &pv_arg.v.matrix);
                result->type = VAL_MATRIX;
                return 1;
            }
            
            /* Correlation matrix: corr(X), corrcoef(X) or corrcoef(X,Y) */
            if (str_eq(name, "corr") || str_eq(name, "corrcoef")) {
                next_token();
                if (!parse_value(&pv_arg)) return 0;
                
                /* Check for two-argument corrcoef(X, Y) */
                if (current_token.type == TOK_COMMA && str_eq(name, "corrcoef")) {
                    value_t arg2;
                    apf mean_x, mean_y, sum_x, sum_y, sum_xy, sum_xx, sum_yy;
                    apf diff_x, diff_y, prod, n_val;
                    int ii, n;
                    
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
                    for (ii = 0; ii < n; ii++) {
                        apf_add(&sum_x, &sum_x, &pv_arg.v.matrix.data[ii].re);
                        apf_add(&sum_y, &sum_y, &arg2.v.matrix.data[ii].re);
                    }
                    apf_from_int(&n_val, n);
                    apf_div(&mean_x, &sum_x, &n_val);
                    apf_div(&mean_y, &sum_y, &n_val);
                    
                    /* Compute correlation */
                    apf_zero(&sum_xy);
                    apf_zero(&sum_xx);
                    apf_zero(&sum_yy);
                    for (ii = 0; ii < n; ii++) {
                        apf_sub(&diff_x, &pv_arg.v.matrix.data[ii].re, &mean_x);
                        apf_sub(&diff_y, &arg2.v.matrix.data[ii].re, &mean_y);
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
                        apf_zero(&result->v.scalar.re);
                    } else {
                        apf_div(&result->v.scalar.re, &sum_xy, &prod);
                    }
                    apf_zero(&result->v.scalar.im);
                    result->type = VAL_SCALAR;
                    return 1;
                }
                
                if (current_token.type != TOK_RPAREN) {
                    printf("Error: expected ')'\n");
                    return 0;
                }
                next_token();
                
                if (pv_arg.type != VAL_MATRIX) {
                    printf("Error: corr requires matrix argument\n");
                    return 0;
                }
                mat_corrcoef(&result->v.matrix, &pv_arg.v.matrix);
                result->type = VAL_MATRIX;
                return 1;
            }
            
            /* K-means: kmeans(X, k) */
            if (str_eq(name, "kmeans")) {
                value_t k_arg;
                int k;
                next_token();
                if (!parse_value(&pv_arg)) return 0;
                if (current_token.type != TOK_COMMA) {
                    printf("Error: expected ',' in kmeans(X, k)\n");
                    return 0;
                }
                next_token();
                if (!parse_value(&k_arg)) return 0;
                if (current_token.type != TOK_RPAREN) {
                    printf("Error: expected ')'\n");
                    return 0;
                }
                next_token();
                
                if (pv_arg.type != VAL_MATRIX) {
                    printf("Error: kmeans requires matrix argument\n");
                    return 0;
                }
                k = (k_arg.type == VAL_SCALAR) ? (int)apf_to_double(&k_arg.v.scalar.re) : 3;
                {
                    matrix_t centroids_tmp = {0, 0, NULL};
                    mat_kmeans(&result->v.matrix, &centroids_tmp, &pv_arg.v.matrix, k);
                }
                result->type = VAL_MATRIX;
                return 1;
            }
            
            /* Cross-tabulation: crosstab(a, b) */
            if (str_eq(name, "crosstab")) {
                value_t b_arg;
                next_token();
                if (!parse_value(&pv_arg)) return 0;
                if (current_token.type != TOK_COMMA) {
                    printf("Error: expected ',' in crosstab(a, b)\n");
                    return 0;
                }
                next_token();
                if (!parse_value(&b_arg)) return 0;
                if (current_token.type != TOK_RPAREN) {
                    printf("Error: expected ')'\n");
                    return 0;
                }
                next_token();
                
                if (pv_arg.type != VAL_MATRIX || b_arg.type != VAL_MATRIX) {
                    printf("Error: crosstab requires matrix arguments\n");
                    return 0;
                }
                mat_crosstab(&result->v.matrix, &pv_arg.v.matrix, &b_arg.v.matrix);
                result->type = VAL_MATRIX;
                return 1;
            }
            
            /* Rand Index: randindex(a, b) */
            if (str_eq(name, "randindex")) {
                value_t b_arg;
                next_token();
                if (!parse_value(&pv_arg)) return 0;
                if (current_token.type != TOK_COMMA) {
                    printf("Error: expected ',' in randindex(a, b)\n");
                    return 0;
                }
                next_token();
                if (!parse_value(&b_arg)) return 0;
                if (current_token.type != TOK_RPAREN) {
                    printf("Error: expected ')'\n");
                    return 0;
                }
                next_token();
                
                if (pv_arg.type != VAL_MATRIX || b_arg.type != VAL_MATRIX) {
                    printf("Error: randindex requires matrix arguments\n");
                    return 0;
                }
                mat_randindex(&result->v.scalar, &pv_arg.v.matrix, &b_arg.v.matrix);
                result->type = VAL_SCALAR;
                return 1;
            }
            
            /* ML: zscore(X) - standardize data */
            if (str_eq(name, "zscore")) {
                next_token();
                if (!parse_value(&pv_arg)) return 0;
                if (current_token.type != TOK_RPAREN) {
                    printf("Error: expected ')'\n");
                    return 0;
                }
                next_token();
                
                if (pv_arg.type != VAL_MATRIX) {
                    printf("Error: zscore requires matrix argument\n");
                    return 0;
                }
                mat_zscore(&result->v.matrix, &pv_arg.v.matrix);
                result->type = VAL_MATRIX;
                return 1;
            }
            
            /* ML: fitcknn(X, Y) or fitcknn(X, Y, k) - train KNN classifier */
            if (str_eq(name, "fitcknn")) {
                value_t y_arg;
                int k = 5;  /* default */
                int model_id;
                
                next_token();
                if (!parse_value(&pv_arg)) return 0;
                if (current_token.type != TOK_COMMA) {
                    printf("Error: expected ',' in fitcknn(X, Y)\n");
                    return 0;
                }
                next_token();
                if (!parse_value(&y_arg)) return 0;
                
                /* Optional k parameter */
                if (current_token.type == TOK_COMMA) {
                    value_t k_val;
                    next_token();
                    if (!parse_value(&k_val)) return 0;
                    k = (int)apf_to_long(&k_val.v.scalar.re);
                }
                
                if (current_token.type != TOK_RPAREN) {
                    printf("Error: expected ')'\n");
                    return 0;
                }
                next_token();
                
                if (pv_arg.type != VAL_MATRIX || y_arg.type != VAL_MATRIX) {
                    printf("Error: fitcknn requires matrix arguments\n");
                    return 0;
                }
                
                if (!ml_fit_knn(&model_id, &pv_arg.v.matrix, &y_arg.v.matrix, k)) {
                    return 0;
                }
                
                /* Return model ID as scalar */
                apf_from_int(&result->v.scalar.re, model_id);
                apf_zero(&result->v.scalar.im);
                result->type = VAL_SCALAR;
                return 1;
            }
            
            /* ML: fitcsvm(X, Y) or fitcsvm(X, Y, standardize) */
            if (str_eq(name, "fitcsvm")) {
                value_t y_arg;
                int standardize = 0;
                int model_id;
                
                next_token();
                if (!parse_value(&pv_arg)) return 0;
                if (current_token.type != TOK_COMMA) {
                    printf("Error: expected ',' in fitcsvm(X, Y)\n");
                    return 0;
                }
                next_token();
                if (!parse_value(&y_arg)) return 0;
                
                /* Optional standardize parameter (1 = true) */
                if (current_token.type == TOK_COMMA) {
                    value_t s_val;
                    next_token();
                    if (!parse_value(&s_val)) return 0;
                    standardize = (int)apf_to_long(&s_val.v.scalar.re);
                }
                
                if (current_token.type != TOK_RPAREN) {
                    printf("Error: expected ')'\n");
                    return 0;
                }
                next_token();
                
                if (pv_arg.type != VAL_MATRIX || y_arg.type != VAL_MATRIX) {
                    printf("Error: fitcsvm requires matrix arguments\n");
                    return 0;
                }
                
                if (!ml_fit_svm(&model_id, &pv_arg.v.matrix, &y_arg.v.matrix, standardize)) {
                    return 0;
                }
                
                apf_from_int(&result->v.scalar.re, model_id);
                apf_zero(&result->v.scalar.im);
                result->type = VAL_SCALAR;
                return 1;
            }
            
            /* ML: fitctree(X, Y) - train decision tree */
            if (str_eq(name, "fitctree")) {
                value_t y_arg;
                int model_id;
                
                next_token();
                if (!parse_value(&pv_arg)) return 0;
                if (current_token.type != TOK_COMMA) {
                    printf("Error: expected ',' in fitctree(X, Y)\n");
                    return 0;
                }
                next_token();
                if (!parse_value(&y_arg)) return 0;
                
                /* Skip any extra parameters */
                while (current_token.type == TOK_COMMA) {
                    value_t dummy;
                    next_token();
                    if (!parse_value(&dummy)) return 0;
                }
                
                if (current_token.type != TOK_RPAREN) {
                    printf("Error: expected ')'\n");
                    return 0;
                }
                next_token();
                
                if (pv_arg.type != VAL_MATRIX || y_arg.type != VAL_MATRIX) {
                    printf("Error: fitctree requires matrix arguments\n");
                    return 0;
                }
                
                if (!ml_fit_tree(&model_id, &pv_arg.v.matrix, &y_arg.v.matrix)) {
                    return 0;
                }
                
                apf_from_int(&result->v.scalar.re, model_id);
                apf_zero(&result->v.scalar.im);
                result->type = VAL_SCALAR;
                return 1;
            }
            
            /* ML: fitcnb(X, Y) - train Naive Bayes */
            if (str_eq(name, "fitcnb")) {
                value_t y_arg;
                int model_id;
                
                next_token();
                if (!parse_value(&pv_arg)) return 0;
                if (current_token.type != TOK_COMMA) {
                    printf("Error: expected ',' in fitcnb(X, Y)\n");
                    return 0;
                }
                next_token();
                if (!parse_value(&y_arg)) return 0;
                
                /* Skip any extra parameters */
                while (current_token.type == TOK_COMMA) {
                    value_t dummy;
                    next_token();
                    if (!parse_value(&dummy)) return 0;
                }
                
                if (current_token.type != TOK_RPAREN) {
                    printf("Error: expected ')'\n");
                    return 0;
                }
                next_token();
                
                if (pv_arg.type != VAL_MATRIX || y_arg.type != VAL_MATRIX) {
                    printf("Error: fitcnb requires matrix arguments\n");
                    return 0;
                }
                
                if (!ml_fit_nb(&model_id, &pv_arg.v.matrix, &y_arg.v.matrix)) {
                    return 0;
                }
                
                apf_from_int(&result->v.scalar.re, model_id);
                apf_zero(&result->v.scalar.im);
                result->type = VAL_SCALAR;
                return 1;
            }
            
            /* ML: predict(model, X) - predict using trained model */
            if (str_eq(name, "predict")) {
                value_t model_arg;
                int model_id;
                
                next_token();
                if (!parse_value(&model_arg)) return 0;
                if (current_token.type != TOK_COMMA) {
                    printf("Error: expected ',' in predict(model, X)\n");
                    return 0;
                }
                next_token();
                if (!parse_value(&pv_arg)) return 0;
                
                if (current_token.type != TOK_RPAREN) {
                    printf("Error: expected ')'\n");
                    return 0;
                }
                next_token();
                
                if (model_arg.type != VAL_SCALAR) {
                    printf("Error: predict requires model ID as first argument\n");
                    return 0;
                }
                if (pv_arg.type != VAL_MATRIX) {
                    printf("Error: predict requires matrix as second argument\n");
                    return 0;
                }
                
                model_id = (int)apf_to_long(&model_arg.v.scalar.re);
                if (!ml_predict(&result->v.matrix, model_id, &pv_arg.v.matrix)) {
                    return 0;
                }
                result->type = VAL_MATRIX;
                return 1;
            }
            
            /* ML: cvpartition(Y, holdout_pct) - create train/test split */
            if (str_eq(name, "cvpartition")) {
                int holdout_pct = 30;  /* default 30% test */
                
                next_token();
                if (!parse_value(&pv_arg)) return 0;
                
                /* Optional holdout percentage (0-100) */
                if (current_token.type == TOK_COMMA) {
                    value_t h_val;
                    next_token();
                    if (!parse_value(&h_val)) return 0;
                    holdout_pct = (int)apf_to_long(&h_val.v.scalar.re);
                }
                
                if (current_token.type != TOK_RPAREN) {
                    printf("Error: expected ')'\n");
                    return 0;
                }
                next_token();
                
                if (pv_arg.type != VAL_MATRIX) {
                    printf("Error: cvpartition requires matrix argument\n");
                    return 0;
                }
                
                /* Return training indices */
                {
                    matrix_t train_idx = {0,0,NULL}, test_idx = {0,0,NULL};
                    int n_samples;
                    n_samples = pv_arg.v.matrix.rows;
                    ml_cvpartition(&train_idx, &test_idx, n_samples, holdout_pct);
                    mat_copy(&result->v.matrix, &train_idx);
                }
                result->type = VAL_MATRIX;
                return 1;
            }
            
            /* ML: confusionmat(Ytrue, Ypred) - confusion matrix */
            if (str_eq(name, "confusionmat")) {
                value_t y_pred;
                
                next_token();
                if (!parse_value(&pv_arg)) return 0;
                if (current_token.type != TOK_COMMA) {
                    printf("Error: expected ',' in confusionmat(Ytrue, Ypred)\n");
                    return 0;
                }
                next_token();
                if (!parse_value(&y_pred)) return 0;
                
                if (current_token.type != TOK_RPAREN) {
                    printf("Error: expected ')'\n");
                    return 0;
                }
                next_token();
                
                if (pv_arg.type != VAL_MATRIX || y_pred.type != VAL_MATRIX) {
                    printf("Error: confusionmat requires matrix arguments\n");
                    return 0;
                }
                
                ml_confusionmat(&result->v.matrix, &pv_arg.v.matrix, &y_pred.v.matrix);
                result->type = VAL_MATRIX;
                return 1;
            }
            
            /* ML: accuracy(Ytrue, Ypred) - classification accuracy */
            if (str_eq(name, "accuracy")) {
                value_t y_pred;
                
                next_token();
                if (!parse_value(&pv_arg)) return 0;
                if (current_token.type != TOK_COMMA) {
                    printf("Error: expected ',' in accuracy(Ytrue, Ypred)\n");
                    return 0;
                }
                next_token();
                if (!parse_value(&y_pred)) return 0;
                
                if (current_token.type != TOK_RPAREN) {
                    printf("Error: expected ')'\n");
                    return 0;
                }
                next_token();
                
                if (pv_arg.type != VAL_MATRIX || y_pred.type != VAL_MATRIX) {
                    printf("Error: accuracy requires matrix arguments\n");
                    return 0;
                }
                
                ml_accuracy(&result->v.scalar, &pv_arg.v.matrix, &y_pred.v.matrix);
                result->type = VAL_SCALAR;
                return 1;
            }
            
            /* Array manipulation functions that return matrices */
            if (str_eq(name, "fliplr") || str_eq(name, "flipud") || str_eq(name, "flip") ||
                str_eq(name, "sort") || str_eq(name, "sortrows") || 
                str_eq(name, "cumsum") || str_eq(name, "cumprod") ||
                str_eq(name, "cummax") || str_eq(name, "cummin") ||
                str_eq(name, "diff") || str_eq(name, "find") || str_eq(name, "unique") ||
                str_eq(name, "rot90") || str_eq(name, "normalize") ||
                str_eq(name, "triu") || str_eq(name, "tril") || str_eq(name, "squeeze") ||
                str_eq(name, "rescale") ||
                str_eq(name, "dateshift") || str_eq(name, "startofmonth_vec")) {
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
                    mat_zero(&result->v.matrix, pv_arg.v.matrix.rows, pv_arg.v.matrix.cols);
                    for (r = 0; r < pv_arg.v.matrix.rows; r++) {
                        for (c = 0; c < pv_arg.v.matrix.cols; c++) {
                            MAT_AT(&result->v.matrix, r, pv_arg.v.matrix.cols - 1 - c) = 
                                MAT_AT(&pv_arg.v.matrix, r, c);
                        }
                    }
                } else if (str_eq(name, "flipud")) {
                    mat_zero(&result->v.matrix, pv_arg.v.matrix.rows, pv_arg.v.matrix.cols);
                    for (r = 0; r < pv_arg.v.matrix.rows; r++) {
                        for (c = 0; c < pv_arg.v.matrix.cols; c++) {
                            MAT_AT(&result->v.matrix, pv_arg.v.matrix.rows - 1 - r, c) = 
                                MAT_AT(&pv_arg.v.matrix, r, c);
                        }
                    }
                } else if (str_eq(name, "flip")) {
                    /* flip: flip vector (like fliplr for row vec, flipud for col vec) */
                    mat_zero(&result->v.matrix, pv_arg.v.matrix.rows, pv_arg.v.matrix.cols);
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
                    mat_zero(&result->v.matrix, pv_arg.v.matrix.cols, pv_arg.v.matrix.rows);
                    for (r = 0; r < pv_arg.v.matrix.rows; r++) {
                        for (c = 0; c < pv_arg.v.matrix.cols; c++) {
                            MAT_AT(&result->v.matrix, pv_arg.v.matrix.cols - 1 - c, r) = 
                                MAT_AT(&pv_arg.v.matrix, r, c);
                        }
                    }
                } else if (str_eq(name, "normalize")) {
                    /* normalize: normalize vector to unit length */
                    apf sum_sq, norm_val, sq;
                    mat_zero(&result->v.matrix, pv_arg.v.matrix.rows, pv_arg.v.matrix.cols);
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
                    /* Then remove duplicates - allocate max possible size first */
                    mat_zero(&result->v.matrix, 1, n > 0 ? n : 1);
                    if (!result->v.matrix.data) return 0;
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
                        if (!is_dup) {
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
                } else if (str_eq(name, "sortrows")) {
                    /* sortrows: sort matrix rows by first column (bubble sort) */
                    int nrows, ncols, i, j, c;
                    mat_copy(&result->v.matrix, &pv_arg.v.matrix);
                    nrows = result->v.matrix.rows;
                    ncols = result->v.matrix.cols;
                    for (i = 0; i < nrows - 1; i++) {
                        for (j = 0; j < nrows - 1 - i; j++) {
                            /* Compare by first column */
                            if (apf_cmp(&MAT_AT(&result->v.matrix, j, 0).re,
                                       &MAT_AT(&result->v.matrix, j+1, 0).re) > 0) {
                                /* Swap entire rows */
                                for (c = 0; c < ncols; c++) {
                                    apfc tmp = MAT_AT(&result->v.matrix, j, c);
                                    MAT_AT(&result->v.matrix, j, c) = MAT_AT(&result->v.matrix, j+1, c);
                                    MAT_AT(&result->v.matrix, j+1, c) = tmp;
                                }
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
                    mat_zero(&result->v.matrix, pv_arg.v.matrix.rows, pv_arg.v.matrix.cols);
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
                    mat_zero(&result->v.matrix, pv_arg.v.matrix.rows, pv_arg.v.matrix.cols);
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
                } else if (str_eq(name, "cummax")) {
                    /* cummax: cumulative maximum */
                    apf current_max;
                    int first = 1;
                    
                    mat_zero(&result->v.matrix, pv_arg.v.matrix.rows, pv_arg.v.matrix.cols);
                    
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
                    
                    mat_zero(&result->v.matrix, pv_arg.v.matrix.rows, pv_arg.v.matrix.cols);
                    
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
                    
                    mat_zero(&result->v.matrix, pv_arg.v.matrix.rows, pv_arg.v.matrix.cols);
                    
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
                } else if (str_eq(name, "dateshift") || str_eq(name, "startofmonth_vec")) {
                    /* dateshift: Apply start-of-month to each element (like MATLAB dateshift(x,'start','month')) */
                    static const int days_in_month[] = {31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31};
                    
                    mat_zero(&result->v.matrix, pv_arg.v.matrix.rows, pv_arg.v.matrix.cols);
                    if (!result->v.matrix.data) return 0;
                    
                    for (r = 0; r < pv_arg.v.matrix.rows; r++) {
                        for (c = 0; c < pv_arg.v.matrix.cols; c++) {
                            long ts = apf_to_long(&MAT_AT(&pv_arg.v.matrix, r, c).re);
                            long days_since_epoch = ts / 86400;
                            long year = 1970, month, i, new_days;
                            
                            /* Find year */
                            while (1) {
                                int year_days = 365;
                                if ((year % 4 == 0 && year % 100 != 0) || (year % 400 == 0)) year_days = 366;
                                if (days_since_epoch < year_days) break;
                                days_since_epoch -= year_days;
                                year++;
                            }
                            
                            /* Find month */
                            month = 1;
                            for (i = 0; i < 12; i++) {
                                int dm = days_in_month[i];
                                if (i == 1 && ((year % 4 == 0 && year % 100 != 0) || (year % 400 == 0))) dm++;
                                if (days_since_epoch < dm) break;
                                days_since_epoch -= dm;
                                month++;
                            }
                            
                            /* Calculate days from 1970 to start of this month */
                            new_days = 0;
                            for (i = 1970; i < year; i++) {
                                new_days += 365;
                                if ((i % 4 == 0 && i % 100 != 0) || (i % 400 == 0)) new_days++;
                            }
                            for (i = 1; i < month; i++) {
                                new_days += days_in_month[i-1];
                                if (i == 2 && ((year % 4 == 0 && year % 100 != 0) || (year % 400 == 0))) new_days++;
                            }
                            
                            apf_from_int(&MAT_AT(&result->v.matrix, r, c).re, new_days * 86400);
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
                ensure_result_matrix(&result->v.matrix, pv_arg.v.matrix.rows, pv_arg.v.matrix.cols);
                
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
                ensure_result_matrix(&result->v.matrix, pv_arg.v.matrix.rows, pv_arg.v.matrix.cols);
                
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
            
            /* xcorr(a, b) - cross-correlation of two vectors */
            if (str_eq(name, "xcorr")) {
                value_t arg2;
                
                next_token();
                if (!parse_value(&pv_arg)) return 0;
                if (current_token.type != TOK_COMMA) {
                    /* xcorr(a) = autocorrelation */
                    arg2 = pv_arg;
                } else {
                    next_token();
                    if (!parse_value(&arg2)) return 0;
                }
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
                
                result->type = VAL_MATRIX;
                mat_xcorr(&result->v.matrix, &pv_arg.v.matrix, &arg2.v.matrix);
                return 1;
            }
            
            /* xcov(a, b) - cross-covariance */
            if (str_eq(name, "xcov")) {
                value_t arg2;
                
                next_token();
                if (!parse_value(&pv_arg)) return 0;
                if (current_token.type != TOK_COMMA) {
                    arg2 = pv_arg;
                } else {
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
                if (arg2.type == VAL_SCALAR) {
                    apfc saved = arg2.v.scalar;
                    mat_zero(&arg2.v.matrix, 1, 1);
                    MAT_AT(&arg2.v.matrix, 0, 0) = saved;
                    arg2.type = VAL_MATRIX;
                }
                
                result->type = VAL_MATRIX;
                mat_xcov(&result->v.matrix, &pv_arg.v.matrix, &arg2.v.matrix);
                return 1;
            }
            
            /* histcounts(data) or histcounts(data, nbins) */
            if (str_eq(name, "histcounts")) {
                int nbins = 10;  /* default */
                
                next_token();
                if (!parse_value(&pv_arg)) return 0;
                if (current_token.type == TOK_COMMA) {
                    apfc nbins_val;
                    next_token();
                    if (!parse_expr(&nbins_val)) return 0;
                    nbins = (int)apf_to_double(&nbins_val.re);
                    if (nbins < 1) nbins = 1;
                    if (nbins > 1000) nbins = 1000;
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
                
                result->type = VAL_MATRIX;
                mat_histcounts(&result->v.matrix, &pv_arg.v.matrix, nbins);
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
            
            /* movavg(v, k) - simple moving average (alias for movmean) */
            if (str_eq(name, "movavg")) {
                apfc k_val;
                int n, k, i, j;
                
                next_token();
                if (!parse_value(&pv_arg)) return 0;
                if (current_token.type != TOK_COMMA) {
                    printf("Error: movavg requires (v, k)\n");
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
                
                result->type = VAL_MATRIX;
                mat_zero(&result->v.matrix, pv_arg.v.matrix.rows, pv_arg.v.matrix.cols);
                
                for (i = 0; i < n; i++) {
                    apf sum, count;
                    int ri = i % pv_arg.v.matrix.rows;
                    int ci = i / pv_arg.v.matrix.rows;
                    int start = i - k + 1;
                    int cnt = 0;
                    
                    if (start < 0) start = 0;
                    
                    apf_zero(&sum);
                    for (j = start; j <= i; j++) {
                        int rj = j % pv_arg.v.matrix.rows;
                        int cj = j / pv_arg.v.matrix.rows;
                        apf_add(&sum, &sum, &MAT_AT(&pv_arg.v.matrix, rj, cj).re);
                        cnt++;
                    }
                    
                    if (cnt > 0) {
                        apf_from_int(&count, cnt);
                        apf_div(&MAT_AT(&result->v.matrix, ri, ci).re, &sum, &count);
                    }
                    apf_zero(&MAT_AT(&result->v.matrix, ri, ci).im);
                }
                return 1;
            }
            
            /* ewma(v, alpha) - exponentially weighted moving average */
            if (str_eq(name, "ewma") || str_eq(name, "expsmooth")) {
                apfc alpha_val;
                int n, i;
                apf alpha, one_minus_alpha, one, prev, curr, term;
                
                next_token();
                if (!parse_value(&pv_arg)) return 0;
                if (current_token.type != TOK_COMMA) {
                    printf("Error: ewma requires (v, alpha)\n");
                    return 0;
                }
                next_token();
                if (!parse_expr(&alpha_val)) return 0;
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
                apf_copy(&alpha, &alpha_val.re);
                apf_from_int(&one, 1);
                apf_sub(&one_minus_alpha, &one, &alpha);
                
                result->type = VAL_MATRIX;
                mat_zero(&result->v.matrix, pv_arg.v.matrix.rows, pv_arg.v.matrix.cols);
                
                /* First value is just the first data point */
                apf_copy(&prev, &MAT_AT(&pv_arg.v.matrix, 0, 0).re);
                apf_copy(&MAT_AT(&result->v.matrix, 0, 0).re, &prev);
                
                for (i = 1; i < n; i++) {
                    int ri = i % pv_arg.v.matrix.rows;
                    int ci = i / pv_arg.v.matrix.rows;
                    
                    /* ewma[i] = alpha * x[i] + (1-alpha) * ewma[i-1] */
                    apf_mul(&curr, &alpha, &MAT_AT(&pv_arg.v.matrix, ri, ci).re);
                    apf_mul(&term, &one_minus_alpha, &prev);
                    apf_add(&prev, &curr, &term);
                    apf_copy(&MAT_AT(&result->v.matrix, ri, ci).re, &prev);
                    apf_zero(&MAT_AT(&result->v.matrix, ri, ci).im);
                }
                return 1;
            }
            
            /* lag(v, k) - lag time series by k periods */
            if (str_eq(name, "lag")) {
                apfc k_val;
                int n, k, i;
                apf nan_val;
                
                next_token();
                if (!parse_value(&pv_arg)) return 0;
                
                k = 1;  /* default lag */
                if (current_token.type == TOK_COMMA) {
                    next_token();
                    if (!parse_expr(&k_val)) return 0;
                    k = apf_to_long(&k_val.re);
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
                
                n = pv_arg.v.matrix.rows * pv_arg.v.matrix.cols;
                if (k < 0) k = 0;
                if (k > n) k = n;
                
                apf_set_nan(&nan_val);
                
                result->type = VAL_MATRIX;
                mat_zero(&result->v.matrix, pv_arg.v.matrix.rows, pv_arg.v.matrix.cols);
                
                for (i = 0; i < n; i++) {
                    int ri = i % pv_arg.v.matrix.rows;
                    int ci = i / pv_arg.v.matrix.rows;
                    
                    if (i < k) {
                        apf_copy(&MAT_AT(&result->v.matrix, ri, ci).re, &nan_val);
                    } else {
                        int src = i - k;
                        int sr = src % pv_arg.v.matrix.rows;
                        int sc = src / pv_arg.v.matrix.rows;
                        apf_copy(&MAT_AT(&result->v.matrix, ri, ci).re, 
                                 &MAT_AT(&pv_arg.v.matrix, sr, sc).re);
                    }
                    apf_zero(&MAT_AT(&result->v.matrix, ri, ci).im);
                }
                return 1;
            }
            
            /* lead(v, k) - lead time series by k periods */
            if (str_eq(name, "lead")) {
                apfc k_val;
                int n, k, i;
                apf nan_val;
                
                next_token();
                if (!parse_value(&pv_arg)) return 0;
                
                k = 1;  /* default lead */
                if (current_token.type == TOK_COMMA) {
                    next_token();
                    if (!parse_expr(&k_val)) return 0;
                    k = apf_to_long(&k_val.re);
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
                
                n = pv_arg.v.matrix.rows * pv_arg.v.matrix.cols;
                if (k < 0) k = 0;
                if (k > n) k = n;
                
                apf_set_nan(&nan_val);
                
                result->type = VAL_MATRIX;
                mat_zero(&result->v.matrix, pv_arg.v.matrix.rows, pv_arg.v.matrix.cols);
                
                for (i = 0; i < n; i++) {
                    int ri = i % pv_arg.v.matrix.rows;
                    int ci = i / pv_arg.v.matrix.rows;
                    
                    if (i + k >= n) {
                        apf_copy(&MAT_AT(&result->v.matrix, ri, ci).re, &nan_val);
                    } else {
                        int src = i + k;
                        int sr = src % pv_arg.v.matrix.rows;
                        int sc = src / pv_arg.v.matrix.rows;
                        apf_copy(&MAT_AT(&result->v.matrix, ri, ci).re, 
                                 &MAT_AT(&pv_arg.v.matrix, sr, sc).re);
                    }
                    apf_zero(&MAT_AT(&result->v.matrix, ri, ci).im);
                }
                return 1;
            }
            
            /* detrend(v) - remove linear trend */
            if (str_eq(name, "detrend")) {
                int n, i;
                apf sumx, sumy, sumxy, sumxx, xbar, ybar, slope, intercept;
                apf x, y, fitted, count, tmp;
                
                next_token();
                if (!parse_value(&pv_arg)) return 0;
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
                
                /* Compute linear regression: y = slope*x + intercept */
                apf_zero(&sumx); apf_zero(&sumy);
                apf_zero(&sumxy); apf_zero(&sumxx);
                
                for (i = 0; i < n; i++) {
                    int ri = i % pv_arg.v.matrix.rows;
                    int ci = i / pv_arg.v.matrix.rows;
                    apf_from_int(&x, i);
                    apf_copy(&y, &MAT_AT(&pv_arg.v.matrix, ri, ci).re);
                    apf_add(&sumx, &sumx, &x);
                    apf_add(&sumy, &sumy, &y);
                    apf_mul(&tmp, &x, &y);
                    apf_add(&sumxy, &sumxy, &tmp);
                    apf_mul(&tmp, &x, &x);
                    apf_add(&sumxx, &sumxx, &tmp);
                }
                
                apf_from_int(&count, n);
                apf_div(&xbar, &sumx, &count);
                apf_div(&ybar, &sumy, &count);
                
                /* slope = (sumxy - n*xbar*ybar) / (sumxx - n*xbar*xbar) */
                {
                    apf num, den;
                    apf_mul(&tmp, &xbar, &ybar);
                    apf_mul(&tmp, &tmp, &count);
                    apf_sub(&num, &sumxy, &tmp);
                    apf_mul(&tmp, &xbar, &xbar);
                    apf_mul(&tmp, &tmp, &count);
                    apf_sub(&den, &sumxx, &tmp);
                    apf_div(&slope, &num, &den);
                }
                
                /* intercept = ybar - slope*xbar */
                apf_mul(&tmp, &slope, &xbar);
                apf_sub(&intercept, &ybar, &tmp);
                
                result->type = VAL_MATRIX;
                mat_zero(&result->v.matrix, pv_arg.v.matrix.rows, pv_arg.v.matrix.cols);
                
                for (i = 0; i < n; i++) {
                    int ri = i % pv_arg.v.matrix.rows;
                    int ci = i / pv_arg.v.matrix.rows;
                    apf_from_int(&x, i);
                    apf_mul(&fitted, &slope, &x);
                    apf_add(&fitted, &fitted, &intercept);
                    apf_sub(&MAT_AT(&result->v.matrix, ri, ci).re, 
                            &MAT_AT(&pv_arg.v.matrix, ri, ci).re, &fitted);
                    apf_zero(&MAT_AT(&result->v.matrix, ri, ci).im);
                }
                return 1;
            }
            
            /* trend(v) - extract linear trend */
            if (str_eq(name, "trend")) {
                int n, i;
                apf sumx, sumy, sumxy, sumxx, xbar, ybar, slope, intercept;
                apf x, y, fitted, count, tmp;
                
                next_token();
                if (!parse_value(&pv_arg)) return 0;
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
                
                /* Compute linear regression: y = slope*x + intercept */
                apf_zero(&sumx); apf_zero(&sumy);
                apf_zero(&sumxy); apf_zero(&sumxx);
                
                for (i = 0; i < n; i++) {
                    int ri = i % pv_arg.v.matrix.rows;
                    int ci = i / pv_arg.v.matrix.rows;
                    apf_from_int(&x, i);
                    apf_copy(&y, &MAT_AT(&pv_arg.v.matrix, ri, ci).re);
                    apf_add(&sumx, &sumx, &x);
                    apf_add(&sumy, &sumy, &y);
                    apf_mul(&tmp, &x, &y);
                    apf_add(&sumxy, &sumxy, &tmp);
                    apf_mul(&tmp, &x, &x);
                    apf_add(&sumxx, &sumxx, &tmp);
                }
                
                apf_from_int(&count, n);
                apf_div(&xbar, &sumx, &count);
                apf_div(&ybar, &sumy, &count);
                
                /* slope = (sumxy - n*xbar*ybar) / (sumxx - n*xbar*xbar) */
                {
                    apf num, den;
                    apf_mul(&tmp, &xbar, &ybar);
                    apf_mul(&tmp, &tmp, &count);
                    apf_sub(&num, &sumxy, &tmp);
                    apf_mul(&tmp, &xbar, &xbar);
                    apf_mul(&tmp, &tmp, &count);
                    apf_sub(&den, &sumxx, &tmp);
                    apf_div(&slope, &num, &den);
                }
                
                /* intercept = ybar - slope*xbar */
                apf_mul(&tmp, &slope, &xbar);
                apf_sub(&intercept, &ybar, &tmp);
                
                result->type = VAL_MATRIX;
                mat_zero(&result->v.matrix, pv_arg.v.matrix.rows, pv_arg.v.matrix.cols);
                
                for (i = 0; i < n; i++) {
                    int ri = i % pv_arg.v.matrix.rows;
                    int ci = i / pv_arg.v.matrix.rows;
                    apf_from_int(&x, i);
                    apf_mul(&fitted, &slope, &x);
                    apf_add(&fitted, &fitted, &intercept);
                    apf_copy(&MAT_AT(&result->v.matrix, ri, ci).re, &fitted);
                    apf_zero(&MAT_AT(&result->v.matrix, ri, ci).im);
                }
                return 1;
            }
            
            /* pctchange(v) - percentage change */
            if (str_eq(name, "pctchange") || str_eq(name, "percentchange")) {
                int n, i;
                apf nan_val, hundred;
                
                next_token();
                if (!parse_value(&pv_arg)) return 0;
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
                apf_set_nan(&nan_val);
                apf_from_int(&hundred, 100);
                
                result->type = VAL_MATRIX;
                mat_zero(&result->v.matrix, pv_arg.v.matrix.rows, pv_arg.v.matrix.cols);
                
                /* First value is NaN */
                apf_copy(&MAT_AT(&result->v.matrix, 0, 0).re, &nan_val);
                
                for (i = 1; i < n; i++) {
                    int ri = i % pv_arg.v.matrix.rows;
                    int ci = i / pv_arg.v.matrix.rows;
                    int pr = (i-1) % pv_arg.v.matrix.rows;
                    int pc = (i-1) / pv_arg.v.matrix.rows;
                    apf curr, prev, diff, pct;
                    
                    apf_copy(&curr, &MAT_AT(&pv_arg.v.matrix, ri, ci).re);
                    apf_copy(&prev, &MAT_AT(&pv_arg.v.matrix, pr, pc).re);
                    
                    if (apf_is_zero(&prev)) {
                        apf_copy(&MAT_AT(&result->v.matrix, ri, ci).re, &nan_val);
                    } else {
                        apf_sub(&diff, &curr, &prev);
                        apf_div(&pct, &diff, &prev);
                        apf_mul(&MAT_AT(&result->v.matrix, ri, ci).re, &pct, &hundred);
                    }
                    apf_zero(&MAT_AT(&result->v.matrix, ri, ci).im);
                }
                return 1;
            }
            
            /* autocorr(v, maxlag) - autocorrelation function */
            if (str_eq(name, "autocorr")) {
                apfc maxlag_val;
                int n, maxlag, lag, i;
                apf mean, var, sum, count, tmp;
                
                next_token();
                if (!parse_value(&pv_arg)) return 0;
                
                maxlag = 20;  /* default */
                if (current_token.type == TOK_COMMA) {
                    next_token();
                    if (!parse_expr(&maxlag_val)) return 0;
                    maxlag = apf_to_long(&maxlag_val.re);
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
                
                n = pv_arg.v.matrix.rows * pv_arg.v.matrix.cols;
                if (maxlag >= n) maxlag = n - 1;
                if (maxlag < 1) maxlag = 1;
                
                /* Compute mean */
                apf_zero(&sum);
                for (i = 0; i < n; i++) {
                    int ri = i % pv_arg.v.matrix.rows;
                    int ci = i / pv_arg.v.matrix.rows;
                    apf_add(&sum, &sum, &MAT_AT(&pv_arg.v.matrix, ri, ci).re);
                }
                apf_from_int(&count, n);
                apf_div(&mean, &sum, &count);
                
                /* Compute variance */
                apf_zero(&var);
                for (i = 0; i < n; i++) {
                    int ri = i % pv_arg.v.matrix.rows;
                    int ci = i / pv_arg.v.matrix.rows;
                    apf_sub(&tmp, &MAT_AT(&pv_arg.v.matrix, ri, ci).re, &mean);
                    apf_mul(&tmp, &tmp, &tmp);
                    apf_add(&var, &var, &tmp);
                }
                
                result->type = VAL_MATRIX;
                mat_zero(&result->v.matrix, maxlag + 1, 1);
                
                /* Compute autocorrelations for each lag */
                for (lag = 0; lag <= maxlag; lag++) {
                    apf acov;
                    apf_zero(&acov);
                    
                    for (i = 0; i < n - lag; i++) {
                        int ri = i % pv_arg.v.matrix.rows;
                        int ci = i / pv_arg.v.matrix.rows;
                        int rj = (i + lag) % pv_arg.v.matrix.rows;
                        int cj = (i + lag) / pv_arg.v.matrix.rows;
                        apf xi, xj;
                        
                        apf_sub(&xi, &MAT_AT(&pv_arg.v.matrix, ri, ci).re, &mean);
                        apf_sub(&xj, &MAT_AT(&pv_arg.v.matrix, rj, cj).re, &mean);
                        apf_mul(&tmp, &xi, &xj);
                        apf_add(&acov, &acov, &tmp);
                    }
                    
                    if (!apf_is_zero(&var)) {
                        apf_div(&MAT_AT(&result->v.matrix, lag, 0).re, &acov, &var);
                    }
                    apf_zero(&MAT_AT(&result->v.matrix, lag, 0).im);
                }
                return 1;
            }
            
            /* forecast(v, periods) - simple linear forecast */
            if (str_eq(name, "forecast")) {
                apfc periods_val;
                int n, periods, i;
                apf sumx, sumy, sumxy, sumxx, xbar, ybar, slope, intercept;
                apf x, y, count, tmp;
                
                next_token();
                if (!parse_value(&pv_arg)) return 0;
                if (current_token.type != TOK_COMMA) {
                    printf("Error: forecast requires (v, periods)\n");
                    return 0;
                }
                next_token();
                if (!parse_expr(&periods_val)) return 0;
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
                periods = apf_to_long(&periods_val.re);
                if (periods < 1) periods = 1;
                
                /* Compute linear regression for trend */
                apf_zero(&sumx); apf_zero(&sumy);
                apf_zero(&sumxy); apf_zero(&sumxx);
                
                for (i = 0; i < n; i++) {
                    int ri = i % pv_arg.v.matrix.rows;
                    int ci = i / pv_arg.v.matrix.rows;
                    apf_from_int(&x, i);
                    apf_copy(&y, &MAT_AT(&pv_arg.v.matrix, ri, ci).re);
                    apf_add(&sumx, &sumx, &x);
                    apf_add(&sumy, &sumy, &y);
                    apf_mul(&tmp, &x, &y);
                    apf_add(&sumxy, &sumxy, &tmp);
                    apf_mul(&tmp, &x, &x);
                    apf_add(&sumxx, &sumxx, &tmp);
                }
                
                apf_from_int(&count, n);
                apf_div(&xbar, &sumx, &count);
                apf_div(&ybar, &sumy, &count);
                
                /* slope = (sumxy - n*xbar*ybar) / (sumxx - n*xbar*xbar) */
                {
                    apf num, den;
                    apf_mul(&tmp, &xbar, &ybar);
                    apf_mul(&tmp, &tmp, &count);
                    apf_sub(&num, &sumxy, &tmp);
                    apf_mul(&tmp, &xbar, &xbar);
                    apf_mul(&tmp, &tmp, &count);
                    apf_sub(&den, &sumxx, &tmp);
                    apf_div(&slope, &num, &den);
                }
                
                /* intercept = ybar - slope*xbar */
                apf_mul(&tmp, &slope, &xbar);
                apf_sub(&intercept, &ybar, &tmp);
                
                /* Generate forecasts */
                result->type = VAL_MATRIX;
                mat_zero(&result->v.matrix, periods, 1);
                
                for (i = 0; i < periods; i++) {
                    apf forecast_x;
                    apf_from_int(&forecast_x, n + i);
                    apf_mul(&tmp, &slope, &forecast_x);
                    apf_add(&MAT_AT(&result->v.matrix, i, 0).re, &tmp, &intercept);
                    apf_zero(&MAT_AT(&result->v.matrix, i, 0).im);
                }
                return 1;
            }
            
            /* polyfit(x, y, n) - polynomial curve fitting */
            if (str_eq(name, "polyfit")) {
                value_t y_arg;
                apfc deg_val;
                int n, deg, i, j, k;
                apf *x_vals, *y_vals;
                apf *A, *b, *coeffs;
                
                next_token();
                if (!parse_value(&pv_arg)) return 0;
                if (current_token.type != TOK_COMMA) {
                    printf("Error: polyfit requires (x, y, degree)\n");
                    return 0;
                }
                next_token();
                if (!parse_value(&y_arg)) return 0;
                if (current_token.type != TOK_COMMA) {
                    printf("Error: polyfit requires (x, y, degree)\n");
                    return 0;
                }
                next_token();
                if (!parse_expr(&deg_val)) return 0;
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
                if (y_arg.type == VAL_SCALAR) {
                    apfc saved = y_arg.v.scalar;
                    mat_zero(&y_arg.v.matrix, 1, 1);
                    MAT_AT(&y_arg.v.matrix, 0, 0) = saved;
                    y_arg.type = VAL_MATRIX;
                }
                
                n = pv_arg.v.matrix.rows * pv_arg.v.matrix.cols;
                deg = apf_to_long(&deg_val.re);
                if (deg < 0) deg = 0;
                if (deg > 10) deg = 10;
                if (deg >= n) deg = n - 1;
                
                /* Allocate temporary arrays */
                x_vals = (apf *)malloc(n * sizeof(apf));
                y_vals = (apf *)malloc(n * sizeof(apf));
                A = (apf *)malloc((deg+1) * (deg+1) * sizeof(apf));
                b = (apf *)malloc((deg+1) * sizeof(apf));
                coeffs = (apf *)malloc((deg+1) * sizeof(apf));
                
                if (!x_vals || !y_vals || !A || !b || !coeffs) {
                    free(x_vals); free(y_vals); free(A); free(b); free(coeffs);
                    printf("Error: memory allocation failed\n");
                    return 0;
                }
                
                /* Extract x and y values */
                for (i = 0; i < n; i++) {
                    int ri = i % pv_arg.v.matrix.rows;
                    int ci = i / pv_arg.v.matrix.rows;
                    apf_copy(&x_vals[i], &MAT_AT(&pv_arg.v.matrix, ri, ci).re);
                    ri = i % y_arg.v.matrix.rows;
                    ci = i / y_arg.v.matrix.rows;
                    apf_copy(&y_vals[i], &MAT_AT(&y_arg.v.matrix, ri, ci).re);
                }
                
                /* Build normal equations: A * coeffs = b */
                /* A[j][k] = sum(x^(j+k)), b[j] = sum(y * x^j) */
                for (j = 0; j <= deg; j++) {
                    apf_zero(&b[j]);
                    for (k = 0; k <= deg; k++) {
                        apf_zero(&A[j*(deg+1)+k]);
                    }
                }
                
                for (i = 0; i < n; i++) {
                    apf xpow;
                    apf_from_int(&xpow, 1);
                    for (j = 0; j <= deg; j++) {
                        apf xpow2, tmp;
                        apf_copy(&xpow2, &xpow);
                        for (k = 0; k <= deg; k++) {
                            apf_add(&A[j*(deg+1)+k], &A[j*(deg+1)+k], &xpow2);
                            apf_mul(&xpow2, &xpow2, &x_vals[i]);
                        }
                        apf_mul(&tmp, &y_vals[i], &xpow);
                        apf_add(&b[j], &b[j], &tmp);
                        apf_mul(&xpow, &xpow, &x_vals[i]);
                    }
                }
                
                /* Solve using Gaussian elimination */
                for (i = 0; i <= deg; i++) {
                    apf_copy(&coeffs[i], &b[i]);
                }
                
                for (i = 0; i <= deg; i++) {
                    /* Find pivot */
                    int max_row = i;
                    apf max_val, abs_val;
                    apf_abs(&max_val, &A[i*(deg+1)+i]);
                    for (k = i + 1; k <= deg; k++) {
                        apf_abs(&abs_val, &A[k*(deg+1)+i]);
                        if (apf_cmp(&abs_val, &max_val) > 0) {
                            apf_copy(&max_val, &abs_val);
                            max_row = k;
                        }
                    }
                    
                    /* Swap rows */
                    if (max_row != i) {
                        apf tmp;
                        for (k = 0; k <= deg; k++) {
                            apf_copy(&tmp, &A[i*(deg+1)+k]);
                            apf_copy(&A[i*(deg+1)+k], &A[max_row*(deg+1)+k]);
                            apf_copy(&A[max_row*(deg+1)+k], &tmp);
                        }
                        apf_copy(&tmp, &coeffs[i]);
                        apf_copy(&coeffs[i], &coeffs[max_row]);
                        apf_copy(&coeffs[max_row], &tmp);
                    }
                    
                    /* Eliminate */
                    for (k = i + 1; k <= deg; k++) {
                        apf factor, tmp;
                        if (apf_is_zero(&A[i*(deg+1)+i])) continue;
                        apf_div(&factor, &A[k*(deg+1)+i], &A[i*(deg+1)+i]);
                        for (j = i; j <= deg; j++) {
                            apf_mul(&tmp, &factor, &A[i*(deg+1)+j]);
                            apf_sub(&A[k*(deg+1)+j], &A[k*(deg+1)+j], &tmp);
                        }
                        apf_mul(&tmp, &factor, &coeffs[i]);
                        apf_sub(&coeffs[k], &coeffs[k], &tmp);
                    }
                }
                
                /* Back substitution */
                for (i = deg; i >= 0; i--) {
                    apf sum;
                    apf_zero(&sum);
                    for (k = i + 1; k <= deg; k++) {
                        apf tmp;
                        apf_mul(&tmp, &A[i*(deg+1)+k], &coeffs[k]);
                        apf_add(&sum, &sum, &tmp);
                    }
                    apf_sub(&coeffs[i], &coeffs[i], &sum);
                    if (!apf_is_zero(&A[i*(deg+1)+i])) {
                        apf_div(&coeffs[i], &coeffs[i], &A[i*(deg+1)+i]);
                    }
                }
                
                /* Return coefficients in descending order */
                result->type = VAL_MATRIX;
                mat_zero(&result->v.matrix, 1, deg + 1);
                for (i = 0; i <= deg; i++) {
                    apf_copy(&MAT_AT(&result->v.matrix, 0, deg - i).re, &coeffs[i]);
                    apf_zero(&MAT_AT(&result->v.matrix, 0, deg - i).im);
                }
                
                free(x_vals); free(y_vals); free(A); free(b); free(coeffs);
                return 1;
            }
            
            /* linreg(x, y) - simple linear regression */
            if (str_eq(name, "linreg") || str_eq(name, "regress")) {
                value_t y_arg;
                int n, i;
                apf sumx, sumy, sumxy, sumxx, sumyy;
                apf xbar, ybar, slope, intercept, r2;
                apf count, tmp, ss_tot, ss_res;
                
                next_token();
                if (!parse_value(&pv_arg)) return 0;
                if (current_token.type != TOK_COMMA) {
                    printf("Error: linreg requires (x, y)\n");
                    return 0;
                }
                next_token();
                if (!parse_value(&y_arg)) return 0;
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
                if (y_arg.type == VAL_SCALAR) {
                    apfc saved = y_arg.v.scalar;
                    mat_zero(&y_arg.v.matrix, 1, 1);
                    MAT_AT(&y_arg.v.matrix, 0, 0) = saved;
                    y_arg.type = VAL_MATRIX;
                }
                
                n = pv_arg.v.matrix.rows * pv_arg.v.matrix.cols;
                
                apf_zero(&sumx); apf_zero(&sumy);
                apf_zero(&sumxy); apf_zero(&sumxx); apf_zero(&sumyy);
                
                for (i = 0; i < n; i++) {
                    int ri = i % pv_arg.v.matrix.rows;
                    int ci = i / pv_arg.v.matrix.rows;
                    int rj = i % y_arg.v.matrix.rows;
                    int cj = i / y_arg.v.matrix.rows;
                    apf x, y;
                    apf_copy(&x, &MAT_AT(&pv_arg.v.matrix, ri, ci).re);
                    apf_copy(&y, &MAT_AT(&y_arg.v.matrix, rj, cj).re);
                    
                    apf_add(&sumx, &sumx, &x);
                    apf_add(&sumy, &sumy, &y);
                    apf_mul(&tmp, &x, &y);
                    apf_add(&sumxy, &sumxy, &tmp);
                    apf_mul(&tmp, &x, &x);
                    apf_add(&sumxx, &sumxx, &tmp);
                    apf_mul(&tmp, &y, &y);
                    apf_add(&sumyy, &sumyy, &tmp);
                }
                
                apf_from_int(&count, n);
                apf_div(&xbar, &sumx, &count);
                apf_div(&ybar, &sumy, &count);
                
                /* slope = (sumxy - n*xbar*ybar) / (sumxx - n*xbar^2) */
                {
                    apf num, den;
                    apf_mul(&tmp, &xbar, &ybar);
                    apf_mul(&tmp, &tmp, &count);
                    apf_sub(&num, &sumxy, &tmp);
                    apf_mul(&tmp, &xbar, &xbar);
                    apf_mul(&tmp, &tmp, &count);
                    apf_sub(&den, &sumxx, &tmp);
                    if (!apf_is_zero(&den)) {
                        apf_div(&slope, &num, &den);
                    } else {
                        apf_zero(&slope);
                    }
                }
                
                /* intercept = ybar - slope*xbar */
                apf_mul(&tmp, &slope, &xbar);
                apf_sub(&intercept, &ybar, &tmp);
                
                /* R-squared: 1 - SS_res/SS_tot */
                apf_zero(&ss_tot);
                apf_zero(&ss_res);
                for (i = 0; i < n; i++) {
                    int ri = i % pv_arg.v.matrix.rows;
                    int ci = i / pv_arg.v.matrix.rows;
                    int rj = i % y_arg.v.matrix.rows;
                    int cj = i / y_arg.v.matrix.rows;
                    apf x, y, y_pred, diff;
                    apf_copy(&x, &MAT_AT(&pv_arg.v.matrix, ri, ci).re);
                    apf_copy(&y, &MAT_AT(&y_arg.v.matrix, rj, cj).re);
                    
                    apf_mul(&y_pred, &slope, &x);
                    apf_add(&y_pred, &y_pred, &intercept);
                    
                    apf_sub(&diff, &y, &ybar);
                    apf_mul(&tmp, &diff, &diff);
                    apf_add(&ss_tot, &ss_tot, &tmp);
                    
                    apf_sub(&diff, &y, &y_pred);
                    apf_mul(&tmp, &diff, &diff);
                    apf_add(&ss_res, &ss_res, &tmp);
                }
                
                if (!apf_is_zero(&ss_tot)) {
                    apf one;
                    apf_from_int(&one, 1);
                    apf_div(&tmp, &ss_res, &ss_tot);
                    apf_sub(&r2, &one, &tmp);
                } else {
                    apf_from_int(&r2, 1);
                }
                
                /* Return [slope, intercept, r2] */
                result->type = VAL_MATRIX;
                mat_zero(&result->v.matrix, 1, 3);
                apf_copy(&MAT_AT(&result->v.matrix, 0, 0).re, &slope);
                apf_copy(&MAT_AT(&result->v.matrix, 0, 1).re, &intercept);
                apf_copy(&MAT_AT(&result->v.matrix, 0, 2).re, &r2);
                return 1;
            }
            
            /* cumtrapz(y) or cumtrapz(x, y) - cumulative trapezoidal integration */
            if (str_eq(name, "cumtrapz")) {
                value_t y_arg;
                int n, i, has_x = 0;
                apf *x_vals, *y_vals;
                apf sum, half, dx, avg;
                
                next_token();
                if (!parse_value(&pv_arg)) return 0;
                
                if (current_token.type == TOK_COMMA) {
                    next_token();
                    if (!parse_value(&y_arg)) return 0;
                    has_x = 1;
                } else {
                    y_arg = pv_arg;
                    has_x = 0;
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
                if (y_arg.type == VAL_SCALAR) {
                    apfc saved = y_arg.v.scalar;
                    mat_zero(&y_arg.v.matrix, 1, 1);
                    MAT_AT(&y_arg.v.matrix, 0, 0) = saved;
                    y_arg.type = VAL_MATRIX;
                }
                
                n = y_arg.v.matrix.rows * y_arg.v.matrix.cols;
                
                x_vals = (apf *)malloc(n * sizeof(apf));
                y_vals = (apf *)malloc(n * sizeof(apf));
                
                if (!x_vals || !y_vals) {
                    free(x_vals); free(y_vals);
                    printf("Error: memory allocation failed\n");
                    return 0;
                }
                
                for (i = 0; i < n; i++) {
                    int ri = i % y_arg.v.matrix.rows;
                    int ci = i / y_arg.v.matrix.rows;
                    apf_copy(&y_vals[i], &MAT_AT(&y_arg.v.matrix, ri, ci).re);
                    if (has_x) {
                        ri = i % pv_arg.v.matrix.rows;
                        ci = i / pv_arg.v.matrix.rows;
                        apf_copy(&x_vals[i], &MAT_AT(&pv_arg.v.matrix, ri, ci).re);
                    } else {
                        apf_from_int(&x_vals[i], i);
                    }
                }
                
                result->type = VAL_MATRIX;
                mat_zero(&result->v.matrix, y_arg.v.matrix.rows, y_arg.v.matrix.cols);
                
                apf_zero(&sum);
                apf_from_int(&half, 1);
                {
                    apf two;
                    apf_from_int(&two, 2);
                    apf_div(&half, &half, &two);
                }
                
                /* First value is 0 */
                apf_zero(&MAT_AT(&result->v.matrix, 0, 0).re);
                
                for (i = 1; i < n; i++) {
                    int ri = i % result->v.matrix.rows;
                    int ci = i / result->v.matrix.rows;
                    
                    apf_sub(&dx, &x_vals[i], &x_vals[i-1]);
                    apf_add(&avg, &y_vals[i], &y_vals[i-1]);
                    apf_mul(&avg, &avg, &half);
                    apf_mul(&avg, &avg, &dx);
                    apf_add(&sum, &sum, &avg);
                    apf_copy(&MAT_AT(&result->v.matrix, ri, ci).re, &sum);
                    apf_zero(&MAT_AT(&result->v.matrix, ri, ci).im);
                }
                
                free(x_vals);
                free(y_vals);
                return 1;
            }
            
            /* quantile(x, p) - alias for prctile but p is 0-1 instead of 0-100 */
            if (str_eq(name, "quantile")) {
                apfc p_val;
                int n, i, idx;
                apf p, p100;
                
                next_token();
                if (!parse_value(&pv_arg)) return 0;
                if (current_token.type != TOK_COMMA) {
                    printf("Error: quantile requires (x, p)\n");
                    return 0;
                }
                next_token();
                if (!parse_expr(&p_val)) return 0;
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
                
                /* Sort the data */
                mat_copy(&pv_tmp_mat, &pv_arg.v.matrix);
                for (i = 0; i < n - 1; i++) {
                    int j;
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
                
                /* Convert p (0-1) to index */
                apf_copy(&p, &p_val.re);
                apf_from_int(&p100, n - 1);
                apf_mul(&p, &p, &p100);
                idx = apf_to_long(&p);
                if (idx < 0) idx = 0;
                if (idx >= n) idx = n - 1;
                
                {
                    int ri = idx % pv_tmp_mat.rows;
                    int ci = idx / pv_tmp_mat.rows;
                    result->type = VAL_SCALAR;
                    result->v.scalar = MAT_AT(&pv_tmp_mat, ri, ci);
                }
                return 1;
            }
            
            /* fillmissing(v, method) - fill NaN values */
            if (str_eq(name, "fillmissing")) {
                int n, i, j, count;
                apf mean_val, sum, cnt, prev;
                int has_method = 0;
                char method[32] = "linear";
                
                next_token();
                if (!parse_value(&pv_arg)) return 0;
                
                if (current_token.type == TOK_COMMA) {
                    next_token();
                    /* Parse method string */
                    if (current_token.type == TOK_STRING) {
                        strncpy(method, current_token.str_value, 31);
                        method[31] = '\0';
                        next_token();
                        has_method = 1;
                    }
                }
                if (current_token.type != TOK_RPAREN) {
                    printf("Error: expected ')'\n");
                    return 0;
                }
                next_token();
                (void)has_method;
                
                if (pv_arg.type == VAL_SCALAR) {
                    apfc saved = pv_arg.v.scalar;
                    mat_zero(&pv_arg.v.matrix, 1, 1);
                    MAT_AT(&pv_arg.v.matrix, 0, 0) = saved;
                    pv_arg.type = VAL_MATRIX;
                }
                
                n = pv_arg.v.matrix.rows * pv_arg.v.matrix.cols;
                result->type = VAL_MATRIX;
                mat_copy(&result->v.matrix, &pv_arg.v.matrix);
                
                if (str_eq(method, "mean") || str_eq(method, "constant")) {
                    /* Fill with mean of non-NaN values */
                    apf_zero(&sum);
                    count = 0;
                    for (i = 0; i < n; i++) {
                        int ri = i % pv_arg.v.matrix.rows;
                        int ci = i / pv_arg.v.matrix.rows;
                        if (MAT_AT(&pv_arg.v.matrix, ri, ci).re.cls != APF_CLASS_NAN) {
                            apf_add(&sum, &sum, &MAT_AT(&pv_arg.v.matrix, ri, ci).re);
                            count++;
                        }
                    }
                    if (count > 0) {
                        apf_from_int(&cnt, count);
                        apf_div(&mean_val, &sum, &cnt);
                    } else {
                        apf_zero(&mean_val);
                    }
                    for (i = 0; i < n; i++) {
                        int ri = i % result->v.matrix.rows;
                        int ci = i / result->v.matrix.rows;
                        if (MAT_AT(&result->v.matrix, ri, ci).re.cls == APF_CLASS_NAN) {
                            apf_copy(&MAT_AT(&result->v.matrix, ri, ci).re, &mean_val);
                        }
                    }
                } else if (str_eq(method, "previous")) {
                    /* Fill with previous non-NaN value */
                    apf_zero(&prev);
                    for (i = 0; i < n; i++) {
                        int ri = i % result->v.matrix.rows;
                        int ci = i / result->v.matrix.rows;
                        if (MAT_AT(&result->v.matrix, ri, ci).re.cls == APF_CLASS_NAN) {
                            apf_copy(&MAT_AT(&result->v.matrix, ri, ci).re, &prev);
                        } else {
                            apf_copy(&prev, &MAT_AT(&result->v.matrix, ri, ci).re);
                        }
                    }
                } else {
                    /* linear interpolation (default) */
                    for (i = 0; i < n; i++) {
                        int ri = i % result->v.matrix.rows;
                        int ci = i / result->v.matrix.rows;
                        if (MAT_AT(&result->v.matrix, ri, ci).re.cls == APF_CLASS_NAN) {
                            /* Find prev and next non-NaN */
                            int prev_idx = -1, next_idx = -1;
                            apf y0, y1, t, one_minus_t, diff;
                            for (j = i - 1; j >= 0; j--) {
                                int rj = j % pv_arg.v.matrix.rows;
                                int cj = j / pv_arg.v.matrix.rows;
                                if (MAT_AT(&pv_arg.v.matrix, rj, cj).re.cls != APF_CLASS_NAN) {
                                    prev_idx = j;
                                    break;
                                }
                            }
                            for (j = i + 1; j < n; j++) {
                                int rj = j % pv_arg.v.matrix.rows;
                                int cj = j / pv_arg.v.matrix.rows;
                                if (MAT_AT(&pv_arg.v.matrix, rj, cj).re.cls != APF_CLASS_NAN) {
                                    next_idx = j;
                                    break;
                                }
                            }
                            if (prev_idx >= 0 && next_idx >= 0) {
                                int rp = prev_idx % pv_arg.v.matrix.rows;
                                int cp = prev_idx / pv_arg.v.matrix.rows;
                                int rn = next_idx % pv_arg.v.matrix.rows;
                                int cn = next_idx / pv_arg.v.matrix.rows;
                                apf_copy(&y0, &MAT_AT(&pv_arg.v.matrix, rp, cp).re);
                                apf_copy(&y1, &MAT_AT(&pv_arg.v.matrix, rn, cn).re);
                                /* t = (i - prev_idx) / (next_idx - prev_idx) */
                                apf_from_int(&t, i - prev_idx);
                                apf_from_int(&diff, next_idx - prev_idx);
                                apf_div(&t, &t, &diff);
                                /* result = y0 * (1-t) + y1 * t */
                                apf_from_int(&one_minus_t, 1);
                                apf_sub(&one_minus_t, &one_minus_t, &t);
                                apf_mul(&y0, &y0, &one_minus_t);
                                apf_mul(&y1, &y1, &t);
                                apf_add(&MAT_AT(&result->v.matrix, ri, ci).re, &y0, &y1);
                            } else if (prev_idx >= 0) {
                                int rp = prev_idx % pv_arg.v.matrix.rows;
                                int cp = prev_idx / pv_arg.v.matrix.rows;
                                apf_copy(&MAT_AT(&result->v.matrix, ri, ci).re, 
                                        &MAT_AT(&pv_arg.v.matrix, rp, cp).re);
                            } else if (next_idx >= 0) {
                                int rn = next_idx % pv_arg.v.matrix.rows;
                                int cn = next_idx / pv_arg.v.matrix.rows;
                                apf_copy(&MAT_AT(&result->v.matrix, ri, ci).re,
                                        &MAT_AT(&pv_arg.v.matrix, rn, cn).re);
                            }
                        }
                    }
                }
                return 1;
            }
            
            /* rmmissing(v) - remove NaN values */
            if (str_eq(name, "rmmissing")) {
                int n, i, count;
                
                next_token();
                if (!parse_value(&pv_arg)) return 0;
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
                result->type = VAL_MATRIX;
                
                /* Count non-NaN values */
                count = 0;
                for (i = 0; i < n; i++) {
                    int ri = i % pv_arg.v.matrix.rows;
                    int ci = i / pv_arg.v.matrix.rows;
                    if (MAT_AT(&pv_arg.v.matrix, ri, ci).re.cls != APF_CLASS_NAN) {
                        count++;
                    }
                }
                
                if (count == 0) {
                    mat_zero(&result->v.matrix, 1, 1);
                    return 1;
                }
                
                /* Copy non-NaN values */
                if (pv_arg.v.matrix.cols == 1) {
                    mat_zero(&result->v.matrix, count, 1);
                } else {
                    mat_zero(&result->v.matrix, 1, count);
                }
                
                count = 0;
                for (i = 0; i < n; i++) {
                    int ri = i % pv_arg.v.matrix.rows;
                    int ci = i / pv_arg.v.matrix.rows;
                    if (MAT_AT(&pv_arg.v.matrix, ri, ci).re.cls != APF_CLASS_NAN) {
                        if (pv_arg.v.matrix.cols == 1) {
                            MAT_AT(&result->v.matrix, count, 0) = MAT_AT(&pv_arg.v.matrix, ri, ci);
                        } else {
                            MAT_AT(&result->v.matrix, 0, count) = MAT_AT(&pv_arg.v.matrix, ri, ci);
                        }
                        count++;
                    }
                }
                return 1;
            }
            
            /* smoothdata(v, method, window) - smooth data */
            if (str_eq(name, "smoothdata")) {
                int n, i, j, window;
                char method[32] = "movmean";
                apfc window_val;
                
                next_token();
                if (!parse_value(&pv_arg)) return 0;
                
                window = 5;  /* default window */
                
                if (current_token.type == TOK_COMMA) {
                    next_token();
                    if (current_token.type == TOK_STRING) {
                        strncpy(method, current_token.str_value, 31);
                        method[31] = '\0';
                        next_token();
                    }
                    if (current_token.type == TOK_COMMA) {
                        next_token();
                        if (!parse_expr(&window_val)) return 0;
                        window = apf_to_long(&window_val.re);
                    }
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
                
                n = pv_arg.v.matrix.rows * pv_arg.v.matrix.cols;
                if (window < 1) window = 1;
                if (window > n) window = n;
                
                result->type = VAL_MATRIX;
                mat_zero(&result->v.matrix, pv_arg.v.matrix.rows, pv_arg.v.matrix.cols);
                
                if (str_eq(method, "gaussian")) {
                    /* Gaussian smoothing - use weights */
                    for (i = 0; i < n; i++) {
                        apf sum, weight_sum, w, diff_sq, sigma, two;
                        int ri = i % pv_arg.v.matrix.rows;
                        int ci = i / pv_arg.v.matrix.rows;
                        int half = window / 2;
                        
                        apf_zero(&sum);
                        apf_zero(&weight_sum);
                        apf_from_int(&sigma, window / 4 + 1);
                        apf_from_int(&two, 2);
                        
                        for (j = -half; j <= half; j++) {
                            int idx = i + j;
                            if (idx >= 0 && idx < n) {
                                int rj = idx % pv_arg.v.matrix.rows;
                                int cj = idx / pv_arg.v.matrix.rows;
                                apf jval, tmp;
                                /* w = exp(-j^2 / (2*sigma^2)) */
                                apf_from_int(&jval, j);
                                apf_mul(&diff_sq, &jval, &jval);
                                apf_mul(&tmp, &sigma, &sigma);
                                apf_mul(&tmp, &tmp, &two);
                                apf_div(&diff_sq, &diff_sq, &tmp);
                                apf_neg(&diff_sq, &diff_sq);
                                apfx_exp(&w, &diff_sq);
                                
                                apf_mul(&tmp, &w, &MAT_AT(&pv_arg.v.matrix, rj, cj).re);
                                apf_add(&sum, &sum, &tmp);
                                apf_add(&weight_sum, &weight_sum, &w);
                            }
                        }
                        apf_div(&MAT_AT(&result->v.matrix, ri, ci).re, &sum, &weight_sum);
                        apf_zero(&MAT_AT(&result->v.matrix, ri, ci).im);
                    }
                } else {
                    /* movmean (default) */
                    for (i = 0; i < n; i++) {
                        apf sum, cnt_apf;
                        int ri = i % pv_arg.v.matrix.rows;
                        int ci = i / pv_arg.v.matrix.rows;
                        int half = window / 2;
                        int start = i - half;
                        int end = i + half;
                        int cnt = 0;
                        
                        if (start < 0) start = 0;
                        if (end >= n) end = n - 1;
                        
                        apf_zero(&sum);
                        for (j = start; j <= end; j++) {
                            int rj = j % pv_arg.v.matrix.rows;
                            int cj = j / pv_arg.v.matrix.rows;
                            apf_add(&sum, &sum, &MAT_AT(&pv_arg.v.matrix, rj, cj).re);
                            cnt++;
                        }
                        apf_from_int(&cnt_apf, cnt);
                        apf_div(&MAT_AT(&result->v.matrix, ri, ci).re, &sum, &cnt_apf);
                        apf_zero(&MAT_AT(&result->v.matrix, ri, ci).im);
                    }
                }
                return 1;
            }
            
            /* winsorize(v, pct) - winsorize data by trimming percentiles */
            if (str_eq(name, "winsorize")) {
                apfc pct_val;
                int n, i, low_idx, high_idx;
                apf low_val, high_val, pct_dec, cnt;
                
                next_token();
                if (!parse_value(&pv_arg)) return 0;
                
                /* Default 5% winsorization */
                apf_from_int(&pct_dec, 5);
                {
                    apf hundred;
                    apf_from_int(&hundred, 100);
                    apf_div(&pct_dec, &pct_dec, &hundred);
                }
                
                if (current_token.type == TOK_COMMA) {
                    next_token();
                    if (!parse_expr(&pct_val)) return 0;
                    apf_copy(&pct_dec, &pct_val.re);
                    {
                        apf hundred;
                        apf_from_int(&hundred, 100);
                        apf_div(&pct_dec, &pct_dec, &hundred);
                    }
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
                
                n = pv_arg.v.matrix.rows * pv_arg.v.matrix.cols;
                
                /* Sort data to find percentiles */
                mat_copy(&pv_tmp_mat, &pv_arg.v.matrix);
                for (i = 0; i < n - 1; i++) {
                    int j;
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
                
                /* Calculate indices for low and high bounds */
                apf_from_int(&cnt, n);
                apf_mul(&cnt, &cnt, &pct_dec);
                low_idx = apf_to_long(&cnt);
                if (low_idx < 0) low_idx = 0;
                if (low_idx >= n) low_idx = n - 1;
                high_idx = n - 1 - low_idx;
                if (high_idx < low_idx) high_idx = low_idx;
                
                {
                    int rl = low_idx % pv_tmp_mat.rows;
                    int cl = low_idx / pv_tmp_mat.rows;
                    int rh = high_idx % pv_tmp_mat.rows;
                    int ch = high_idx / pv_tmp_mat.rows;
                    apf_copy(&low_val, &MAT_AT(&pv_tmp_mat, rl, cl).re);
                    apf_copy(&high_val, &MAT_AT(&pv_tmp_mat, rh, ch).re);
                }
                
                result->type = VAL_MATRIX;
                mat_zero(&result->v.matrix, pv_arg.v.matrix.rows, pv_arg.v.matrix.cols);
                
                for (i = 0; i < n; i++) {
                    int ri = i % pv_arg.v.matrix.rows;
                    int ci = i / pv_arg.v.matrix.rows;
                    
                    if (apf_cmp(&MAT_AT(&pv_arg.v.matrix, ri, ci).re, &low_val) < 0) {
                        apf_copy(&MAT_AT(&result->v.matrix, ri, ci).re, &low_val);
                    } else if (apf_cmp(&MAT_AT(&pv_arg.v.matrix, ri, ci).re, &high_val) > 0) {
                        apf_copy(&MAT_AT(&result->v.matrix, ri, ci).re, &high_val);
                    } else {
                        apf_copy(&MAT_AT(&result->v.matrix, ri, ci).re, 
                                &MAT_AT(&pv_arg.v.matrix, ri, ci).re);
                    }
                    apf_zero(&MAT_AT(&result->v.matrix, ri, ci).im);
                }
                return 1;
            }
            
            /* bounds(v) - return [min, max] of vector */
            if (str_eq(name, "bounds") || str_eq(name, "minmax")) {
                int n, i;
                apf min_val, max_val;
                int first = 1;
                
                next_token();
                if (!parse_value(&pv_arg)) return 0;
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
                
                for (i = 0; i < n; i++) {
                    int ri = i % pv_arg.v.matrix.rows;
                    int ci = i / pv_arg.v.matrix.rows;
                    if (first) {
                        apf_copy(&min_val, &MAT_AT(&pv_arg.v.matrix, ri, ci).re);
                        apf_copy(&max_val, &MAT_AT(&pv_arg.v.matrix, ri, ci).re);
                        first = 0;
                    } else {
                        if (apf_cmp(&MAT_AT(&pv_arg.v.matrix, ri, ci).re, &min_val) < 0) {
                            apf_copy(&min_val, &MAT_AT(&pv_arg.v.matrix, ri, ci).re);
                        }
                        if (apf_cmp(&MAT_AT(&pv_arg.v.matrix, ri, ci).re, &max_val) > 0) {
                            apf_copy(&max_val, &MAT_AT(&pv_arg.v.matrix, ri, ci).re);
                        }
                    }
                }
                
                result->type = VAL_MATRIX;
                mat_zero(&result->v.matrix, 1, 2);
                apf_copy(&MAT_AT(&result->v.matrix, 0, 0).re, &min_val);
                apf_zero(&MAT_AT(&result->v.matrix, 0, 0).im);
                apf_copy(&MAT_AT(&result->v.matrix, 0, 1).re, &max_val);
                apf_zero(&MAT_AT(&result->v.matrix, 0, 1).im);
                return 1;
            }
            
            /* quantile(v, q) - compute quantile(s) */
            if (str_eq(name, "quantile")) {
                apfc q_val;
                int n, i, idx;
                apf q, cnt, pos, frac, one_minus_frac, low_v, high_v;
                
                next_token();
                if (!parse_value(&pv_arg)) return 0;
                if (current_token.type != TOK_COMMA) {
                    printf("Error: quantile requires (v, q)\n");
                    return 0;
                }
                next_token();
                if (!parse_expr(&q_val)) return 0;
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
                
                apf_copy(&q, &q_val.re);
                n = pv_arg.v.matrix.rows * pv_arg.v.matrix.cols;
                
                /* Sort data */
                mat_copy(&pv_tmp_mat, &pv_arg.v.matrix);
                for (i = 0; i < n - 1; i++) {
                    int j;
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
                
                /* Calculate position */
                apf_from_int(&cnt, n - 1);
                apf_mul(&pos, &q, &cnt);
                idx = apf_to_long(&pos);
                if (idx < 0) idx = 0;
                if (idx >= n - 1) idx = n - 2;
                
                /* Linear interpolation */
                {
                    apf idx_apf;
                    apf_from_int(&idx_apf, idx);
                    apf_sub(&frac, &pos, &idx_apf);
                }
                apf_from_int(&one_minus_frac, 1);
                apf_sub(&one_minus_frac, &one_minus_frac, &frac);
                
                {
                    int rl = idx % pv_tmp_mat.rows;
                    int cl = idx / pv_tmp_mat.rows;
                    int rh = (idx + 1) % pv_tmp_mat.rows;
                    int ch = (idx + 1) / pv_tmp_mat.rows;
                    apf_copy(&low_v, &MAT_AT(&pv_tmp_mat, rl, cl).re);
                    apf_copy(&high_v, &MAT_AT(&pv_tmp_mat, rh, ch).re);
                }
                
                apf_mul(&low_v, &low_v, &one_minus_frac);
                apf_mul(&high_v, &high_v, &frac);
                
                result->type = VAL_SCALAR;
                apf_add(&result->v.scalar.re, &low_v, &high_v);
                apf_zero(&result->v.scalar.im);
                return 1;
            }
            
            /* clip(v, lo, hi) - clip values to range [lo, hi] */
            if (str_eq(name, "clip")) {
                apfc lo_val, hi_val;
                int n, i;
                
                next_token();
                if (!parse_value(&pv_arg)) return 0;
                if (current_token.type != TOK_COMMA) {
                    printf("Error: clip requires (v, lo, hi)\n");
                    return 0;
                }
                next_token();
                if (!parse_expr(&lo_val)) return 0;
                if (current_token.type != TOK_COMMA) {
                    printf("Error: clip requires (v, lo, hi)\n");
                    return 0;
                }
                next_token();
                if (!parse_expr(&hi_val)) return 0;
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
                result->type = VAL_MATRIX;
                mat_zero(&result->v.matrix, pv_arg.v.matrix.rows, pv_arg.v.matrix.cols);
                
                for (i = 0; i < n; i++) {
                    int ri = i % pv_arg.v.matrix.rows;
                    int ci = i / pv_arg.v.matrix.rows;
                    
                    if (apf_cmp(&MAT_AT(&pv_arg.v.matrix, ri, ci).re, &lo_val.re) < 0) {
                        apf_copy(&MAT_AT(&result->v.matrix, ri, ci).re, &lo_val.re);
                    } else if (apf_cmp(&MAT_AT(&pv_arg.v.matrix, ri, ci).re, &hi_val.re) > 0) {
                        apf_copy(&MAT_AT(&result->v.matrix, ri, ci).re, &hi_val.re);
                    } else {
                        apf_copy(&MAT_AT(&result->v.matrix, ri, ci).re, 
                                &MAT_AT(&pv_arg.v.matrix, ri, ci).re);
                    }
                    apf_zero(&MAT_AT(&result->v.matrix, ri, ci).im);
                }
                return 1;
            }
            
            /* npv(rate, cashflows) - Net Present Value */
            if (str_eq(name, "npv")) {
                apfc rate_val;
                int n, i;
                apf rate, one, discount, pv, cf, sum;
                
                next_token();
                if (!parse_expr(&rate_val)) return 0;
                if (current_token.type != TOK_COMMA) {
                    printf("Error: npv requires (rate, cashflows)\n");
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
                    apfc saved = pv_arg.v.scalar;
                    mat_zero(&pv_arg.v.matrix, 1, 1);
                    MAT_AT(&pv_arg.v.matrix, 0, 0) = saved;
                    pv_arg.type = VAL_MATRIX;
                }
                
                n = pv_arg.v.matrix.rows * pv_arg.v.matrix.cols;
                apf_copy(&rate, &rate_val.re);
                apf_from_int(&one, 1);
                apf_add(&discount, &one, &rate);  /* 1 + rate */
                
                apf_zero(&sum);
                apf_copy(&pv, &one);  /* discount factor starts at 1 */
                
                for (i = 0; i < n; i++) {
                    int ri = i % pv_arg.v.matrix.rows;
                    int ci = i / pv_arg.v.matrix.rows;
                    apf term;
                    
                    apf_div(&pv, &pv, &discount);  /* pv = pv / (1+rate) */
                    apf_copy(&cf, &MAT_AT(&pv_arg.v.matrix, ri, ci).re);
                    apf_mul(&term, &cf, &pv);
                    apf_add(&sum, &sum, &term);
                }
                
                result->type = VAL_SCALAR;
                apf_copy(&result->v.scalar.re, &sum);
                apf_zero(&result->v.scalar.im);
                return 1;
            }
            
            /* irr(cashflows) - Internal Rate of Return using Newton-Raphson */
            if (str_eq(name, "irr")) {
                int n, i, iter;
                apf rate, one, npv_val, dnpv, discount, cf, term, dterm;
                apf delta, tol, neg_one;
                
                next_token();
                if (!parse_value(&pv_arg)) return 0;
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
                
                /* Initial guess: 10% */
                apf_from_double(&rate, 0.1);
                apf_from_int(&one, 1);
                apf_from_int(&neg_one, -1);
                apf_from_double(&tol, 1e-10);
                
                /* Newton-Raphson iteration */
                for (iter = 0; iter < 100; iter++) {
                    apf_zero(&npv_val);
                    apf_zero(&dnpv);
                    apf_add(&discount, &one, &rate);
                    
                    for (i = 0; i < n; i++) {
                        int ri = i % pv_arg.v.matrix.rows;
                        int ci = i / pv_arg.v.matrix.rows;
                        apf power, iapf;
                        
                        apf_copy(&cf, &MAT_AT(&pv_arg.v.matrix, ri, ci).re);
                        
                        /* term = cf / (1+rate)^(i+1) */
                        apf_from_int(&iapf, i + 1);
                        apfx_pow(&power, &discount, &iapf);
                        apf_div(&term, &cf, &power);
                        apf_add(&npv_val, &npv_val, &term);
                        
                        /* dterm = -(i+1) * cf / (1+rate)^(i+2) */
                        apf_from_int(&iapf, i + 2);
                        apfx_pow(&power, &discount, &iapf);
                        apf_div(&dterm, &cf, &power);
                        apf_from_int(&iapf, -(i + 1));
                        apf_mul(&dterm, &dterm, &iapf);
                        apf_add(&dnpv, &dnpv, &dterm);
                    }
                    
                    /* Check convergence */
                    {
                        apf abs_npv;
                        apf_abs(&abs_npv, &npv_val);
                        if (apf_cmp(&abs_npv, &tol) < 0) break;
                    }
                    
                    /* rate = rate - npv/dnpv */
                    if (!apf_is_zero(&dnpv)) {
                        apf_div(&delta, &npv_val, &dnpv);
                        apf_sub(&rate, &rate, &delta);
                    }
                }
                
                result->type = VAL_SCALAR;
                apf_copy(&result->v.scalar.re, &rate);
                apf_zero(&result->v.scalar.im);
                return 1;
            }
            
            /* payback(cashflows) - Payback period */
            if (str_eq(name, "payback")) {
                int n, i;
                apf cumsum, zero_apf;
                
                next_token();
                if (!parse_value(&pv_arg)) return 0;
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
                apf_zero(&cumsum);
                apf_zero(&zero_apf);
                
                for (i = 0; i < n; i++) {
                    int ri = i % pv_arg.v.matrix.rows;
                    int ci = i / pv_arg.v.matrix.rows;
                    apf prev_cumsum;
                    
                    apf_copy(&prev_cumsum, &cumsum);
                    apf_add(&cumsum, &cumsum, &MAT_AT(&pv_arg.v.matrix, ri, ci).re);
                    
                    /* Check if we crossed zero */
                    if (apf_cmp(&prev_cumsum, &zero_apf) < 0 && apf_cmp(&cumsum, &zero_apf) >= 0) {
                        /* Linear interpolation for fractional period */
                        apf frac, period, neg_prev;
                        apf_neg(&neg_prev, &prev_cumsum);
                        apf_sub(&frac, &cumsum, &prev_cumsum);
                        apf_div(&frac, &neg_prev, &frac);
                        apf_from_int(&period, i);
                        apf_add(&result->v.scalar.re, &period, &frac);
                        apf_zero(&result->v.scalar.im);
                        result->type = VAL_SCALAR;
                        return 1;
                    }
                }
                
                /* Never paid back */
                result->type = VAL_SCALAR;
                apf_set_inf(&result->v.scalar.re, 0);
                apf_zero(&result->v.scalar.im);
                return 1;
            }
            
            /* cagr(start, end, years) - Compound Annual Growth Rate */
            if (str_eq(name, "cagr")) {
                apfc start_val, end_val, years_val;
                apf ratio, one, inv_years, growth;
                
                next_token();
                if (!parse_expr(&start_val)) return 0;
                if (current_token.type != TOK_COMMA) {
                    printf("Error: cagr requires (start, end, years)\n");
                    return 0;
                }
                next_token();
                if (!parse_expr(&end_val)) return 0;
                if (current_token.type != TOK_COMMA) {
                    printf("Error: cagr requires (start, end, years)\n");
                    return 0;
                }
                next_token();
                if (!parse_expr(&years_val)) return 0;
                if (current_token.type != TOK_RPAREN) {
                    printf("Error: expected ')'\n");
                    return 0;
                }
                next_token();
                
                /* CAGR = (end/start)^(1/years) - 1 */
                apf_div(&ratio, &end_val.re, &start_val.re);
                apf_from_int(&one, 1);
                apf_div(&inv_years, &one, &years_val.re);
                apfx_pow(&growth, &ratio, &inv_years);
                apf_sub(&result->v.scalar.re, &growth, &one);
                apf_zero(&result->v.scalar.im);
                result->type = VAL_SCALAR;
                return 1;
            }
            
            /* compound(principal, rate, n, t) - Compound interest */
            if (str_eq(name, "compound")) {
                apfc principal, rate, n_periods, time_val;
                apf base, exponent, result_val, one;
                
                next_token();
                if (!parse_expr(&principal)) return 0;
                if (current_token.type != TOK_COMMA) {
                    printf("Error: compound requires (principal, rate, n, t)\n");
                    return 0;
                }
                next_token();
                if (!parse_expr(&rate)) return 0;
                if (current_token.type != TOK_COMMA) {
                    printf("Error: compound requires (principal, rate, n, t)\n");
                    return 0;
                }
                next_token();
                if (!parse_expr(&n_periods)) return 0;
                if (current_token.type != TOK_COMMA) {
                    printf("Error: compound requires (principal, rate, n, t)\n");
                    return 0;
                }
                next_token();
                if (!parse_expr(&time_val)) return 0;
                if (current_token.type != TOK_RPAREN) {
                    printf("Error: expected ')'\n");
                    return 0;
                }
                next_token();
                
                /* A = P * (1 + r/n)^(n*t) */
                apf_from_int(&one, 1);
                apf_div(&base, &rate.re, &n_periods.re);
                apf_add(&base, &one, &base);
                apf_mul(&exponent, &n_periods.re, &time_val.re);
                apfx_pow(&result_val, &base, &exponent);
                apf_mul(&result->v.scalar.re, &principal.re, &result_val);
                apf_zero(&result->v.scalar.im);
                result->type = VAL_SCALAR;
                return 1;
            }
            
            /* roi(gain, cost) - Return on Investment */
            if (str_eq(name, "roi")) {
                apfc gain, cost;
                apf profit, hundred;
                
                next_token();
                if (!parse_expr(&gain)) return 0;
                if (current_token.type != TOK_COMMA) {
                    printf("Error: roi requires (gain, cost)\n");
                    return 0;
                }
                next_token();
                if (!parse_expr(&cost)) return 0;
                if (current_token.type != TOK_RPAREN) {
                    printf("Error: expected ')'\n");
                    return 0;
                }
                next_token();
                
                /* ROI = (gain - cost) / cost * 100 */
                apf_sub(&profit, &gain.re, &cost.re);
                apf_div(&result->v.scalar.re, &profit, &cost.re);
                apf_from_int(&hundred, 100);
                apf_mul(&result->v.scalar.re, &result->v.scalar.re, &hundred);
                apf_zero(&result->v.scalar.im);
                result->type = VAL_SCALAR;
                return 1;
            }
            
            /* crossover(x, y) - Find crossover points of two series */
            if (str_eq(name, "crossover")) {
                value_t arg2;
                int n, i, count;
                
                next_token();
                if (!parse_value(&pv_arg)) return 0;
                if (current_token.type != TOK_COMMA) {
                    printf("Error: crossover requires (x, y)\n");
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
                
                n = pv_arg.v.matrix.rows * pv_arg.v.matrix.cols;
                if (n > arg2.v.matrix.rows * arg2.v.matrix.cols) {
                    n = arg2.v.matrix.rows * arg2.v.matrix.cols;
                }
                
                /* Count crossovers first */
                count = 0;
                for (i = 1; i < n; i++) {
                    int r1 = (i-1) % pv_arg.v.matrix.rows;
                    int c1 = (i-1) / pv_arg.v.matrix.rows;
                    int r2 = i % pv_arg.v.matrix.rows;
                    int c2 = i / pv_arg.v.matrix.rows;
                    apf diff_prev, diff_curr;
                    
                    apf_sub(&diff_prev, &MAT_AT(&pv_arg.v.matrix, r1, c1).re,
                            &MAT_AT(&arg2.v.matrix, r1, c1).re);
                    apf_sub(&diff_curr, &MAT_AT(&pv_arg.v.matrix, r2, c2).re,
                            &MAT_AT(&arg2.v.matrix, r2, c2).re);
                    
                    /* Check for sign change */
                    if ((diff_prev.sign != diff_curr.sign) || 
                        (apf_is_zero(&diff_prev) != apf_is_zero(&diff_curr))) {
                        count++;
                    }
                }
                
                if (count == 0) {
                    result->type = VAL_MATRIX;
                    mat_zero(&result->v.matrix, 1, 1);
                    return 1;
                }
                
                result->type = VAL_MATRIX;
                mat_zero(&result->v.matrix, count, 1);
                
                count = 0;
                for (i = 1; i < n; i++) {
                    int r1 = (i-1) % pv_arg.v.matrix.rows;
                    int c1 = (i-1) / pv_arg.v.matrix.rows;
                    int r2 = i % pv_arg.v.matrix.rows;
                    int c2 = i / pv_arg.v.matrix.rows;
                    apf diff_prev, diff_curr;
                    
                    apf_sub(&diff_prev, &MAT_AT(&pv_arg.v.matrix, r1, c1).re,
                            &MAT_AT(&arg2.v.matrix, r1, c1).re);
                    apf_sub(&diff_curr, &MAT_AT(&pv_arg.v.matrix, r2, c2).re,
                            &MAT_AT(&arg2.v.matrix, r2, c2).re);
                    
                    if ((diff_prev.sign != diff_curr.sign) || 
                        (apf_is_zero(&diff_prev) != apf_is_zero(&diff_curr))) {
                        apf_from_int(&MAT_AT(&result->v.matrix, count, 0).re, i);
                        apf_zero(&MAT_AT(&result->v.matrix, count, 0).im);
                        count++;
                    }
                }
                return 1;
            }
            
            /* peaks(x) - Find local maxima indices */
            if (str_eq(name, "peaks") || str_eq(name, "findpeaks")) {
                int n, i, count;
                
                next_token();
                if (!parse_value(&pv_arg)) return 0;
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
                
                /* Count peaks first */
                count = 0;
                for (i = 1; i < n - 1; i++) {
                    int rp = (i-1) % pv_arg.v.matrix.rows;
                    int cp = (i-1) / pv_arg.v.matrix.rows;
                    int rc = i % pv_arg.v.matrix.rows;
                    int cc = i / pv_arg.v.matrix.rows;
                    int rn = (i+1) % pv_arg.v.matrix.rows;
                    int cn = (i+1) / pv_arg.v.matrix.rows;
                    
                    if (apf_cmp(&MAT_AT(&pv_arg.v.matrix, rc, cc).re,
                               &MAT_AT(&pv_arg.v.matrix, rp, cp).re) > 0 &&
                        apf_cmp(&MAT_AT(&pv_arg.v.matrix, rc, cc).re,
                               &MAT_AT(&pv_arg.v.matrix, rn, cn).re) > 0) {
                        count++;
                    }
                }
                
                if (count == 0) {
                    result->type = VAL_MATRIX;
                    mat_zero(&result->v.matrix, 1, 1);
                    return 1;
                }
                
                result->type = VAL_MATRIX;
                mat_zero(&result->v.matrix, count, 1);
                
                count = 0;
                for (i = 1; i < n - 1; i++) {
                    int rp = (i-1) % pv_arg.v.matrix.rows;
                    int cp = (i-1) / pv_arg.v.matrix.rows;
                    int rc = i % pv_arg.v.matrix.rows;
                    int cc = i / pv_arg.v.matrix.rows;
                    int rn = (i+1) % pv_arg.v.matrix.rows;
                    int cn = (i+1) / pv_arg.v.matrix.rows;
                    
                    if (apf_cmp(&MAT_AT(&pv_arg.v.matrix, rc, cc).re,
                               &MAT_AT(&pv_arg.v.matrix, rp, cp).re) > 0 &&
                        apf_cmp(&MAT_AT(&pv_arg.v.matrix, rc, cc).re,
                               &MAT_AT(&pv_arg.v.matrix, rn, cn).re) > 0) {
                        apf_from_int(&MAT_AT(&result->v.matrix, count, 0).re, i + 1);
                        apf_zero(&MAT_AT(&result->v.matrix, count, 0).im);
                        count++;
                    }
                }
                return 1;
            }
            
            /* valleys(x) - Find local minima indices */
            if (str_eq(name, "valleys") || str_eq(name, "findvalleys")) {
                int n, i, count;
                
                next_token();
                if (!parse_value(&pv_arg)) return 0;
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
                
                /* Count valleys first */
                count = 0;
                for (i = 1; i < n - 1; i++) {
                    int rp = (i-1) % pv_arg.v.matrix.rows;
                    int cp = (i-1) / pv_arg.v.matrix.rows;
                    int rc = i % pv_arg.v.matrix.rows;
                    int cc = i / pv_arg.v.matrix.rows;
                    int rn = (i+1) % pv_arg.v.matrix.rows;
                    int cn = (i+1) / pv_arg.v.matrix.rows;
                    
                    if (apf_cmp(&MAT_AT(&pv_arg.v.matrix, rc, cc).re,
                               &MAT_AT(&pv_arg.v.matrix, rp, cp).re) < 0 &&
                        apf_cmp(&MAT_AT(&pv_arg.v.matrix, rc, cc).re,
                               &MAT_AT(&pv_arg.v.matrix, rn, cn).re) < 0) {
                        count++;
                    }
                }
                
                if (count == 0) {
                    result->type = VAL_MATRIX;
                    mat_zero(&result->v.matrix, 1, 1);
                    return 1;
                }
                
                result->type = VAL_MATRIX;
                mat_zero(&result->v.matrix, count, 1);
                
                count = 0;
                for (i = 1; i < n - 1; i++) {
                    int rp = (i-1) % pv_arg.v.matrix.rows;
                    int cp = (i-1) / pv_arg.v.matrix.rows;
                    int rc = i % pv_arg.v.matrix.rows;
                    int cc = i / pv_arg.v.matrix.rows;
                    int rn = (i+1) % pv_arg.v.matrix.rows;
                    int cn = (i+1) / pv_arg.v.matrix.rows;
                    
                    if (apf_cmp(&MAT_AT(&pv_arg.v.matrix, rc, cc).re,
                               &MAT_AT(&pv_arg.v.matrix, rp, cp).re) < 0 &&
                        apf_cmp(&MAT_AT(&pv_arg.v.matrix, rc, cc).re,
                               &MAT_AT(&pv_arg.v.matrix, rn, cn).re) < 0) {
                        apf_from_int(&MAT_AT(&result->v.matrix, count, 0).re, i + 1);
                        apf_zero(&MAT_AT(&result->v.matrix, count, 0).im);
                        count++;
                    }
                }
                return 1;
            }
            
            /* zscore_vec(x) - Z-score normalization for vectors */
            if (str_eq(name, "zscore_vec") || str_eq(name, "standardize")) {
                int n, i;
                apf mean_val, std_val, sum, sumsq, count, tmp, variance;
                
                next_token();
                if (!parse_value(&pv_arg)) return 0;
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
                
                /* Compute mean */
                apf_zero(&sum);
                for (i = 0; i < n; i++) {
                    int ri = i % pv_arg.v.matrix.rows;
                    int ci = i / pv_arg.v.matrix.rows;
                    apf_add(&sum, &sum, &MAT_AT(&pv_arg.v.matrix, ri, ci).re);
                }
                apf_from_int(&count, n);
                apf_div(&mean_val, &sum, &count);
                
                /* Compute std dev */
                apf_zero(&sumsq);
                for (i = 0; i < n; i++) {
                    int ri = i % pv_arg.v.matrix.rows;
                    int ci = i / pv_arg.v.matrix.rows;
                    apf_sub(&tmp, &MAT_AT(&pv_arg.v.matrix, ri, ci).re, &mean_val);
                    apf_mul(&tmp, &tmp, &tmp);
                    apf_add(&sumsq, &sumsq, &tmp);
                }
                apf_div(&variance, &sumsq, &count);
                apf_sqrt(&std_val, &variance);
                
                result->type = VAL_MATRIX;
                mat_zero(&result->v.matrix, pv_arg.v.matrix.rows, pv_arg.v.matrix.cols);
                
                for (i = 0; i < n; i++) {
                    int ri = i % pv_arg.v.matrix.rows;
                    int ci = i / pv_arg.v.matrix.rows;
                    apf_sub(&tmp, &MAT_AT(&pv_arg.v.matrix, ri, ci).re, &mean_val);
                    if (!apf_is_zero(&std_val)) {
                        apf_div(&MAT_AT(&result->v.matrix, ri, ci).re, &tmp, &std_val);
                    } else {
                        apf_zero(&MAT_AT(&result->v.matrix, ri, ci).re);
                    }
                    apf_zero(&MAT_AT(&result->v.matrix, ri, ci).im);
                }
                return 1;
            }
            
            /* minmax_scale(x, a, b) - Scale to [a, b] range */
            if (str_eq(name, "minmax_scale") || str_eq(name, "minmaxscale")) {
                apfc a_val, b_val;
                int n, i;
                apf min_val, max_val, range_in, range_out, a, b;
                int first = 1;
                
                next_token();
                if (!parse_value(&pv_arg)) return 0;
                
                /* Default to [0, 1] */
                apf_zero(&a);
                apf_from_int(&b, 1);
                
                if (current_token.type == TOK_COMMA) {
                    next_token();
                    if (!parse_expr(&a_val)) return 0;
                    apf_copy(&a, &a_val.re);
                    if (current_token.type == TOK_COMMA) {
                        next_token();
                        if (!parse_expr(&b_val)) return 0;
                        apf_copy(&b, &b_val.re);
                    }
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
                
                n = pv_arg.v.matrix.rows * pv_arg.v.matrix.cols;
                
                /* Find min and max */
                for (i = 0; i < n; i++) {
                    int ri = i % pv_arg.v.matrix.rows;
                    int ci = i / pv_arg.v.matrix.rows;
                    if (first) {
                        apf_copy(&min_val, &MAT_AT(&pv_arg.v.matrix, ri, ci).re);
                        apf_copy(&max_val, &MAT_AT(&pv_arg.v.matrix, ri, ci).re);
                        first = 0;
                    } else {
                        if (apf_cmp(&MAT_AT(&pv_arg.v.matrix, ri, ci).re, &min_val) < 0) {
                            apf_copy(&min_val, &MAT_AT(&pv_arg.v.matrix, ri, ci).re);
                        }
                        if (apf_cmp(&MAT_AT(&pv_arg.v.matrix, ri, ci).re, &max_val) > 0) {
                            apf_copy(&max_val, &MAT_AT(&pv_arg.v.matrix, ri, ci).re);
                        }
                    }
                }
                
                apf_sub(&range_in, &max_val, &min_val);
                apf_sub(&range_out, &b, &a);
                
                result->type = VAL_MATRIX;
                mat_zero(&result->v.matrix, pv_arg.v.matrix.rows, pv_arg.v.matrix.cols);
                
                for (i = 0; i < n; i++) {
                    int ri = i % pv_arg.v.matrix.rows;
                    int ci = i / pv_arg.v.matrix.rows;
                    apf normalized, scaled;
                    
                    if (!apf_is_zero(&range_in)) {
                        apf_sub(&normalized, &MAT_AT(&pv_arg.v.matrix, ri, ci).re, &min_val);
                        apf_div(&normalized, &normalized, &range_in);
                        apf_mul(&scaled, &normalized, &range_out);
                        apf_add(&MAT_AT(&result->v.matrix, ri, ci).re, &scaled, &a);
                    } else {
                        apf_copy(&MAT_AT(&result->v.matrix, ri, ci).re, &a);
                    }
                    apf_zero(&MAT_AT(&result->v.matrix, ri, ci).im);
                }
                return 1;
            }
            
            /* movstd(x, window) - Moving standard deviation */
            if (str_eq(name, "movstd")) {
                apfc window_val;
                int n, window, i, j;
                
                next_token();
                if (!parse_value(&pv_arg)) return 0;
                if (current_token.type != TOK_COMMA) {
                    printf("Error: movstd requires (x, window)\n");
                    return 0;
                }
                next_token();
                if (!parse_expr(&window_val)) return 0;
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
                window = apf_to_long(&window_val.re);
                if (window < 1) window = 1;
                if (window > n) window = n;
                
                result->type = VAL_MATRIX;
                mat_zero(&result->v.matrix, pv_arg.v.matrix.rows, pv_arg.v.matrix.cols);
                
                for (i = 0; i < n; i++) {
                    int ri = i % pv_arg.v.matrix.rows;
                    int ci = i / pv_arg.v.matrix.rows;
                    int half = window / 2;
                    int start = i - half;
                    int end = start + window;
                    int cnt = 0;
                    apf sum, sumsq, mean_val, var_val, count, tmp;
                    
                    if (start < 0) start = 0;
                    if (end > n) end = n;
                    
                    apf_zero(&sum);
                    for (j = start; j < end; j++) {
                        int rj = j % pv_arg.v.matrix.rows;
                        int cj = j / pv_arg.v.matrix.rows;
                        apf_add(&sum, &sum, &MAT_AT(&pv_arg.v.matrix, rj, cj).re);
                        cnt++;
                    }
                    apf_from_int(&count, cnt);
                    apf_div(&mean_val, &sum, &count);
                    
                    apf_zero(&sumsq);
                    for (j = start; j < end; j++) {
                        int rj = j % pv_arg.v.matrix.rows;
                        int cj = j / pv_arg.v.matrix.rows;
                        apf_sub(&tmp, &MAT_AT(&pv_arg.v.matrix, rj, cj).re, &mean_val);
                        apf_mul(&tmp, &tmp, &tmp);
                        apf_add(&sumsq, &sumsq, &tmp);
                    }
                    
                    if (cnt > 1) {
                        apf_from_int(&count, cnt - 1);
                        apf_div(&var_val, &sumsq, &count);
                        apf_sqrt(&MAT_AT(&result->v.matrix, ri, ci).re, &var_val);
                    } else {
                        apf_zero(&MAT_AT(&result->v.matrix, ri, ci).re);
                    }
                    apf_zero(&MAT_AT(&result->v.matrix, ri, ci).im);
                }
                return 1;
            }
            
            /* movmax(x, window) - Moving maximum */
            if (str_eq(name, "movmax")) {
                apfc window_val;
                int n, window, i, j;
                
                next_token();
                if (!parse_value(&pv_arg)) return 0;
                if (current_token.type != TOK_COMMA) {
                    printf("Error: movmax requires (x, window)\n");
                    return 0;
                }
                next_token();
                if (!parse_expr(&window_val)) return 0;
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
                window = apf_to_long(&window_val.re);
                if (window < 1) window = 1;
                if (window > n) window = n;
                
                result->type = VAL_MATRIX;
                mat_zero(&result->v.matrix, pv_arg.v.matrix.rows, pv_arg.v.matrix.cols);
                
                for (i = 0; i < n; i++) {
                    int ri = i % pv_arg.v.matrix.rows;
                    int ci = i / pv_arg.v.matrix.rows;
                    int half = window / 2;
                    int start = i - half;
                    int end = start + window;
                    int first = 1;
                    apf max_val;
                    
                    if (start < 0) start = 0;
                    if (end > n) end = n;
                    
                    for (j = start; j < end; j++) {
                        int rj = j % pv_arg.v.matrix.rows;
                        int cj = j / pv_arg.v.matrix.rows;
                        if (first || apf_cmp(&MAT_AT(&pv_arg.v.matrix, rj, cj).re, &max_val) > 0) {
                            apf_copy(&max_val, &MAT_AT(&pv_arg.v.matrix, rj, cj).re);
                            first = 0;
                        }
                    }
                    apf_copy(&MAT_AT(&result->v.matrix, ri, ci).re, &max_val);
                    apf_zero(&MAT_AT(&result->v.matrix, ri, ci).im);
                }
                return 1;
            }
            
            /* movmin(x, window) - Moving minimum */
            if (str_eq(name, "movmin")) {
                apfc window_val;
                int n, window, i, j;
                
                next_token();
                if (!parse_value(&pv_arg)) return 0;
                if (current_token.type != TOK_COMMA) {
                    printf("Error: movmin requires (x, window)\n");
                    return 0;
                }
                next_token();
                if (!parse_expr(&window_val)) return 0;
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
                window = apf_to_long(&window_val.re);
                if (window < 1) window = 1;
                if (window > n) window = n;
                
                result->type = VAL_MATRIX;
                mat_zero(&result->v.matrix, pv_arg.v.matrix.rows, pv_arg.v.matrix.cols);
                
                for (i = 0; i < n; i++) {
                    int ri = i % pv_arg.v.matrix.rows;
                    int ci = i / pv_arg.v.matrix.rows;
                    int half = window / 2;
                    int start = i - half;
                    int end = start + window;
                    int first = 1;
                    apf min_val;
                    
                    if (start < 0) start = 0;
                    if (end > n) end = n;
                    
                    for (j = start; j < end; j++) {
                        int rj = j % pv_arg.v.matrix.rows;
                        int cj = j / pv_arg.v.matrix.rows;
                        if (first || apf_cmp(&MAT_AT(&pv_arg.v.matrix, rj, cj).re, &min_val) < 0) {
                            apf_copy(&min_val, &MAT_AT(&pv_arg.v.matrix, rj, cj).re);
                            first = 0;
                        }
                    }
                    apf_copy(&MAT_AT(&result->v.matrix, ri, ci).re, &min_val);
                    apf_zero(&MAT_AT(&result->v.matrix, ri, ci).im);
                }
                return 1;
            }
            
            /* cumreturn(x) - Cumulative return from percentage returns */
            if (str_eq(name, "cumreturn") || str_eq(name, "cumret")) {
                int n, i;
                apf product, one, hundred, pct;
                
                next_token();
                if (!parse_value(&pv_arg)) return 0;
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
                apf_from_int(&one, 1);
                apf_from_int(&hundred, 100);
                apf_copy(&product, &one);
                
                result->type = VAL_MATRIX;
                mat_zero(&result->v.matrix, pv_arg.v.matrix.rows, pv_arg.v.matrix.cols);
                
                for (i = 0; i < n; i++) {
                    int ri = i % pv_arg.v.matrix.rows;
                    int ci = i / pv_arg.v.matrix.rows;
                    apf factor;
                    
                    /* Convert percentage to factor: 1 + pct/100 */
                    apf_div(&pct, &MAT_AT(&pv_arg.v.matrix, ri, ci).re, &hundred);
                    apf_add(&factor, &one, &pct);
                    apf_mul(&product, &product, &factor);
                    
                    /* Store cumulative return as percentage */
                    apf_sub(&factor, &product, &one);
                    apf_mul(&MAT_AT(&result->v.matrix, ri, ci).re, &factor, &hundred);
                    apf_zero(&MAT_AT(&result->v.matrix, ri, ci).im);
                }
                return 1;
            }
            
            /* drawdown(x) - Maximum drawdown from peak */
            if (str_eq(name, "drawdown")) {
                int n, i;
                apf peak, dd, hundred;
                
                next_token();
                if (!parse_value(&pv_arg)) return 0;
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
                apf_from_int(&hundred, 100);
                
                result->type = VAL_MATRIX;
                mat_zero(&result->v.matrix, pv_arg.v.matrix.rows, pv_arg.v.matrix.cols);
                
                apf_copy(&peak, &MAT_AT(&pv_arg.v.matrix, 0, 0).re);
                
                for (i = 0; i < n; i++) {
                    int ri = i % pv_arg.v.matrix.rows;
                    int ci = i / pv_arg.v.matrix.rows;
                    
                    /* Update peak */
                    if (apf_cmp(&MAT_AT(&pv_arg.v.matrix, ri, ci).re, &peak) > 0) {
                        apf_copy(&peak, &MAT_AT(&pv_arg.v.matrix, ri, ci).re);
                    }
                    
                    /* Calculate drawdown as percentage from peak */
                    if (!apf_is_zero(&peak)) {
                        apf_sub(&dd, &MAT_AT(&pv_arg.v.matrix, ri, ci).re, &peak);
                        apf_div(&dd, &dd, &peak);
                        apf_mul(&MAT_AT(&result->v.matrix, ri, ci).re, &dd, &hundred);
                    } else {
                        apf_zero(&MAT_AT(&result->v.matrix, ri, ci).re);
                    }
                    apf_zero(&MAT_AT(&result->v.matrix, ri, ci).im);
                }
                return 1;
            }
            
            /* sharpe(returns, riskfree) - Sharpe ratio */
            if (str_eq(name, "sharpe") || str_eq(name, "sharperatio")) {
                apfc rf_val;
                int n, i;
                apf mean_val, std_val, sum, sumsq, count, tmp, rf, excess;
                
                next_token();
                if (!parse_value(&pv_arg)) return 0;
                
                apf_zero(&rf);  /* Default risk-free rate = 0 */
                if (current_token.type == TOK_COMMA) {
                    next_token();
                    if (!parse_expr(&rf_val)) return 0;
                    apf_copy(&rf, &rf_val.re);
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
                
                n = pv_arg.v.matrix.rows * pv_arg.v.matrix.cols;
                
                /* Compute mean */
                apf_zero(&sum);
                for (i = 0; i < n; i++) {
                    int ri = i % pv_arg.v.matrix.rows;
                    int ci = i / pv_arg.v.matrix.rows;
                    apf_add(&sum, &sum, &MAT_AT(&pv_arg.v.matrix, ri, ci).re);
                }
                apf_from_int(&count, n);
                apf_div(&mean_val, &sum, &count);
                
                /* Compute std dev */
                apf_zero(&sumsq);
                for (i = 0; i < n; i++) {
                    int ri = i % pv_arg.v.matrix.rows;
                    int ci = i / pv_arg.v.matrix.rows;
                    apf_sub(&tmp, &MAT_AT(&pv_arg.v.matrix, ri, ci).re, &mean_val);
                    apf_mul(&tmp, &tmp, &tmp);
                    apf_add(&sumsq, &sumsq, &tmp);
                }
                apf_from_int(&count, n - 1);
                apf_div(&tmp, &sumsq, &count);
                apf_sqrt(&std_val, &tmp);
                
                /* Sharpe = (mean - rf) / std */
                apf_sub(&excess, &mean_val, &rf);
                if (!apf_is_zero(&std_val)) {
                    apf_div(&result->v.scalar.re, &excess, &std_val);
                } else {
                    apf_zero(&result->v.scalar.re);
                }
                apf_zero(&result->v.scalar.im);
                result->type = VAL_SCALAR;
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
            
            /* meshgrid(x, y) - create 2D grid matrices
             * Returns X matrix where each row is a copy of x
             */
            if (str_eq(name, "meshgrid")) {
                value_t arg2;
                int nx, ny, i, j;
                
                next_token();
                if (!parse_value(&pv_arg)) return 0;
                if (current_token.type != TOK_COMMA) {
                    printf("Error: meshgrid requires two arguments\n");
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
                
                nx = pv_arg.v.matrix.rows * pv_arg.v.matrix.cols;
                ny = arg2.v.matrix.rows * arg2.v.matrix.cols;
                
                result->type = VAL_MATRIX;
                mat_zero(&result->v.matrix, ny, nx);
                
                /* X: repeat x vector for each row */
                for (i = 0; i < ny; i++) {
                    for (j = 0; j < nx; j++) {
                        MAT_AT(&result->v.matrix, i, j) = pv_arg.v.matrix.data[j];
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
            
            /* binopdf(k, n, p) - Binomial PMF */
            if (str_eq(name, "binopdf")) {
                apfc k_arg, n_arg, p_arg;
                double k, n, p, res;
                
                next_token();
                if (!parse_expr(&k_arg)) return 0;
                if (current_token.type != TOK_COMMA) {
                    printf("Error: binopdf requires 3 arguments (k, n, p)\n");
                    return 0;
                }
                next_token();
                if (!parse_expr(&n_arg)) return 0;
                if (current_token.type != TOK_COMMA) {
                    printf("Error: binopdf requires 3 arguments (k, n, p)\n");
                    return 0;
                }
                next_token();
                if (!parse_expr(&p_arg)) return 0;
                if (current_token.type != TOK_RPAREN) {
                    printf("Error: expected ')'\n");
                    return 0;
                }
                next_token();
                
                k = apf_to_double(&k_arg.re);
                n = apf_to_double(&n_arg.re);
                p = apf_to_double(&p_arg.re);
                
                /* binopdf = C(n,k) * p^k * (1-p)^(n-k) */
                {
                    int ki = (int)k, ni = (int)n;
                    if (ki < 0 || ki > ni || p < 0 || p > 1) {
                        res = 0;
                    } else if (p == 0) {
                        res = (ki == 0) ? 1.0 : 0.0;
                    } else if (p == 1) {
                        res = (ki == ni) ? 1.0 : 0.0;
                    } else {
                        res = exp(lgamma(ni + 1) - lgamma(ki + 1) - lgamma(ni - ki + 1) +
                                  ki * log(p) + (ni - ki) * log(1 - p));
                    }
                }
                
                apf_from_double(&result->v.scalar.re, res);
                apf_zero(&result->v.scalar.im);
                result->type = VAL_SCALAR;
                return 1;
            }
            
            /* binocdf(k, n, p) - Binomial CDF */
            if (str_eq(name, "binocdf")) {
                apfc k_arg, n_arg, p_arg;
                double k, n, p, res;
                int i, ki, ni;
                
                next_token();
                if (!parse_expr(&k_arg)) return 0;
                if (current_token.type != TOK_COMMA) {
                    printf("Error: binocdf requires 3 arguments (k, n, p)\n");
                    return 0;
                }
                next_token();
                if (!parse_expr(&n_arg)) return 0;
                if (current_token.type != TOK_COMMA) {
                    printf("Error: binocdf requires 3 arguments (k, n, p)\n");
                    return 0;
                }
                next_token();
                if (!parse_expr(&p_arg)) return 0;
                if (current_token.type != TOK_RPAREN) {
                    printf("Error: expected ')'\n");
                    return 0;
                }
                next_token();
                
                k = apf_to_double(&k_arg.re);
                n = apf_to_double(&n_arg.re);
                p = apf_to_double(&p_arg.re);
                ki = (int)k;
                ni = (int)n;
                
                res = 0;
                for (i = 0; i <= ki && i <= ni; i++) {
                    double pmf = exp(lgamma(ni + 1) - lgamma(i + 1) - lgamma(ni - i + 1) +
                                     i * log(p) + (ni - i) * log(1 - p));
                    res += pmf;
                }
                if (res > 1) res = 1;
                
                apf_from_double(&result->v.scalar.re, res);
                apf_zero(&result->v.scalar.im);
                result->type = VAL_SCALAR;
                return 1;
            }
            
            /* unifpdf(x, a, b) - Uniform PDF on [a, b] */
            if (str_eq(name, "unifpdf")) {
                apfc x_arg, a_arg, b_arg;
                double x, a, b, res;
                
                next_token();
                if (!parse_expr(&x_arg)) return 0;
                a = 0; b = 1;  /* defaults */
                if (current_token.type == TOK_COMMA) {
                    next_token();
                    if (!parse_expr(&a_arg)) return 0;
                    a = apf_to_double(&a_arg.re);
                    if (current_token.type == TOK_COMMA) {
                        next_token();
                        if (!parse_expr(&b_arg)) return 0;
                        b = apf_to_double(&b_arg.re);
                    }
                }
                if (current_token.type != TOK_RPAREN) {
                    printf("Error: expected ')'\n");
                    return 0;
                }
                next_token();
                
                x = apf_to_double(&x_arg.re);
                if (b <= a) res = 0;
                else if (x < a || x > b) res = 0;
                else res = 1.0 / (b - a);
                
                apf_from_double(&result->v.scalar.re, res);
                apf_zero(&result->v.scalar.im);
                result->type = VAL_SCALAR;
                return 1;
            }
            
            /* unifcdf(x, a, b) - Uniform CDF on [a, b] */
            if (str_eq(name, "unifcdf")) {
                apfc x_arg, a_arg, b_arg;
                double x, a, b, res;
                
                next_token();
                if (!parse_expr(&x_arg)) return 0;
                a = 0; b = 1;
                if (current_token.type == TOK_COMMA) {
                    next_token();
                    if (!parse_expr(&a_arg)) return 0;
                    a = apf_to_double(&a_arg.re);
                    if (current_token.type == TOK_COMMA) {
                        next_token();
                        if (!parse_expr(&b_arg)) return 0;
                        b = apf_to_double(&b_arg.re);
                    }
                }
                if (current_token.type != TOK_RPAREN) {
                    printf("Error: expected ')'\n");
                    return 0;
                }
                next_token();
                
                x = apf_to_double(&x_arg.re);
                if (b <= a) res = 0;
                else if (x <= a) res = 0;
                else if (x >= b) res = 1;
                else res = (x - a) / (b - a);
                
                apf_from_double(&result->v.scalar.re, res);
                apf_zero(&result->v.scalar.im);
                result->type = VAL_SCALAR;
                return 1;
            }
            
            /* fpdf(x, d1, d2) - F-distribution PDF */
            if (str_eq(name, "fpdf")) {
                apfc x_arg, d1_arg, d2_arg;
                double x, d1, d2, res;
                
                next_token();
                if (!parse_expr(&x_arg)) return 0;
                if (current_token.type != TOK_COMMA) {
                    printf("Error: fpdf requires 3 arguments (x, d1, d2)\n");
                    return 0;
                }
                next_token();
                if (!parse_expr(&d1_arg)) return 0;
                if (current_token.type != TOK_COMMA) {
                    printf("Error: fpdf requires 3 arguments (x, d1, d2)\n");
                    return 0;
                }
                next_token();
                if (!parse_expr(&d2_arg)) return 0;
                if (current_token.type != TOK_RPAREN) {
                    printf("Error: expected ')'\n");
                    return 0;
                }
                next_token();
                
                x = apf_to_double(&x_arg.re);
                d1 = apf_to_double(&d1_arg.re);
                d2 = apf_to_double(&d2_arg.re);
                
                if (x < 0 || d1 <= 0 || d2 <= 0) {
                    res = 0;
                } else if (x == 0) {
                    res = (d1 < 2) ? INFINITY : ((d1 == 2) ? 1.0 : 0.0);
                } else {
                    double num = pow(d1 * x, d1 / 2) * pow(d2, d2 / 2);
                    double den = pow(d1 * x + d2, (d1 + d2) / 2);
                    res = (num / den) / (x * exp(lgamma(d1/2) + lgamma(d2/2) - lgamma((d1+d2)/2)));
                }
                
                apf_from_double(&result->v.scalar.re, res);
                apf_zero(&result->v.scalar.im);
                result->type = VAL_SCALAR;
                return 1;
            }
            
            /* fcdf(x, d1, d2) - F-distribution CDF using incomplete beta */
            if (str_eq(name, "fcdf")) {
                apfc x_arg, d1_arg, d2_arg;
                double x, d1, d2, res;
                extern double stat_betainc(double, double, double);
                
                next_token();
                if (!parse_expr(&x_arg)) return 0;
                if (current_token.type != TOK_COMMA) {
                    printf("Error: fcdf requires 3 arguments (x, d1, d2)\n");
                    return 0;
                }
                next_token();
                if (!parse_expr(&d1_arg)) return 0;
                if (current_token.type != TOK_COMMA) {
                    printf("Error: fcdf requires 3 arguments (x, d1, d2)\n");
                    return 0;
                }
                next_token();
                if (!parse_expr(&d2_arg)) return 0;
                if (current_token.type != TOK_RPAREN) {
                    printf("Error: expected ')'\n");
                    return 0;
                }
                next_token();
                
                x = apf_to_double(&x_arg.re);
                d1 = apf_to_double(&d1_arg.re);
                d2 = apf_to_double(&d2_arg.re);
                
                if (x <= 0) res = 0;
                else if (d1 <= 0 || d2 <= 0) res = 0;
                else res = stat_betainc(d1 * x / (d1 * x + d2), d1 / 2, d2 / 2);
                
                apf_from_double(&result->v.scalar.re, res);
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
            
            /* lineplot(y) or lineplot(x, y) - ASCII line plot */
            if (str_eq(name, "lineplot")) {
                value_t x_val, y_val;
                matrix_t *x_ptr = NULL;
                int i, px, py;
                double xmin, xmax, ymin, ymax;
                char plot[20][65];
                int has_x = 0;
                
                next_token();
                if (!parse_value(&y_val)) return 0;
                if (y_val.type != VAL_MATRIX) {
                    printf("Error: lineplot requires vector\n");
                    return 0;
                }
                
                /* Check for optional x argument */
                if (current_token.type == TOK_COMMA) {
                    x_val = y_val;  /* First arg is x */
                    x_ptr = &x_val.v.matrix;
                    has_x = 1;
                    next_token();
                    if (!parse_value(&y_val)) return 0;
                    if (y_val.type != VAL_MATRIX) {
                        printf("Error: lineplot y must be vector\n");
                        return 0;
                    }
                }
                
                if (current_token.type != TOK_RPAREN) {
                    printf("Error: expected ')'\n");
                    return 0;
                }
                next_token();
                
                /* Find bounds */
                ymin = ymax = apf_to_double(&MAT_AT(&y_val.v.matrix, 0, 0).re);
                if (has_x) {
                    xmin = xmax = apf_to_double(&MAT_AT(x_ptr, 0, 0).re);
                } else {
                    xmin = 1; xmax = y_val.v.matrix.rows;
                }
                
                for (i = 0; i < y_val.v.matrix.rows; i++) {
                    double y = apf_to_double(&MAT_AT(&y_val.v.matrix, i, 0).re);
                    if (y < ymin) ymin = y;
                    if (y > ymax) ymax = y;
                    if (has_x) {
                        double x = apf_to_double(&MAT_AT(x_ptr, i, 0).re);
                        if (x < xmin) xmin = x;
                        if (x > xmax) xmax = x;
                    }
                }
                
                /* Add margins */
                if (ymax == ymin) { ymax += 1; ymin -= 1; }
                { double m = (ymax - ymin) * 0.05; ymin -= m; ymax += m; }
                if (xmax == xmin) xmax = xmin + 1;
                
                /* Initialize plot */
                for (py = 0; py < 20; py++) {
                    for (px = 0; px < 64; px++) {
                        plot[py][px] = ' ';
                    }
                    plot[py][64] = '\0';
                }
                
                /* Draw axes */
                for (px = 0; px < 60; px++) plot[19][px] = '-';
                for (py = 0; py < 20; py++) plot[py][0] = '|';
                plot[19][0] = '+';
                
                /* Plot points and connect with lines */
                for (i = 0; i < y_val.v.matrix.rows; i++) {
                    double x = has_x ? apf_to_double(&MAT_AT(x_ptr, i, 0).re) : (i + 1);
                    double y = apf_to_double(&MAT_AT(&y_val.v.matrix, i, 0).re);
                    
                    px = 1 + (int)((x - xmin) / (xmax - xmin) * 58);
                    py = 18 - (int)((y - ymin) / (ymax - ymin) * 18);
                    
                    if (px >= 1 && px < 60 && py >= 0 && py < 19) {
                        plot[py][px] = '*';
                    }
                }
                
                /* Print plot */
                printf("\n");
                printf("  %8.2f |%s\n", ymax, plot[0] + 1);
                for (py = 1; py < 19; py++) {
                    if (py == 9) {
                        printf("  %8.2f |%s\n", (ymax + ymin) / 2, plot[py] + 1);
                    } else {
                        printf("           |%s\n", plot[py] + 1);
                    }
                }
                printf("  %8.2f +%s\n", ymin, plot[19] + 1);
                printf("           %-8.2f%50.2f\n\n", xmin, xmax);
                
                result->type = VAL_SCALAR;
                apf_zero(&result->v.scalar.re);
                apf_zero(&result->v.scalar.im);
                return 1;
            }
            
            /* scatter(x, y) or scatter(x, y, groups) or scatter(x, y, groups, "chars") */
            if (str_eq(name, "scatter")) {
                extern void mat_scatter(const matrix_t *x, const matrix_t *y, const matrix_t *groups, const char *chars);
                value_t x_val, y_val;
                matrix_t *groups_ptr = NULL;
                matrix_t groups_mat;
                char chars[64] = "";
                
                next_token();
                if (!parse_value(&x_val)) return 0;
                if (x_val.type != VAL_MATRIX) {
                    printf("Error: scatter requires vector x\n");
                    return 0;
                }
                
                if (current_token.type != TOK_COMMA) {
                    printf("Error: scatter requires (x, y) or (x, y, groups)\n");
                    return 0;
                }
                next_token();
                if (!parse_value(&y_val)) return 0;
                if (y_val.type != VAL_MATRIX) {
                    printf("Error: scatter requires vector y\n");
                    return 0;
                }
                
                /* Optional groups argument */
                if (current_token.type == TOK_COMMA) {
                    next_token();
                    if (!parse_value(&pv_arg)) return 0;
                    if (pv_arg.type == VAL_MATRIX) {
                        groups_mat = pv_arg.v.matrix;
                        groups_ptr = &groups_mat;
                    }
                }
                
                /* Optional chars argument */
                if (current_token.type == TOK_COMMA) {
                    next_token();
                    if (current_token.type == TOK_STRING) {
                        strncpy(chars, current_token.str_value, 63);
                        chars[63] = '\0';
                        next_token();
                    }
                }
                
                if (current_token.type != TOK_RPAREN) {
                    printf("Error: expected ')'\n");
                    return 0;
                }
                next_token();
                
                mat_scatter(&x_val.v.matrix, &y_val.v.matrix, groups_ptr, chars);
                
                /* Return empty result (scatter is a display function) */
                result->type = VAL_SCALAR;
                apf_zero(&result->v.scalar.re);
                apf_zero(&result->v.scalar.im);
                return 1;
            }
            
            /* toprevenue(data, n) - top N customers by revenue */
            if (str_eq(name, "toprevenue")) {
                extern void mat_toprevenue(matrix_t *result, const matrix_t *data, int top_n);
                int top_n = 10;  /* Default to top 10 */
                apfc n_val;
                
                next_token();
                if (!parse_value(&pv_arg)) return 0;
                if (pv_arg.type != VAL_MATRIX) {
                    printf("Error: toprevenue requires matrix input\n");
                    return 0;
                }
                
                /* Optional second argument */
                if (current_token.type == TOK_COMMA) {
                    next_token();
                    if (!parse_expr(&n_val)) return 0;
                    top_n = apf_to_long(&n_val.re);
                    if (top_n < 1) top_n = 1;
                }
                
                if (current_token.type != TOK_RPAREN) {
                    printf("Error: expected ')'\n");
                    return 0;
                }
                next_token();
                
                result->type = VAL_MATRIX;
                mat_toprevenue(&result->v.matrix, &pv_arg.v.matrix, top_n);
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
            
            /* polyfit(x, y, n) - polynomial curve fitting */
            if (str_eq(name, "polyfit")) {
                value_t y_arg;
                apfc n_val;
                int nx, ny, degree, i, j, k;
                apf *X_data, *y_data, *XtX, *Xty, *coeffs;
                
                next_token();
                if (!parse_value(&pv_arg)) return 0;
                if (current_token.type != TOK_COMMA) {
                    printf("Error: polyfit requires (x, y, n)\n");
                    return 0;
                }
                next_token();
                if (!parse_value(&y_arg)) return 0;
                if (current_token.type != TOK_COMMA) {
                    printf("Error: polyfit requires (x, y, n)\n");
                    return 0;
                }
                next_token();
                if (!parse_expr(&n_val)) return 0;
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
                if (y_arg.type == VAL_SCALAR) {
                    apfc saved = y_arg.v.scalar;
                    mat_zero(&y_arg.v.matrix, 1, 1);
                    MAT_AT(&y_arg.v.matrix, 0, 0) = saved;
                    y_arg.type = VAL_MATRIX;
                }
                
                nx = pv_arg.v.matrix.rows * pv_arg.v.matrix.cols;
                ny = y_arg.v.matrix.rows * y_arg.v.matrix.cols;
                degree = apf_to_long(&n_val.re);
                
                if (nx != ny || nx < degree + 1) {
                    printf("Error: polyfit needs matching x,y with enough points\n");
                    return 0;
                }
                if (degree < 0) degree = 0;
                if (degree > 10) degree = 10;
                
                /* Simple least squares: solve (X'X)c = X'y */
                /* X is Vandermonde matrix */
                X_data = (apf*)malloc(nx * (degree + 1) * sizeof(apf));
                y_data = (apf*)malloc(nx * sizeof(apf));
                XtX = (apf*)malloc((degree + 1) * (degree + 1) * sizeof(apf));
                Xty = (apf*)malloc((degree + 1) * sizeof(apf));
                coeffs = (apf*)malloc((degree + 1) * sizeof(apf));
                
                if (!X_data || !y_data || !XtX || !Xty || !coeffs) {
                    free(X_data); free(y_data); free(XtX); free(Xty); free(coeffs);
                    printf("Error: memory allocation failed\n");
                    return 0;
                }
                
                /* Build Vandermonde matrix and copy y */
                for (i = 0; i < nx; i++) {
                    int ri = i % pv_arg.v.matrix.rows;
                    int ci = i / pv_arg.v.matrix.rows;
                    apf xi, xpow;
                    apf_copy(&xi, &MAT_AT(&pv_arg.v.matrix, ri, ci).re);
                    apf_from_int(&xpow, 1);
                    for (j = degree; j >= 0; j--) {
                        apf_copy(&X_data[i * (degree + 1) + (degree - j)], &xpow);
                        apf_mul(&xpow, &xpow, &xi);
                    }
                    ri = i % y_arg.v.matrix.rows;
                    ci = i / y_arg.v.matrix.rows;
                    apf_copy(&y_data[i], &MAT_AT(&y_arg.v.matrix, ri, ci).re);
                }
                
                /* Compute X'X */
                for (i = 0; i <= degree; i++) {
                    for (j = 0; j <= degree; j++) {
                        apf_zero(&XtX[i * (degree + 1) + j]);
                        for (k = 0; k < nx; k++) {
                            apf prod;
                            apf_mul(&prod, &X_data[k * (degree + 1) + i], &X_data[k * (degree + 1) + j]);
                            apf_add(&XtX[i * (degree + 1) + j], &XtX[i * (degree + 1) + j], &prod);
                        }
                    }
                }
                
                /* Compute X'y */
                for (i = 0; i <= degree; i++) {
                    apf_zero(&Xty[i]);
                    for (k = 0; k < nx; k++) {
                        apf prod;
                        apf_mul(&prod, &X_data[k * (degree + 1) + i], &y_data[k]);
                        apf_add(&Xty[i], &Xty[i], &prod);
                    }
                }
                
                /* Solve using Gaussian elimination */
                for (i = 0; i <= degree; i++) {
                    apf_copy(&coeffs[i], &Xty[i]);
                }
                
                for (i = 0; i <= degree; i++) {
                    /* Find pivot */
                    int pivot = i;
                    for (j = i + 1; j <= degree; j++) {
                        apf ai, aj;
                        apf_abs(&ai, &XtX[pivot * (degree + 1) + i]);
                        apf_abs(&aj, &XtX[j * (degree + 1) + i]);
                        if (apf_cmp(&aj, &ai) > 0) pivot = j;
                    }
                    /* Swap rows */
                    if (pivot != i) {
                        for (k = 0; k <= degree; k++) {
                            apf tmp;
                            apf_copy(&tmp, &XtX[i * (degree + 1) + k]);
                            apf_copy(&XtX[i * (degree + 1) + k], &XtX[pivot * (degree + 1) + k]);
                            apf_copy(&XtX[pivot * (degree + 1) + k], &tmp);
                        }
                        {
                            apf tmp;
                            apf_copy(&tmp, &coeffs[i]);
                            apf_copy(&coeffs[i], &coeffs[pivot]);
                            apf_copy(&coeffs[pivot], &tmp);
                        }
                    }
                    /* Eliminate */
                    for (j = i + 1; j <= degree; j++) {
                        apf factor;
                        if (!apf_is_zero(&XtX[i * (degree + 1) + i])) {
                            apf_div(&factor, &XtX[j * (degree + 1) + i], &XtX[i * (degree + 1) + i]);
                            for (k = i; k <= degree; k++) {
                                apf prod;
                                apf_mul(&prod, &factor, &XtX[i * (degree + 1) + k]);
                                apf_sub(&XtX[j * (degree + 1) + k], &XtX[j * (degree + 1) + k], &prod);
                            }
                            {
                                apf prod;
                                apf_mul(&prod, &factor, &coeffs[i]);
                                apf_sub(&coeffs[j], &coeffs[j], &prod);
                            }
                        }
                    }
                }
                
                /* Back substitution */
                for (i = degree; i >= 0; i--) {
                    for (j = i + 1; j <= degree; j++) {
                        apf prod;
                        apf_mul(&prod, &XtX[i * (degree + 1) + j], &coeffs[j]);
                        apf_sub(&coeffs[i], &coeffs[i], &prod);
                    }
                    if (!apf_is_zero(&XtX[i * (degree + 1) + i])) {
                        apf_div(&coeffs[i], &coeffs[i], &XtX[i * (degree + 1) + i]);
                    }
                }
                
                result->type = VAL_MATRIX;
                mat_zero(&result->v.matrix, 1, degree + 1);
                for (i = 0; i <= degree; i++) {
                    apf_copy(&MAT_AT(&result->v.matrix, 0, i).re, &coeffs[i]);
                    apf_zero(&MAT_AT(&result->v.matrix, 0, i).im);
                }
                
                free(X_data); free(y_data); free(XtX); free(Xty); free(coeffs);
                return 1;
            }
            
            /* cumtrapz(y) or cumtrapz(x, y) - cumulative trapezoidal integration */
            if (str_eq(name, "cumtrapz")) {
                value_t y_arg;
                int n, i, has_x = 0;
                apf sum, half, dx, avg;
                
                next_token();
                if (!parse_value(&pv_arg)) return 0;
                
                if (current_token.type == TOK_COMMA) {
                    next_token();
                    if (!parse_value(&y_arg)) return 0;
                    has_x = 1;
                } else {
                    y_arg = pv_arg;
                    has_x = 0;
                }
                
                if (current_token.type != TOK_RPAREN) {
                    printf("Error: expected ')'\n");
                    return 0;
                }
                next_token();
                
                if (y_arg.type == VAL_SCALAR) {
                    apfc saved = y_arg.v.scalar;
                    mat_zero(&y_arg.v.matrix, 1, 1);
                    MAT_AT(&y_arg.v.matrix, 0, 0) = saved;
                    y_arg.type = VAL_MATRIX;
                }
                if (has_x && pv_arg.type == VAL_SCALAR) {
                    apfc saved = pv_arg.v.scalar;
                    mat_zero(&pv_arg.v.matrix, 1, 1);
                    MAT_AT(&pv_arg.v.matrix, 0, 0) = saved;
                    pv_arg.type = VAL_MATRIX;
                }
                
                n = y_arg.v.matrix.rows * y_arg.v.matrix.cols;
                result->type = VAL_MATRIX;
                mat_zero(&result->v.matrix, y_arg.v.matrix.rows, y_arg.v.matrix.cols);
                
                apf_zero(&sum);
                apf_from_int(&half, 1);
                {
                    apf two;
                    apf_from_int(&two, 2);
                    apf_div(&half, &half, &two);
                }
                
                /* First element is 0 */
                apf_zero(&MAT_AT(&result->v.matrix, 0, 0).re);
                apf_zero(&MAT_AT(&result->v.matrix, 0, 0).im);
                
                for (i = 1; i < n; i++) {
                    int ri = i % y_arg.v.matrix.rows;
                    int ci = i / y_arg.v.matrix.rows;
                    int ri_prev = (i-1) % y_arg.v.matrix.rows;
                    int ci_prev = (i-1) / y_arg.v.matrix.rows;
                    
                    if (has_x) {
                        int xi = i % pv_arg.v.matrix.rows;
                        int xci = i / pv_arg.v.matrix.rows;
                        int xi_prev = (i-1) % pv_arg.v.matrix.rows;
                        int xci_prev = (i-1) / pv_arg.v.matrix.rows;
                        apf_sub(&dx, &MAT_AT(&pv_arg.v.matrix, xi, xci).re,
                               &MAT_AT(&pv_arg.v.matrix, xi_prev, xci_prev).re);
                    } else {
                        apf_from_int(&dx, 1);
                    }
                    
                    apf_add(&avg, &MAT_AT(&y_arg.v.matrix, ri, ci).re,
                           &MAT_AT(&y_arg.v.matrix, ri_prev, ci_prev).re);
                    apf_mul(&avg, &avg, &half);
                    apf_mul(&avg, &avg, &dx);
                    apf_add(&sum, &sum, &avg);
                    
                    apf_copy(&MAT_AT(&result->v.matrix, ri, ci).re, &sum);
                    apf_zero(&MAT_AT(&result->v.matrix, ri, ci).im);
                }
                return 1;
            }
            
            /* diff2(v) - second difference */
            if (str_eq(name, "diff2")) {
                int n, i;
                
                next_token();
                if (!parse_value(&pv_arg)) return 0;
                if (current_token.type != TOK_RPAREN) {
                    printf("Error: expected ')'\n");
                    return 0;
                }
                next_token();
                
                if (pv_arg.type == VAL_SCALAR) {
                    result->type = VAL_MATRIX;
                    mat_zero(&result->v.matrix, 1, 1);
                    return 1;
                }
                
                n = pv_arg.v.matrix.rows * pv_arg.v.matrix.cols;
                if (n < 3) {
                    result->type = VAL_MATRIX;
                    mat_zero(&result->v.matrix, 1, 1);
                    return 1;
                }
                
                result->type = VAL_MATRIX;
                if (pv_arg.v.matrix.cols == 1) {
                    mat_zero(&result->v.matrix, n - 2, 1);
                } else {
                    mat_zero(&result->v.matrix, 1, n - 2);
                }
                
                for (i = 0; i < n - 2; i++) {
                    int ri0 = i % pv_arg.v.matrix.rows;
                    int ci0 = i / pv_arg.v.matrix.rows;
                    int ri1 = (i+1) % pv_arg.v.matrix.rows;
                    int ci1 = (i+1) / pv_arg.v.matrix.rows;
                    int ri2 = (i+2) % pv_arg.v.matrix.rows;
                    int ci2 = (i+2) / pv_arg.v.matrix.rows;
                    int ro = i % result->v.matrix.rows;
                    int co = i / result->v.matrix.rows;
                    apf d1, d2;
                    
                    apf_sub(&d1, &MAT_AT(&pv_arg.v.matrix, ri1, ci1).re,
                           &MAT_AT(&pv_arg.v.matrix, ri0, ci0).re);
                    apf_sub(&d2, &MAT_AT(&pv_arg.v.matrix, ri2, ci2).re,
                           &MAT_AT(&pv_arg.v.matrix, ri1, ci1).re);
                    apf_sub(&MAT_AT(&result->v.matrix, ro, co).re, &d2, &d1);
                    apf_zero(&MAT_AT(&result->v.matrix, ro, co).im);
                }
                return 1;
            }
            
            /* zscore_vec(v) - z-score standardization (vectorized) */
            if (str_eq(name, "zscore_vec") || str_eq(name, "standardize")) {
                int n, i;
                apf mean_val, std_val, sum, sum_sq, count, tmp, var;
                
                next_token();
                if (!parse_value(&pv_arg)) return 0;
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
                
                /* Compute mean */
                apf_zero(&sum);
                for (i = 0; i < n; i++) {
                    int ri = i % pv_arg.v.matrix.rows;
                    int ci = i / pv_arg.v.matrix.rows;
                    apf_add(&sum, &sum, &MAT_AT(&pv_arg.v.matrix, ri, ci).re);
                }
                apf_from_int(&count, n);
                apf_div(&mean_val, &sum, &count);
                
                /* Compute std */
                apf_zero(&sum_sq);
                for (i = 0; i < n; i++) {
                    int ri = i % pv_arg.v.matrix.rows;
                    int ci = i / pv_arg.v.matrix.rows;
                    apf_sub(&tmp, &MAT_AT(&pv_arg.v.matrix, ri, ci).re, &mean_val);
                    apf_mul(&tmp, &tmp, &tmp);
                    apf_add(&sum_sq, &sum_sq, &tmp);
                }
                apf_div(&var, &sum_sq, &count);
                apf_sqrt(&std_val, &var);
                
                result->type = VAL_MATRIX;
                mat_zero(&result->v.matrix, pv_arg.v.matrix.rows, pv_arg.v.matrix.cols);
                
                for (i = 0; i < n; i++) {
                    int ri = i % pv_arg.v.matrix.rows;
                    int ci = i / pv_arg.v.matrix.rows;
                    apf_sub(&tmp, &MAT_AT(&pv_arg.v.matrix, ri, ci).re, &mean_val);
                    if (!apf_is_zero(&std_val)) {
                        apf_div(&MAT_AT(&result->v.matrix, ri, ci).re, &tmp, &std_val);
                    } else {
                        apf_zero(&MAT_AT(&result->v.matrix, ri, ci).re);
                    }
                    apf_zero(&MAT_AT(&result->v.matrix, ri, ci).im);
                }
                return 1;
            }
            
            /* center(v) - subtract mean (center data) */
            if (str_eq(name, "center")) {
                int n, i;
                apf mean_val, sum, count;
                
                next_token();
                if (!parse_value(&pv_arg)) return 0;
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
                
                /* Compute mean */
                apf_zero(&sum);
                for (i = 0; i < n; i++) {
                    int ri = i % pv_arg.v.matrix.rows;
                    int ci = i / pv_arg.v.matrix.rows;
                    apf_add(&sum, &sum, &MAT_AT(&pv_arg.v.matrix, ri, ci).re);
                }
                apf_from_int(&count, n);
                apf_div(&mean_val, &sum, &count);
                
                result->type = VAL_MATRIX;
                mat_zero(&result->v.matrix, pv_arg.v.matrix.rows, pv_arg.v.matrix.cols);
                
                for (i = 0; i < n; i++) {
                    int ri = i % pv_arg.v.matrix.rows;
                    int ci = i / pv_arg.v.matrix.rows;
                    apf_sub(&MAT_AT(&result->v.matrix, ri, ci).re, 
                           &MAT_AT(&pv_arg.v.matrix, ri, ci).re, &mean_val);
                    apf_zero(&MAT_AT(&result->v.matrix, ri, ci).im);
                }
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
                mat_zero(&result->v.matrix, 1, n);
                if (!result->v.matrix.data) return 0;
                
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
            matrix_t tmp = {0, 0, NULL};
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
        mat_zero(&result->v.matrix, 1, n);
        if (!result->v.matrix.data) return 0;
        
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
