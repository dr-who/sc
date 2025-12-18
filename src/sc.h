/* scalc.h - Stu's Calculator
 * Main header with common definitions
 * C89 compliant for Watcom C / DOS
 */
#ifndef SCALC_H
#define SCALC_H

#include "config.h"  /* Must be first - sets feature flags */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <time.h>
#include "apf.h"
#include "apfx.h"
#include "apfc.h"
#include "matrix.h"   /* Types needed by parser.c even if functions not linked */
#include "mathx.h"

/* Optional modules - only include if feature is enabled and not DOS */
#if !defined(SCALC_MEDIUM) && !defined(SCALC_TINY) && !defined(SCALC_MINIMAL) && !defined(SCALC_VIC20)
#include "apf_native.h"  /* Native float conversion (uses double) */
#include "orbital.h"
#endif

/* Platform detection */
#ifdef __WATCOMC__
#include <conio.h>
#include <dos.h>
#define HAVE_CONIO 1
#else
#ifdef _WIN32
#include <conio.h>
#define HAVE_CONIO 1
#else
/* Unix/Linux */
#include <termios.h>
#include <unistd.h>
#define HAVE_TERMIOS 1
#endif
#endif

/* Configuration */
#ifndef MAX_INPUT
#define MAX_INPUT 256
#endif
#ifndef MAX_HISTORY
#define MAX_HISTORY 50
#endif
#ifndef MAX_FUNCTIONS
#define MAX_FUNCTIONS 26
#endif
#ifndef MAX_FUNC_BODY
#define MAX_FUNC_BODY 128
#endif
/* Matrix and user function limits */

/* ========== Token Types ========== */
typedef enum {
    TOK_END, TOK_NUM, TOK_IMAG, TOK_FUNC, TOK_VAR, TOK_ANS,
    TOK_PLUS, TOK_MINUS, TOK_MUL, TOK_DIV,
    TOK_POW, TOK_FACT, TOK_LPAREN, TOK_RPAREN,
    TOK_ASSIGN, TOK_EQUAL, TOK_NE, TOK_LT, TOK_LE, TOK_GT, TOK_GE, TOK_APPROX,
    TOK_AND, TOK_OR, TOK_NOT, TOK_XOR,
    TOK_SEMI, TOK_ERROR,
    TOK_STRING,  /* String literal in single quotes */
    TOK_IDENT,   /* Identifier (for field access) */
    TOK_DOT,     /* . for field access */
    /* Matrix tokens */
    TOK_LBRACKET, TOK_RBRACKET, TOK_COLON, TOK_COMMA,
    /* Element-wise operators */
    TOK_DOT_MUL, TOK_DOT_DIV, TOK_DOT_POW, TOK_DOT_BACKSLASH,
    /* Transpose */
    TOK_TRANSPOSE, TOK_DOT_TRANSPOSE,
    /* Matrix division */
    TOK_BACKSLASH,
    /* Keywords */
    TOK_FOR, TOK_IN, TOK_ENDFOR
} token_type_t;

typedef struct {
    token_type_t type;
    apf value;
    char func_name[32];
    char str_value[64];  /* For string literals */
} token_t;

/* ========== Value Type (scalar or matrix) ========== */
typedef enum { VAL_SCALAR, VAL_MATRIX } value_type_t;

typedef struct {
    value_type_t type;
    union {
        apfc scalar;
        matrix_t matrix;
    } v;
} value_t;

/* ========== User Functions ========== */
typedef struct {
    int defined;
    char param;
    char body[MAX_FUNC_BODY];
} user_func_t;

/* ========== Variables ========== */
/* Separate storage for scalars and matrices to save memory:
 * - 26 scalar variables (a-z): ~2KB total
 * - 4 matrix variables: ~6KB total
 * Total: ~8KB instead of ~39KB if we stored value_t for each
 */
#ifndef MAX_SCALAR_VARS
#define MAX_SCALAR_VARS 26
#endif
#ifndef MAX_MATRIX_VARS
#define MAX_MATRIX_VARS 4
#endif

typedef struct {
    int defined;
    apfc val;
} scalar_var_t;

typedef struct {
    int defined;
    char name;      /* which letter a-z */
    matrix_t val;
} matrix_var_t;

/* Named variables for MATLAB compatibility */
#ifdef HAVE_NAMED_VARS
#ifndef MAX_NAMED_VARS
#define MAX_NAMED_VARS 32
#endif
#ifndef MAX_VAR_NAME
#define MAX_VAR_NAME 32
#endif

typedef struct {
    int defined;
    char name[MAX_VAR_NAME];
    value_type_t type;
    union {
        apfc scalar;
        matrix_t matrix;
    } val;
} named_var_t;

extern named_var_t named_vars[MAX_NAMED_VARS];

/* Named variable functions */
int find_named_var(const char *name);
int create_named_var(const char *name);
int set_named_scalar(const char *name, const apfc *val);
int set_named_matrix(const char *name, const matrix_t *val);
int get_named_var(const char *name, value_t *result);
#endif /* HAVE_NAMED_VARS */

/* ========== Mode ========== */
typedef enum { 
    MODE_DECIMAL,   /* Show all significant digits */
    MODE_FRACTION,  /* Show as fraction if possible */
    MODE_HEX,       /* Hexadecimal output */
    MODE_BIN,       /* Binary output */
    MODE_IEEE,      /* IEEE 754 bit representation */
    MODE_FIX,       /* Fixed decimal places */
    MODE_SCI,       /* Scientific notation */
    MODE_ENG        /* Engineering notation (exp % 3 == 0) */
} calc_mode_t;

/* Maximum display digits - calculated from precision
 * decimal_digits ≈ AP_BITS * log10(2) ≈ AP_BITS * 0.301
 * Use integer math: (AP_BITS * 301) / 1000
 */
#define INTERNAL_DIGITS ((AP_BITS * 301) / 1000)
#define MAX_DISPLAY_DIGITS (INTERNAL_DIGITS - 1)  /* One less for guard digit */
#define DEFAULT_DIGITS 16

/* ========== Global State ========== */
extern calc_mode_t current_mode;
extern int display_digits;       /* 0 = show all, else round to N significant figures */
extern int display_fixed_places; /* For FIX/SCI/ENG: number of decimal places */
extern value_t last_ans;    /* Previous answer (ans variable) */
extern int last_ans_valid;  /* Is last_ans set? */
extern int result_is_boolean; /* Result should be displayed as true/false */
extern int angle_mode;      /* 0=radians, 1=degrees, 2=gradians */
#define ANGLE_RAD  0
#define ANGLE_DEG  1
#define ANGLE_GRAD 2
extern user_func_t user_funcs[MAX_FUNCTIONS];
extern scalar_var_t scalar_vars[MAX_SCALAR_VARS];
extern matrix_var_t SC_FAR matrix_vars[MAX_MATRIX_VARS];
extern const char *input_ptr;
extern token_t current_token;

/* ========== Utility Functions ========== */
int str_eq(const char *a, const char *b);
int str_starts(const char *str, const char *prefix);

/* ========== Lexer (lexer.c) ========== */
void skip_whitespace(void);
void next_token(void);

/* ========== Parser (parser.c) ========== */
int parse_expr(apfc *result);
int parse_value(value_t *result);
int parse_matrix(matrix_t *result);
void parser_init_temps(void);
int ensure_result_matrix(matrix_t *m, int rows, int cols);

/* ========== Runtime (runtime.c) ========== */
void init_user_funcs(void);
void init_variables(void);
#ifdef HAVE_NAMED_VARS
void init_named_vars(void);
#endif
int get_func_index(const char *name);
int get_var_index(const char *name);
int get_var_value(int idx, value_t *val);
int set_var_value(int idx, const value_t *val);
int is_var_defined(int idx);
int is_var_matrix(int idx);
int eval_user_func(apfc *result, int func_idx, const apfc *arg);

/* Evaluate any expression with x substituted for a given value */
int eval_expr_with_x(apfc *result, const char *expr, const apfc *x_val);
int eval_expr_with_var(apfc *result, const char *expr, char var_name, const apfc *val);

/* ========== Format (format.c) ========== */
void print_result(const apfc *val);
void print_value(const value_t *val);
void print_value_ex(const value_t *val, int quiet);
void value_to_scalar(apfc *r, const value_t *val);
void value_from_scalar(value_t *r, const apfc *val);
void value_from_matrix(value_t *r, const matrix_t *m);
void apf_to_hex_str(const apf *val);
void apf_to_bin_str(const apf *val);
void apf_ieee_display(const apf *val);
void apfc_ieee_display(const apfc *val);

/* ========== Commands (commands.c) ========== */
void print_help(void);
void print_constants(void);
void print_features(void);
void print_commands(void);
void run_demo(void);
void print_date_time(void);
void handle_mode(const char *input);
void run_tests(void);
void run_bench(void);

/* ========== REPL (repl.c) ========== */
int read_line(char *buf, int max_len);
void dos_init_screen(void);
void print_prompt(void);
int do_plot(const char *expr, char var, apf *xmin, apf *xmax);
int do_textplot(const char *expr, char var, apf *xmin, apf *xmax);
int do_lorenz_text(long sigma, long rho, long beta_num, long beta_den, int steps);
int do_rossler_text(long a_num, long a_den, long b_num, long b_den, 
                    long c_num, long c_den, int steps);
int do_parametric_text(const char *xfunc, const char *yfunc,
                       apf *tmin, apf *tmax, int steps);

/* Evaluate a single expression line (for demos) */
int eval_expr_line(const char *line, int quiet);

#endif /* SCALC_H */
