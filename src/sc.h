/* scalc.h - Stu's Calculator
 * Main header with common definitions
 * C89 compliant for Watcom C / DOS
 */
#ifndef SCALC_H
#define SCALC_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <time.h>
#include "apf.h"
#include "apfx.h"
#include "apfc.h"
#include "matrix.h"
#include "mathx.h"
#include "orbital.h"

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
    TOK_ASSIGN, TOK_SEMI, TOK_ERROR,
    /* Matrix tokens */
    TOK_LBRACKET, TOK_RBRACKET, TOK_COLON, TOK_COMMA,
    /* Element-wise operators */
    TOK_DOT_MUL, TOK_DOT_DIV, TOK_DOT_POW,
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
    char func_name[16];
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

/* ========== Mode ========== */
typedef enum { MODE_DECIMAL, MODE_FRACTION } calc_mode_t;

/* Maximum display digits (one less than internal for guard digit) */
#define MAX_DISPLAY_DIGITS 37
#define INTERNAL_DIGITS 38
#define DEFAULT_DIGITS 16

/* ========== Global State ========== */
extern calc_mode_t current_mode;
extern int display_digits;  /* 0 = show all, else round to N significant figures */
extern value_t last_ans;    /* Previous answer (ans variable) */
extern int last_ans_valid;  /* Is last_ans set? */
extern user_func_t user_funcs[MAX_FUNCTIONS];
extern scalar_var_t scalar_vars[MAX_SCALAR_VARS];
extern matrix_var_t matrix_vars[MAX_MATRIX_VARS];
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

/* ========== Runtime (runtime.c) ========== */
void init_user_funcs(void);
void init_variables(void);
int get_func_index(const char *name);
int get_var_index(const char *name);
int get_var_value(int idx, value_t *val);
int set_var_value(int idx, const value_t *val);
int is_var_defined(int idx);
int is_var_matrix(int idx);
int eval_user_func(apfc *result, int func_idx, const apfc *arg);

/* ========== Format (format.c) ========== */
void print_result(const apfc *val);
void print_value(const value_t *val);
void value_to_scalar(apfc *r, const value_t *val);
void value_from_scalar(value_t *r, const apfc *val);
void value_from_matrix(value_t *r, const matrix_t *m);
void apf_to_hex_str(const apf *val);
void apf_to_bin_str(const apf *val);

/* ========== Commands (commands.c) ========== */
void print_help(void);
void print_date_time(void);
void handle_mode(const char *input);
void run_tests(void);
void run_bench(void);

/* ========== REPL (repl.c) ========== */
int read_line(char *buf, int max_len);
void dos_init_screen(void);
void print_prompt(void);
int do_plot(const char *func_name, apf *xmin, apf *xmax);
int do_textplot(const char *func_name, apf *xmin, apf *xmax);
int do_lorenz_text(double sigma, double rho, double beta, int steps);
int do_rossler_text(double a, double b, double c, int steps);
int do_parametric_text(const char *xfunc, const char *yfunc,
                       double tmin, double tmax, int steps);

#endif /* SCALC_H */
