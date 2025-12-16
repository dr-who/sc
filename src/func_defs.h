/*
 * func_defs.h - Function documentation and dispatch system
 * 
 * All built-in functions are documented with doxygen-style comments
 * that are used to generate help, demos, tests, and man pages.
 *
 * C89 compliant
 */
#ifndef FUNC_DEFS_H
#define FUNC_DEFS_H

#include <stddef.h>  /* For NULL */

/* Maximum examples per function */
#define MAX_EXAMPLES 6

/* Function documentation entry */
typedef struct {
    const char *name;           /* Function name */
    const char *category;       /* Category for grouping */
    const char *syntax;         /* Usage syntax */
    const char *description;    /* Full description */
    const char *examples[MAX_EXAMPLES]; /* Example expressions (NULL-terminated) */
    const char *see_also;       /* Related functions */
} FuncDoc;

/* Category definitions */
#define CAT_TRIG        "Trigonometry"
#define CAT_EXP         "Exponential"
#define CAT_LOG         "Logarithmic"
#define CAT_COMPLEX     "Complex"
#define CAT_ROUND       "Rounding"
#define CAT_NUMBER      "Number Theory"
#define CAT_SPECIAL     "Special Functions"
#define CAT_STATS       "Statistics"
#define CAT_PROB        "Probability"
#define CAT_LINALG      "Linear Algebra"
#define CAT_MATRIX      "Matrix"
#define CAT_MATQUERY    "Matrix Query"
#define CAT_SIGNAL      "Signal Processing"
#define CAT_POLY        "Polynomials"
#define CAT_ANGLE       "Angle Conversion"
#define CAT_SET         "Set Operations"
#define CAT_LOGIC       "Logical"
#define CAT_BITWISE     "Bitwise"
#define CAT_COORD       "Coordinates"
#define CAT_UTIL        "Utility"
#define CAT_CONST       "Constants"
#define CAT_DATASCI     "Data Science"
#define CAT_TEXT        "Text Analysis"

/* ========== Function Tables by Category ========== */

extern const FuncDoc func_trig_docs[];
extern const FuncDoc func_exp_docs[];
extern const FuncDoc func_log_docs[];
extern const FuncDoc func_complex_docs[];
extern const FuncDoc func_round_docs[];
extern const FuncDoc func_number_docs[];
extern const FuncDoc func_special_docs[];
extern const FuncDoc func_stats_docs[];
extern const FuncDoc func_prob_docs[];
extern const FuncDoc func_linalg_docs[];
extern const FuncDoc func_matrix_docs[];
extern const FuncDoc func_matquery_docs[];
extern const FuncDoc func_signal_docs[];
extern const FuncDoc func_poly_docs[];
extern const FuncDoc func_angle_docs[];
extern const FuncDoc func_set_docs[];
extern const FuncDoc func_logic_docs[];
extern const FuncDoc func_bitwise_docs[];
extern const FuncDoc func_coord_docs[];
extern const FuncDoc func_util_docs[];
extern const FuncDoc func_const_docs[];
extern const FuncDoc func_datasci_docs[];
extern const FuncDoc func_text_docs[];

/* ========== API Functions ========== */

/* Find function documentation by name */
const FuncDoc *func_find_doc(const char *name);

/* Check if function exists */
int func_doc_exists(const char *name);

/* Count total documented functions */
int func_doc_count(void);

/* Show help for a function */
void func_show_help(const char *name);

/* Run demo for a function (execute examples) */
void func_show_demo(const char *name);

/* List all functions by category */
void func_list_all(void);

/* Run all function tests from examples */
int func_run_tests(void);

/* Generate test script to stdout */
void func_gen_test_script(void);

/* Generate man page section to stdout */
void func_gen_manpage(void);

#endif /* FUNC_DEFS_H */
