/*
 * func_registry.h - Function registration and documentation system
 * 
 * Each function is documented with doxygen-style comments that are
 * extracted by tools/gen_docs.sh to generate:
 *   - help.c entries for runtime help
 *   - man page sections
 *   - test cases
 *
 * C89 compliant
 */
#ifndef FUNC_REGISTRY_H
#define FUNC_REGISTRY_H

#include "sc.h"

/*
 * Function documentation comment format (parsed by tools/gen_docs.sh):
 *
 * /**
 *  * @func function_name
 *  * @category Category Name
 *  * @syntax syntax_string
 *  * @desc Description of the function
 *  * @example expression -> expected_result
 *  * @example another_expression -> another_result
 *  * /
 *
 * The @example lines are used for:
 *   - Interactive demos (demo funcname)
 *   - Test cases (make test)
 *   - Man page examples
 */

/* Function types */
typedef enum {
    FUNC_SCALAR,      /* f(x) -> scalar */
    FUNC_SCALAR2,     /* f(x,y) -> scalar */
    FUNC_SCALAR3,     /* f(x,y,z) -> scalar */
    FUNC_MATRIX,      /* f(A) -> matrix */
    FUNC_MATRIX2,     /* f(A,B) -> matrix */
    FUNC_MULTI_RET,   /* [a,b] = f(x) multi-value return */
    FUNC_SPECIAL      /* custom parsing */
} func_type_t;

/* Function entry in registry */
typedef struct {
    const char *name;
    const char *category;
    const char *syntax;
    const char *description;
    const char *examples[4];  /* NULL-terminated */
    func_type_t type;
    void *impl;               /* function pointer */
} func_entry_t;

/* Category registry - each func_*.c exports one of these */
typedef struct {
    const char *name;
    const func_entry_t *functions;
    int count;
} func_category_t;

/* Global function list (built from all categories) */
extern func_category_t *func_categories[];
extern int func_category_count;

/* Lookup functions */
const func_entry_t *func_find(const char *name);
const func_entry_t *func_find_in_category(const char *category, const char *name);

/* List functions */
void func_list_all(void);
void func_list_category(const char *category);
int func_count_all(void);

/* Help/demo functions */
void func_show_help(const char *name);
void func_show_demo(const char *name);

/* Test/bench all registered functions */
int func_run_tests(void);
void func_run_bench(void);

/* Registration (called at startup) */
void func_init_all(void);

#endif /* FUNC_REGISTRY_H */
