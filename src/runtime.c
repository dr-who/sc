/* runtime.c - Variables, user functions, and built-in functions
 * C89 compliant for Watcom C / DOS
 */
#include "sc.h"

/* Global state */
calc_mode_t current_mode = MODE_DECIMAL;
int display_digits = DEFAULT_DIGITS;  /* Default 15 (format long) */
int display_fixed_places = 6;         /* Default 6 decimal places for FIX/SCI/ENG */
value_t last_ans;      /* Previous answer */
int last_ans_valid = 0;  /* Is last_ans set? */
int result_is_boolean = 0; /* Result should be displayed as true/false */
int angle_mode = ANGLE_RAD; /* 0=radians (default), 1=degrees, 2=gradians */
user_func_t user_funcs[MAX_FUNCTIONS];
scalar_var_t scalar_vars[MAX_SCALAR_VARS];
matrix_var_t SC_FAR matrix_vars[MAX_MATRIX_VARS];

/* ========== Utility Functions ========== */

int str_eq(const char *a, const char *b)
{
    while (*a && *b) {
        char ca = *a, cb = *b;
        if (ca >= 'A' && ca <= 'Z') ca += 32;
        if (cb >= 'A' && cb <= 'Z') cb += 32;
        if (ca != cb) return 0;
        a++; b++;
    }
    return *a == *b;
}

int str_starts(const char *str, const char *prefix)
{
    while (*prefix) {
        char cs = *str, cp = *prefix;
        if (cs >= 'A' && cs <= 'Z') cs += 32;
        if (cp >= 'A' && cp <= 'Z') cp += 32;
        if (cs != cp) return 0;
        str++; prefix++;
    }
    return 1;
}

/* ========== User Functions ========== */

void init_user_funcs(void)
{
    int i;
    for (i = 0; i < MAX_FUNCTIONS; i++) {
        user_funcs[i].defined = 0;
    }
}

int get_func_index(const char *name)
{
    char c = name[0];
    /* Convert uppercase to lowercase */
    if (c >= 'A' && c <= 'Z') c = c - 'A' + 'a';
    if (c >= 'a' && c <= 'z' && name[1] == '\0') {
        return c - 'a';
    }
    return -1;
}

/* ========== Variables ========== */

void init_variables(void)
{
    int i;
    for (i = 0; i < MAX_SCALAR_VARS; i++) {
        scalar_vars[i].defined = 0;
    }
    for (i = 0; i < MAX_MATRIX_VARS; i++) {
        matrix_vars[i].defined = 0;
    }
}

int get_var_index(const char *name)
{
    char c = name[0];
    /* Convert uppercase to lowercase */
    if (c >= 'A' && c <= 'Z') c = c - 'A' + 'a';
    if (c >= 'a' && c <= 'z' && name[1] == '\0') {
        return c - 'a';  /* 0-25 for a-z */
    }
    return -1;
}

/* Get variable value - check scalar first, then matrix */
int get_var_value(int idx, value_t *val)
{
    int i;
    if (idx < 0 || idx >= MAX_SCALAR_VARS) return 0;
    
    /* Check scalar first */
    if (scalar_vars[idx].defined) {
        val->type = VAL_SCALAR;
        val->v.scalar = scalar_vars[idx].val;
        return 1;
    }
    
    /* Check matrix slots */
    for (i = 0; i < MAX_MATRIX_VARS; i++) {
        if (matrix_vars[i].defined && matrix_vars[i].name == 'a' + idx) {
            val->type = VAL_MATRIX;
            mat_copy(&val->v.matrix, &matrix_vars[i].val);
            return 1;
        }
    }
    
    return 0;  /* Not defined */
}

/* Set variable value */
int set_var_value(int idx, const value_t *val)
{
    int i;
    char name;
    
    if (idx < 0 || idx >= MAX_SCALAR_VARS) return 0;
    name = 'a' + idx;
    
    if (val->type == VAL_SCALAR) {
        /* Clear any matrix with this name */
        for (i = 0; i < MAX_MATRIX_VARS; i++) {
            if (matrix_vars[i].defined && matrix_vars[i].name == name) {
                matrix_vars[i].defined = 0;
                break;
            }
        }
        /* Set scalar */
        scalar_vars[idx].defined = 1;
        scalar_vars[idx].val = val->v.scalar;
        return 1;
    } else {
        /* Matrix - clear scalar first */
        scalar_vars[idx].defined = 0;
        
        /* Find existing slot or empty slot */
        for (i = 0; i < MAX_MATRIX_VARS; i++) {
            if (matrix_vars[i].defined && matrix_vars[i].name == name) {
                /* Reuse existing slot - use persistent storage */
                mat_copy_persist(&matrix_vars[i].val, &val->v.matrix);
                return 1;
            }
        }
        for (i = 0; i < MAX_MATRIX_VARS; i++) {
            if (!matrix_vars[i].defined) {
                /* Use empty slot - use persistent storage */
                matrix_vars[i].defined = 1;
                matrix_vars[i].name = name;
                mat_copy_persist(&matrix_vars[i].val, &val->v.matrix);
                return 1;
            }
        }
        /* No slots available */
        printf("Error: max %d matrix variables (clear one with: x=0)\n", MAX_MATRIX_VARS);
        return 0;
    }
}

/* Check if variable is defined */
int is_var_defined(int idx)
{
    int i;
    if (idx < 0 || idx >= MAX_SCALAR_VARS) return 0;
    
    if (scalar_vars[idx].defined) return 1;
    
    for (i = 0; i < MAX_MATRIX_VARS; i++) {
        if (matrix_vars[i].defined && matrix_vars[i].name == 'a' + idx) {
            return 1;
        }
    }
    return 0;
}

/* Check if variable is a matrix */
int is_var_matrix(int idx)
{
    int i;
    if (idx < 0 || idx >= MAX_SCALAR_VARS) return 0;
    
    for (i = 0; i < MAX_MATRIX_VARS; i++) {
        if (matrix_vars[i].defined && matrix_vars[i].name == 'a' + idx) {
            return 1;
        }
    }
    return 0;
}

/* ========== User Function Evaluation ========== */

/* Forward declarations for parser */
static apfc param_value_c;
static char param_name = 0;
static int in_user_func = 0;

int eval_user_func(apfc *result, int func_idx, const apfc *arg)
{
    const char *saved_input;
    token_t saved_token;
    int ok;
    
    if (!user_funcs[func_idx].defined) {
        printf("Error: function not defined\n");
        return 0;
    }
    
    saved_input = input_ptr;
    saved_token = current_token;
    param_value_c = *arg;
    param_name = user_funcs[func_idx].param;
    in_user_func = 1;
    
    input_ptr = user_funcs[func_idx].body;
    next_token();
    ok = parse_expr(result);
    
    in_user_func = 0;
    input_ptr = saved_input;
    current_token = saved_token;
    
    return ok;
}

/* Evaluate an arbitrary expression with x substituted */
int eval_expr_with_x(apfc *result, const char *expr, const apfc *x_val)
{
    const char *saved_input;
    token_t saved_token;
    int ok;
    
    saved_input = input_ptr;
    saved_token = current_token;
    param_value_c = *x_val;
    param_name = 'x';
    in_user_func = 1;
    
    input_ptr = expr;
    next_token();
    ok = parse_expr(result);
    
    in_user_func = 0;
    input_ptr = saved_input;
    current_token = saved_token;
    
    return ok;
}

/* Evaluate an arbitrary expression with any variable substituted */
int eval_expr_with_var(apfc *result, const char *expr, char var_name, const apfc *val)
{
    const char *saved_input;
    token_t saved_token;
    int ok;
    
    saved_input = input_ptr;
    saved_token = current_token;
    param_value_c = *val;
    param_name = var_name;
    in_user_func = 1;
    
    input_ptr = expr;
    next_token();
    ok = parse_expr(result);
    
    in_user_func = 0;
    input_ptr = saved_input;
    current_token = saved_token;
    
    return ok;
}

/* Accessors for parser module */
int is_in_user_func(void) { return in_user_func; }
char get_param_name(void) { return param_name; }
apfc *get_param_value(void) { return &param_value_c; }

/* ========== Named Variables (MATLAB compatibility) ========== */
#ifdef HAVE_NAMED_VARS

named_var_t named_vars[MAX_NAMED_VARS];

void init_named_vars(void)
{
    int i;
    for (i = 0; i < MAX_NAMED_VARS; i++) {
        named_vars[i].defined = 0;
        named_vars[i].name[0] = '\0';
    }
}

int find_named_var(const char *name)
{
    int i;
    for (i = 0; i < MAX_NAMED_VARS; i++) {
        if (named_vars[i].defined && str_eq(named_vars[i].name, name)) {
            return i;
        }
    }
    return -1;
}

int create_named_var(const char *name)
{
    int i, idx;
    size_t len;
    
    /* Check if already exists */
    idx = find_named_var(name);
    if (idx >= 0) return idx;
    
    /* Find empty slot */
    for (i = 0; i < MAX_NAMED_VARS; i++) {
        if (!named_vars[i].defined) {
            len = strlen(name);
            if (len >= MAX_VAR_NAME) len = MAX_VAR_NAME - 1;
            memcpy(named_vars[i].name, name, len);
            named_vars[i].name[len] = '\0';
            named_vars[i].defined = 1;
            return i;
        }
    }
    
    printf("Error: too many named variables (max %d)\n", MAX_NAMED_VARS);
    return -1;
}

int set_named_scalar(const char *name, const apfc *val)
{
    int idx = create_named_var(name);
    if (idx < 0) return 0;
    
    named_vars[idx].type = VAL_SCALAR;
    named_vars[idx].val.scalar = *val;
    return 1;
}

int set_named_matrix(const char *name, const matrix_t *val)
{
    int idx;
    
    /* Validate input */
    if (!val || !val->data || val->rows <= 0 || val->cols <= 0) {
        printf("Error: Invalid matrix data\n");
        return 0;
    }
    
    /* Create or get variable slot */
    idx = create_named_var(name);
    if (idx < 0) return 0;
    
    named_vars[idx].type = VAL_MATRIX;
    /* Copy matrix data to persistent storage */
    if (!mat_copy_persist(&named_vars[idx].val.matrix, val)) {
        return 0;
    }
    return 1;
}

int get_named_var(const char *name, value_t *result)
{
    int idx = find_named_var(name);
    if (idx < 0) return 0;
    
    result->type = named_vars[idx].type;
    if (named_vars[idx].type == VAL_SCALAR) {
        result->v.scalar = named_vars[idx].val.scalar;
    } else {
        /* Copy to arena for use in current command */
        mat_copy(&result->v.matrix, &named_vars[idx].val.matrix);
    }
    return 1;
}

#endif /* HAVE_NAMED_VARS */
