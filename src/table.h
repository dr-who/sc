/* table.h - MATLAB-style tables and categorical arrays */
#ifndef TABLE_H
#define TABLE_H

#include "matrix.h"

/* ========== Categorical Array ========== */

#define MAX_CATEGORIES 64
#define MAX_CAT_LABEL 32
#define MAX_CAT_ELEMENTS 256

typedef struct {
    char labels[MAX_CATEGORIES][MAX_CAT_LABEL];
    int num_categories;
    int indices[MAX_CAT_ELEMENTS];
    int num_elements;
} categorical_t;

void cat_init(categorical_t *c);
int cat_find_or_add(categorical_t *c, const char *label);
int cat_add_element(categorical_t *c, const char *label);
const char *cat_get_label(const categorical_t *c, int i);
void cat_print(const categorical_t *c);
int cat_num_elements(const categorical_t *c);
int cat_max_width(const categorical_t *c);

/* Global categorical storage */
void cat_init_all(void);
int cat_store(const char *name, const categorical_t *c);
categorical_t *cat_get(const char *name);

/* ========== Table Structure ========== */

#define MAX_TABLE_COLS 16
#define MAX_COL_NAME 32

typedef enum {
    COL_NUMERIC,
    COL_CATEGORICAL
} col_type_t;

typedef struct {
    char name[MAX_COL_NAME];
    col_type_t type;
    union {
        matrix_t numeric;
        categorical_t cat;
    } data;
} table_col_t;

typedef struct {
    table_col_t cols[MAX_TABLE_COLS];
    int num_cols;
    int num_rows;
} table_t;

void table_init_all(void);
void table_init(table_t *t);
int table_add_numeric_col(table_t *t, const char *name, const matrix_t *data);
int table_add_categorical_col(table_t *t, const char *name, const categorical_t *data);
void table_print(const table_t *t);
int table_store(const char *name, const table_t *t);
table_t *table_get(const char *name);
matrix_t *table_get_column(table_t *t, const char *col_name);
int table_get_cat_indices(table_t *t, const char *col_name, matrix_t *result);
int table_groupsummary(table_t *result, const table_t *src, const char *group_col, const char *method, const char *data_col);
int table_sortrows(table_t *t);
int table_height(const table_t *t);

#endif /* TABLE_H */
int table_groupsummary_binned(table_t *result, const table_t *src,
                              const char *group_col, const char *bin_type,
                              const char *method, const char *data_col);
