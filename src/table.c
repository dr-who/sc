/* table.c - MATLAB-style tables and categorical arrays
 * C89 compliant for Watcom C / DOS
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "apf.h"
#include "apfc.h"
#include "matrix.h"

/* ========== Categorical Array ========== */

#define MAX_CATEGORIES 64
#define MAX_CAT_LABEL 32
#define MAX_CAT_ELEMENTS 256

typedef struct {
    char labels[MAX_CATEGORIES][MAX_CAT_LABEL];  /* Category names */
    int num_categories;                           /* Number of unique categories */
    int indices[MAX_CAT_ELEMENTS];               /* Index into labels for each element */
    int num_elements;                            /* Number of elements */
} categorical_t;

/* Initialize categorical array */
void cat_init(categorical_t *c)
{
    c->num_categories = 0;
    c->num_elements = 0;
}

/* Find or add a category label, returns index */
int cat_find_or_add(categorical_t *c, const char *label)
{
    int i;
    size_t len;
    
    /* Search existing */
    for (i = 0; i < c->num_categories; i++) {
        if (strcmp(c->labels[i], label) == 0) {
            return i;
        }
    }
    
    /* Add new */
    if (c->num_categories >= MAX_CATEGORIES) {
        return -1;
    }
    
    len = strlen(label);
    if (len >= MAX_CAT_LABEL) len = MAX_CAT_LABEL - 1;
    memcpy(c->labels[c->num_categories], label, len);
    c->labels[c->num_categories][len] = '\0';
    
    return c->num_categories++;
}

/* Add element to categorical array */
int cat_add_element(categorical_t *c, const char *label)
{
    int idx;
    
    if (c->num_elements >= MAX_CAT_ELEMENTS) {
        return 0;
    }
    
    idx = cat_find_or_add(c, label);
    if (idx < 0) return 0;
    
    c->indices[c->num_elements++] = idx;
    return 1;
}

/* Get label for element at position i */
const char *cat_get_label(const categorical_t *c, int i)
{
    if (i < 0 || i >= c->num_elements) return "";
    return c->labels[c->indices[i]];
}

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
        matrix_t numeric;      /* For numeric columns */
        categorical_t cat;     /* For categorical columns */
    } data;
} table_col_t;

typedef struct {
    table_col_t cols[MAX_TABLE_COLS];
    int num_cols;
    int num_rows;
} table_t;

/* Global table storage */
#define MAX_TABLES 16
static table_t tables[MAX_TABLES];
static char table_names[MAX_TABLES][MAX_COL_NAME];
static int table_defined[MAX_TABLES];

void table_init_all(void)
{
    int i;
    for (i = 0; i < MAX_TABLES; i++) {
        table_defined[i] = 0;
    }
}

/* Initialize a table */
void table_init(table_t *t)
{
    t->num_cols = 0;
    t->num_rows = 0;
}

/* Add numeric column */
int table_add_numeric_col(table_t *t, const char *name, const matrix_t *data)
{
    extern int mat_copy_persist(matrix_t *dst, const matrix_t *src);
    table_col_t *col;
    size_t len;
    
    if (t->num_cols >= MAX_TABLE_COLS) return 0;
    
    col = &t->cols[t->num_cols];
    len = strlen(name);
    if (len >= MAX_COL_NAME) len = MAX_COL_NAME - 1;
    memcpy(col->name, name, len);
    col->name[len] = '\0';
    
    col->type = COL_NUMERIC;
    /* Use persistent storage so data survives between commands */
    mat_copy_persist(&col->data.numeric, data);
    
    /* Update row count */
    if (t->num_rows == 0) {
        t->num_rows = data->rows;
    } else if (t->num_rows != data->rows) {
        printf("Warning: column '%s' has different row count\n", name);
    }
    
    t->num_cols++;
    return 1;
}

/* Add categorical column */
int table_add_categorical_col(table_t *t, const char *name, const categorical_t *data)
{
    table_col_t *col;
    size_t len;
    
    if (t->num_cols >= MAX_TABLE_COLS) return 0;
    
    col = &t->cols[t->num_cols];
    len = strlen(name);
    if (len >= MAX_COL_NAME) len = MAX_COL_NAME - 1;
    memcpy(col->name, name, len);
    col->name[len] = '\0';
    
    col->type = COL_CATEGORICAL;
    memcpy(&col->data.cat, data, sizeof(categorical_t));
    
    /* Update row count */
    if (t->num_rows == 0) {
        t->num_rows = data->num_elements;
    } else if (t->num_rows != data->num_elements) {
        printf("Warning: column '%s' has different row count\n", name);
    }
    
    t->num_cols++;
    return 1;
}

/* Print table */
void table_print(const table_t *t)
{
    int i, j;
    int col_widths[MAX_TABLE_COLS];
    char buf[64];
    
    if (t->num_cols == 0 || t->num_rows == 0) {
        printf("Empty table\n");
        return;
    }
    
    /* Calculate column widths */
    for (j = 0; j < t->num_cols; j++) {
        col_widths[j] = (int)strlen(t->cols[j].name);
        
        if (t->cols[j].type == COL_CATEGORICAL) {
            /* Check category label widths */
            for (i = 0; i < t->cols[j].data.cat.num_categories; i++) {
                int len = (int)strlen(t->cols[j].data.cat.labels[i]);
                if (len > col_widths[j]) col_widths[j] = len;
            }
        } else {
            /* Check numeric widths */
            for (i = 0; i < t->num_rows; i++) {
                int len;
                apf_to_str(buf, sizeof(buf), &MAT_AT(&t->cols[j].data.numeric, i, 0).re, 6);
                len = (int)strlen(buf);
                if (len > col_widths[j]) col_widths[j] = len;
            }
        }
        
        /* Minimum width */
        if (col_widths[j] < 4) col_widths[j] = 4;
        /* Add padding */
        col_widths[j] += 2;
    }
    
    /* Print header */
    printf("    ");
    for (j = 0; j < t->num_cols; j++) {
        printf("%*s", col_widths[j], t->cols[j].name);
    }
    printf("\n");
    
    /* Print separator */
    printf("    ");
    for (j = 0; j < t->num_cols; j++) {
        int k;
        for (k = 0; k < col_widths[j]; k++) {
            printf("_");
        }
    }
    printf("\n\n");
    
    /* Print rows */
    for (i = 0; i < t->num_rows; i++) {
        printf("    ");
        for (j = 0; j < t->num_cols; j++) {
            if (t->cols[j].type == COL_CATEGORICAL) {
                printf("%*s", col_widths[j], cat_get_label(&t->cols[j].data.cat, i));
            } else {
                apf_to_str(buf, sizeof(buf), &MAT_AT(&t->cols[j].data.numeric, i, 0).re, 6);
                printf("%*s", col_widths[j], buf);
            }
        }
        printf("\n");
    }
}

/* Print specific rows of table (for head/tail) */
void table_print_rows(const table_t *t, int start, int count)
{
    int i, j;
    int col_widths[MAX_TABLE_COLS];
    char buf[64];
    int end;
    
    if (t->num_cols == 0 || t->num_rows == 0) {
        printf("Empty table\n");
        return;
    }
    
    end = start + count;
    if (end > t->num_rows) end = t->num_rows;
    if (start < 0) start = 0;
    
    /* Calculate column widths */
    for (j = 0; j < t->num_cols; j++) {
        col_widths[j] = (int)strlen(t->cols[j].name);
        
        if (t->cols[j].type == COL_CATEGORICAL) {
            for (i = 0; i < t->cols[j].data.cat.num_categories; i++) {
                int len = (int)strlen(t->cols[j].data.cat.labels[i]);
                if (len > col_widths[j]) col_widths[j] = len;
            }
        } else {
            for (i = start; i < end; i++) {
                int len;
                apf_to_str(buf, sizeof(buf), &MAT_AT(&t->cols[j].data.numeric, i, 0).re, 6);
                len = (int)strlen(buf);
                if (len > col_widths[j]) col_widths[j] = len;
            }
        }
        
        if (col_widths[j] < 4) col_widths[j] = 4;
        col_widths[j] += 2;
    }
    
    /* Print header */
    printf("    ");
    for (j = 0; j < t->num_cols; j++) {
        printf("%*s", col_widths[j], t->cols[j].name);
    }
    printf("\n");
    
    /* Print separator */
    printf("    ");
    for (j = 0; j < t->num_cols; j++) {
        int k;
        for (k = 0; k < col_widths[j]; k++) {
            printf("_");
        }
    }
    printf("\n\n");
    
    /* Print rows */
    for (i = start; i < end; i++) {
        printf("    ");
        for (j = 0; j < t->num_cols; j++) {
            if (t->cols[j].type == COL_CATEGORICAL) {
                printf("%*s", col_widths[j], cat_get_label(&t->cols[j].data.cat, i));
            } else {
                apf_to_str(buf, sizeof(buf), &MAT_AT(&t->cols[j].data.numeric, i, 0).re, 6);
                printf("%*s", col_widths[j], buf);
            }
        }
        printf("\n");
    }
    
    printf("\n(%d rows shown, %d total)\n", end - start, t->num_rows);
}

/* Store table with name */
int table_store(const char *name, const table_t *t)
{
    int i;
    size_t len;
    
    
    
    /* Find existing or empty slot */
    for (i = 0; i < MAX_TABLES; i++) {
        if (table_defined[i] && strcmp(table_names[i], name) == 0) {
            /* Replace existing */
            memcpy(&tables[i], t, sizeof(table_t));
            return 1;
        }
    }
    
    /* Find empty slot */
    for (i = 0; i < MAX_TABLES; i++) {
        if (!table_defined[i]) {
            len = strlen(name);
            if (len >= MAX_COL_NAME) len = MAX_COL_NAME - 1;
            memcpy(table_names[i], name, len);
            table_names[i][len] = '\0';
            memcpy(&tables[i], t, sizeof(table_t));
            table_defined[i] = 1;
            return 1;
        }
    }
    
    printf("Error: maximum tables reached\n");
    return 0;
}

/* Retrieve table by name */
table_t *table_get(const char *name)
{
    int i;
    for (i = 0; i < MAX_TABLES; i++) {
        if (table_defined[i] && strcmp(table_names[i], name) == 0) {
            return &tables[i];
        }
    }
    return NULL;
}

/* ========== Global categorical storage ========== */

#define MAX_CATEGORICALS 16
static categorical_t categoricals[MAX_CATEGORICALS];
static char cat_names[MAX_CATEGORICALS][MAX_COL_NAME];
static int cat_defined[MAX_CATEGORICALS];

void cat_init_all(void)
{
    int i;
    for (i = 0; i < MAX_CATEGORICALS; i++) {
        cat_defined[i] = 0;
    }
}

/* Store categorical with name */
int cat_store(const char *name, const categorical_t *c)
{
    int i;
    size_t len;
    
    /* Find existing or empty slot */
    for (i = 0; i < MAX_CATEGORICALS; i++) {
        if (cat_defined[i] && strcmp(cat_names[i], name) == 0) {
            memcpy(&categoricals[i], c, sizeof(categorical_t));
            return 1;
        }
    }
    
    for (i = 0; i < MAX_CATEGORICALS; i++) {
        if (!cat_defined[i]) {
            len = strlen(name);
            if (len >= MAX_COL_NAME) len = MAX_COL_NAME - 1;
            memcpy(cat_names[i], name, len);
            cat_names[i][len] = '\0';
            memcpy(&categoricals[i], c, sizeof(categorical_t));
            cat_defined[i] = 1;
            return 1;
        }
    }
    
    return 0;
}

/* Retrieve categorical by name */
categorical_t *cat_get(const char *name)
{
    int i;
    for (i = 0; i < MAX_CATEGORICALS; i++) {
        if (cat_defined[i] && strcmp(cat_names[i], name) == 0) {
            return &categoricals[i];
        }
    }
    return NULL;
}

/* Print categorical */
void cat_print(const categorical_t *c)
{
    int i;
    for (i = 0; i < c->num_elements; i++) {
        printf("     %s\n", cat_get_label(c, i));
    }
}

/* Get number of elements in categorical */
int cat_num_elements(const categorical_t *c)
{
    if (!c) return 0;
    return c->num_elements;
}

/* Get max label width in categorical */
int cat_max_width(const categorical_t *c)
{
    int i, max_w = 0;
    if (!c) return 0;
    for (i = 0; i < c->num_categories; i++) {
        int len = (int)strlen(c->labels[i]);
        if (len > max_w) max_w = len;
    }
    return max_w;
}

/* Get numeric column from table by name, returns pointer to matrix or NULL */
matrix_t *table_get_column(table_t *t, const char *col_name)
{
    int i;
    if (!t) return NULL;
    
    for (i = 0; i < t->num_cols; i++) {
        if (strcmp(t->cols[i].name, col_name) == 0) {
            if (t->cols[i].type == COL_NUMERIC) {
                return &t->cols[i].data.numeric;
            }
            return NULL;  /* Not a numeric column */
        }
    }
    return NULL;
}

/* Get categorical column indices as numeric vector */
int table_get_cat_indices(table_t *t, const char *col_name, matrix_t *result)
{
    int i, j;
    if (!t) return 0;
    
    for (i = 0; i < t->num_cols; i++) {
        if (strcmp(t->cols[i].name, col_name) == 0) {
            if (t->cols[i].type == COL_CATEGORICAL) {
                categorical_t *c = &t->cols[i].data.cat;
                mat_zero(result, c->num_elements, 1);
                if (!result->data) return 0;
                
                for (j = 0; j < c->num_elements; j++) {
                    apf_from_int(&MAT_AT(result, j, 0).re, c->indices[j] + 1);
                    apf_zero(&MAT_AT(result, j, 0).im);
                }
                return 1;
            }
            return 0;
        }
    }
    return 0;
}

/* Get column index by name, returns -1 if not found */
int table_find_col(table_t *t, const char *col_name)
{
    int i;
    if (!t) return -1;
    
    for (i = 0; i < t->num_cols; i++) {
        if (strcmp(t->cols[i].name, col_name) == 0) {
            return i;
        }
    }
    return -1;
}

/* groupsummary: Group table by one column and apply aggregation to another
 * Creates a new table with unique values of group_col and aggregated values
 * method: "min", "max", "mean", "sum", "count", "numel"
 */
int table_groupsummary(table_t *result, const table_t *src, 
                       const char *group_col, const char *method, const char *data_col)
{
    extern int mat_copy_persist(matrix_t *dst, const matrix_t *src);
    matrix_t *grp_data, *val_data;
    int grp_idx, val_idx;
    int i, j, n_groups;
    long *group_ids;   /* unique group values */
    int *group_counts; /* count per group */
    apf *group_sums;   /* sum per group (for mean/sum) */
    apf *group_mins;   /* min per group */
    apf *group_maxs;   /* max per group */
    int max_groups = 256;
    
    table_init(result);
    
    /* Find columns */
    grp_idx = table_find_col((table_t*)src, group_col);
    val_idx = table_find_col((table_t*)src, data_col);
    
    if (grp_idx < 0 || val_idx < 0) {
        printf("Error: column not found\n");
        return 0;
    }
    
    if (src->cols[grp_idx].type != COL_NUMERIC || src->cols[val_idx].type != COL_NUMERIC) {
        printf("Error: groupsummary requires numeric columns\n");
        return 0;
    }
    
    grp_data = &((table_t*)src)->cols[grp_idx].data.numeric;
    val_data = &((table_t*)src)->cols[val_idx].data.numeric;
    
    /* Allocate temporary arrays */
    group_ids = (long*)malloc(max_groups * sizeof(long));
    group_counts = (int*)malloc(max_groups * sizeof(int));
    group_sums = (apf*)malloc(max_groups * sizeof(apf));
    group_mins = (apf*)malloc(max_groups * sizeof(apf));
    group_maxs = (apf*)malloc(max_groups * sizeof(apf));
    
    if (!group_ids || !group_counts || !group_sums || !group_mins || !group_maxs) {
        free(group_ids); free(group_counts); free(group_sums);
        free(group_mins); free(group_maxs);
        return 0;
    }
    
    n_groups = 0;
    
    /* Process each row */
    for (i = 0; i < src->num_rows; i++) {
        long grp_val = apf_to_long(&MAT_AT(grp_data, i, 0).re);
        apf *val = &MAT_AT(val_data, i, 0).re;
        int found = -1;
        
        /* Find existing group */
        for (j = 0; j < n_groups; j++) {
            if (group_ids[j] == grp_val) {
                found = j;
                break;
            }
        }
        
        if (found < 0) {
            /* New group */
            if (n_groups >= max_groups) continue;
            found = n_groups++;
            group_ids[found] = grp_val;
            group_counts[found] = 0;
            apf_zero(&group_sums[found]);
            apf_copy(&group_mins[found], val);
            apf_copy(&group_maxs[found], val);
        }
        
        /* Update aggregates */
        group_counts[found]++;
        {
            apf tmp;
            apf_add(&tmp, &group_sums[found], val);
            apf_copy(&group_sums[found], &tmp);
        }
        if (apf_cmp(val, &group_mins[found]) < 0) {
            apf_copy(&group_mins[found], val);
        }
        if (apf_cmp(val, &group_maxs[found]) > 0) {
            apf_copy(&group_maxs[found], val);
        }
    }
    
    /* Create result columns */
    {
        matrix_t grp_result, val_result;
        extern int mat_zero_persist(matrix_t *m, int rows, int cols);
        char result_col_name[64];
        
        mat_zero_persist(&grp_result, n_groups, 1);
        mat_zero_persist(&val_result, n_groups, 1);
        
        for (i = 0; i < n_groups; i++) {
            apf_from_int(&MAT_AT(&grp_result, i, 0).re, group_ids[i]);
            apf_zero(&MAT_AT(&grp_result, i, 0).im);
            
            if (strcmp(method, "min") == 0) {
                apf_copy(&MAT_AT(&val_result, i, 0).re, &group_mins[i]);
            } else if (strcmp(method, "max") == 0) {
                apf_copy(&MAT_AT(&val_result, i, 0).re, &group_maxs[i]);
            } else if (strcmp(method, "sum") == 0) {
                apf_copy(&MAT_AT(&val_result, i, 0).re, &group_sums[i]);
            } else if (strcmp(method, "mean") == 0) {
                apf cnt;
                apf_from_int(&cnt, group_counts[i]);
                apf_div(&MAT_AT(&val_result, i, 0).re, &group_sums[i], &cnt);
            } else if (strcmp(method, "count") == 0 || strcmp(method, "numel") == 0) {
                apf_from_int(&MAT_AT(&val_result, i, 0).re, group_counts[i]);
            }
            apf_zero(&MAT_AT(&val_result, i, 0).im);
        }
        
        /* Add columns to result table */
        table_add_numeric_col(result, group_col, &grp_result);
        
        /* Create result column name like "min_time" */
        sprintf(result_col_name, "%s_%s", method, data_col);
        table_add_numeric_col(result, result_col_name, &val_result);
    }
    
    free(group_ids);
    free(group_counts);
    free(group_sums);
    free(group_mins);
    free(group_maxs);
    
    return 1;
}

/* Sort table by first column (ascending) */
int table_sortrows(table_t *t)
{
    int i, j, k;
    
    if (!t || t->num_rows <= 1 || t->num_cols == 0) return 1;
    
    /* Simple bubble sort by first column */
    for (i = 0; i < t->num_rows - 1; i++) {
        for (j = 0; j < t->num_rows - 1 - i; j++) {
            int swap = 0;
            
            if (t->cols[0].type == COL_NUMERIC) {
                if (apf_cmp(&MAT_AT(&t->cols[0].data.numeric, j, 0).re,
                           &MAT_AT(&t->cols[0].data.numeric, j+1, 0).re) > 0) {
                    swap = 1;
                }
            }
            
            if (swap) {
                /* Swap all columns for rows j and j+1 */
                for (k = 0; k < t->num_cols; k++) {
                    if (t->cols[k].type == COL_NUMERIC) {
                        apfc tmp = MAT_AT(&t->cols[k].data.numeric, j, 0);
                        MAT_AT(&t->cols[k].data.numeric, j, 0) = MAT_AT(&t->cols[k].data.numeric, j+1, 0);
                        MAT_AT(&t->cols[k].data.numeric, j+1, 0) = tmp;
                    } else if (t->cols[k].type == COL_CATEGORICAL) {
                        int tmp = t->cols[k].data.cat.indices[j];
                        t->cols[k].data.cat.indices[j] = t->cols[k].data.cat.indices[j+1];
                        t->cols[k].data.cat.indices[j+1] = tmp;
                    }
                }
            }
        }
    }
    return 1;
}

/* Get table row count */
int table_height(const table_t *t)
{
    return t ? t->num_rows : 0;
}


/* Helper: bin timestamp by type
 * Returns binned value (start of period as Unix timestamp, or day-of-week 1-7)
 */
static long bin_timestamp(long ts, const char *bin_type)
{
    static const int days_in_month[] = {31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31};
    long days_since_epoch = ts / 86400;
    long year = 1970, month = 1, day = 1;
    long hour, minute, dow;
    int i, is_leap, year_days, dm;
    
    /* Find year */
    while (1) {
        is_leap = ((year % 4 == 0 && year % 100 != 0) || (year % 400 == 0));
        year_days = 365 + is_leap;
        if (days_since_epoch < year_days) break;
        days_since_epoch -= year_days;
        year++;
    }
    
    /* Find month */
    is_leap = ((year % 4 == 0 && year % 100 != 0) || (year % 400 == 0));
    for (i = 0; i < 12; i++) {
        dm = days_in_month[i] + (i == 1 && is_leap ? 1 : 0);
        if (days_since_epoch < dm) {
            month = i + 1;
            day = days_since_epoch + 1;
            break;
        }
        days_since_epoch -= dm;
    }
    
    hour = (ts % 86400) / 3600;
    minute = (ts % 3600) / 60;
    
    /* Suppress unused variable warnings */
    (void)day;
    (void)hour;
    (void)minute;
    
    /* Calculate day of week (0=Thu for Jan 1, 1970) */
    dow = ((ts / 86400) + 4) % 7;  /* 0=Sun, 1=Mon, ..., 6=Sat */
    if (dow == 0) dow = 7;  /* Make 1=Mon, 7=Sun */
    
    if (strcmp(bin_type, "year") == 0) {
        /* Return start of year as timestamp */
        long new_days = 0;
        for (i = 1970; i < year; i++) {
            new_days += 365;
            if ((i % 4 == 0 && i % 100 != 0) || (i % 400 == 0)) new_days++;
        }
        return new_days * 86400;
    }
    else if (strcmp(bin_type, "month") == 0) {
        /* Return start of month as timestamp */
        long new_days = 0;
        for (i = 1970; i < year; i++) {
            new_days += 365;
            if ((i % 4 == 0 && i % 100 != 0) || (i % 400 == 0)) new_days++;
        }
        for (i = 1; i < month; i++) {
            new_days += days_in_month[i-1];
            if (i == 2 && is_leap) new_days++;
        }
        return new_days * 86400;
    }
    else if (strcmp(bin_type, "day") == 0) {
        return (ts / 86400) * 86400;
    }
    else if (strcmp(bin_type, "hour") == 0) {
        return (ts / 3600) * 3600;
    }
    else if (strcmp(bin_type, "minute") == 0) {
        return (ts / 60) * 60;
    }
    else if (strcmp(bin_type, "dayname") == 0 || strcmp(bin_type, "weekday") == 0) {
        return dow;  /* 1=Mon, 7=Sun */
    }
    else if (strcmp(bin_type, "week") == 0) {
        /* Return start of week (Monday) */
        long day_ts = (ts / 86400) * 86400;
        dow = ((ts / 86400) + 4) % 7;  /* 0=Sun */
        if (dow == 0) dow = 7;
        return day_ts - (dow - 1) * 86400;
    }
    else if (strcmp(bin_type, "quarter") == 0) {
        /* Return start of quarter */
        int q_month = ((month - 1) / 3) * 3 + 1;
        long new_days = 0;
        for (i = 1970; i < year; i++) {
            new_days += 365;
            if ((i % 4 == 0 && i % 100 != 0) || (i % 400 == 0)) new_days++;
        }
        for (i = 1; i < q_month; i++) {
            new_days += days_in_month[i-1];
            if (i == 2 && is_leap) new_days++;
        }
        return new_days * 86400;
    }
    
    return ts;  /* Default: no binning */
}

/* groupsummary with time-based binning */
int table_groupsummary_binned(table_t *result, const table_t *src,
                              const char *group_col, const char *bin_type,
                              const char *method, const char *data_col)
{
    matrix_t *grp_data, *val_data;
    int grp_idx, val_idx;
    int i, j, n_groups;
    long *group_keys;
    apf *group_sums, *group_mins, *group_maxs;
    int *group_counts;
    int is_sum, is_mean, is_min, is_max, is_numel;
    char result_col_name[64];
    
    if (!src || !result) return 0;
    
    is_sum = (strcmp(method, "sum") == 0);
    is_mean = (strcmp(method, "mean") == 0);
    is_min = (strcmp(method, "min") == 0);
    is_max = (strcmp(method, "max") == 0);
    is_numel = (strcmp(method, "numel") == 0 || strcmp(method, "count") == 0);
    
    grp_idx = table_find_col((table_t*)src, group_col);
    val_idx = table_find_col((table_t*)src, data_col);
    
    if (grp_idx < 0) {
        printf("Error: group column '%s' not found\n", group_col);
        return 0;
    }
    if (val_idx < 0 && !is_numel) {
        printf("Error: data column '%s' not found\n", data_col);
        return 0;
    }
    
    grp_data = &((table_t*)src)->cols[grp_idx].data.numeric;
    if (!is_numel) {
        val_data = &((table_t*)src)->cols[val_idx].data.numeric;
    }
    
    /* Allocate working arrays */
    group_keys = (long *)malloc(src->num_rows * sizeof(long));
    group_sums = (apf *)malloc(src->num_rows * sizeof(apf));
    group_mins = (apf *)malloc(src->num_rows * sizeof(apf));
    group_maxs = (apf *)malloc(src->num_rows * sizeof(apf));
    group_counts = (int *)malloc(src->num_rows * sizeof(int));
    
    if (!group_keys || !group_sums || !group_counts) {
        free(group_keys); free(group_sums); free(group_mins); 
        free(group_maxs); free(group_counts);
        return 0;
    }
    
    /* Initialize */
    n_groups = 0;
    
    /* Process each row */
    for (i = 0; i < src->num_rows; i++) {
        long raw_key = apf_to_long(&MAT_AT(grp_data, i, 0).re);
        long binned_key = bin_timestamp(raw_key, bin_type);
        int found = -1;
        apf val;
        
        if (!is_numel) {
            apf_copy(&val, &MAT_AT(val_data, i, 0).re);
        }
        
        /* Find existing group */
        for (j = 0; j < n_groups; j++) {
            if (group_keys[j] == binned_key) {
                found = j;
                break;
            }
        }
        
        if (found < 0) {
            /* New group */
            found = n_groups++;
            group_keys[found] = binned_key;
            apf_zero(&group_sums[found]);
            if (!is_numel) {
                apf_copy(&group_mins[found], &val);
                apf_copy(&group_maxs[found], &val);
            }
            group_counts[found] = 0;
        }
        
        /* Aggregate */
        if ((is_sum || is_mean) && !is_numel) {
            apf tmp;
            apf_add(&tmp, &group_sums[found], &val);
            apf_copy(&group_sums[found], &tmp);
        }
        if (is_min && !is_numel) {
            if (apf_cmp(&val, &group_mins[found]) < 0) {
                apf_copy(&group_mins[found], &val);
            }
        }
        if (is_max && !is_numel) {
            if (apf_cmp(&val, &group_maxs[found]) > 0) {
                apf_copy(&group_maxs[found], &val);
            }
        }
        group_counts[found]++;
    }
    
    /* Build result table */
    table_init(result);
    
    /* Group key column */
    {
        matrix_t key_mat;
        char col_name[64];
        mat_zero(&key_mat, n_groups, 1);
        for (j = 0; j < n_groups; j++) {
            apf_from_int(&MAT_AT(&key_mat, j, 0).re, group_keys[j]);
            apf_zero(&MAT_AT(&key_mat, j, 0).im);
        }
        sprintf(col_name, "%s_%s", bin_type, group_col);
        table_add_numeric_col(result, col_name, &key_mat);
    }
    
    /* Result column */
    {
        matrix_t result_mat;
        mat_zero(&result_mat, n_groups, 1);
        
        for (j = 0; j < n_groups; j++) {
            if (is_numel) {
                apf_from_int(&MAT_AT(&result_mat, j, 0).re, group_counts[j]);
            } else if (is_mean) {
                apf cnt;
                apf_from_int(&cnt, group_counts[j]);
                apf_div(&MAT_AT(&result_mat, j, 0).re, &group_sums[j], &cnt);
            } else if (is_min) {
                apf_copy(&MAT_AT(&result_mat, j, 0).re, &group_mins[j]);
            } else if (is_max) {
                apf_copy(&MAT_AT(&result_mat, j, 0).re, &group_maxs[j]);
            } else {
                apf_copy(&MAT_AT(&result_mat, j, 0).re, &group_sums[j]);
            }
            apf_zero(&MAT_AT(&result_mat, j, 0).im);
        }
        
        sprintf(result_col_name, "%s_%s", method, is_numel ? group_col : data_col);
        table_add_numeric_col(result, result_col_name, &result_mat);
    }
    
    free(group_keys);
    free(group_sums);
    free(group_mins);
    free(group_maxs);
    free(group_counts);
    
    return 1;
}
