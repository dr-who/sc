/* matrix_ops.c - Element-wise operations, concatenation, and indexing
 * C89 compliant for Watcom C / DOS
 * 16-bit clean: int is 16-bit, long is 32-bit
 */
#include <stdio.h>
#include <string.h>
#include "matrix.h"
#include "apfx.h"

/* ========== Element-wise Operations ========== */

void mat_add(matrix_t *r, const matrix_t *a, const matrix_t *b)
{
    int i;
    if (!mat_same_size(a, b)) {
        printf("Error: matrix dimensions must match for addition\n");
        mat_zero(r, 1, 1);
        return;
    }
    r->rows = a->rows;
    r->cols = a->cols;
    for (i = 0; i < a->rows * a->cols; i++) {
        apfc_add(&r->data[i], &a->data[i], &b->data[i]);
    }
}

void mat_sub(matrix_t *r, const matrix_t *a, const matrix_t *b)
{
    int i;
    if (!mat_same_size(a, b)) {
        printf("Error: matrix dimensions must match for subtraction\n");
        mat_zero(r, 1, 1);
        return;
    }
    r->rows = a->rows;
    r->cols = a->cols;
    for (i = 0; i < a->rows * a->cols; i++) {
        apfc_sub(&r->data[i], &a->data[i], &b->data[i]);
    }
}

void mat_neg(matrix_t *r, const matrix_t *a)
{
    int i;
    r->rows = a->rows;
    r->cols = a->cols;
    for (i = 0; i < a->rows * a->cols; i++) {
        apfc_neg(&r->data[i], &a->data[i]);
    }
}

void mat_scale(matrix_t *r, const matrix_t *a, const apfc *s)
{
    int i;
    r->rows = a->rows;
    r->cols = a->cols;
    for (i = 0; i < a->rows * a->cols; i++) {
        apfc_mul(&r->data[i], &a->data[i], s);
    }
}

void mat_add_scalar(matrix_t *r, const matrix_t *a, const apfc *s)
{
    int i;
    r->rows = a->rows;
    r->cols = a->cols;
    for (i = 0; i < a->rows * a->cols; i++) {
        apfc_add(&r->data[i], &a->data[i], s);
    }
}

/* Element-wise multiplication: A .* B */
void mat_elemwise_mul(matrix_t *r, const matrix_t *a, const matrix_t *b)
{
    int i;
    if (!mat_same_size(a, b)) {
        printf("Error: matrix dimensions must match for .* operation\n");
        mat_zero(r, 1, 1);
        return;
    }
    r->rows = a->rows;
    r->cols = a->cols;
    for (i = 0; i < a->rows * a->cols; i++) {
        apfc_mul(&r->data[i], &a->data[i], &b->data[i]);
    }
}

/* Element-wise division: A ./ B */
void mat_elemwise_div(matrix_t *r, const matrix_t *a, const matrix_t *b)
{
    int i;
    if (!mat_same_size(a, b)) {
        printf("Error: matrix dimensions must match for ./ operation\n");
        mat_zero(r, 1, 1);
        return;
    }
    r->rows = a->rows;
    r->cols = a->cols;
    for (i = 0; i < a->rows * a->cols; i++) {
        apfc_div(&r->data[i], &a->data[i], &b->data[i]);
    }
}

/* Element-wise power: A .^ scalar */
void mat_elemwise_pow(matrix_t *r, const matrix_t *a, const apfc *p)
{
    int i;
    r->rows = a->rows;
    r->cols = a->cols;
    for (i = 0; i < a->rows * a->cols; i++) {
        apfc_pow(&r->data[i], &a->data[i], p);
    }
}

/* Element-wise power: A .^ B (matrix) */
void mat_elemwise_pow_mat(matrix_t *r, const matrix_t *a, const matrix_t *b)
{
    int i;
    if (!mat_same_size(a, b)) {
        printf("Error: matrix dimensions must match for .^ operation\n");
        mat_zero(r, 1, 1);
        return;
    }
    r->rows = a->rows;
    r->cols = a->cols;
    for (i = 0; i < a->rows * a->cols; i++) {
        apfc_pow(&r->data[i], &a->data[i], &b->data[i]);
    }
}

/* ========== Concatenation ========== */

/* Horizontal concatenation: [A B] */
int mat_hcat(matrix_t *r, const matrix_t *a, const matrix_t *b)
{
    int i, j;
    int new_cols;
    
    if (a->rows != b->rows) {
        printf("Error: row count must match for horizontal concatenation\n");
        return 0;
    }
    
    new_cols = a->cols + b->cols;
    if (new_cols > MAT_MAX_COLS) {
        printf("Error: result matrix too wide\n");
        return 0;
    }
    
    r->rows = a->rows;
    r->cols = new_cols;
    
    for (i = 0; i < a->rows; i++) {
        /* Copy columns from A */
        for (j = 0; j < a->cols; j++) {
            MAT_AT(r, i, j) = MAT_AT(a, i, j);
        }
        /* Copy columns from B */
        for (j = 0; j < b->cols; j++) {
            MAT_AT(r, i, a->cols + j) = MAT_AT(b, i, j);
        }
    }
    
    return 1;
}

/* Vertical concatenation: [A; B] */
int mat_vcat(matrix_t *r, const matrix_t *a, const matrix_t *b)
{
    int i, j;
    int new_rows;
    
    if (a->cols != b->cols) {
        printf("Error: column count must match for vertical concatenation\n");
        return 0;
    }
    
    new_rows = a->rows + b->rows;
    if (new_rows > MAT_MAX_ROWS) {
        printf("Error: result matrix too tall\n");
        return 0;
    }
    
    r->rows = new_rows;
    r->cols = a->cols;
    
    /* Copy rows from A */
    for (i = 0; i < a->rows; i++) {
        for (j = 0; j < a->cols; j++) {
            MAT_AT(r, i, j) = MAT_AT(a, i, j);
        }
    }
    /* Copy rows from B */
    for (i = 0; i < b->rows; i++) {
        for (j = 0; j < b->cols; j++) {
            MAT_AT(r, a->rows + i, j) = MAT_AT(b, i, j);
        }
    }
    
    return 1;
}

/* ========== Indexing/Slicing ========== */

/* Get a single row as a 1xN matrix */
void mat_get_row(matrix_t *r, const matrix_t *m, int row)
{
    int j;
    
    if (row < 0 || row >= m->rows) {
        printf("Error: row index out of bounds\n");
        mat_zero(r, 1, 1);
        return;
    }
    
    r->rows = 1;
    r->cols = m->cols;
    
    for (j = 0; j < m->cols; j++) {
        MAT_AT(r, 0, j) = MAT_AT(m, row, j);
    }
}

/* Get a single column as an Mx1 matrix */
void mat_get_col(matrix_t *r, const matrix_t *m, int col)
{
    int i;
    
    if (col < 0 || col >= m->cols) {
        printf("Error: column index out of bounds\n");
        mat_zero(r, 1, 1);
        return;
    }
    
    r->rows = m->rows;
    r->cols = 1;
    
    for (i = 0; i < m->rows; i++) {
        MAT_AT(r, i, 0) = MAT_AT(m, i, col);
    }
}

/* Set a single row from a vector */
void mat_set_row(matrix_t *m, int row, const matrix_t *v)
{
    int j;
    
    if (row < 0 || row >= m->rows) {
        printf("Error: row index out of bounds\n");
        return;
    }
    if (v->rows * v->cols != m->cols) {
        printf("Error: vector length must match matrix column count\n");
        return;
    }
    
    for (j = 0; j < m->cols; j++) {
        MAT_AT(m, row, j) = v->data[j];
    }
}

/* Set a single column from a vector */
void mat_set_col(matrix_t *m, int col, const matrix_t *v)
{
    int i;
    
    if (col < 0 || col >= m->cols) {
        printf("Error: column index out of bounds\n");
        return;
    }
    if (v->rows * v->cols != m->rows) {
        printf("Error: vector length must match matrix row count\n");
        return;
    }
    
    for (i = 0; i < m->rows; i++) {
        MAT_AT(m, i, col) = v->data[i];
    }
}

/* Extract submatrix: rows [r1,r2], cols [c1,c2] (inclusive, 0-indexed) */
void mat_submat(matrix_t *r, const matrix_t *m, int r1, int r2, int c1, int c2)
{
    int i, j;
    int new_rows, new_cols;
    
    /* Clamp to valid range */
    if (r1 < 0) r1 = 0;
    if (c1 < 0) c1 = 0;
    if (r2 >= m->rows) r2 = m->rows - 1;
    if (c2 >= m->cols) c2 = m->cols - 1;
    
    if (r1 > r2 || c1 > c2) {
        printf("Error: invalid submatrix range\n");
        mat_zero(r, 1, 1);
        return;
    }
    
    new_rows = r2 - r1 + 1;
    new_cols = c2 - c1 + 1;
    
    r->rows = new_rows;
    r->cols = new_cols;
    
    for (i = 0; i < new_rows; i++) {
        for (j = 0; j < new_cols; j++) {
            MAT_AT(r, i, j) = MAT_AT(m, r1 + i, c1 + j);
        }
    }
}
