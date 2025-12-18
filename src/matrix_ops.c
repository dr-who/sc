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
    mat_zero(r, a->rows, a->cols);
    if (!r->data) return;
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
    mat_zero(r, a->rows, a->cols);
    if (!r->data) return;
    for (i = 0; i < a->rows * a->cols; i++) {
        apfc_sub(&r->data[i], &a->data[i], &b->data[i]);
    }
}

void mat_neg(matrix_t *r, const matrix_t *a)
{
    int i;
    mat_zero(r, a->rows, a->cols);
    if (!r->data) return;
    for (i = 0; i < a->rows * a->cols; i++) {
        apfc_neg(&r->data[i], &a->data[i]);
    }
}

void mat_scale(matrix_t *r, const matrix_t *a, const apfc *s)
{
    int i;
    mat_zero(r, a->rows, a->cols);
    if (!r->data) return;
    for (i = 0; i < a->rows * a->cols; i++) {
        apfc_mul(&r->data[i], &a->data[i], s);
    }
}

void mat_add_scalar(matrix_t *r, const matrix_t *a, const apfc *s)
{
    int i;
    mat_zero(r, a->rows, a->cols);
    if (!r->data) return;
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
    mat_zero(r, a->rows, a->cols);
    if (!r->data) return;
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
    mat_zero(r, a->rows, a->cols);
    if (!r->data) return;
    for (i = 0; i < a->rows * a->cols; i++) {
        apfc_div(&r->data[i], &a->data[i], &b->data[i]);
    }
}

/* Element-wise power: A .^ scalar */
void mat_elemwise_pow(matrix_t *r, const matrix_t *a, const apfc *p)
{
    int i;
    mat_zero(r, a->rows, a->cols);
    if (!r->data) return;
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
    mat_zero(r, a->rows, a->cols);
    if (!r->data) return;
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
    
    mat_zero(r, a->rows, new_cols);
    if (!r->data) return 0;
    
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
    
    mat_zero(r, new_rows, a->cols);
    if (!r->data) return 0;
    
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
    
    mat_zero(r, 1, m->cols);
    if (!r->data) return;
    
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
    
    mat_zero(r, m->rows, 1);
    if (!r->data) return;
    
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
    
    mat_zero(r, new_rows, new_cols);
    if (!r->data) return;
    
    for (i = 0; i < new_rows; i++) {
        for (j = 0; j < new_cols; j++) {
            MAT_AT(r, i, j) = MAT_AT(m, r1 + i, c1 + j);
        }
    }
}

/* ============================================================
 * Array manipulation functions
 * ============================================================ */

/* Flip matrix left-right */
void mat_fliplr(matrix_t *r, const matrix_t *m)
{
    int i, j;
    mat_zero(r, m->rows, m->cols);
    if (!r->data) return;
    for (i = 0; i < m->rows; i++) {
        for (j = 0; j < m->cols; j++) {
            MAT_AT(r, i, m->cols - 1 - j) = MAT_AT(m, i, j);
        }
    }
}

/* Flip matrix up-down */
void mat_flipud(matrix_t *r, const matrix_t *m)
{
    int i, j;
    mat_zero(r, m->rows, m->cols);
    if (!r->data) return;
    for (i = 0; i < m->rows; i++) {
        for (j = 0; j < m->cols; j++) {
            MAT_AT(r, m->rows - 1 - i, j) = MAT_AT(m, i, j);
        }
    }
}

/* Flip vector (horizontal for row, vertical for column/matrix) */
void mat_flip(matrix_t *r, const matrix_t *m)
{
    if (m->rows == 1) {
        mat_fliplr(r, m);
    } else {
        mat_flipud(r, m);
    }
}

/* Rotate matrix 90 degrees counterclockwise */
void mat_rot90(matrix_t *r, const matrix_t *m)
{
    int i, j;
    mat_zero(r, m->cols, m->rows);
    if (!r->data) return;
    for (i = 0; i < m->rows; i++) {
        for (j = 0; j < m->cols; j++) {
            MAT_AT(r, m->cols - 1 - j, i) = MAT_AT(m, i, j);
        }
    }
}

/* Upper triangular part */
void mat_triu(matrix_t *r, const matrix_t *m)
{
    int i, j;
    mat_zero(r, m->rows, m->cols);
    if (!r->data) return;
    for (i = 0; i < m->rows; i++) {
        for (j = i; j < m->cols; j++) {
            MAT_AT(r, i, j) = MAT_AT(m, i, j);
        }
    }
}

/* Lower triangular part */
void mat_tril(matrix_t *r, const matrix_t *m)
{
    int i, j;
    mat_zero(r, m->rows, m->cols);
    if (!r->data) return;
    for (i = 0; i < m->rows; i++) {
        for (j = 0; j <= i && j < m->cols; j++) {
            MAT_AT(r, i, j) = MAT_AT(m, i, j);
        }
    }
}

/* Cumulative sum */
void mat_cumsum(matrix_t *r, const matrix_t *m)
{
    int i, j, n;
    apfc sum;
    n = m->rows * m->cols;
    mat_zero(r, m->rows, m->cols);
    if (!r->data) return;
    apf_zero(&sum.re);
    apf_zero(&sum.im);
    for (i = 0; i < n; i++) {
        int ri = i % m->rows;
        int ci = i / m->rows;
        apfc_add(&sum, &sum, &MAT_AT(m, ri, ci));
        MAT_AT(r, ri, ci) = sum;
    }
    (void)j;
}

/* Cumulative product */
void mat_cumprod(matrix_t *r, const matrix_t *m)
{
    int i, j, n;
    apfc prod;
    n = m->rows * m->cols;
    mat_zero(r, m->rows, m->cols);
    if (!r->data) return;
    apf_from_int(&prod.re, 1);
    apf_zero(&prod.im);
    for (i = 0; i < n; i++) {
        int ri = i % m->rows;
        int ci = i / m->rows;
        apfc_mul(&prod, &prod, &MAT_AT(m, ri, ci));
        MAT_AT(r, ri, ci) = prod;
    }
    (void)j;
}

/* Differences (diff) */
void mat_diff(matrix_t *r, const matrix_t *m)
{
    int i, n;
    n = m->rows * m->cols;
    if (n <= 1) {
        mat_zero(r, 1, 1);
        return;
    }
    mat_zero(r, m->rows == 1 ? 1 : m->rows, m->cols == 1 || m->rows > 1 ? m->cols : m->cols - 1);
    if (!r->data) return;
    
    if (m->rows == 1) {
        /* Row vector */
        for (i = 0; i < m->cols - 1; i++) {
            apfc_sub(&MAT_AT(r, 0, i), &MAT_AT(m, 0, i + 1), &MAT_AT(m, 0, i));
        }
        r->cols = m->cols - 1;
    } else {
        /* Column vector or matrix - diff along columns */
        int j;
        mat_zero(r, m->rows - 1, m->cols);
        for (j = 0; j < m->cols; j++) {
            for (i = 0; i < m->rows - 1; i++) {
                apfc_sub(&MAT_AT(r, i, j), &MAT_AT(m, i + 1, j), &MAT_AT(m, i, j));
            }
        }
    }
}

/* Normalize to unit length */
void mat_normalize(matrix_t *r, const matrix_t *m)
{
    int i, j;
    apf sum_sq, norm_val, sq;
    mat_zero(r, m->rows, m->cols);
    if (!r->data) return;
    
    /* Compute Frobenius norm */
    apf_zero(&sum_sq);
    for (i = 0; i < m->rows; i++) {
        for (j = 0; j < m->cols; j++) {
            apfc *v = &MAT_AT(m, i, j);
            apf_mul(&sq, &v->re, &v->re);
            apf_add(&sum_sq, &sum_sq, &sq);
            apf_mul(&sq, &v->im, &v->im);
            apf_add(&sum_sq, &sum_sq, &sq);
        }
    }
    apf_sqrt(&norm_val, &sum_sq);
    
    if (apf_is_zero(&norm_val)) {
        mat_copy(r, m);
        return;
    }
    
    /* Divide each element */
    for (i = 0; i < m->rows; i++) {
        for (j = 0; j < m->cols; j++) {
            apf_div(&MAT_AT(r, i, j).re, &MAT_AT(m, i, j).re, &norm_val);
            apf_div(&MAT_AT(r, i, j).im, &MAT_AT(m, i, j).im, &norm_val);
        }
    }
}

/* Sort elements */
void mat_sort(matrix_t *r, const matrix_t *m)
{
    int i, j, n;
    mat_copy(r, m);
    if (!r->data) return;
    n = r->rows * r->cols;
    
    /* Bubble sort - simple for small matrices */
    for (i = 0; i < n - 1; i++) {
        for (j = 0; j < n - 1 - i; j++) {
            int r1 = j % r->rows;
            int c1 = j / r->rows;
            int r2 = (j + 1) % r->rows;
            int c2 = (j + 1) / r->rows;
            if (apf_cmp(&MAT_AT(r, r1, c1).re, &MAT_AT(r, r2, c2).re) > 0) {
                apfc tmp = MAT_AT(r, r1, c1);
                MAT_AT(r, r1, c1) = MAT_AT(r, r2, c2);
                MAT_AT(r, r2, c2) = tmp;
            }
        }
    }
}

/* Squeeze - remove singleton dimensions */
void mat_squeeze(matrix_t *r, const matrix_t *m)
{
    mat_copy(r, m);
    if (r->rows == 1 && r->cols > 1) {
        /* Already a row vector */
    } else if (r->cols == 1 && r->rows > 1) {
        /* Already a column vector */
    } else if (r->rows == 1 && r->cols == 1) {
        /* Scalar - keep as 1x1 */
    }
    /* For now just copy - full squeeze would need dimension info */
}

/* Cumulative maximum */
void mat_cummax(matrix_t *r, const matrix_t *m)
{
    int i, n;
    apfc maxval;
    n = m->rows * m->cols;
    mat_zero(r, m->rows, m->cols);
    if (!r->data || n == 0) return;
    
    /* Initialize with first element */
    maxval = MAT_AT(m, 0, 0);
    MAT_AT(r, 0, 0) = maxval;
    
    for (i = 1; i < n; i++) {
        int ri = i % m->rows;
        int ci = i / m->rows;
        if (apf_cmp(&MAT_AT(m, ri, ci).re, &maxval.re) > 0) {
            maxval = MAT_AT(m, ri, ci);
        }
        MAT_AT(r, ri, ci) = maxval;
    }
}

/* Cumulative minimum */
void mat_cummin(matrix_t *r, const matrix_t *m)
{
    int i, n;
    apfc minval;
    n = m->rows * m->cols;
    mat_zero(r, m->rows, m->cols);
    if (!r->data || n == 0) return;
    
    /* Initialize with first element */
    minval = MAT_AT(m, 0, 0);
    MAT_AT(r, 0, 0) = minval;
    
    for (i = 1; i < n; i++) {
        int ri = i % m->rows;
        int ci = i / m->rows;
        if (apf_cmp(&MAT_AT(m, ri, ci).re, &minval.re) < 0) {
            minval = MAT_AT(m, ri, ci);
        }
        MAT_AT(r, ri, ci) = minval;
    }
}

/* Unique sorted elements */
void mat_unique(matrix_t *r, const matrix_t *m)
{
    int i, j, n, count;
    matrix_t sorted;
    n = m->rows * m->cols;
    
    if (n == 0) {
        mat_zero(r, 1, 1);
        return;
    }
    
    /* First sort */
    mat_sort(&sorted, m);
    
    /* Allocate max possible size */
    mat_zero(r, 1, n);
    if (!r->data) return;
    
    /* Remove duplicates */
    count = 0;
    for (i = 0; i < n; i++) {
        int ri = i % sorted.rows;
        int ci = i / sorted.rows;
        int is_dup = 0;
        for (j = 0; j < count; j++) {
            if (apf_eq(&MAT_AT(&sorted, ri, ci).re, &MAT_AT(r, 0, j).re) &&
                apf_eq(&MAT_AT(&sorted, ri, ci).im, &MAT_AT(r, 0, j).im)) {
                is_dup = 1;
                break;
            }
        }
        if (!is_dup) {
            MAT_AT(r, 0, count) = MAT_AT(&sorted, ri, ci);
            count++;
        }
    }
    r->cols = count > 0 ? count : 1;
}

/* Find nonzero elements - returns indices */
void mat_find(matrix_t *r, const matrix_t *m)
{
    int i, n, count;
    n = m->rows * m->cols;
    
    /* Allocate max possible size */
    mat_zero(r, n > 0 ? n : 1, 1);
    if (!r->data) return;
    
    count = 0;
    for (i = 0; i < n; i++) {
        int ri = i % m->rows;
        int ci = i / m->rows;
        if (!apf_is_zero(&MAT_AT(m, ri, ci).re) || !apf_is_zero(&MAT_AT(m, ri, ci).im)) {
            apf_from_int(&MAT_AT(r, count, 0).re, i + 1); /* 1-indexed */
            apf_zero(&MAT_AT(r, count, 0).im);
            count++;
        }
    }
    r->rows = count > 0 ? count : 1;
    if (count == 0) {
        /* Empty result - return empty column vector */
        apf_zero(&MAT_AT(r, 0, 0).re);
    }
}

/* Sort rows by first column */
void mat_sortrows(matrix_t *r, const matrix_t *m)
{
    int nrows, ncols, i, j, c;
    mat_copy(r, m);
    if (!r->data) return;
    
    nrows = r->rows;
    ncols = r->cols;
    
    /* Bubble sort rows by first column */
    for (i = 0; i < nrows - 1; i++) {
        for (j = 0; j < nrows - 1 - i; j++) {
            if (apf_cmp(&MAT_AT(r, j, 0).re, &MAT_AT(r, j + 1, 0).re) > 0) {
                /* Swap entire rows */
                for (c = 0; c < ncols; c++) {
                    apfc tmp = MAT_AT(r, j, c);
                    MAT_AT(r, j, c) = MAT_AT(r, j + 1, c);
                    MAT_AT(r, j + 1, c) = tmp;
                }
            }
        }
    }
}

/* Rescale to [0,1] range */
void mat_rescale(matrix_t *r, const matrix_t *m)
{
    int i, j, n;
    apf minval, maxval, range;
    
    mat_zero(r, m->rows, m->cols);
    if (!r->data) return;
    
    n = m->rows * m->cols;
    if (n == 0) return;
    
    /* Find min and max */
    minval = MAT_AT(m, 0, 0).re;
    maxval = MAT_AT(m, 0, 0).re;
    for (i = 0; i < m->rows; i++) {
        for (j = 0; j < m->cols; j++) {
            if (apf_cmp(&MAT_AT(m, i, j).re, &minval) < 0) {
                minval = MAT_AT(m, i, j).re;
            }
            if (apf_cmp(&MAT_AT(m, i, j).re, &maxval) > 0) {
                maxval = MAT_AT(m, i, j).re;
            }
        }
    }
    
    /* Compute range */
    apf_sub(&range, &maxval, &minval);
    
    if (apf_is_zero(&range)) {
        /* All elements are the same - return zeros */
        return;
    }
    
    /* Rescale: (x - min) / range */
    for (i = 0; i < m->rows; i++) {
        for (j = 0; j < m->cols; j++) {
            apf tmp;
            apf_sub(&tmp, &MAT_AT(m, i, j).re, &minval);
            apf_div(&MAT_AT(r, i, j).re, &tmp, &range);
            apf_zero(&MAT_AT(r, i, j).im);
        }
    }
}

/* ============================================================
 * Additional statistical reduce functions
 * ============================================================ */

/* Range: max - min */
void mat_range(apfc *r, const matrix_t *m)
{
    apfc minval, maxval;
    mat_min(&minval, m);
    mat_max(&maxval, m);
    apf_sub(&r->re, &maxval.re, &minval.re);
    apf_zero(&r->im);
}

/* RMS: sqrt(mean(x^2)) */
void mat_rms(apfc *r, const matrix_t *m)
{
    int i, j, n;
    apf sum_sq, sq, n_val, mean_sq;
    n = m->rows * m->cols;
    if (n == 0) {
        apf_zero(&r->re);
        apf_zero(&r->im);
        return;
    }
    apf_zero(&sum_sq);
    for (i = 0; i < m->rows; i++) {
        for (j = 0; j < m->cols; j++) {
            apf re_sq, im_sq;
            apf_mul(&re_sq, &MAT_AT(m, i, j).re, &MAT_AT(m, i, j).re);
            apf_mul(&im_sq, &MAT_AT(m, i, j).im, &MAT_AT(m, i, j).im);
            apf_add(&sq, &re_sq, &im_sq);
            apf_add(&sum_sq, &sum_sq, &sq);
        }
    }
    apf_from_int(&n_val, n);
    apf_div(&mean_sq, &sum_sq, &n_val);
    apf_sqrt(&r->re, &mean_sq);
    apf_zero(&r->im);
}

/* Geometric mean: exp(mean(log(x))) */
void mat_geomean(apfc *r, const matrix_t *m)
{
    int i, j, n;
    apf sum_log, log_val, n_val;
    n = m->rows * m->cols;
    if (n == 0) {
        apf_zero(&r->re);
        apf_zero(&r->im);
        return;
    }
    apf_zero(&sum_log);
    for (i = 0; i < m->rows; i++) {
        for (j = 0; j < m->cols; j++) {
            apfx_log(&log_val, &MAT_AT(m, i, j).re);
            apf_add(&sum_log, &sum_log, &log_val);
        }
    }
    apf_from_int(&n_val, n);
    apf_div(&sum_log, &sum_log, &n_val);
    apfx_exp(&r->re, &sum_log);
    apf_zero(&r->im);
}

/* Harmonic mean: n / sum(1/x) */
void mat_harmmean(apfc *r, const matrix_t *m)
{
    int i, j, n;
    apf sum_inv, inv_val, n_val, one;
    n = m->rows * m->cols;
    if (n == 0) {
        apf_zero(&r->re);
        apf_zero(&r->im);
        return;
    }
    apf_zero(&sum_inv);
    apf_from_int(&one, 1);
    for (i = 0; i < m->rows; i++) {
        for (j = 0; j < m->cols; j++) {
            apf_div(&inv_val, &one, &MAT_AT(m, i, j).re);
            apf_add(&sum_inv, &sum_inv, &inv_val);
        }
    }
    apf_from_int(&n_val, n);
    apf_div(&r->re, &n_val, &sum_inv);
    apf_zero(&r->im);
}

/* Sum of squares */
void mat_sumsq(apfc *r, const matrix_t *m)
{
    int i, j;
    apf sum_sq, sq;
    apf_zero(&sum_sq);
    for (i = 0; i < m->rows; i++) {
        for (j = 0; j < m->cols; j++) {
            apf re_sq, im_sq;
            apf_mul(&re_sq, &MAT_AT(m, i, j).re, &MAT_AT(m, i, j).re);
            apf_mul(&im_sq, &MAT_AT(m, i, j).im, &MAT_AT(m, i, j).im);
            apf_add(&sq, &re_sq, &im_sq);
            apf_add(&sum_sq, &sum_sq, &sq);
        }
    }
    apf_copy(&r->re, &sum_sq);
    apf_zero(&r->im);
}

/* Mean of squares */
void mat_meansq(apfc *r, const matrix_t *m)
{
    int n;
    apf n_val;
    mat_sumsq(r, m);
    n = m->rows * m->cols;
    if (n > 0) {
        apf_from_int(&n_val, n);
        apf_div(&r->re, &r->re, &n_val);
    }
}

/* Number of elements */
void mat_numel_r(apfc *r, const matrix_t *m)
{
    apf_from_int(&r->re, m->rows * m->cols);
    apf_zero(&r->im);
}

/* Length: max(rows, cols) */
void mat_length_r(apfc *r, const matrix_t *m)
{
    int len = (m->rows > m->cols) ? m->rows : m->cols;
    apf_from_int(&r->re, len);
    apf_zero(&r->im);
}

/* Number of dimensions */
void mat_ndims_r(apfc *r, const matrix_t *m)
{
    int nd;
    if (m->rows == 1 && m->cols == 1) {
        nd = 0;  /* scalar */
    } else if (m->rows == 1 || m->cols == 1) {
        nd = 1;  /* vector */
    } else {
        nd = 2;  /* matrix */
    }
    apf_from_int(&r->re, nd);
    apf_zero(&r->im);
}

/* Number of rows */
void mat_nrows_r(apfc *r, const matrix_t *m)
{
    apf_from_int(&r->re, m->rows);
    apf_zero(&r->im);
}

/* Number of columns */
void mat_ncols_r(apfc *r, const matrix_t *m)
{
    apf_from_int(&r->re, m->cols);
    apf_zero(&r->im);
}

/* Width (same as cols for matrix) */
void mat_width_r(apfc *r, const matrix_t *m)
{
    apf_from_int(&r->re, m->cols);
    apf_zero(&r->im);
}

/* Height (same as rows for matrix) */
void mat_height_r(apfc *r, const matrix_t *m)
{
    apf_from_int(&r->re, m->rows);
    apf_zero(&r->im);
}

/* Skewness: E[(X-μ)³] / σ³ */
void mat_skewness(apfc *r, const matrix_t *m)
{
    int i, j, n;
    apf mean_val, sum_val, std_val, m3, tmp, n_val, var_val;
    n = m->rows * m->cols;
    if (n == 0) { apf_zero(&r->re); apf_zero(&r->im); return; }
    
    /* Calculate mean */
    apf_zero(&sum_val);
    for (i = 0; i < m->rows; i++) {
        for (j = 0; j < m->cols; j++) {
            apf_add(&sum_val, &sum_val, &MAT_AT(m, i, j).re);
        }
    }
    apf_from_int(&n_val, n);
    apf_div(&mean_val, &sum_val, &n_val);
    
    /* Calculate variance and third moment */
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
        apf_div(&r->re, &m3, &tmp);
    } else {
        apf_zero(&r->re);
    }
    apf_zero(&r->im);
}

/* Excess kurtosis: E[(X-μ)⁴] / σ⁴ - 3 */
void mat_kurtosis(apfc *r, const matrix_t *m)
{
    int i, j, n;
    apf mean_val, sum_val, m4, tmp, n_val, var_val, three;
    n = m->rows * m->cols;
    if (n == 0) { apf_zero(&r->re); apf_zero(&r->im); return; }
    
    /* Calculate mean */
    apf_zero(&sum_val);
    for (i = 0; i < m->rows; i++) {
        for (j = 0; j < m->cols; j++) {
            apf_add(&sum_val, &sum_val, &MAT_AT(m, i, j).re);
        }
    }
    apf_from_int(&n_val, n);
    apf_div(&mean_val, &sum_val, &n_val);
    
    /* Calculate variance and fourth moment */
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
        apf_div(&r->re, &m4, &tmp);
        apf_from_int(&three, 3);
        apf_sub(&r->re, &r->re, &three);
    } else {
        apf_zero(&r->re);
    }
    apf_zero(&r->im);
}

/* Mean absolute deviation: mean(|x - mean(x)|) */
void mat_mad(apfc *r, const matrix_t *m)
{
    int i, j, n;
    apf mean_val, sum_val, n_val, mad_sum;
    n = m->rows * m->cols;
    if (n == 0) { apf_zero(&r->re); apf_zero(&r->im); return; }
    
    /* Calculate mean */
    apf_zero(&sum_val);
    for (i = 0; i < m->rows; i++) {
        for (j = 0; j < m->cols; j++) {
            apf_add(&sum_val, &sum_val, &MAT_AT(m, i, j).re);
        }
    }
    apf_from_int(&n_val, n);
    apf_div(&mean_val, &sum_val, &n_val);
    
    /* Calculate mean absolute deviation */
    apf_zero(&mad_sum);
    for (i = 0; i < m->rows; i++) {
        for (j = 0; j < m->cols; j++) {
            apf diff, abs_diff;
            apf_sub(&diff, &MAT_AT(m, i, j).re, &mean_val);
            apf_abs(&abs_diff, &diff);
            apf_add(&mad_sum, &mad_sum, &abs_diff);
        }
    }
    apf_div(&r->re, &mad_sum, &n_val);
    apf_zero(&r->im);
}

/* any: true if any element is nonzero */
void mat_any(apfc *r, const matrix_t *m)
{
    int i, j, found = 0;
    for (i = 0; i < m->rows && !found; i++) {
        for (j = 0; j < m->cols && !found; j++) {
            if (!apf_is_zero(&MAT_AT(m, i, j).re) || 
                !apf_is_zero(&MAT_AT(m, i, j).im)) {
                found = 1;
            }
        }
    }
    apf_from_int(&r->re, found);
    apf_zero(&r->im);
}

/* all: true if all elements are nonzero */
void mat_all(apfc *r, const matrix_t *m)
{
    int i, j, all_nonzero = 1;
    for (i = 0; i < m->rows && all_nonzero; i++) {
        for (j = 0; j < m->cols && all_nonzero; j++) {
            if (apf_is_zero(&MAT_AT(m, i, j).re) && 
                apf_is_zero(&MAT_AT(m, i, j).im)) {
                all_nonzero = 0;
            }
        }
    }
    apf_from_int(&r->re, all_nonzero);
    apf_zero(&r->im);
}

/* isempty: true if matrix has zero dimensions */
void mat_isempty(apfc *r, const matrix_t *m)
{
    int is_empty = (m->rows == 0 || m->cols == 0);
    apf_from_int(&r->re, is_empty);
    apf_zero(&r->im);
}

/* numel: number of elements */
void mat_numel(apfc *r, const matrix_t *m)
{
    apf_from_int(&r->re, m->rows * m->cols);
    apf_zero(&r->im);
}

/* length: max dimension */
/* Note: mat_length_r is the reduce version, already defined above */

/* ndims: number of dimensions (always 2 for our matrices) */
void mat_ndims(apfc *r, const matrix_t *m)
{
    (void)m;
    apf_from_int(&r->re, 2);
    apf_zero(&r->im);
}

/* isscalar: true if 1x1 matrix */
void mat_isscalar(apfc *r, const matrix_t *m)
{
    apf_from_int(&r->re, (m->rows == 1 && m->cols == 1) ? 1 : 0);
    apf_zero(&r->im);
}

/* isvector: true if row or column vector */
void mat_isvector(apfc *r, const matrix_t *m)
{
    apf_from_int(&r->re, (m->rows == 1 || m->cols == 1) ? 1 : 0);
    apf_zero(&r->im);
}

/* isrow: true if row vector */
void mat_isrow(apfc *r, const matrix_t *m)
{
    apf_from_int(&r->re, (m->rows == 1) ? 1 : 0);
    apf_zero(&r->im);
}

/* iscolumn: true if column vector */
void mat_iscolumn(apfc *r, const matrix_t *m)
{
    apf_from_int(&r->re, (m->cols == 1) ? 1 : 0);
    apf_zero(&r->im);
}

/* ismatrix: true if 2D (always true for us) */
void mat_ismatrix(apfc *r, const matrix_t *m)
{
    (void)m;
    apf_from_int(&r->re, 1);
    apf_zero(&r->im);
}

/* issquare: true if square matrix */
void mat_issquare(apfc *r, const matrix_t *m)
{
    apf_from_int(&r->re, (m->rows == m->cols) ? 1 : 0);
    apf_zero(&r->im);
}

/* skewness: E[(X-μ)³] / σ³ */

/* center - subtract mean from each element */
void mat_center(matrix_t *r, const matrix_t *m)
{
    int i, j, n;
    apf mean_val, sum, count;
    n = m->rows * m->cols;
    if (n == 0) {
        mat_zero(r, m->rows, m->cols);
        return;
    }
    apf_zero(&sum);
    for (i = 0; i < m->rows; i++) {
        for (j = 0; j < m->cols; j++) {
            apf_add(&sum, &sum, &MAT_AT(m, i, j).re);
        }
    }
    apf_from_int(&count, n);
    apf_div(&mean_val, &sum, &count);
    mat_zero(r, m->rows, m->cols);
    for (i = 0; i < m->rows; i++) {
        for (j = 0; j < m->cols; j++) {
            apf_sub(&MAT_AT(r, i, j).re, &MAT_AT(m, i, j).re, &mean_val);
        }
    }
}

/* gradient - numerical gradient (forward difference) */
void mat_gradient(matrix_t *r, const matrix_t *m)
{
    int n = m->rows * m->cols;
    int i;
    if (n == 0) {
        mat_zero(r, m->rows, m->cols);
        return;
    }
    if (m->rows == 1 || m->cols == 1) {
        /* Vector case */
        mat_zero(r, m->rows, m->cols);
        for (i = 0; i < n - 1; i++) {
            int ri = (m->rows == 1) ? 0 : i;
            int ci = (m->cols == 1) ? 0 : i;
            int ri1 = (m->rows == 1) ? 0 : i + 1;
            int ci1 = (m->cols == 1) ? 0 : i + 1;
            apf_sub(&MAT_AT(r, ri, ci).re, &MAT_AT(m, ri1, ci1).re, &MAT_AT(m, ri, ci).re);
        }
        if (n > 0) {
            int ri = (m->rows == 1) ? 0 : n - 1;
            int ci = (m->cols == 1) ? 0 : n - 1;
            apf_zero(&MAT_AT(r, ri, ci).re);
        }
    } else {
        mat_copy(r, m);  /* For 2D matrices, just copy for now */
    }
}

/* diff2 - second order difference */
void mat_diff2(matrix_t *r, const matrix_t *m)
{
    int n = m->rows * m->cols;
    int i;
    if (n < 3) {
        mat_zero(r, 1, n > 2 ? n - 2 : 0);
        return;
    }
    if (m->rows == 1 || m->cols == 1) {
        mat_zero(r, m->rows == 1 ? 1 : m->rows - 2, m->cols == 1 ? 1 : m->cols - 2);
        for (i = 0; i < n - 2; i++) {
            int ri = (m->rows == 1) ? 0 : i;
            int ci = (m->cols == 1) ? 0 : i;
            int ri1 = (m->rows == 1) ? 0 : i + 1;
            int ci1 = (m->cols == 1) ? 0 : i + 1;
            int ri2 = (m->rows == 1) ? 0 : i + 2;
            int ci2 = (m->cols == 1) ? 0 : i + 2;
            apf d1, d2;
            apf_sub(&d1, &MAT_AT(m, ri1, ci1).re, &MAT_AT(m, ri, ci).re);
            apf_sub(&d2, &MAT_AT(m, ri2, ci2).re, &MAT_AT(m, ri1, ci1).re);
            apf_sub(&MAT_AT(r, ri, ci).re, &d2, &d1);
        }
    } else {
        mat_zero(r, 0, 0);  /* For 2D matrices, return empty */
    }
}

/* detrend - remove linear trend */
void mat_detrend(matrix_t *r, const matrix_t *m)
{
    int n = m->rows * m->cols;
    int i;
    apf sumx, sumy, sumxy, sumxx, xbar, ybar, slope, intercept;
    apf x, y, fitted, count, tmp;
    
    if (n == 0) {
        mat_zero(r, m->rows, m->cols);
        return;
    }
    
    /* Calculate linear regression coefficients */
    apf_zero(&sumx); apf_zero(&sumy); apf_zero(&sumxy); apf_zero(&sumxx);
    for (i = 0; i < n; i++) {
        int ri = (m->rows == 1) ? 0 : (m->cols == 1 ? i : i / m->cols);
        int ci = (m->cols == 1) ? 0 : (m->rows == 1 ? i : i % m->cols);
        apf_from_int(&x, i);
        apf_copy(&y, &MAT_AT(m, ri, ci).re);
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
    apf_mul(&tmp, &xbar, &ybar);
    apf_mul(&tmp, &tmp, &count);
    apf_sub(&slope, &sumxy, &tmp);
    apf_mul(&tmp, &xbar, &xbar);
    apf_mul(&tmp, &tmp, &count);
    apf_sub(&intercept, &sumxx, &tmp);
    if (!apf_is_zero(&intercept)) {
        apf_div(&slope, &slope, &intercept);
    }
    
    /* intercept = ybar - slope * xbar */
    apf_mul(&tmp, &slope, &xbar);
    apf_sub(&intercept, &ybar, &tmp);
    
    /* Subtract fitted values */
    mat_zero(r, m->rows, m->cols);
    for (i = 0; i < n; i++) {
        int ri = (m->rows == 1) ? 0 : (m->cols == 1 ? i : i / m->cols);
        int ci = (m->cols == 1) ? 0 : (m->rows == 1 ? i : i % m->cols);
        apf_from_int(&x, i);
        apf_mul(&fitted, &slope, &x);
        apf_add(&fitted, &fitted, &intercept);
        apf_sub(&MAT_AT(r, ri, ci).re, &MAT_AT(m, ri, ci).re, &fitted);
    }
}

/* nnz: count of non-zero elements */
void mat_nnz(apfc *r, const matrix_t *m)
{
    int i, j, count = 0;
    for (i = 0; i < m->rows; i++) {
        for (j = 0; j < m->cols; j++) {
            if (!apf_is_zero(&MAT_AT(m, i, j).re) ||
                !apf_is_zero(&MAT_AT(m, i, j).im)) {
                count++;
            }
        }
    }
    apf_from_int(&r->re, count);
    apf_zero(&r->im);
}

/* mode: most frequent element */
void mat_mode(apfc *r, const matrix_t *m)
{
    int n, i, j, max_count, curr_count;
    apfc mode_val;
    
    n = m->rows * m->cols;
    if (n == 0) {
        apf_zero(&r->re);
        apf_zero(&r->im);
        return;
    }
    
    max_count = 0;
    apf_zero(&mode_val.re);
    apf_zero(&mode_val.im);
    
    for (i = 0; i < n; i++) {
        int ri = i % m->rows;
        int ci = i / m->rows;
        apfc val = MAT_AT(m, ri, ci);
        curr_count = 0;
        
        for (j = 0; j < n; j++) {
            int rj = j % m->rows;
            int cj = j / m->rows;
            if (apf_eq(&val.re, &MAT_AT(m, rj, cj).re)) {
                curr_count++;
            }
        }
        
        if (curr_count > max_count || 
            (curr_count == max_count && apf_cmp(&val.re, &mode_val.re) < 0)) {
            max_count = curr_count;
            mode_val = val;
        }
    }
    
    r->re = mode_val.re;
    r->im = mode_val.im;
}

/* iqr: interquartile range (Q3 - Q1) */
void mat_iqr(apfc *r, const matrix_t *m)
{
    int n, i, j;
    matrix_t sorted;
    apf q1, q3, pos_val, n_minus_1, quarter, three_quarters;
    int idx;
    
    n = m->rows * m->cols;
    if (n == 0) {
        apf_zero(&r->re);
        apf_zero(&r->im);
        return;
    }
    
    /* Copy and sort */
    mat_copy(&sorted, m);
    for (i = 0; i < n - 1; i++) {
        for (j = 0; j < n - 1 - i; j++) {
            int r1 = j % sorted.rows;
            int c1 = j / sorted.rows;
            int r2 = (j+1) % sorted.rows;
            int c2 = (j+1) / sorted.rows;
            if (apf_cmp(&MAT_AT(&sorted, r1, c1).re,
                       &MAT_AT(&sorted, r2, c2).re) > 0) {
                apfc tmp = MAT_AT(&sorted, r1, c1);
                MAT_AT(&sorted, r1, c1) = MAT_AT(&sorted, r2, c2);
                MAT_AT(&sorted, r2, c2) = tmp;
            }
        }
    }
    
    apf_from_int(&n_minus_1, n - 1);
    
    /* Q1 (25th percentile) */
    apf_from_str(&quarter, "0.25");
    apf_mul(&pos_val, &quarter, &n_minus_1);
    idx = apf_to_long(&pos_val);
    if (idx < 0) idx = 0;
    if (idx >= n - 1) {
        int ri = (n-1) % sorted.rows;
        int ci = (n-1) / sorted.rows;
        q1 = MAT_AT(&sorted, ri, ci).re;
    } else {
        int r1 = idx % sorted.rows;
        int c1 = idx / sorted.rows;
        int r2 = (idx+1) % sorted.rows;
        int c2 = (idx+1) / sorted.rows;
        apf frac, one_minus, tmp, idx_apf;
        apf_from_int(&idx_apf, idx);
        apf_sub(&frac, &pos_val, &idx_apf);
        apf_from_int(&one_minus, 1);
        apf_sub(&one_minus, &one_minus, &frac);
        apf_mul(&q1, &MAT_AT(&sorted, r1, c1).re, &one_minus);
        apf_mul(&tmp, &MAT_AT(&sorted, r2, c2).re, &frac);
        apf_add(&q1, &q1, &tmp);
    }
    
    /* Q3 (75th percentile) */
    apf_from_str(&three_quarters, "0.75");
    apf_mul(&pos_val, &three_quarters, &n_minus_1);
    idx = apf_to_long(&pos_val);
    if (idx < 0) idx = 0;
    if (idx >= n - 1) {
        int ri = (n-1) % sorted.rows;
        int ci = (n-1) / sorted.rows;
        q3 = MAT_AT(&sorted, ri, ci).re;
    } else {
        int r1 = idx % sorted.rows;
        int c1 = idx / sorted.rows;
        int r2 = (idx+1) % sorted.rows;
        int c2 = (idx+1) / sorted.rows;
        apf frac, one_minus, tmp, idx_apf;
        apf_from_int(&idx_apf, idx);
        apf_sub(&frac, &pos_val, &idx_apf);
        apf_from_int(&one_minus, 1);
        apf_sub(&one_minus, &one_minus, &frac);
        apf_mul(&q3, &MAT_AT(&sorted, r1, c1).re, &one_minus);
        apf_mul(&tmp, &MAT_AT(&sorted, r2, c2).re, &frac);
        apf_add(&q3, &q3, &tmp);
    }
    
    /* IQR = Q3 - Q1 */
    apf_sub(&r->re, &q3, &q1);
    apf_zero(&r->im);
}

