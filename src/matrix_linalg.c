/* matrix_linalg.c - Linear algebra operations for scalc
 * C89 compliant for Watcom C / DOS
 * 16-bit clean: int is 16-bit, long is 32-bit
 */
#include <stdio.h>
#include <string.h>
#include "matrix.h"
#include "apfx.h"

/* ========== LU Decomposition ========== */

int mat_lu(matrix_t *L, matrix_t *U, int *perm, const matrix_t *m)
{
    int i, j, k, n, pivot_row;
    apfc tmp, factor;
    apf abs_pivot, abs_test, zero_thresh;
    
    if (!mat_is_square(m)) {
        printf("Error: LU decomposition requires square matrix\n");
        return 0;
    }
    
    n = m->rows;
    mat_copy(U, m);
    mat_identity(L, n);
    for (i = 0; i < n; i++) perm[i] = i;
    
    apf_from_int(&zero_thresh, 1);
    zero_thresh.exp = -200;
    
    for (k = 0; k < n - 1; k++) {
        pivot_row = k;
        apfc_abs(&abs_pivot, &MAT_AT(U, k, k));
        
        for (i = k + 1; i < n; i++) {
            apfc_abs(&abs_test, &MAT_AT(U, i, k));
            if (apf_lt(&abs_pivot, &abs_test)) {
                abs_pivot = abs_test;
                pivot_row = i;
            }
        }
        
        if (apf_lt(&abs_pivot, &zero_thresh)) {
            return 0;  /* Singular */
        }
        
        if (pivot_row != k) {
            int p_tmp = perm[k];
            perm[k] = perm[pivot_row];
            perm[pivot_row] = p_tmp;
            
            for (j = 0; j < n; j++) {
                tmp = MAT_AT(U, k, j);
                MAT_AT(U, k, j) = MAT_AT(U, pivot_row, j);
                MAT_AT(U, pivot_row, j) = tmp;
            }
            for (j = 0; j < k; j++) {
                tmp = MAT_AT(L, k, j);
                MAT_AT(L, k, j) = MAT_AT(L, pivot_row, j);
                MAT_AT(L, pivot_row, j) = tmp;
            }
        }
        
        for (i = k + 1; i < n; i++) {
            apfc prod;
            apfc_div(&factor, &MAT_AT(U, i, k), &MAT_AT(U, k, k));
            MAT_AT(L, i, k) = factor;
            for (j = k; j < n; j++) {
                apfc_mul(&prod, &factor, &MAT_AT(U, k, j));
                apfc_sub(&MAT_AT(U, i, j), &MAT_AT(U, i, j), &prod);
            }
        }
    }
    
    return 1;
}

/* ========== Solve Linear System ========== */

int mat_solve(matrix_t *x, const matrix_t *A, const matrix_t *b)
{
    matrix_t L, U;
    int perm[MAT_MAX_ROWS];
    int i, j, n;
    matrix_t y, pb;
    apfc sum, prod;
    
    if (!mat_is_square(A)) {
        printf("Error: system matrix must be square\n");
        return 0;
    }
    if (A->rows != b->rows) {
        printf("Error: dimension mismatch in A\\b\n");
        return 0;
    }
    
    n = A->rows;
    
    if (!mat_lu(&L, &U, perm, A)) {
        printf("Error: matrix is singular\n");
        return 0;
    }
    
    /* Apply permutation to b */
    mat_zero(&pb, n, b->cols);
    for (i = 0; i < n; i++) {
        for (j = 0; j < b->cols; j++) {
            MAT_AT(&pb, i, j) = MAT_AT(b, perm[i], j);
        }
    }
    
    /* Forward substitution: L*y = pb */
    mat_zero(&y, n, b->cols);
    for (j = 0; j < b->cols; j++) {
        for (i = 0; i < n; i++) {
            int k;
            sum = MAT_AT(&pb, i, j);
            for (k = 0; k < i; k++) {
                apfc_mul(&prod, &MAT_AT(&L, i, k), &MAT_AT(&y, k, j));
                apfc_sub(&sum, &sum, &prod);
            }
            MAT_AT(&y, i, j) = sum;
        }
    }
    
    /* Back substitution: U*x = y */
    mat_zero(x, n, b->cols);
    for (j = 0; j < b->cols; j++) {
        for (i = n - 1; i >= 0; i--) {
            int k;
            sum = MAT_AT(&y, i, j);
            for (k = i + 1; k < n; k++) {
                apfc_mul(&prod, &MAT_AT(&U, i, k), &MAT_AT(x, k, j));
                apfc_sub(&sum, &sum, &prod);
            }
            apfc_div(&MAT_AT(x, i, j), &sum, &MAT_AT(&U, i, i));
        }
    }
    
    return 1;
}

/* A\B - left division (solve A*X = B) */
int mat_mldivide(matrix_t *r, const matrix_t *A, const matrix_t *B)
{
    return mat_solve(r, A, B);
}

/* A/B - right division (A * inv(B)) */
int mat_mrdivide(matrix_t *r, const matrix_t *A, const matrix_t *B)
{
    matrix_t B_inv;
    
    if (!mat_inv(&B_inv, B)) {
        return 0;
    }
    mat_mul(r, A, &B_inv);
    return 1;
}

/* ========== SVD (2x2 only) ========== */

int mat_svd2(apfc *s1, apfc *s2, const matrix_t *m)
{
    /* For 2x2 matrix, singular values are sqrt of eigenvalues of A'*A */
    matrix_t ata;
    apfc eig1, eig2;
    
    if (m->rows != 2 || m->cols != 2) {
        printf("Error: svd only works for 2x2 matrices\n");
        return 0;
    }
    
    mat_conj_transpose(&ata, m);
    {
        matrix_t tmp;
        mat_mul(&tmp, &ata, m);
        mat_copy(&ata, &tmp);
    }
    
    if (!mat_eig2(&eig1, &eig2, &ata)) {
        return 0;
    }
    
    apfc_sqrt(s1, &eig1);
    apfc_sqrt(s2, &eig2);
    
    return 2;
}

/* ========== Additional Norms (1-norm, inf-norm) ========== */

void mat_norm_1(apfc *r, const matrix_t *m)
{
    int i, j;
    apf max_col, col_sum, abs_val;
    
    apf_zero(&max_col);
    
    for (j = 0; j < m->cols; j++) {
        apf_zero(&col_sum);
        for (i = 0; i < m->rows; i++) {
            apfc_abs(&abs_val, &MAT_AT(m, i, j));
            apf_add(&col_sum, &col_sum, &abs_val);
        }
        if (apf_gt(&col_sum, &max_col)) {
            max_col = col_sum;
        }
    }
    
    r->re = max_col;
    apf_zero(&r->im);
}

void mat_norm_inf(apfc *r, const matrix_t *m)
{
    int i, j;
    apf max_row, row_sum, abs_val;
    
    apf_zero(&max_row);
    
    for (i = 0; i < m->rows; i++) {
        apf_zero(&row_sum);
        for (j = 0; j < m->cols; j++) {
            apfc_abs(&abs_val, &MAT_AT(m, i, j));
            apf_add(&row_sum, &row_sum, &abs_val);
        }
        if (apf_gt(&row_sum, &max_row)) {
            max_row = row_sum;
        }
    }
    
    r->re = max_row;
    apf_zero(&r->im);
}

/* ========== Aggregations ========== */

void mat_sum(apfc *r, const matrix_t *m)
{
    int i;
    apfc sum;
    
    apf_zero(&sum.re);
    apf_zero(&sum.im);
    
    for (i = 0; i < m->rows * m->cols; i++) {
        apfc tmp;
        apfc_add(&tmp, &sum, &m->data[i]);
        sum = tmp;
    }
    
    *r = sum;
}

void mat_sum_rows(matrix_t *r, const matrix_t *m)
{
    int i, j;
    
    mat_zero(r, 1, m->cols);
    
    for (j = 0; j < m->cols; j++) {
        apfc sum;
        apf_zero(&sum.re);
        apf_zero(&sum.im);
        for (i = 0; i < m->rows; i++) {
            apfc tmp;
            apfc_add(&tmp, &sum, &MAT_AT(m, i, j));
            sum = tmp;
        }
        MAT_AT(r, 0, j) = sum;
    }
}

void mat_sum_cols(matrix_t *r, const matrix_t *m)
{
    int i, j;
    
    mat_zero(r, m->rows, 1);
    
    for (i = 0; i < m->rows; i++) {
        apfc sum;
        apf_zero(&sum.re);
        apf_zero(&sum.im);
        for (j = 0; j < m->cols; j++) {
            apfc tmp;
            apfc_add(&tmp, &sum, &MAT_AT(m, i, j));
            sum = tmp;
        }
        MAT_AT(r, i, 0) = sum;
    }
}

void mat_mean(apfc *r, const matrix_t *m)
{
    apfc sum, count;
    int n = m->rows * m->cols;
    
    mat_sum(&sum, m);
    apf_from_int(&count.re, n);
    apf_zero(&count.im);
    apfc_div(r, &sum, &count);
}

void mat_min(apfc *r, const matrix_t *m)
{
    int i;
    apf min_abs, cur_abs;
    int min_idx = 0;
    
    if (m->rows * m->cols == 0) {
        apf_zero(&r->re);
        apf_zero(&r->im);
        return;
    }
    
    apfc_abs(&min_abs, &m->data[0]);
    
    for (i = 1; i < m->rows * m->cols; i++) {
        apfc_abs(&cur_abs, &m->data[i]);
        if (apf_lt(&cur_abs, &min_abs)) {
            min_abs = cur_abs;
            min_idx = i;
        }
    }
    
    *r = m->data[min_idx];
}

void mat_max(apfc *r, const matrix_t *m)
{
    int i;
    apf max_abs, cur_abs;
    int max_idx = 0;
    
    if (m->rows * m->cols == 0) {
        apf_zero(&r->re);
        apf_zero(&r->im);
        return;
    }
    
    apfc_abs(&max_abs, &m->data[0]);
    
    for (i = 1; i < m->rows * m->cols; i++) {
        apfc_abs(&cur_abs, &m->data[i]);
        if (apf_gt(&cur_abs, &max_abs)) {
            max_abs = cur_abs;
            max_idx = i;
        }
    }
    
    *r = m->data[max_idx];
}

/* Product of all elements */
void mat_prod(apfc *r, const matrix_t *m)
{
    int i;
    apfc prod;
    apfc tmp;
    
    if (m->rows * m->cols == 0) {
        apf_from_int(&r->re, 1);
        apf_zero(&r->im);
        return;
    }
    
    prod = m->data[0];
    
    for (i = 1; i < m->rows * m->cols; i++) {
        apfc_mul(&tmp, &prod, &m->data[i]);
        prod = tmp;
    }
    
    *r = prod;
}

/* ========== Standard Deviation ========== */

void mat_std(apfc *r, const matrix_t *m)
{
    int i;
    int n;
    apfc mean_val;
    apfc sum_sq;
    apfc diff;
    apfc sq;
    apfc count;
    apfc variance;
    apfc tmp;
    
    n = m->rows * m->cols;
    if (n <= 1) {
        apf_zero(&r->re);
        apf_zero(&r->im);
        return;
    }
    
    /* Calculate mean */
    mat_mean(&mean_val, m);
    
    /* Calculate sum of squared differences */
    apf_zero(&sum_sq.re);
    apf_zero(&sum_sq.im);
    
    for (i = 0; i < n; i++) {
        apfc_sub(&diff, &m->data[i], &mean_val);
        apfc_mul(&sq, &diff, &diff);
        apfc_add(&tmp, &sum_sq, &sq);
        sum_sq = tmp;
    }
    
    /* Divide by n-1 (sample standard deviation) */
    apf_from_int(&count.re, n - 1);
    apf_zero(&count.im);
    apfc_div(&variance, &sum_sq, &count);
    
    /* Square root */
    apfc_sqrt(r, &variance);
}

/* ========== Variance ========== */

void mat_var(apfc *r, const matrix_t *m)
{
    int i;
    int n;
    apfc mean_val;
    apfc sum_sq;
    apfc diff;
    apfc sq;
    apfc count;
    apfc tmp;
    
    n = m->rows * m->cols;
    if (n <= 1) {
        apf_zero(&r->re);
        apf_zero(&r->im);
        return;
    }
    
    /* Calculate mean */
    mat_mean(&mean_val, m);
    
    /* Calculate sum of squared differences */
    apf_zero(&sum_sq.re);
    apf_zero(&sum_sq.im);
    
    for (i = 0; i < n; i++) {
        apfc_sub(&diff, &m->data[i], &mean_val);
        apfc_mul(&sq, &diff, &diff);
        apfc_add(&tmp, &sum_sq, &sq);
        sum_sq = tmp;
    }
    
    /* Divide by n-1 (sample variance) */
    apf_from_int(&count.re, n - 1);
    apf_zero(&count.im);
    apfc_div(r, &sum_sq, &count);
}

/* ========== Median ========== */

/* Simple bubble sort for small arrays - sorts by real part */
static void sort_apfc_array(apfc *arr, int n)
{
    int i;
    int j;
    apfc tmp;
    
    for (i = 0; i < n - 1; i++) {
        for (j = 0; j < n - i - 1; j++) {
            /* Compare by real part */
            if (apf_gt(&arr[j].re, &arr[j+1].re)) {
                tmp = arr[j];
                arr[j] = arr[j+1];
                arr[j+1] = tmp;
            }
        }
    }
}

void mat_median(apfc *r, const matrix_t *m)
{
    apfc sorted[MAT_MAX_ELEM];
    int i;
    int n;
    apfc sum;
    apfc two;
    
    n = m->rows * m->cols;
    if (n == 0) {
        apf_zero(&r->re);
        apf_zero(&r->im);
        return;
    }
    
    /* Copy to temp array */
    for (i = 0; i < n; i++) {
        sorted[i] = m->data[i];
    }
    
    /* Sort */
    sort_apfc_array(sorted, n);
    
    /* Find median */
    if (n % 2 == 1) {
        /* Odd count: middle element */
        *r = sorted[n / 2];
    } else {
        /* Even count: average of two middle elements */
        apfc_add(&sum, &sorted[n/2 - 1], &sorted[n/2]);
        apf_from_int(&two.re, 2);
        apf_zero(&two.im);
        apfc_div(r, &sum, &two);
    }
}

/* ========== QR Decomposition (Gram-Schmidt) ========== */

int mat_qr(matrix_t *Q, matrix_t *R, const matrix_t *A)
{
    int i;
    int j;
    int k;
    int n;
    matrix_t V;
    apfc dot;
    apfc norm_val;
    apfc tmp;
    apfc conj_q;
    apfc prod;
    apfc diff;
    apf norm_sq;
    apf zero_thresh;
    apf abs_val;
    apf sq;
    
    if (!mat_is_square(A)) {
        printf("Error: QR decomposition requires square matrix\n");
        return 0;
    }
    
    n = A->rows;
    mat_copy(&V, A);
    mat_zero(Q, n, n);
    mat_zero(R, n, n);
    
    apf_from_int(&zero_thresh, 1);
    zero_thresh.exp = -200;
    
    /* Modified Gram-Schmidt */
    for (j = 0; j < n; j++) {
        /* Copy column j of A to V */
        for (i = 0; i < n; i++) {
            MAT_AT(&V, i, j) = MAT_AT(A, i, j);
        }
        
        /* Subtract projections onto previous Q columns */
        for (k = 0; k < j; k++) {
            /* R[k,j] = Q[:,k]' * V[:,j] */
            apf_zero(&dot.re);
            apf_zero(&dot.im);
            for (i = 0; i < n; i++) {
                apfc_conj(&conj_q, &MAT_AT(Q, i, k));
                apfc_mul(&prod, &conj_q, &MAT_AT(&V, i, j));
                apfc_add(&tmp, &dot, &prod);
                dot = tmp;
            }
            MAT_AT(R, k, j) = dot;
            
            /* V[:,j] = V[:,j] - R[k,j] * Q[:,k] */
            for (i = 0; i < n; i++) {
                apfc_mul(&prod, &dot, &MAT_AT(Q, i, k));
                apfc_sub(&diff, &MAT_AT(&V, i, j), &prod);
                MAT_AT(&V, i, j) = diff;
            }
        }
        
        /* R[j,j] = ||V[:,j]|| */
        apf_zero(&norm_sq);
        for (i = 0; i < n; i++) {
            apfc_abs(&abs_val, &MAT_AT(&V, i, j));
            apf_mul(&sq, &abs_val, &abs_val);
            apf_add(&norm_sq, &norm_sq, &sq);
        }
        apf_sqrt(&norm_val.re, &norm_sq);
        apf_zero(&norm_val.im);
        MAT_AT(R, j, j) = norm_val;
        
        /* Check for zero column (rank deficient) */
        if (apf_lt(&norm_val.re, &zero_thresh)) {
            /* Set Q column to zero */
            for (i = 0; i < n; i++) {
                apf_zero(&MAT_AT(Q, i, j).re);
                apf_zero(&MAT_AT(Q, i, j).im);
            }
        } else {
            /* Q[:,j] = V[:,j] / R[j,j] */
            for (i = 0; i < n; i++) {
                apfc_div(&MAT_AT(Q, i, j), &MAT_AT(&V, i, j), &norm_val);
            }
        }
    }
    
    return 1;
}

/* ========== Eigenvalues (QR iteration) ========== */

int mat_eig(matrix_t *eigenvalues, const matrix_t *A)
{
    matrix_t T;
    matrix_t Q;
    matrix_t R;
    matrix_t tmp;
    int i;
    int j;
    int n;
    int iter;
    int max_iter;
    apf off_diag_norm;
    apf tolerance;
    apf prev_norm;
    apf abs_val;
    apf sum;
    apfc eig1;
    apfc eig2;
    
    if (!mat_is_square(A)) {
        printf("Error: eigenvalues require square matrix\n");
        return 0;
    }
    
    n = A->rows;
    max_iter = 200;
    
    /* Special case: 1x1 */
    if (n == 1) {
        mat_zero(eigenvalues, 1, 1);
        MAT_AT(eigenvalues, 0, 0) = MAT_AT(A, 0, 0);
        return 1;
    }
    
    /* Special case: 2x2 - use closed form */
    if (n == 2) {
        mat_eig2(&eig1, &eig2, A);
        mat_zero(eigenvalues, 2, 1);
        MAT_AT(eigenvalues, 0, 0) = eig1;
        MAT_AT(eigenvalues, 1, 0) = eig2;
        return 2;
    }
    
    /* QR iteration */
    mat_copy(&T, A);
    
    apf_from_int(&tolerance, 1);
    tolerance.exp = -80;  /* Convergence threshold */
    
    apf_from_int(&prev_norm, 1);
    prev_norm.exp = 100;  /* Start with large value */
    
    for (iter = 0; iter < max_iter; iter++) {
        /* QR decomposition */
        if (!mat_qr(&Q, &R, &T)) {
            return 0;
        }
        
        /* T = R * Q */
        mat_mul(&tmp, &R, &Q);
        mat_copy(&T, &tmp);
        
        /* Check convergence: sum of absolute values of off-diagonal elements */
        apf_zero(&off_diag_norm);
        for (i = 0; i < n; i++) {
            for (j = 0; j < n; j++) {
                if (i != j) {
                    apfc_abs(&abs_val, &MAT_AT(&T, i, j));
                    apf_add(&sum, &off_diag_norm, &abs_val);
                    off_diag_norm = sum;
                }
            }
        }
        
        /* Check if converged */
        if (apf_lt(&off_diag_norm, &tolerance)) {
            break;
        }
        
        /* Check if making progress */
        if (iter > 50 && !apf_lt(&off_diag_norm, &prev_norm)) {
            /* Not converging, use current approximation */
            break;
        }
        prev_norm = off_diag_norm;
    }
    
    /* Extract eigenvalues from diagonal */
    mat_zero(eigenvalues, n, 1);
    for (i = 0; i < n; i++) {
        MAT_AT(eigenvalues, i, 0) = MAT_AT(&T, i, i);
    }
    
    return n;
}
