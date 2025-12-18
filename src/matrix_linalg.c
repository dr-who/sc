/* matrix_linalg.c - Linear algebra operations for scalc
 * C89 compliant for Watcom C / DOS
 * 16-bit clean: int is 16-bit, long is 32-bit
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "matrix.h"
#include "apfx.h"
#include "apf_native.h"

#ifndef NAN
#define NAN (0.0/0.0)
#endif

#ifndef INFINITY
#define INFINITY HUGE_VAL
#endif

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
    
    /* Create very small threshold: 1e-30 */
    {
        apf ten;
        int count;
        apf_from_int(&zero_thresh, 1);
        apf_from_int(&ten, 10);
        for (count = 0; count < 30; count++) {
            apf temp;
            apf_div(&temp, &zero_thresh, &ten);
            apf_copy(&zero_thresh, &temp);
        }
    }
    
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
    matrix_t L = {0, 0, NULL}, U = {0, 0, NULL};
    int perm[MAT_MAX_ROWS];
    int i, j, n;
    matrix_t y = {0, 0, NULL}, pb = {0, 0, NULL};
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
    matrix_t B_inv = {0, 0, NULL};
    
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
    matrix_t ata = {0, 0, NULL};
    apfc eig1, eig2;
    
    if (m->rows != 2 || m->cols != 2) {
        printf("Error: svd only works for 2x2 matrices\n");
        return 0;
    }
    
    mat_conj_transpose(&ata, m);
    {
        matrix_t tmp = {0, 0, NULL};
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
    matrix_t V = {0, 0, NULL};
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
    if (!V.data) return 0;  /* Allocation failed */
    mat_zero(Q, n, n);
    if (!Q->data) return 0;  /* Allocation failed */
    mat_zero(R, n, n);
    if (!R->data) return 0;  /* Allocation failed */
    
    /* Create very small threshold: 1e-30 */
    {
        apf ten;
        int count;
        apf_from_int(&zero_thresh, 1);
        apf_from_int(&ten, 10);
        for (count = 0; count < 30; count++) {
            apf temp;
            apf_div(&temp, &zero_thresh, &ten);
            apf_copy(&zero_thresh, &temp);
        }
    }
    
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
    matrix_t T = {0, 0, NULL};
    matrix_t Q = {0, 0, NULL};
    matrix_t R = {0, 0, NULL};
    matrix_t tmp = {0, 0, NULL};
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
    
    /* Create small tolerance: 1e-30 */
    {
        apf ten;
        int cnt;
        apf_from_int(&tolerance, 1);
        apf_from_int(&ten, 10);
        for (cnt = 0; cnt < 30; cnt++) {
            apf tm;
            apf_div(&tm, &tolerance, &ten);
            apf_copy(&tolerance, &tm);
        }
    }
    
    /* Create large prev_norm: 1e100 */
    {
        apf ten;
        int cnt;
        apf_from_int(&prev_norm, 1);
        apf_from_int(&ten, 10);
        for (cnt = 0; cnt < 100; cnt++) {
            apf tm;
            apf_mul(&tm, &prev_norm, &ten);
            apf_copy(&prev_norm, &tm);
        }
    }
    
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

/* ========== Cholesky Decomposition ========== */
/* For symmetric positive-definite matrices: A = L * L^T */

int mat_chol(matrix_t *L, const matrix_t *A)
{
    int i, j, k, n;
    apfc sum, tmp, diag;
    apf zero_val, test_val;
    
    if (!mat_is_square(A)) {
        printf("Error: Cholesky requires square matrix\n");
        return 0;
    }
    
    n = A->rows;
    mat_zero(L, n, n);
    
    apf_zero(&zero_val);
    
    for (i = 0; i < n; i++) {
        for (j = 0; j <= i; j++) {
            apfc_zero(&sum);
            
            if (j == i) {
                /* Diagonal element: L[i][i] = sqrt(A[i][i] - sum(L[i][k]^2)) */
                for (k = 0; k < j; k++) {
                    apfc sq;
                    apfc_mul(&sq, &MAT_AT(L, i, k), &MAT_AT(L, i, k));
                    apfc_add(&sum, &sum, &sq);
                }
                apfc_sub(&diag, &MAT_AT(A, i, i), &sum);
                
                /* Check if positive (matrix must be positive definite) */
                apfc_abs(&test_val, &diag);
                if (diag.re.sign || apf_lt(&test_val, &zero_val)) {
                    printf("Error: matrix is not positive definite\n");
                    return 0;
                }
                
                apfc_sqrt(&MAT_AT(L, i, j), &diag);
            } else {
                /* Off-diagonal: L[i][j] = (A[i][j] - sum(L[i][k]*L[j][k])) / L[j][j] */
                for (k = 0; k < j; k++) {
                    apfc prod;
                    apfc_mul(&prod, &MAT_AT(L, i, k), &MAT_AT(L, j, k));
                    apfc_add(&sum, &sum, &prod);
                }
                apfc_sub(&tmp, &MAT_AT(A, i, j), &sum);
                apfc_div(&MAT_AT(L, i, j), &tmp, &MAT_AT(L, j, j));
            }
        }
    }
    
    return 1;
}

/* ========== Singular Value Decomposition (SVD) ========== */
/* A = U * S * V^T where U, V are orthogonal, S is diagonal */
/* Uses iterative method for general matrices */

int mat_svd(matrix_t *U, matrix_t *S, matrix_t *V, const matrix_t *A)
{
    int m, n, min_mn;
    int i, j, iter;
    int max_iter;
    matrix_t AtA = {0, 0, NULL};
    matrix_t Q = {0, 0, NULL}, R = {0, 0, NULL};
    matrix_t tmp = {0, 0, NULL}, tmp2 = {0, 0, NULL};
    matrix_t At = {0, 0, NULL};
    apf tolerance, off_norm, sum, abs_val;
    int converged = 0;
    static int svd_debug = -1;  /* -1 = unset, 0 = off, 1 = on */
    
    /* Check for debug flag on first call */
    if (svd_debug < 0) {
        svd_debug = (getenv("SVD_DEBUG")) ? 1 : 0;
    }
    
    m = A->rows;
    n = A->cols;
    min_mn = (m < n) ? m : n;
    max_iter = 100;  /* Increased iterations */
    
    if (svd_debug) printf("SVD: input %dx%d\n", m, n);
    
    /* Validate input */
    if (!A || !A->data || m <= 0 || n <= 0) {
        if (svd_debug) printf("SVD: invalid input\n");
        return 0;
    }
    
    /* Special case: 2x2 - compute full SVD analytically */
    if (m == 2 && n == 2) {
        apfc a, b, c, d;
        apf aa, bb, cc, dd, ab, cd, ac, bd;
        apf sum1, sum2, diff1, diff2, prod1, prod2;
        apf s1_sq, s2_sq, s1, s2;
        apf theta, ct, st;
        apf tmp_r, two;
        
        /* A = [[a,b],[c,d]] */
        a = MAT_AT(A, 0, 0);
        b = MAT_AT(A, 0, 1);
        c = MAT_AT(A, 1, 0);
        d = MAT_AT(A, 1, 1);
        
        /* For real matrices, use explicit formulas */
        /* A^T*A = [[a²+c², ab+cd], [ab+cd, b²+d²]] */
        apf_mul(&aa, &a.re, &a.re);
        apf_mul(&bb, &b.re, &b.re);
        apf_mul(&cc, &c.re, &c.re);
        apf_mul(&dd, &d.re, &d.re);
        apf_mul(&ab, &a.re, &b.re);
        apf_mul(&cd, &c.re, &d.re);
        apf_mul(&ac, &a.re, &c.re);
        apf_mul(&bd, &b.re, &d.re);
        
        /* sum1 = a²+c², sum2 = b²+d² */
        apf_add(&sum1, &aa, &cc);
        apf_add(&sum2, &bb, &dd);
        
        /* prod1 = ab+cd */
        apf_add(&prod1, &ab, &cd);
        
        /* Eigenvalues of A^T*A: ((sum1+sum2) ± sqrt((sum1-sum2)²+4*prod1²))/2 */
        apf_add(&tmp_r, &sum1, &sum2);  /* sum1+sum2 */
        apf_sub(&diff1, &sum1, &sum2);  /* sum1-sum2 */
        apf_mul(&diff2, &diff1, &diff1); /* (sum1-sum2)² */
        apf_mul(&prod2, &prod1, &prod1); /* prod1² */
        apf_from_int(&two, 4);
        apf_mul(&prod2, &prod2, &two);   /* 4*prod1² */
        apf_add(&diff2, &diff2, &prod2); /* (sum1-sum2)²+4*prod1² */
        apf_sqrt(&diff2, &diff2);        /* sqrt(...) */
        
        apf_from_int(&two, 2);
        apf_add(&s1_sq, &tmp_r, &diff2);
        apf_div(&s1_sq, &s1_sq, &two);   /* larger eigenvalue */
        apf_sub(&s2_sq, &tmp_r, &diff2);
        apf_div(&s2_sq, &s2_sq, &two);   /* smaller eigenvalue */
        
        /* Singular values */
        if (apf_is_neg(&s1_sq)) apf_zero(&s1_sq);
        if (apf_is_neg(&s2_sq)) apf_zero(&s2_sq);
        apf_sqrt(&s1, &s1_sq);
        apf_sqrt(&s2, &s2_sq);
        
        /* Set S */
        mat_zero(S, 2, 2);
        MAT_AT(S, 0, 0).re = s1;
        MAT_AT(S, 1, 1).re = s2;
        
        /* Compute angle for V (eigenvectors of A^T*A) */
        /* theta = 0.5 * atan2(2*prod1, sum1-sum2) */
        apf_from_int(&two, 2);
        apf_mul(&tmp_r, &prod1, &two);  /* 2*prod1 */
        apfx_atan2(&theta, &tmp_r, &diff1);
        apf_div(&theta, &theta, &two);
        
        apfx_cos(&ct, &theta);
        apfx_sin(&st, &theta);
        
        /* V = [[cos,-sin],[sin,cos]] */
        mat_zero(V, 2, 2);
        MAT_AT(V, 0, 0).re = ct;
        apf_neg(&MAT_AT(V, 0, 1).re, &st);
        MAT_AT(V, 1, 0).re = st;
        MAT_AT(V, 1, 1).re = ct;
        
        /* Compute U = A*V*S^(-1) */
        /* u1 = (1/s1) * A * v1, u2 = (1/s2) * A * v2 */
        mat_zero(U, 2, 2);
        
        if (!apf_iszero(&s1)) {
            /* U[:,0] = (A * V[:,0]) / s1 */
            apf u00, u10, tmp1, tmp2r;
            apf_mul(&tmp1, &a.re, &ct);
            apf_mul(&tmp2r, &b.re, &st);
            apf_add(&u00, &tmp1, &tmp2r);
            apf_div(&MAT_AT(U, 0, 0).re, &u00, &s1);
            
            apf_mul(&tmp1, &c.re, &ct);
            apf_mul(&tmp2r, &d.re, &st);
            apf_add(&u10, &tmp1, &tmp2r);
            apf_div(&MAT_AT(U, 1, 0).re, &u10, &s1);
        }
        
        if (!apf_iszero(&s2)) {
            /* U[:,1] = (A * V[:,1]) / s2 */
            apf u01, u11, tmp1, tmp2r, neg_st;
            apf_neg(&neg_st, &st);
            
            apf_mul(&tmp1, &a.re, &neg_st);
            apf_mul(&tmp2r, &b.re, &ct);
            apf_add(&u01, &tmp1, &tmp2r);
            apf_div(&MAT_AT(U, 0, 1).re, &u01, &s2);
            
            apf_mul(&tmp1, &c.re, &neg_st);
            apf_mul(&tmp2r, &d.re, &ct);
            apf_add(&u11, &tmp1, &tmp2r);
            apf_div(&MAT_AT(U, 1, 1).re, &u11, &s2);
        } else {
            /* If s2=0, set U[:,1] to orthogonal complement of U[:,0] */
            MAT_AT(U, 0, 1).re = MAT_AT(U, 1, 0).re;
            apf_neg(&MAT_AT(U, 1, 1).re, &MAT_AT(U, 0, 0).re);
        }
        
        return 1;
    }
    
    /* Compute A^T * A for eigenvalues -> singular values */
    mat_transpose(&At, A);
    if (!At.data) {
        if (svd_debug) printf("SVD: transpose failed\n");
        return 0;  /* Allocation failed */
    }
    if (svd_debug) printf("SVD: transpose ok %dx%d\n", At.rows, At.cols);
    
    mat_mul(&AtA, &At, A);  /* n x n */
    if (!AtA.data) {
        if (svd_debug) printf("SVD: multiply failed\n");
        return 0;  /* Allocation failed */
    }
    if (svd_debug) printf("SVD: AtA ok %dx%d\n", AtA.rows, AtA.cols);
    
    /* Check if AtA is already diagonal (e.g., for identity or orthogonal matrices) */
    {
        int is_diagonal = 1;
        apf diag_thresh;
        apf ten;
        int cnt;
        
        /* Create small threshold: 1e-20 */
        apf_from_int(&diag_thresh, 1);
        apf_from_int(&ten, 10);
        for (cnt = 0; cnt < 20; cnt++) {
            apf tm;
            apf_div(&tm, &diag_thresh, &ten);
            apf_copy(&diag_thresh, &tm);
        }
        
        for (i = 0; i < n && is_diagonal; i++) {
            for (j = 0; j < n && is_diagonal; j++) {
                if (i != j) {
                    apf abs_off;
                    apfc_abs(&abs_off, &MAT_AT(&AtA, i, j));
                    if (apf_gt(&abs_off, &diag_thresh)) {
                        is_diagonal = 0;
                    }
                }
            }
        }
        
        if (svd_debug) printf("SVD: is_diagonal=%d\n", is_diagonal);
        
        if (is_diagonal) {
            /* AtA is diagonal - SVD is trivial */
            /* V = I, S[i,i] = sqrt(AtA[i,i]), U = A * V * S^(-1) */
            mat_identity(V, n);
            if (!V->data) {
                if (svd_debug) printf("SVD: V identity alloc failed\n");
                return 0;
            }
            
            mat_zero(S, m, n);
            if (!S->data) {
                if (svd_debug) printf("SVD: S alloc failed\n");
                return 0;
            }
            for (i = 0; i < min_mn; i++) {
                apf sv;
                apfc_abs(&sv, &MAT_AT(&AtA, i, i));
                apf_sqrt(&MAT_AT(S, i, i).re, &sv);
                apf_zero(&MAT_AT(S, i, i).im);
            }
            
            /* U = A * V * S^(-1) = A * S^(-1) for diagonal S */
            mat_zero(U, m, m);
            if (!U->data) {
                if (svd_debug) printf("SVD: U alloc failed\n");
                return 0;
            }
            for (i = 0; i < m; i++) {
                for (j = 0; j < min_mn; j++) {
                    if (!apf_iszero(&MAT_AT(S, j, j).re)) {
                        apfc_div(&MAT_AT(U, i, j), &MAT_AT(A, i, j), &MAT_AT(S, j, j));
                    }
                }
            }
            /* Set remaining diagonal of U to 1 */
            for (i = min_mn; i < m; i++) {
                apf_from_int(&MAT_AT(U, i, i).re, 1);
                apf_zero(&MAT_AT(U, i, i).im);
            }
            
            if (svd_debug) printf("SVD: diagonal path success\n");
            return 1;
        }
    }
    
    if (svd_debug) printf("SVD: using general QR path\n");
    
    /* Initialize V with identity */
    mat_identity(V, n);
    if (!V->data) {
        if (svd_debug) printf("SVD: V identity alloc failed (general)\n");
        return 0;  /* Allocation failed */
    }
    
    /* Create small tolerance: 1e-30 */
    {
        apf ten;
        int cnt;
        apf_from_int(&tolerance, 1);
        apf_from_int(&ten, 10);
        for (cnt = 0; cnt < 30; cnt++) {
            apf tm;
            apf_div(&tm, &tolerance, &ten);
            apf_copy(&tolerance, &tm);
        }
    }
    
    /* Power iteration on A^T * A to get V */
    mat_copy(&tmp, &AtA);
    if (!tmp.data) return 0;  /* Allocation failed */
    
    for (iter = 0; iter < max_iter; iter++) {
        /* Check convergence FIRST before allocating more memory */
        apf_zero(&off_norm);
        for (i = 0; i < n; i++) {
            for (j = 0; j < n; j++) {
                if (i != j) {
                    apfc_abs(&abs_val, &MAT_AT(&tmp, i, j));
                    apf_add(&sum, &off_norm, &abs_val);
                    off_norm = sum;
                }
            }
        }
        
        /* Check if converged (off-diagonal is small) */
        if (apf_lt(&off_norm, &tolerance) || apf_iszero(&off_norm)) {
            converged = 1;
            break;
        }
        
        /* Not converged - do QR iteration */
        if (!mat_qr(&Q, &R, &tmp)) {
            /* QR failed - try to use current state */
            break;
        }
        if (!Q.data || !R.data) break;  /* Allocation failed */
        
        /* Update V = V * Q */
        mat_mul(&tmp2, V, &Q);
        if (!tmp2.data) break;  /* Allocation failed */
        mat_copy(V, &tmp2);
        if (!V->data) break;  /* Allocation failed */
        
        /* tmp = R * Q (similar transformation) */
        mat_mul(&tmp2, &R, &Q);
        if (!tmp2.data) break;  /* Allocation failed */
        mat_copy(&tmp, &tmp2);
        if (!tmp.data) break;  /* Allocation failed */
    }
    
    /* If we didn't converge via QR, try to use AtA directly for diagonal dominant cases */
    if (!converged && iter >= max_iter) {
        /* Revert to AtA - might work for nearly diagonal matrices */
        mat_copy(&tmp, &AtA);
    }
    
    /* Singular values are square roots of diagonal of converged matrix */
    mat_zero(S, m, n);
    if (!S->data) return 0;  /* Allocation failed */
    for (i = 0; i < min_mn; i++) {
        apfc sv;
        apfc_abs(&sv.re, &MAT_AT(&tmp, i, i));
        apf_sqrt(&sv.re, &sv.re);
        apf_zero(&sv.im);
        MAT_AT(S, i, i) = sv;
    }
    
    /* Compute U = A * V * S^(-1) */
    mat_mul(&tmp, A, V);
    if (!tmp.data) return 0;  /* Allocation failed */
    mat_zero(U, m, m);
    if (!U->data) return 0;  /* Allocation failed */
    for (i = 0; i < m; i++) {
        for (j = 0; j < min_mn; j++) {
            apf sv_val;
            apfc_abs(&sv_val, &MAT_AT(S, j, j));
            if (!apf_iszero(&sv_val)) {
                apfc_div(&MAT_AT(U, i, j), &MAT_AT(&tmp, i, j), &MAT_AT(S, j, j));
            }
        }
    }
    
    /* Complete U to orthonormal basis (simplified - just set remaining to identity) */
    for (i = min_mn; i < m; i++) {
        apf_from_int(&MAT_AT(U, i, i).re, 1);
        apf_zero(&MAT_AT(U, i, i).im);
    }
    
    return 1;
}

/* ========== Null Space ========== */
/* Returns basis vectors for null space of A */

int mat_null(matrix_t *N, const matrix_t *A)
{
    matrix_t R = {0, 0, NULL};
    int m, n, i, j, k;
    int rank_val, null_dim;
    int pivot_cols[MAT_MAX_COLS];
    int free_cols[MAT_MAX_COLS];
    int num_pivots, num_free;
    int pivot_row;
    apf thresh, abs_pivot, abs_test;
    apfc factor, prod;
    
    m = A->rows;
    n = A->cols;
    
    /* Copy A to R for row reduction */
    mat_copy(&R, A);
    
    /* Very small threshold for zero detection: 1e-30 */
    {
        apf ten;
        int cnt;
        apf_from_int(&thresh, 1);
        apf_from_int(&ten, 10);
        for (cnt = 0; cnt < 30; cnt++) {
            apf tm;
            apf_div(&tm, &thresh, &ten);
            apf_copy(&thresh, &tm);
        }
    }
    
    /* Perform Gaussian elimination to get row echelon form */
    pivot_row = 0;
    num_pivots = 0;
    
    for (j = 0; j < n && pivot_row < m; j++) {
        /* Find pivot in column j */
        int best_row = -1;
        apf_zero(&abs_pivot);
        
        for (i = pivot_row; i < m; i++) {
            apfc_abs(&abs_test, &MAT_AT(&R, i, j));
            if (apf_gt(&abs_test, &abs_pivot)) {
                abs_pivot = abs_test;
                best_row = i;
            }
        }
        
        if (best_row < 0 || apf_lt(&abs_pivot, &thresh)) {
            /* No pivot in this column - it's a free variable */
            continue;
        }
        
        /* Swap rows if needed */
        if (best_row != pivot_row) {
            for (k = 0; k < n; k++) {
                apfc tmp = MAT_AT(&R, pivot_row, k);
                MAT_AT(&R, pivot_row, k) = MAT_AT(&R, best_row, k);
                MAT_AT(&R, best_row, k) = tmp;
            }
        }
        
        /* Scale pivot row to make pivot = 1 */
        {
            apfc pivot_val = MAT_AT(&R, pivot_row, j);
            for (k = j; k < n; k++) {
                apfc_div(&MAT_AT(&R, pivot_row, k), &MAT_AT(&R, pivot_row, k), &pivot_val);
            }
        }
        
        /* Eliminate other rows */
        for (i = 0; i < m; i++) {
            if (i != pivot_row) {
                factor = MAT_AT(&R, i, j);
                for (k = j; k < n; k++) {
                    apfc_mul(&prod, &factor, &MAT_AT(&R, pivot_row, k));
                    apfc_sub(&MAT_AT(&R, i, k), &MAT_AT(&R, i, k), &prod);
                }
            }
        }
        
        pivot_cols[num_pivots++] = j;
        pivot_row++;
    }
    
    /* Identify free columns */
    num_free = 0;
    {
        int p = 0;
        for (j = 0; j < n; j++) {
            if (p < num_pivots && pivot_cols[p] == j) {
                p++;
            } else {
                free_cols[num_free++] = j;
            }
        }
    }
    
    rank_val = num_pivots;
    null_dim = n - rank_val;
    
    if (null_dim == 0) {
        /* Null space is trivial (only zero vector) */
        mat_zero(N, n, 1);
        return 0;
    }
    
    /* Build null space basis: one vector per free variable */
    mat_zero(N, n, null_dim);
    
    for (k = 0; k < null_dim; k++) {
        int fc = free_cols[k];
        
        /* Set free variable to 1 */
        apf_from_int(&MAT_AT(N, fc, k).re, 1);
        apf_zero(&MAT_AT(N, fc, k).im);
        
        /* Set pivot variables from RREF */
        for (i = 0; i < num_pivots; i++) {
            int pc = pivot_cols[i];
            /* N[pc, k] = -R[i, fc] */
            apfc_neg(&MAT_AT(N, pc, k), &MAT_AT(&R, i, fc));
        }
    }
    
    return null_dim;
}

/* ========== Schur Decomposition ========== */
/* A = Q * T * Q^T where Q is orthogonal, T is upper triangular */

int mat_schur(matrix_t *Q, matrix_t *T, const matrix_t *A)
{
    int n, i, iter, max_iter;
    matrix_t Qk = {0,0,NULL}, R = {0,0,NULL};
    matrix_t tmp = {0, 0, NULL};
    apf off_norm, sum, abs_val, tolerance;
    
    if (!mat_is_square(A)) {
        printf("Error: Schur decomposition requires square matrix\n");
        return 0;
    }
    
    n = A->rows;
    max_iter = 300;
    
    mat_copy(T, A);
    mat_identity(Q, n);
    
    /* Create small tolerance: 1e-30 */
    {
        apf ten;
        int cnt;
        apf_from_int(&tolerance, 1);
        apf_from_int(&ten, 10);
        for (cnt = 0; cnt < 30; cnt++) {
            apf tm;
            apf_div(&tm, &tolerance, &ten);
            apf_copy(&tolerance, &tm);
        }
    }
    
    /* QR iteration */
    for (iter = 0; iter < max_iter; iter++) {
        if (!mat_qr(&Qk, &R, T)) break;
        
        /* T = R * Qk */
        mat_mul(&tmp, &R, &Qk);
        mat_copy(T, &tmp);
        
        /* Q = Q * Qk */
        mat_mul(&tmp, Q, &Qk);
        mat_copy(Q, &tmp);
        
        /* Check convergence (sub-diagonal elements) */
        apf_zero(&off_norm);
        for (i = 1; i < n; i++) {
            apfc_abs(&abs_val, &MAT_AT(T, i, i-1));
            apf_add(&sum, &off_norm, &abs_val);
            off_norm = sum;
        }
        
        if (apf_lt(&off_norm, &tolerance)) break;
    }
    
    return 1;
}

/* ========== PCA (Principal Component Analysis) ========== */

int mat_pca(matrix_t *coeff, matrix_t *score, matrix_t *latent, const matrix_t *X)
{
    int m, n, i, j;
    matrix_t centered = {0, 0, NULL};
    apfc col_mean, diff;
    apf sum_r;
    
    m = X->rows;
    n = X->cols;
    
    if (m < 2) {
        printf("Error: PCA requires at least 2 observations\n");
        return 0;
    }
    
    /* Center the data (subtract column means) */
    mat_zero(&centered, m, n);
    for (j = 0; j < n; j++) {
        apf_zero(&sum_r);
        for (i = 0; i < m; i++) {
            apf_add(&sum_r, &sum_r, &MAT_AT(X, i, j).re);
        }
        apf_from_int(&col_mean.re, m);
        apf_div(&col_mean.re, &sum_r, &col_mean.re);
        apf_zero(&col_mean.im);
        
        for (i = 0; i < m; i++) {
            apfc_sub(&diff, &MAT_AT(X, i, j), &col_mean);
            MAT_AT(&centered, i, j) = diff;
        }
    }
    
    /* SVD on centered data: X = U*S*V' */
    /* coeff = V, score = U*S, latent = diag(S)^2/(m-1) */
    {
        matrix_t U = {0,0,NULL}, S = {0,0,NULL}, V = {0,0,NULL};
        int k;
        
        if (!mat_svd(&U, &S, &V, &centered)) {
            printf("Error: PCA SVD failed\n");
            return 0;
        }
        
        mat_copy(coeff, &V);
        
        /* Apply MATLAB sign convention: largest magnitude element in each column is positive */
        for (k = 0; k < coeff->cols; k++) {
            int max_row = 0;
            apf max_abs, abs_val;
            apf_zero(&max_abs);
            
            /* Find row with largest magnitude in this column */
            for (i = 0; i < coeff->rows; i++) {
                apfc_abs(&abs_val, &MAT_AT(coeff, i, k));
                if (apf_gt(&abs_val, &max_abs)) {
                    max_abs = abs_val;
                    max_row = i;
                }
            }
            
            /* If that element is negative, flip the entire column */
            if (apf_is_neg(&MAT_AT(coeff, max_row, k).re)) {
                for (i = 0; i < coeff->rows; i++) {
                    apf_neg(&MAT_AT(coeff, i, k).re, &MAT_AT(coeff, i, k).re);
                    apf_neg(&MAT_AT(coeff, i, k).im, &MAT_AT(coeff, i, k).im);
                }
            }
        }
        
        /* latent = singular values squared / (m-1) = variances */
        mat_zero(latent, n, 1);
        for (i = 0; i < n && i < m; i++) {
            apf var_i, m_minus_1;
            apf_mul(&var_i, &MAT_AT(&S, i, i).re, &MAT_AT(&S, i, i).re);
            apf_from_int(&m_minus_1, m - 1);
            apf_div(&MAT_AT(latent, i, 0).re, &var_i, &m_minus_1);
        }
        
        /* score = centered * coeff */
        mat_mul(score, &centered, coeff);
    }
    
    return 1;
}

int mat_pca_reduce(matrix_t *result, const matrix_t *X, int k)
{
    matrix_t coeff = {0,0,NULL}, score = {0,0,NULL}, latent = {0,0,NULL};
    int i, j, n;
    
    if (k < 1) k = 1;
    
    if (!mat_pca(&coeff, &score, &latent, X)) {
        return 0;
    }
    
    n = score.rows;
    if (k > score.cols) k = score.cols;
    
    mat_zero(result, n, k);
    for (i = 0; i < n; i++) {
        for (j = 0; j < k; j++) {
            MAT_AT(result, i, j) = MAT_AT(&score, i, j);
        }
    }
    
    return 1;
}

/* ========== K-Means Clustering ========== */

int mat_kmeans(matrix_t *idx, matrix_t *centroids, const matrix_t *X, int k)
{
    int m, n, i, j, c, iter;
    int max_iter = 100;
    int changed;
    int counts[MAT_MAX_ROWS];
    apf min_dist, dist, diff_r, diff_sq, sum_sq;
    int best_cluster;
    
    m = X->rows;
    n = X->cols;
    
    if (k < 1) k = 1;
    if (k > m) k = m;
    
    mat_zero(idx, m, 1);
    mat_zero(centroids, k, n);
    
    /* Initialize centroids spread across data */
    {
        int step = (m > k) ? m / k : 1;
        for (i = 0; i < k; i++) {
            for (j = 0; j < n; j++) {
                MAT_AT(centroids, i, j) = MAT_AT(X, (i * step) % m, j);
            }
        }
    }
    
    for (iter = 0; iter < max_iter; iter++) {
        changed = 0;
        
        /* Assignment step */
        for (i = 0; i < m; i++) {
            apf_from_int(&min_dist, 1);
            min_dist.exp += 100;
            best_cluster = 0;
            
            for (c = 0; c < k; c++) {
                apf_zero(&sum_sq);
                for (j = 0; j < n; j++) {
                    apf_sub(&diff_r, &MAT_AT(X, i, j).re, &MAT_AT(centroids, c, j).re);
                    apf_mul(&diff_sq, &diff_r, &diff_r);
                    apf_add(&sum_sq, &sum_sq, &diff_sq);
                }
                dist = sum_sq;
                
                if (apf_lt(&dist, &min_dist)) {
                    min_dist = dist;
                    best_cluster = c;
                }
            }
            
            {
                int old_cluster = apf_to_long(&MAT_AT(idx, i, 0).re);
                if (old_cluster != best_cluster + 1) changed = 1;
                apf_from_int(&MAT_AT(idx, i, 0).re, best_cluster + 1);  /* 1-indexed */
                apf_zero(&MAT_AT(idx, i, 0).im);
            }
        }
        
        /* Update step */
        for (c = 0; c < k; c++) {
            counts[c] = 0;
            for (j = 0; j < n; j++) {
                apf_zero(&MAT_AT(centroids, c, j).re);
                apf_zero(&MAT_AT(centroids, c, j).im);
            }
        }
        
        for (i = 0; i < m; i++) {
            c = apf_to_long(&MAT_AT(idx, i, 0).re) - 1;  /* Convert back to 0-indexed */
            if (c >= 0 && c < k) {
                counts[c]++;
                for (j = 0; j < n; j++) {
                    apf_add(&MAT_AT(centroids, c, j).re, 
                            &MAT_AT(centroids, c, j).re, 
                            &MAT_AT(X, i, j).re);
                }
            }
        }
        
        for (c = 0; c < k; c++) {
            if (counts[c] > 0) {
                apf cnt;
                apf_from_int(&cnt, counts[c]);
                for (j = 0; j < n; j++) {
                    apf_div(&MAT_AT(centroids, c, j).re, 
                            &MAT_AT(centroids, c, j).re, &cnt);
                }
            }
        }
        
        if (!changed) break;
    }
    
    return 1;
}

/* ========== Pairwise Distance ========== */

int mat_pdist(matrix_t *D, const matrix_t *X)
{
    int m, n, i, j, d, pos;
    int num_pairs;
    apf diff_r, diff_sq, sum_sq, dist;
    
    m = X->rows;
    n = X->cols;
    num_pairs = (m * (m - 1)) / 2;
    
    mat_zero(D, 1, num_pairs);
    
    pos = 0;
    for (i = 0; i < m - 1; i++) {
        for (j = i + 1; j < m; j++) {
            apf_zero(&sum_sq);
            for (d = 0; d < n; d++) {
                apf_sub(&diff_r, &MAT_AT(X, i, d).re, &MAT_AT(X, j, d).re);
                apf_mul(&diff_sq, &diff_r, &diff_r);
                apf_add(&sum_sq, &sum_sq, &diff_sq);
            }
            apf_sqrt(&dist, &sum_sq);
            MAT_AT(D, 0, pos).re = dist;
            pos++;
        }
    }
    
    return 1;
}

/* ========== Silhouette Score ========== */

void mat_silhouette(apfc *score, const matrix_t *X, const matrix_t *idx)
{
    int m, n, i, j, d, ci, cj;
    apf sum_sil, sil_i, a_i, b_i;
    apf dist, diff_r, diff_sq, sum_sq;
    apf min_b, tmp, max_ab;
    apf cluster_dist[MAT_MAX_ROWS];
    int cluster_count[MAT_MAX_ROWS];
    int k_clusters;
    
    m = X->rows;
    n = X->cols;
    
    k_clusters = 0;
    for (i = 0; i < m; i++) {
        int c = apf_to_long(&MAT_AT(idx, i, 0).re);
        if (c >= k_clusters) k_clusters = c + 1;
    }
    
    apf_zero(&sum_sil);
    
    for (i = 0; i < m; i++) {
        ci = apf_to_long(&MAT_AT(idx, i, 0).re);
        
        for (j = 0; j < k_clusters; j++) {
            apf_zero(&cluster_dist[j]);
            cluster_count[j] = 0;
        }
        
        for (j = 0; j < m; j++) {
            if (i == j) continue;
            cj = apf_to_long(&MAT_AT(idx, j, 0).re);
            
            apf_zero(&sum_sq);
            for (d = 0; d < n; d++) {
                apf_sub(&diff_r, &MAT_AT(X, i, d).re, &MAT_AT(X, j, d).re);
                apf_mul(&diff_sq, &diff_r, &diff_r);
                apf_add(&sum_sq, &sum_sq, &diff_sq);
            }
            apf_sqrt(&dist, &sum_sq);
            
            apf_add(&cluster_dist[cj], &cluster_dist[cj], &dist);
            cluster_count[cj]++;
        }
        
        if (cluster_count[ci] > 0) {
            apf cnt;
            apf_from_int(&cnt, cluster_count[ci]);
            apf_div(&a_i, &cluster_dist[ci], &cnt);
        } else {
            apf_zero(&a_i);
        }
        
        apf_from_int(&min_b, 1);
        min_b.exp += 100;
        for (j = 0; j < k_clusters; j++) {
            if (j == ci || cluster_count[j] == 0) continue;
            apf_from_int(&tmp, cluster_count[j]);
            apf_div(&tmp, &cluster_dist[j], &tmp);
            if (apf_lt(&tmp, &min_b)) min_b = tmp;
        }
        b_i = min_b;
        
        apf_sub(&sil_i, &b_i, &a_i);
        max_ab = apf_gt(&a_i, &b_i) ? a_i : b_i;
        if (!apf_iszero(&max_ab)) {
            apf_div(&sil_i, &sil_i, &max_ab);
        }
        
        apf_add(&sum_sil, &sum_sil, &sil_i);
    }
    
    if (m > 0) {
        apf_from_int(&tmp, m);
        apf_div(&score->re, &sum_sil, &tmp);
    } else {
        apf_zero(&score->re);
    }
    apf_zero(&score->im);
}


/* ========== Column-wise statistics (MATLAB compatible) ========== */

/* Column means: returns 1xN row vector - pure APF */
void mat_mean_cols(matrix_t *r, const matrix_t *m)
{
    int i, j;
    apf sum, divisor;
    
    mat_zero(r, 1, m->cols);
    if (!r->data) return;
    
    apf_from_int(&divisor, m->rows);
    
    for (j = 0; j < m->cols; j++) {
        apf_zero(&sum);
        for (i = 0; i < m->rows; i++) {
            apf_add(&sum, &sum, &MAT_AT(m, i, j).re);
        }
        apf_div(&MAT_AT(r, 0, j).re, &sum, &divisor);
        apf_zero(&MAT_AT(r, 0, j).im);
    }
}

/* Column std devs: returns 1xN row vector - pure APF */
void mat_std_cols(matrix_t *r, const matrix_t *m)
{
    int i, j;
    apf sum, sum2, mean, var, v_sq, divisor;
    
    mat_zero(r, 1, m->cols);
    if (!r->data) return;
    
    apf_from_int(&divisor, m->rows);
    
    for (j = 0; j < m->cols; j++) {
        apf_zero(&sum);
        apf_zero(&sum2);
        for (i = 0; i < m->rows; i++) {
            apf_add(&sum, &sum, &MAT_AT(m, i, j).re);
            apf_mul(&v_sq, &MAT_AT(m, i, j).re, &MAT_AT(m, i, j).re);
            apf_add(&sum2, &sum2, &v_sq);
        }
        /* mean = sum / n */
        apf_div(&mean, &sum, &divisor);
        /* var = sum2/n - mean^2 */
        apf_div(&var, &sum2, &divisor);
        apf_mul(&v_sq, &mean, &mean);
        apf_sub(&var, &var, &v_sq);
        /* std = sqrt(max(var, 0)) */
        if (apf_cmp_int(&var, 0) > 0) {
            apf_sqrt(&MAT_AT(r, 0, j).re, &var);
        } else {
            apf_zero(&MAT_AT(r, 0, j).re);
        }
        apf_zero(&MAT_AT(r, 0, j).im);
    }
}

/* Covariance matrix: returns NxN matrix for MxN input - pure APF */
void mat_cov(matrix_t *r, const matrix_t *m)
{
    int i, j, k;
    apf *means;
    apf sum, divisor, diff_i, diff_j, prod, cov_ij;
    
    means = (apf *)mat_arena_alloc(m->cols * sizeof(apf));
    if (!means) return;
    
    /* Initialize and compute column means */
    apf_from_int(&divisor, m->rows);
    for (j = 0; j < m->cols; j++) {
        apf_zero(&sum);
        for (i = 0; i < m->rows; i++) {
            apf_add(&sum, &sum, &MAT_AT(m, i, j).re);
        }
        apf_div(&means[j], &sum, &divisor);
    }
    
    /* Compute covariance matrix */
    mat_zero(r, m->cols, m->cols);
    if (!r->data) return;
    
    apf_from_int(&divisor, (m->rows > 1) ? (m->rows - 1) : 1);
    
    for (i = 0; i < m->cols; i++) {
        for (j = 0; j < m->cols; j++) {
            apf_zero(&cov_ij);
            for (k = 0; k < m->rows; k++) {
                apf_sub(&diff_i, &MAT_AT(m, k, i).re, &means[i]);
                apf_sub(&diff_j, &MAT_AT(m, k, j).re, &means[j]);
                apf_mul(&prod, &diff_i, &diff_j);
                apf_add(&cov_ij, &cov_ij, &prod);
            }
            /* Use N-1 for sample covariance */
            apf_div(&MAT_AT(r, i, j).re, &cov_ij, &divisor);
            apf_zero(&MAT_AT(r, i, j).im);
        }
    }
}

/* Correlation matrix: returns NxN matrix */
void mat_corrcoef(matrix_t *r, const matrix_t *m)
{
    int i, j, k;
    apf *means, *stds;
    apf sum, sum2, v_sq, var, divisor, diff_i, diff_j, prod, corr_ij, eps;
    
    means = (apf *)mat_arena_alloc(m->cols * sizeof(apf));
    stds = (apf *)mat_arena_alloc(m->cols * sizeof(apf));
    if (!means || !stds) return;
    
    /* Set epsilon for avoiding division by zero */
    apf_from_str(&eps, "1e-30");
    apf_from_int(&divisor, m->rows);
    
    /* Compute column means and stds */
    for (j = 0; j < m->cols; j++) {
        apf_zero(&sum);
        apf_zero(&sum2);
        for (i = 0; i < m->rows; i++) {
            apf_add(&sum, &sum, &MAT_AT(m, i, j).re);
            apf_mul(&v_sq, &MAT_AT(m, i, j).re, &MAT_AT(m, i, j).re);
            apf_add(&sum2, &sum2, &v_sq);
        }
        apf_div(&means[j], &sum, &divisor);
        /* var = sum2/n - mean^2 */
        apf_div(&var, &sum2, &divisor);
        apf_mul(&v_sq, &means[j], &means[j]);
        apf_sub(&var, &var, &v_sq);
        /* std = sqrt(max(var, eps)) */
        if (apf_cmp(&var, &eps) > 0) {
            apf_sqrt(&stds[j], &var);
        } else {
            apf_copy(&stds[j], &eps);
        }
    }
    
    /* Compute correlation matrix */
    mat_zero(r, m->cols, m->cols);
    if (!r->data) return;
    
    for (i = 0; i < m->cols; i++) {
        for (j = 0; j < m->cols; j++) {
            if (i == j) {
                apf_from_int(&MAT_AT(r, i, j).re, 1);
            } else {
                apf_zero(&corr_ij);
                for (k = 0; k < m->rows; k++) {
                    /* diff_i = (x[k,i] - mean[i]) / std[i] */
                    apf_sub(&diff_i, &MAT_AT(m, k, i).re, &means[i]);
                    apf_div(&diff_i, &diff_i, &stds[i]);
                    /* diff_j = (x[k,j] - mean[j]) / std[j] */
                    apf_sub(&diff_j, &MAT_AT(m, k, j).re, &means[j]);
                    apf_div(&diff_j, &diff_j, &stds[j]);
                    /* corr_ij += diff_i * diff_j */
                    apf_mul(&prod, &diff_i, &diff_j);
                    apf_add(&corr_ij, &corr_ij, &prod);
                }
                apf_div(&MAT_AT(r, i, j).re, &corr_ij, &divisor);
            }
            apf_zero(&MAT_AT(r, i, j).im);
        }
    }
}

/* Cross-tabulation: crosstab(idx1, idx2) returns contingency table */
void mat_crosstab(matrix_t *r, const matrix_t *a, const matrix_t *b)
{
    int i;
    int n = a->rows * a->cols;
    int max_a = 0, max_b = 0;
    int va, vb;
    
    if (n != b->rows * b->cols) {
        printf("Error: crosstab requires same-length vectors\n");
        mat_zero(r, 1, 1);
        return;
    }
    
    /* Find max values to determine table size */
    for (i = 0; i < n; i++) {
        va = (int)apf_to_double(&a->data[i].re);
        vb = (int)apf_to_double(&b->data[i].re);
        if (va > max_a) max_a = va;
        if (vb > max_b) max_b = vb;
    }
    
    if (max_a <= 0 || max_b <= 0) {
        mat_zero(r, 1, 1);
        return;
    }
    
    mat_zero(r, max_a, max_b);
    if (!r->data) return;
    
    /* Count co-occurrences */
    for (i = 0; i < n; i++) {
        va = (int)apf_to_double(&a->data[i].re) - 1; /* 0-based index */
        vb = (int)apf_to_double(&b->data[i].re) - 1;
        if (va >= 0 && va < max_a && vb >= 0 && vb < max_b) {
            double cur = apf_to_double(&MAT_AT(r, va, vb).re);
            apf_from_double(&MAT_AT(r, va, vb).re, cur + 1);
        }
    }
}

/* Rand Index: measure clustering similarity */
void mat_randindex(apfc *r, const matrix_t *a, const matrix_t *b)
{
    int i, j;
    int n = a->rows * a->cols;
    long a11 = 0, a00 = 0;
    long total_pairs;
    double ri;
    
    if (n != b->rows * b->cols || n < 2) {
        apf_zero(&r->re);
        apf_zero(&r->im);
        return;
    }
    
    /* Count pairs */
    for (i = 0; i < n; i++) {
        int ai = (int)apf_to_double(&a->data[i].re);
        int bi = (int)apf_to_double(&b->data[i].re);
        for (j = i + 1; j < n; j++) {
            int aj = (int)apf_to_double(&a->data[j].re);
            int bj = (int)apf_to_double(&b->data[j].re);
            int same_a = (ai == aj);
            int same_b = (bi == bj);
            if (same_a && same_b) a11++;
            else if (!same_a && !same_b) a00++;
        }
    }
    
    total_pairs = (long)n * (n - 1) / 2;
    ri = (double)(a11 + a00) / total_pairs;
    
    apf_from_double(&r->re, ri);
    apf_zero(&r->im);
}

/* ============================================================
 * MATLAB COMPATIBILITY FUNCTIONS
 * ============================================================ */

/* rref - Row Reduced Echelon Form (Gauss-Jordan elimination) */
void mat_rref(matrix_t *r, const matrix_t *a)
{
    int i, j, k, pivot_row;
    int rows = a->rows, cols = a->cols;
    int lead = 0;
    apfc pivot, factor, tmp;
    
    mat_copy(r, a);
    if (!r->data) return;
    
    for (i = 0; i < rows && lead < cols; i++) {
        /* Find pivot - largest absolute value in column */
        pivot_row = i;
        {
            double max_abs = apf_to_double(&MAT_AT(r, i, lead).re);
            if (max_abs < 0) max_abs = -max_abs;
            for (j = i + 1; j < rows; j++) {
                double cur_abs = apf_to_double(&MAT_AT(r, j, lead).re);
                if (cur_abs < 0) cur_abs = -cur_abs;
                if (cur_abs > max_abs) {
                    max_abs = cur_abs;
                    pivot_row = j;
                }
            }
        }
        
        /* Swap rows if needed */
        if (pivot_row != i) {
            for (k = 0; k < cols; k++) {
                tmp = MAT_AT(r, i, k);
                MAT_AT(r, i, k) = MAT_AT(r, pivot_row, k);
                MAT_AT(r, pivot_row, k) = tmp;
            }
        }
        
        pivot = MAT_AT(r, i, lead);
        if (apf_is_zero(&pivot.re) && apf_is_zero(&pivot.im)) {
            lead++;
            i--;
            continue;
        }
        
        /* Scale pivot row */
        for (k = 0; k < cols; k++) {
            apfc_div(&MAT_AT(r, i, k), &MAT_AT(r, i, k), &pivot);
        }
        
        /* Eliminate column in all rows */
        for (j = 0; j < rows; j++) {
            if (j != i) {
                factor = MAT_AT(r, j, lead);
                for (k = 0; k < cols; k++) {
                    apfc_mul(&tmp, &factor, &MAT_AT(r, i, k));
                    apfc_sub(&MAT_AT(r, j, k), &MAT_AT(r, j, k), &tmp);
                }
            }
        }
        lead++;
    }
}

/* orth - Orthonormal basis for column space (via SVD) */
void mat_orth(matrix_t *r, const matrix_t *a)
{
    matrix_t U, S, V;
    int i, j, rank_count = 0;
    double tol;
    
    mat_svd(&U, &S, &V, a);
    
    /* Find numerical rank */
    tol = 1e-10 * (a->rows > a->cols ? a->rows : a->cols);
    for (i = 0; i < S.rows && i < S.cols; i++) {
        if (apf_to_double(&MAT_AT(&S, i, i).re) > tol) {
            rank_count++;
        }
    }
    
    /* Extract first rank_count columns of U */
    if (rank_count == 0) {
        mat_zero(r, a->rows, 1);
        return;
    }
    
    mat_zero(r, U.rows, rank_count);
    for (i = 0; i < U.rows; i++) {
        for (j = 0; j < rank_count; j++) {
            MAT_AT(r, i, j) = MAT_AT(&U, i, j);
        }
    }
}

/* poly - Characteristic polynomial from square matrix eigenvalues
 * Returns coefficients [1, c_{n-1}, ..., c_1, c_0] where
 * p(x) = x^n + c_{n-1}*x^{n-1} + ... + c_1*x + c_0
 */
void mat_poly(matrix_t *r, const matrix_t *a)
{
    matrix_t eig_vals;
    int n, i, j;
    apfc *coeffs;
    
    if (a->rows != a->cols) {
        mat_zero(r, 1, 1);
        return;
    }
    
    n = a->rows;
    mat_eig(&eig_vals, a);
    
    /* Build polynomial from eigenvalues: (x - e1)(x - e2)... */
    mat_zero(r, 1, n + 1);
    coeffs = r->data;
    
    /* Start with [1] */
    apf_from_int(&coeffs[0].re, 1);
    apf_zero(&coeffs[0].im);
    
    /* Multiply by (x - eigenvalue[i]) for each eigenvalue */
    for (i = 0; i < n; i++) {
        apfc ei = eig_vals.data[i];
        /* Shift coefficients and subtract eigenvalue*previous */
        for (j = i + 1; j > 0; j--) {
            apfc tmp;
            apfc_mul(&tmp, &ei, &coeffs[j - 1]);
            apfc_sub(&coeffs[j], &coeffs[j], &tmp);
        }
    }
}

/* roots - Find polynomial roots via companion matrix eigenvalues
 * Input: coefficient vector [a_n, a_{n-1}, ..., a_1, a_0]
 * Output: column vector of roots
 */
void mat_roots(matrix_t *r, const matrix_t *coeffs)
{
    matrix_t companion;
    int n, i;
    apfc lead;
    
    n = coeffs->rows * coeffs->cols - 1;  /* Degree of polynomial */
    if (n < 1) {
        mat_zero(r, 0, 0);
        return;
    }
    
    /* Build companion matrix (second form)
     * For p(x) = a_n*x^n + a_{n-1}*x^{n-1} + ... + a_1*x + a_0
     * 
     * C = [0  0  ... 0  -a_0/a_n ]
     *     [1  0  ... 0  -a_1/a_n ]
     *     [0  1  ... 0  -a_2/a_n ]
     *     [       ...            ]
     *     [0  0  ... 1  -a_{n-1}/a_n]
     *
     * The eigenvalues of C are the roots of p(x).
     */
    mat_zero(&companion, n, n);
    
    /* Leading coefficient a_n */
    lead = coeffs->data[0];
    
    /* Subdiagonal of ones (i,i-1) for i=1..n-1 */
    for (i = 1; i < n; i++) {
        apf_from_int(&MAT_AT(&companion, i, i - 1).re, 1);
    }
    
    /* Last column: -a_i/a_n for i = 0..n-1 */
    for (i = 0; i < n; i++) {
        apfc coef, neg_coef;
        /* coeffs->data[n-i] contains a_i (because coeffs is [a_n, a_{n-1}, ..., a_1, a_0]) */
        apfc_div(&coef, &coeffs->data[n - i], &lead);
        apfc_neg(&neg_coef, &coef);
        MAT_AT(&companion, i, n - 1) = neg_coef;
    }
    
    /* Eigenvalues of companion matrix are the roots */
    mat_eig(r, &companion);
}

/* expm - Matrix exponential via Padé approximation */
void mat_expm(matrix_t *r, const matrix_t *a)
{
    matrix_t X, N, D, tmp, inv_D;
    int n, k;
    int q = 6;  /* Padé order */
    apfc c;
    double norm_est;
    int s = 0;
    
    if (a->rows != a->cols) {
        mat_zero(r, 0, 0);
        return;
    }
    
    n = a->rows;
    mat_copy(&X, a);
    
    /* Scale matrix to reduce norm */
    norm_est = 0;
    for (k = 0; k < n * n; k++) {
        double v = apf_to_double(&a->data[k].re);
        if (v < 0) v = -v;
        norm_est += v;
    }
    norm_est /= n;
    
    while (norm_est > 0.5 && s < 20) {
        apfc half;
        apf_from_double(&half.re, 0.5);
        apf_zero(&half.im);
        mat_scale(&X, &X, &half);
        norm_est *= 0.5;
        s++;
    }
    
    /* Initialize N = I, D = I */
    mat_identity(&N, n);
    mat_identity(&D, n);
    
    /* Padé approximation */
    mat_identity(&tmp, n);
    apf_from_int(&c.re, 1);
    apf_zero(&c.im);
    
    for (k = 1; k <= q; k++) {
        apfc coef;
        double ck = (double)(q - k + 1) / (k * (2 * q - k + 1));
        apf_from_double(&coef.re, ck);
        apf_zero(&coef.im);
        apfc_mul(&c, &c, &coef);
        
        mat_mul(&tmp, &tmp, &X);
        
        {
            matrix_t scaled;
            mat_scale(&scaled, &tmp, &c);
            mat_add(&N, &N, &scaled);
        }
        
        if (k % 2 == 0) {
            matrix_t scaled;
            mat_scale(&scaled, &tmp, &c);
            mat_add(&D, &D, &scaled);
        } else {
            matrix_t scaled;
            apfc neg_c;
            apfc_neg(&neg_c, &c);
            mat_scale(&scaled, &tmp, &neg_c);
            mat_add(&D, &D, &scaled);
        }
    }
    
    /* r = D^(-1) * N */
    mat_inv(&inv_D, &D);
    mat_mul(r, &inv_D, &N);
    
    /* Square back */
    for (k = 0; k < s; k++) {
        mat_mul(&tmp, r, r);
        mat_copy(r, &tmp);
    }
}

/* logm - Matrix logarithm via inverse scaling and squaring */
void mat_logm(matrix_t *r, const matrix_t *a)
{
    matrix_t X, Y, Z, I_n;
    int n, k, max_iter = 50;
    apfc two;
    int s = 0;
    
    if (a->rows != a->cols) {
        mat_zero(r, 0, 0);
        return;
    }
    
    n = a->rows;
    mat_copy(&X, a);
    mat_identity(&I_n, n);
    
    apf_from_int(&two.re, 2);
    apf_zero(&two.im);
    
    /* Scale: repeatedly take square root until close to I */
    for (s = 0; s < 20; s++) {
        double max_diff = 0;
        int i;
        for (i = 0; i < n; i++) {
            double diff = apf_to_double(&MAT_AT(&X, i, i).re) - 1.0;
            if (diff < 0) diff = -diff;
            if (diff > max_diff) max_diff = diff;
        }
        if (max_diff < 0.5) break;
        
        /* X = sqrtm(X) - approximate via Denman-Beavers iteration */
        {
            matrix_t Yk, Zk, Yk1, Zk1, inv_Z;
            int iter;
            mat_copy(&Yk, &X);
            mat_identity(&Zk, n);
            
            for (iter = 0; iter < 10; iter++) {
                mat_inv(&inv_Z, &Zk);
                mat_add(&Yk1, &Yk, &inv_Z);
                {
                    apfc half;
                    apf_from_double(&half.re, 0.5);
                    apf_zero(&half.im);
                    mat_scale(&Yk1, &Yk1, &half);
                }
                
                {
                    matrix_t inv_Y;
                    mat_inv(&inv_Y, &Yk);
                    mat_add(&Zk1, &Zk, &inv_Y);
                    {
                        apfc half;
                        apf_from_double(&half.re, 0.5);
                        apf_zero(&half.im);
                        mat_scale(&Zk1, &Zk1, &half);
                    }
                }
                
                mat_copy(&Yk, &Yk1);
                mat_copy(&Zk, &Zk1);
            }
            mat_copy(&X, &Yk);
        }
    }
    
    /* Now compute log(X) where X is close to I using series */
    /* log(X) = (X-I) - (X-I)^2/2 + (X-I)^3/3 - ... */
    mat_sub(&Y, &X, &I_n);  /* Y = X - I */
    mat_copy(r, &Y);        /* r = Y */
    mat_copy(&Z, &Y);       /* Z = Y (power accumulator) */
    
    for (k = 2; k <= max_iter; k++) {
        apfc coef;
        matrix_t term;
        double sign = (k % 2 == 0) ? -1.0 : 1.0;
        
        mat_mul(&Z, &Z, &Y);  /* Z = Z * Y = Y^k */
        apf_from_double(&coef.re, sign / k);
        apf_zero(&coef.im);
        mat_scale(&term, &Z, &coef);
        mat_add(r, r, &term);
    }
    
    /* Unscale: multiply by 2^s */
    {
        apfc scale;
        apf_from_double(&scale.re, (double)(1 << s));
        apf_zero(&scale.im);
        mat_scale(r, r, &scale);
    }
}

/* ============================================================
 * FFT/IFFT - Fast Fourier Transform (Cooley-Tukey)
 * ============================================================ */

/* Helper: Compute DFT for small sizes or non-power-of-2 */
static void dft_direct(apfc *out, const apfc *in, int n, int inverse)
{
    int k, m;
    double sign = inverse ? 1.0 : -1.0;
    
    for (k = 0; k < n; k++) {
        double sum_re = 0.0, sum_im = 0.0;
        
        for (m = 0; m < n; m++) {
            double angle = sign * 2.0 * 3.14159265358979323846 * k * m / n;
            double cos_a = cos(angle);
            double sin_a = sin(angle);
            double in_re = apf_to_double(&in[m].re);
            double in_im = apf_to_double(&in[m].im);
            
            sum_re += in_re * cos_a - in_im * sin_a;
            sum_im += in_re * sin_a + in_im * cos_a;
        }
        
        if (inverse) {
            sum_re /= n;
            sum_im /= n;
        }
        
        apf_from_double(&out[k].re, sum_re);
        apf_from_double(&out[k].im, sum_im);
    }
}

/* fft - Fast Fourier Transform */
void mat_fft(matrix_t *r, const matrix_t *a)
{
    int n = a->rows * a->cols;
    int i;
    apfc *temp;
    
    mat_zero(r, a->rows, a->cols);
    if (n == 0 || !a->data) return;
    
    temp = (apfc*)malloc(n * sizeof(apfc));
    if (!temp) return;
    
    /* Copy input */
    for (i = 0; i < n; i++) {
        temp[i] = a->data[i];
    }
    
    /* Perform DFT */
    dft_direct(r->data, temp, n, 0);
    
    free(temp);
}

/* ifft - Inverse Fast Fourier Transform */
void mat_ifft(matrix_t *r, const matrix_t *a)
{
    int n = a->rows * a->cols;
    int i;
    apfc *temp;
    
    mat_zero(r, a->rows, a->cols);
    if (n == 0 || !a->data) return;
    
    temp = (apfc*)malloc(n * sizeof(apfc));
    if (!temp) return;
    
    /* Copy input */
    for (i = 0; i < n; i++) {
        temp[i] = a->data[i];
    }
    
    /* Perform inverse DFT */
    dft_direct(r->data, temp, n, 1);
    
    free(temp);
}

/* ============================================================
 * Additional moving window functions
 * ============================================================ */

/* movprod - Moving product */
void mat_movprod(matrix_t *r, const matrix_t *a, int window)
{
    int n = a->rows * a->cols;
    int i, j;
    
    mat_zero(r, a->rows, a->cols);
    if (n == 0 || window < 1) return;
    
    for (i = 0; i < n; i++) {
        double prod = 1.0;
        int start = (i - window / 2 < 0) ? 0 : i - window / 2;
        int end = (i + (window - window/2) > n) ? n : i + (window - window/2);
        
        for (j = start; j < end; j++) {
            prod *= apf_to_double(&a->data[j].re);
        }
        apf_from_double(&r->data[i].re, prod);
    }
}

/* movmedian - Moving median */
void mat_movmedian(matrix_t *r, const matrix_t *a, int window)
{
    int n = a->rows * a->cols;
    int i, j, k;
    double *buf;
    
    mat_zero(r, a->rows, a->cols);
    if (n == 0 || window < 1) return;
    
    buf = (double*)malloc(window * sizeof(double));
    if (!buf) return;
    
    for (i = 0; i < n; i++) {
        int start = (i - window / 2 < 0) ? 0 : i - window / 2;
        int end = (i + (window - window/2) > n) ? n : i + (window - window/2);
        int count = end - start;
        double median;
        
        /* Copy to buffer */
        for (j = start; j < end; j++) {
            buf[j - start] = apf_to_double(&a->data[j].re);
        }
        
        /* Sort buffer (simple bubble sort for small windows) */
        for (j = 0; j < count - 1; j++) {
            for (k = 0; k < count - j - 1; k++) {
                if (buf[k] > buf[k+1]) {
                    double tmp = buf[k];
                    buf[k] = buf[k+1];
                    buf[k+1] = tmp;
                }
            }
        }
        
        /* Get median */
        if (count % 2 == 0) {
            median = (buf[count/2 - 1] + buf[count/2]) / 2.0;
        } else {
            median = buf[count/2];
        }
        
        apf_from_double(&r->data[i].re, median);
    }
    
    free(buf);
}

/* movvar - Moving variance */
void mat_movvar(matrix_t *r, const matrix_t *a, int window)
{
    int n = a->rows * a->cols;
    int i, j;
    
    mat_zero(r, a->rows, a->cols);
    if (n == 0 || window < 1) return;
    
    for (i = 0; i < n; i++) {
        double sum = 0.0, sum_sq = 0.0;
        int start = (i - window / 2 < 0) ? 0 : i - window / 2;
        int end = (i + (window - window/2) > n) ? n : i + (window - window/2);
        int count = end - start;
        double mean, var;
        
        for (j = start; j < end; j++) {
            double v = apf_to_double(&a->data[j].re);
            sum += v;
            sum_sq += v * v;
        }
        
        mean = sum / count;
        var = (sum_sq / count) - (mean * mean);
        if (var < 0) var = 0;  /* Numerical stability */
        
        apf_from_double(&r->data[i].re, var);
    }
}

/* meshgrid - Create 2D grid matrices */
void mat_meshgrid(matrix_t *X, matrix_t *Y, const matrix_t *x, const matrix_t *y)
{
    int nx = x->rows * x->cols;
    int ny = y->rows * y->cols;
    int i, j;
    
    mat_zero(X, ny, nx);
    mat_zero(Y, ny, nx);
    
    /* X: repeat x vector for each row */
    for (i = 0; i < ny; i++) {
        for (j = 0; j < nx; j++) {
            MAT_AT(X, i, j) = x->data[j];
        }
    }
    
    /* Y: repeat y value for each column */
    for (i = 0; i < ny; i++) {
        for (j = 0; j < nx; j++) {
            MAT_AT(Y, i, j) = y->data[i];
        }
    }
}

/* ============================================================
 * DATA ANALYSIS FUNCTIONS
 * ============================================================ */

/* histcounts - Histogram bin counts - pure APF
 * Returns counts for each bin (default 10 bins from min to max)
 */
void mat_histcounts(matrix_t *r, const matrix_t *data, int nbins)
{
    int n = data->rows * data->cols;
    int i, bin;
    apf min_val, max_val, range, bin_width, v, diff, one;
    
    if (n == 0 || nbins < 1) {
        mat_zero(r, 1, nbins > 0 ? nbins : 10);
        return;
    }
    
    /* Find min/max */
    apf_copy(&min_val, &data->data[0].re);
    apf_copy(&max_val, &data->data[0].re);
    for (i = 1; i < n; i++) {
        if (apf_cmp(&data->data[i].re, &min_val) < 0) {
            apf_copy(&min_val, &data->data[i].re);
        }
        if (apf_cmp(&data->data[i].re, &max_val) > 0) {
            apf_copy(&max_val, &data->data[i].re);
        }
    }
    
    apf_sub(&range, &max_val, &min_val);
    if (apf_cmp_int(&range, 0) == 0) {
        apf_from_int(&range, 1);  /* Handle constant data */
    }
    apf_from_int(&v, nbins);
    apf_div(&bin_width, &range, &v);
    
    mat_zero(r, 1, nbins);
    apf_from_int(&one, 1);
    
    /* Count values in each bin */
    for (i = 0; i < n; i++) {
        apf_sub(&diff, &data->data[i].re, &min_val);
        apf_div(&v, &diff, &bin_width);
        bin = (int)apf_to_long(&v);
        if (bin >= nbins) bin = nbins - 1;
        if (bin < 0) bin = 0;
        apf_add(&r->data[bin].re, &r->data[bin].re, &one);
    }
}

/* xcorr - Cross-correlation of two signals - pure APF */
void mat_xcorr(matrix_t *r, const matrix_t *x, const matrix_t *y)
{
    int nx = x->rows * x->cols;
    int ny = y->rows * y->cols;
    int len = nx + ny - 1;
    int i, j, lag;
    apf sum, prod;
    
    mat_zero(r, 1, len);
    
    /* Compute cross-correlation for all lags */
    for (lag = -(ny - 1); lag < nx; lag++) {
        int idx = lag + (ny - 1);  /* Output index */
        apf_zero(&sum);
        
        for (j = 0; j < ny; j++) {
            i = lag + j;
            if (i >= 0 && i < nx) {
                apf_mul(&prod, &x->data[i].re, &y->data[j].re);
                apf_add(&sum, &sum, &prod);
            }
        }
        apf_copy(&r->data[idx].re, &sum);
    }
}

/* xcov - Cross-covariance (mean-removed cross-correlation) - pure APF */
void mat_xcov(matrix_t *r, const matrix_t *x, const matrix_t *y)
{
    int nx = x->rows * x->cols;
    int ny = y->rows * y->cols;
    int i;
    apf mean_x, mean_y, divisor;
    matrix_t x_centered, y_centered;
    
    /* Compute means */
    apf_zero(&mean_x);
    apf_zero(&mean_y);
    for (i = 0; i < nx; i++) apf_add(&mean_x, &mean_x, &x->data[i].re);
    for (i = 0; i < ny; i++) apf_add(&mean_y, &mean_y, &y->data[i].re);
    apf_from_int(&divisor, nx);
    apf_div(&mean_x, &mean_x, &divisor);
    apf_from_int(&divisor, ny);
    apf_div(&mean_y, &mean_y, &divisor);
    
    /* Center the data */
    mat_zero(&x_centered, x->rows, x->cols);
    mat_zero(&y_centered, y->rows, y->cols);
    
    for (i = 0; i < nx; i++) {
        apf_sub(&x_centered.data[i].re, &x->data[i].re, &mean_x);
    }
    for (i = 0; i < ny; i++) {
        apf_sub(&y_centered.data[i].re, &y->data[i].re, &mean_y);
    }
    
    /* Compute cross-correlation of centered data */
    mat_xcorr(r, &x_centered, &y_centered);
}

/* isoutlier - Detect outliers using median absolute deviation (MAD) - pure APF
 * Returns 1 for outliers, 0 otherwise
 */
void mat_isoutlier(matrix_t *r, const matrix_t *data)
{
    int n = data->rows * data->cols;
    int i, j;
    apf *sorted, median, mad, threshold, diff, scale, three, one, tmp;
    
    mat_zero(r, data->rows, data->cols);
    if (n < 3) return;  /* Need at least 3 points */
    
    sorted = (apf*)malloc(n * sizeof(apf));
    if (!sorted) return;
    
    /* Copy and sort for median */
    for (i = 0; i < n; i++) {
        apf_copy(&sorted[i], &data->data[i].re);
    }
    for (i = 0; i < n - 1; i++) {
        for (j = 0; j < n - i - 1; j++) {
            if (apf_cmp(&sorted[j], &sorted[j + 1]) > 0) {
                apf_copy(&tmp, &sorted[j]);
                apf_copy(&sorted[j], &sorted[j + 1]);
                apf_copy(&sorted[j + 1], &tmp);
            }
        }
    }
    
    /* Median */
    if (n % 2 == 0) {
        apf_add(&median, &sorted[n/2 - 1], &sorted[n/2]);
        apf_from_int(&tmp, 2);
        apf_div(&median, &median, &tmp);
    } else {
        apf_copy(&median, &sorted[n/2]);
    }
    
    /* Compute MAD (median absolute deviation) */
    for (i = 0; i < n; i++) {
        apf_sub(&sorted[i], &data->data[i].re, &median);
        if (apf_cmp_int(&sorted[i], 0) < 0) {
            apf_neg(&sorted[i], &sorted[i]);
        }
    }
    for (i = 0; i < n - 1; i++) {
        for (j = 0; j < n - i - 1; j++) {
            if (apf_cmp(&sorted[j], &sorted[j + 1]) > 0) {
                apf_copy(&tmp, &sorted[j]);
                apf_copy(&sorted[j], &sorted[j + 1]);
                apf_copy(&sorted[j + 1], &tmp);
            }
        }
    }
    if (n % 2 == 0) {
        apf_add(&mad, &sorted[n/2 - 1], &sorted[n/2]);
        apf_from_int(&tmp, 2);
        apf_div(&mad, &mad, &tmp);
    } else {
        apf_copy(&mad, &sorted[n/2]);
    }
    
    /* Scale MAD to be consistent with std dev for normal distribution */
    apf_from_str(&scale, "1.4826");
    apf_mul(&mad, &mad, &scale);
    
    /* Threshold: 3 scaled MADs from median */
    apf_from_int(&three, 3);
    apf_mul(&threshold, &three, &mad);
    if (apf_cmp_int(&threshold, 0) == 0) {
        apf_from_str(&threshold, "1e-30");  /* Handle constant data */
    }
    
    /* Mark outliers */
    apf_from_int(&one, 1);
    for (i = 0; i < n; i++) {
        apf_sub(&diff, &data->data[i].re, &median);
        if (apf_cmp_int(&diff, 0) < 0) {
            apf_neg(&diff, &diff);
        }
        if (apf_cmp(&diff, &threshold) > 0) {
            apf_copy(&r->data[i].re, &one);
        }
        /* else already zero from mat_zero */
    }
    
    free(sorted);
}

/* ============================================================
 * NUMERICAL INTEGRATION
 * ============================================================ */

/* integral_trapz - Trapezoidal integration (already have trapz, this is alias) */
/* Note: mat_trapz already exists */

/* integral_simpson - Simpson's rule integration of a vector
 * Assumes uniform spacing
 */
void mat_simpson(matrix_t *r, const matrix_t *y, double h)
{
    int n = y->rows * y->cols;
    int i;
    double sum = 0;
    
    mat_zero(r, 1, 1);
    if (n < 3) {
        /* Fall back to trapezoid for too few points */
        for (i = 0; i < n - 1; i++) {
            double y0 = apf_to_double(&y->data[i].re);
            double y1 = apf_to_double(&y->data[i + 1].re);
            sum += (y0 + y1) * h / 2.0;
        }
        apf_from_double(&r->data[0].re, sum);
        return;
    }
    
    /* Simpson's 1/3 rule */
    sum = apf_to_double(&y->data[0].re) + apf_to_double(&y->data[n - 1].re);
    
    for (i = 1; i < n - 1; i++) {
        double yi = apf_to_double(&y->data[i].re);
        if (i % 2 == 1) {
            sum += 4.0 * yi;
        } else {
            sum += 2.0 * yi;
        }
    }
    
    sum *= h / 3.0;
    apf_from_double(&r->data[0].re, sum);
}

/* ============================================================
 * HESSENBERG FORM
 * ============================================================ */

/* hess - Hessenberg form via Householder reflections */
void mat_hess(matrix_t *H, const matrix_t *A)
{
    int n = A->rows;
    int i, j, k;
    
    if (A->rows != A->cols) {
        mat_zero(H, 0, 0);
        return;
    }
    
    mat_copy(H, A);
    
    for (k = 0; k < n - 2; k++) {
        /* Compute Householder vector for column k below diagonal */
        double norm = 0;
        double alpha, r;
        
        for (i = k + 1; i < n; i++) {
            double v = apf_to_double(&MAT_AT(H, i, k).re);
            norm += v * v;
        }
        norm = sqrt(norm);
        
        if (norm < 1e-15) continue;
        
        alpha = apf_to_double(&MAT_AT(H, k + 1, k).re);
        if (alpha >= 0) {
            alpha = -norm;
        } else {
            alpha = norm;
        }
        
        r = sqrt(0.5 * (alpha * alpha - alpha * apf_to_double(&MAT_AT(H, k + 1, k).re)));
        if (r < 1e-15) continue;
        
        {
            double *v = (double*)malloc(n * sizeof(double));
            if (!v) continue;
            
            for (i = 0; i < n; i++) v[i] = 0;
            v[k + 1] = (apf_to_double(&MAT_AT(H, k + 1, k).re) - alpha) / (2.0 * r);
            for (i = k + 2; i < n; i++) {
                v[i] = apf_to_double(&MAT_AT(H, i, k).re) / (2.0 * r);
            }
            
            /* H = (I - 2*v*v') * H * (I - 2*v*v') */
            /* Apply from left: H = H - 2*v*(v'*H) */
            for (j = k; j < n; j++) {
                double dot = 0;
                for (i = k + 1; i < n; i++) {
                    dot += v[i] * apf_to_double(&MAT_AT(H, i, j).re);
                }
                for (i = k + 1; i < n; i++) {
                    double val = apf_to_double(&MAT_AT(H, i, j).re);
                    apf_from_double(&MAT_AT(H, i, j).re, val - 2.0 * v[i] * dot);
                }
            }
            
            /* Apply from right: H = H - 2*(H*v)*v' */
            for (i = 0; i < n; i++) {
                double dot = 0;
                for (j = k + 1; j < n; j++) {
                    dot += apf_to_double(&MAT_AT(H, i, j).re) * v[j];
                }
                for (j = k + 1; j < n; j++) {
                    double val = apf_to_double(&MAT_AT(H, i, j).re);
                    apf_from_double(&MAT_AT(H, i, j).re, val - 2.0 * dot * v[j]);
                }
            }
            
            free(v);
        }
    }
    
    /* Zero out below subdiagonal */
    for (i = 2; i < n; i++) {
        for (j = 0; j < i - 1; j++) {
            apf_zero(&MAT_AT(H, i, j).re);
            apf_zero(&MAT_AT(H, i, j).im);
        }
    }
}

/* balance - Balance matrix for eigenvalue computation
 * Scales rows and columns to make them roughly equal in norm
 */
void mat_balance(matrix_t *B, const matrix_t *A)
{
    int n = A->rows;
    int i, j, iter;
    int converged;
    double radix = 2.0;
    
    if (A->rows != A->cols) {
        mat_zero(B, 0, 0);
        return;
    }
    
    mat_copy(B, A);
    
    /* Iterate until converged */
    for (iter = 0; iter < 100; iter++) {
        converged = 1;
        
        for (i = 0; i < n; i++) {
            double c = 0, r = 0;
            double f, g, s;
            
            /* Column and row sums of absolute values */
            for (j = 0; j < n; j++) {
                if (j != i) {
                    double val = apf_to_double(&MAT_AT(B, j, i).re);
                    c += (val < 0 ? -val : val);
                    val = apf_to_double(&MAT_AT(B, i, j).re);
                    r += (val < 0 ? -val : val);
                }
            }
            
            if (c == 0 || r == 0) continue;
            
            g = r / radix;
            f = 1.0;
            s = c + r;
            
            while (c < g) {
                f *= radix;
                c *= radix * radix;
            }
            
            g = r * radix;
            while (c > g) {
                f /= radix;
                c /= radix * radix;
            }
            
            if ((c + r) / f < 0.95 * s) {
                converged = 0;
                g = 1.0 / f;
                
                /* Scale column i by g, row i by f */
                for (j = 0; j < n; j++) {
                    double val = apf_to_double(&MAT_AT(B, j, i).re);
                    apf_from_double(&MAT_AT(B, j, i).re, val * g);
                    val = apf_to_double(&MAT_AT(B, i, j).re);
                    apf_from_double(&MAT_AT(B, i, j).re, val * f);
                }
            }
        }
        
        if (converged) break;
    }
}

/* ============================================================
 * NUMERICAL METHODS - ROOT FINDING
 * ============================================================ */

/* bisect - Bisection method for root finding
 * Finds root of f(x) = y_data where x_data are evaluation points
 * This is a simplified version for discrete data
 */
void mat_bisect_interp(matrix_t *r, const matrix_t *x_data, const matrix_t *y_data, double target)
{
    int n = x_data->rows * x_data->cols;
    int i;
    double x0 = 0, x1 = 0, y0 = 0, y1 = 0;
    int found = 0;
    
    mat_zero(r, 1, 1);
    if (n < 2) return;
    
    /* Find interval containing target */
    for (i = 0; i < n - 1; i++) {
        y0 = apf_to_double(&y_data->data[i].re) - target;
        y1 = apf_to_double(&y_data->data[i + 1].re) - target;
        if (y0 * y1 <= 0) {
            x0 = apf_to_double(&x_data->data[i].re);
            x1 = apf_to_double(&x_data->data[i + 1].re);
            found = 1;
            break;
        }
    }
    
    if (!found) {
        apf_from_double(&r->data[0].re, 0.0/0.0);  /* NaN */
        return;
    }
    
    /* Linear interpolation for root */
    if (fabs(y1 - y0) < 1e-15) {
        apf_from_double(&r->data[0].re, (x0 + x1) / 2);
    } else {
        apf_from_double(&r->data[0].re, x0 - y0 * (x1 - x0) / (y1 - y0));
    }
}

/* ============================================================
 * NUMERICAL METHODS - ITERATIVE LINEAR SOLVERS
 * ============================================================ */

/* jacobi - Jacobi iterative method for Ax = b
 * Returns solution vector after specified iterations
 */
void mat_jacobi(matrix_t *r, const matrix_t *A, const matrix_t *b, int max_iter)
{
    int n = A->rows;
    int i, j, iter;
    double *x, *x_new;
    
    if (A->rows != A->cols || A->rows != b->rows * b->cols) {
        mat_zero(r, 0, 0);
        return;
    }
    
    x = (double*)malloc(n * sizeof(double));
    x_new = (double*)malloc(n * sizeof(double));
    if (!x || !x_new) {
        if (x) free(x);
        if (x_new) free(x_new);
        mat_zero(r, 0, 0);
        return;
    }
    
    /* Initialize x to zeros */
    for (i = 0; i < n; i++) x[i] = 0;
    
    for (iter = 0; iter < max_iter; iter++) {
        for (i = 0; i < n; i++) {
            double sum = apf_to_double(&b->data[i].re);
            double diag = apf_to_double(&MAT_AT(A, i, i).re);
            
            for (j = 0; j < n; j++) {
                if (j != i) {
                    sum -= apf_to_double(&MAT_AT(A, i, j).re) * x[j];
                }
            }
            
            x_new[i] = (fabs(diag) > 1e-15) ? sum / diag : 0;
        }
        
        /* Check convergence */
        {
            double max_diff = 0;
            for (i = 0; i < n; i++) {
                double diff = fabs(x_new[i] - x[i]);
                if (diff > max_diff) max_diff = diff;
                x[i] = x_new[i];
            }
            if (max_diff < 1e-10) break;
        }
    }
    
    mat_zero(r, n, 1);
    for (i = 0; i < n; i++) {
        apf_from_double(&r->data[i].re, x[i]);
    }
    
    free(x);
    free(x_new);
}

/* gaussseidel - Gauss-Seidel iterative method for Ax = b */
void mat_gaussseidel(matrix_t *r, const matrix_t *A, const matrix_t *b, int max_iter)
{
    int n = A->rows;
    int i, j, iter;
    double *x;
    
    if (A->rows != A->cols || A->rows != b->rows * b->cols) {
        mat_zero(r, 0, 0);
        return;
    }
    
    x = (double*)malloc(n * sizeof(double));
    if (!x) {
        mat_zero(r, 0, 0);
        return;
    }
    
    /* Initialize x to zeros */
    for (i = 0; i < n; i++) x[i] = 0;
    
    for (iter = 0; iter < max_iter; iter++) {
        double max_diff = 0;
        
        for (i = 0; i < n; i++) {
            double sum = apf_to_double(&b->data[i].re);
            double diag = apf_to_double(&MAT_AT(A, i, i).re);
            double x_old = x[i];
            
            for (j = 0; j < n; j++) {
                if (j != i) {
                    sum -= apf_to_double(&MAT_AT(A, i, j).re) * x[j];
                }
            }
            
            x[i] = (fabs(diag) > 1e-15) ? sum / diag : 0;
            
            {
                double diff = fabs(x[i] - x_old);
                if (diff > max_diff) max_diff = diff;
            }
        }
        
        if (max_diff < 1e-10) break;
    }
    
    mat_zero(r, n, 1);
    for (i = 0; i < n; i++) {
        apf_from_double(&r->data[i].re, x[i]);
    }
    
    free(x);
}

/* ============================================================
 * NUMERICAL METHODS - ODE SOLVERS
 * ============================================================ */

/* euler_step - Single Euler step for dy/dt = f(t, y)
 * y_new = y + h * f(t, y)
 * For demonstration, we assume f is stored as slope data
 */
void mat_euler(matrix_t *r, const matrix_t *y0, const matrix_t *slopes, double h, int steps)
{
    int n = y0->rows * y0->cols;
    int i, step;
    double *y;
    
    if (steps < 1) steps = 1;
    
    y = (double*)malloc(n * sizeof(double));
    if (!y) {
        mat_zero(r, 0, 0);
        return;
    }
    
    /* Initialize */
    for (i = 0; i < n; i++) {
        y[i] = apf_to_double(&y0->data[i].re);
    }
    
    /* Euler steps using average slope (simplified) */
    for (step = 0; step < steps && step < slopes->cols; step++) {
        for (i = 0; i < n && i < slopes->rows; i++) {
            double slope = apf_to_double(&MAT_AT(slopes, i, step).re);
            y[i] += h * slope;
        }
    }
    
    mat_zero(r, n, 1);
    for (i = 0; i < n; i++) {
        apf_from_double(&r->data[i].re, y[i]);
    }
    
    free(y);
}

/* rk4_demo - Runge-Kutta 4th order demonstration
 * For a simple ODE dy/dt = f(y) where f is given as a coefficient
 * y' = coef * y  =>  y(t) = y0 * exp(coef * t)
 */
void mat_rk4_exp(matrix_t *r, double y0, double coef, double h, int steps)
{
    int i;
    double y = y0;
    double k1, k2, k3, k4;
    
    mat_zero(r, 1, steps + 1);
    apf_from_double(&r->data[0].re, y);
    
    for (i = 0; i < steps; i++) {
        k1 = coef * y;
        k2 = coef * (y + h * k1 / 2);
        k3 = coef * (y + h * k2 / 2);
        k4 = coef * (y + h * k3);
        
        y += h * (k1 + 2*k2 + 2*k3 + k4) / 6;
        apf_from_double(&r->data[i + 1].re, y);
    }
}

/* ============================================================
 * NUMERICAL METHODS - INTEGRATION
 * ============================================================ */

/* simpson - Simpson's rule integration of equally-spaced data */
void mat_simpson(matrix_t *r, const matrix_t *y, double h);  /* Already declared */

/* quad_trapz - Adaptive trapezoidal integration
 * Simple adaptive algorithm that doubles intervals
 */
void mat_quad_trapz(matrix_t *r, const matrix_t *y, double h, double tol)
{
    int n = y->rows * y->cols;
    int i;
    double coarse = 0, fine = 0;
    
    (void)tol;  /* For future adaptive refinement */
    
    mat_zero(r, 1, 1);
    if (n < 2) return;
    
    /* Coarse estimate */
    coarse = apf_to_double(&y->data[0].re) + apf_to_double(&y->data[n-1].re);
    for (i = 1; i < n - 1; i++) {
        coarse += 2 * apf_to_double(&y->data[i].re);
    }
    coarse *= h / 2;
    
    /* If we had midpoint values, we could refine. For now, use Simpson's */
    if (n >= 3) {
        fine = apf_to_double(&y->data[0].re) + apf_to_double(&y->data[n-1].re);
        for (i = 1; i < n - 1; i++) {
            if (i % 2 == 1) {
                fine += 4 * apf_to_double(&y->data[i].re);
            } else {
                fine += 2 * apf_to_double(&y->data[i].re);
            }
        }
        fine *= h / 3;
        apf_from_double(&r->data[0].re, fine);
    } else {
        apf_from_double(&r->data[0].re, coarse);
    }
}

/* ============================================================
 * NUMERICAL METHODS - INTERPOLATION
 * ============================================================ */

/* lagrange - Lagrange polynomial interpolation at point xi */
void mat_lagrange(matrix_t *r, const matrix_t *x, const matrix_t *y, double xi)
{
    int n = x->rows * x->cols;
    int i, j;
    double result = 0;
    
    mat_zero(r, 1, 1);
    if (n < 1) return;
    
    for (i = 0; i < n; i++) {
        double term = apf_to_double(&y->data[i].re);
        double xi_val = apf_to_double(&x->data[i].re);
        
        for (j = 0; j < n; j++) {
            if (j != i) {
                double xj = apf_to_double(&x->data[j].re);
                term *= (xi - xj) / (xi_val - xj);
            }
        }
        result += term;
    }
    
    apf_from_double(&r->data[0].re, result);
}

/* pchip - Piecewise Cubic Hermite Interpolating Polynomial
 * Shape-preserving interpolation
 */
void mat_pchip(matrix_t *r, const matrix_t *x, const matrix_t *y, const matrix_t *xi)
{
    int n = x->rows * x->cols;
    int ni = xi->rows * xi->cols;
    int i, j;
    double *h, *delta, *d;
    
    if (n < 2) {
        mat_zero(r, 0, 0);
        return;
    }
    
    h = (double*)malloc((n - 1) * sizeof(double));
    delta = (double*)malloc((n - 1) * sizeof(double));
    d = (double*)malloc(n * sizeof(double));
    
    if (!h || !delta || !d) {
        if (h) free(h);
        if (delta) free(delta);
        if (d) free(d);
        mat_zero(r, 0, 0);
        return;
    }
    
    /* Compute intervals and slopes */
    for (i = 0; i < n - 1; i++) {
        h[i] = apf_to_double(&x->data[i + 1].re) - apf_to_double(&x->data[i].re);
        delta[i] = (apf_to_double(&y->data[i + 1].re) - apf_to_double(&y->data[i].re)) / h[i];
    }
    
    /* Compute derivatives using shape-preserving formula */
    d[0] = delta[0];
    d[n - 1] = delta[n - 2];
    
    for (i = 1; i < n - 1; i++) {
        if (delta[i - 1] * delta[i] > 0) {
            /* Harmonic mean weighted by interval lengths */
            double w1 = 2 * h[i] + h[i - 1];
            double w2 = h[i] + 2 * h[i - 1];
            d[i] = (w1 + w2) / (w1 / delta[i - 1] + w2 / delta[i]);
        } else {
            d[i] = 0;  /* Extremum point */
        }
    }
    
    /* Interpolate at xi points */
    mat_zero(r, 1, ni);
    
    for (j = 0; j < ni; j++) {
        double xij = apf_to_double(&xi->data[j].re);
        double result;
        int k;
        
        /* Find interval */
        k = 0;
        for (i = 0; i < n - 1; i++) {
            if (xij >= apf_to_double(&x->data[i].re) && 
                xij <= apf_to_double(&x->data[i + 1].re)) {
                k = i;
                break;
            }
        }
        
        /* Hermite cubic */
        {
            double x0 = apf_to_double(&x->data[k].re);
            double x1 = apf_to_double(&x->data[k + 1].re);
            double y0 = apf_to_double(&y->data[k].re);
            double y1 = apf_to_double(&y->data[k + 1].re);
            double t = (xij - x0) / (x1 - x0);
            double t2 = t * t;
            double t3 = t2 * t;
            double h00 = 2*t3 - 3*t2 + 1;
            double h10 = t3 - 2*t2 + t;
            double h01 = -2*t3 + 3*t2;
            double h11 = t3 - t2;
            
            result = h00 * y0 + h10 * (x1 - x0) * d[k] + 
                     h01 * y1 + h11 * (x1 - x0) * d[k + 1];
        }
        
        apf_from_double(&r->data[j].re, result);
    }
    
    free(h);
    free(delta);
    free(d);
}
