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
    matrix_t AtA;
    matrix_t Q, R;
    matrix_t tmp, tmp2;
    matrix_t At;
    apf tolerance, off_norm, sum, abs_val;
    
    m = A->rows;
    n = A->cols;
    min_mn = (m < n) ? m : n;
    max_iter = 200;
    
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
    mat_mul(&AtA, &At, A);  /* n x n */
    
    /* Initialize V with identity */
    mat_identity(V, n);
    
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
    for (iter = 0; iter < max_iter; iter++) {
        if (!mat_qr(&Q, &R, &tmp)) break;
        
        /* Update V */
        mat_mul(&tmp2, V, &Q);
        mat_copy(V, &tmp2);
        
        /* tmp = R * Q (similar transformation) */
        mat_mul(&tmp2, &R, &Q);
        mat_copy(&tmp, &tmp2);
        
        /* Check convergence */
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
        if (apf_lt(&off_norm, &tolerance)) break;
    }
    
    /* Singular values are square roots of diagonal of converged matrix */
    mat_zero(S, m, n);
    for (i = 0; i < min_mn; i++) {
        apfc sv;
        apfc_abs(&sv.re, &MAT_AT(&tmp, i, i));
        apf_sqrt(&sv.re, &sv.re);
        apf_zero(&sv.im);
        MAT_AT(S, i, i) = sv;
    }
    
    /* Compute U = A * V * S^(-1) */
    mat_mul(&tmp, A, V);
    mat_zero(U, m, m);
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
    matrix_t R;
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
    matrix_t Qk, R;
    matrix_t tmp;
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

