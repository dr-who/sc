/* matrix.c - Matrix operations for scalc
 * C89 compliant for Watcom C / DOS
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "matrix.h"
#include "apfx.h"

/* Global display digits setting */
extern int display_digits;

/* ========== Matrix Arena Allocator (resets each command) ========== */
#if defined(PLATFORM_DOS) || defined(SCALC_MEDIUM)
/* DOS: Use static arrays */
static char mat_arena_static[MAT_ARENA_DEFAULT];
static char mat_persist_static[MAT_PERSIST_DEFAULT];
static char *mat_arena = mat_arena_static;
static char *mat_persist = mat_persist_static;
static size_t mat_arena_total = MAT_ARENA_DEFAULT;
static size_t mat_persist_total = MAT_PERSIST_DEFAULT;
#else
/* Linux/modern: Use dynamically allocated memory */
static char *mat_arena = NULL;
static char *mat_persist = NULL;
static size_t mat_arena_total = 0;
static size_t mat_persist_total = 0;
#endif

static size_t mat_arena_pos = 0;
static size_t mat_persist_pos = 0;

/* Initialize memory pools */
int mat_memory_init(size_t arena_bytes, size_t persist_bytes)
{
#if defined(PLATFORM_DOS) || defined(SCALC_MEDIUM)
    /* DOS uses static allocation, ignore parameters */
    (void)arena_bytes;
    (void)persist_bytes;
    return 1;
#else
    /* Use defaults if 0 passed */
    if (arena_bytes == 0) arena_bytes = MAT_ARENA_DEFAULT;
    if (persist_bytes == 0) persist_bytes = MAT_PERSIST_DEFAULT;
    
    /* Free any existing allocation */
    mat_memory_free();
    
    /* Allocate arena */
    mat_arena = (char *)malloc(arena_bytes);
    if (!mat_arena) {
        fprintf(stderr, "Error: Failed to allocate %lu MB arena\n",
                (unsigned long)(arena_bytes / (1024*1024)));
        return 0;
    }
    mat_arena_total = arena_bytes;
    
    /* Allocate persistent storage */
    mat_persist = (char *)malloc(persist_bytes);
    if (!mat_persist) {
        fprintf(stderr, "Error: Failed to allocate %lu MB persistent storage\n",
                (unsigned long)(persist_bytes / (1024*1024)));
        free(mat_arena);
        mat_arena = NULL;
        mat_arena_total = 0;
        return 0;
    }
    mat_persist_total = persist_bytes;
    
    mat_arena_pos = 0;
    mat_persist_pos = 0;
    
    return 1;
#endif
}

void mat_memory_free(void)
{
#if !defined(PLATFORM_DOS) && !defined(SCALC_MEDIUM)
    if (mat_arena) {
        free(mat_arena);
        mat_arena = NULL;
    }
    if (mat_persist) {
        free(mat_persist);
        mat_persist = NULL;
    }
    mat_arena_total = 0;
    mat_persist_total = 0;
#endif
}

void mat_arena_reset(void)
{
    if (mat_arena && mat_arena_total > 0) {
        /* Only poison a reasonable amount to avoid slow startup on huge arenas */
        size_t poison_size = mat_arena_pos < 1024*1024 ? mat_arena_pos : 1024*1024;
        if (poison_size > 0) {
            memset(mat_arena, 0xCD, poison_size);
        }
    }
    mat_arena_pos = 0;
}

void *mat_arena_alloc(size_t bytes)
{
    size_t aligned_pos;
    void *ptr;
    
    if (!mat_arena || mat_arena_total == 0) {
        printf("Error: Matrix arena not initialized\n");
        return NULL;
    }
    
    /* Align to 8 bytes */
    aligned_pos = (mat_arena_pos + 7) & ~(size_t)7;
    
    if (aligned_pos + bytes > mat_arena_total) {
        printf("Error: Matrix arena exhausted (%lu bytes requested, %lu available)\n",
               (unsigned long)bytes, (unsigned long)(mat_arena_total - aligned_pos));
        return NULL;
    }
    
    ptr = mat_arena + aligned_pos;
    mat_arena_pos = aligned_pos + bytes;
    
    return ptr;
}

size_t mat_arena_remaining(void)
{
    return mat_arena_total - mat_arena_pos;
}

size_t mat_arena_size(void)
{
    return mat_arena_total;
}

/* ========== Persistent Storage (for named variables) ========== */

void mat_persist_reset(void)
{
    mat_persist_pos = 0;
}

void *mat_persist_alloc(size_t bytes)
{
    size_t aligned_pos;
    void *ptr;
    
    if (!mat_persist || mat_persist_total == 0) {
        printf("Error: Persistent storage not initialized\n");
        return NULL;
    }
    
    aligned_pos = (mat_persist_pos + 7) & ~(size_t)7;
    
    if (aligned_pos + bytes > mat_persist_total) {
        printf("Error: Persistent storage exhausted (%lu MB used of %lu MB)\n",
               (unsigned long)(aligned_pos / (1024*1024)),
               (unsigned long)(mat_persist_total / (1024*1024)));
        return NULL;
    }
    
    ptr = mat_persist + aligned_pos;
    mat_persist_pos = aligned_pos + bytes;
    
    return ptr;
}

void mat_persist_free(void *ptr)
{
    /* Simple allocator - no individual free, just reset all */
    (void)ptr;
}

size_t mat_persist_remaining(void)
{
    return mat_persist_total - mat_persist_pos;
}

size_t mat_persist_size(void)
{
    return mat_persist_total;
}

/* Copy matrix to persistent storage */
int mat_copy_persist(matrix_t *dst, const matrix_t *src)
{
    size_t bytes;
    int i, n;
    
    if (!src || !src->data || src->rows <= 0 || src->cols <= 0) {
        dst->rows = 0;
        dst->cols = 0;
        dst->data = NULL;
        return 0;
    }
    
    bytes = (size_t)src->rows * (size_t)src->cols * sizeof(apfc);
    dst->data = (apfc *)mat_persist_alloc(bytes);
    if (!dst->data) return 0;
    
    dst->rows = src->rows;
    dst->cols = src->cols;
    n = src->rows * src->cols;
    for (i = 0; i < n; i++) {
        dst->data[i] = src->data[i];
    }
    return 1;
}

/* Allocate zero matrix in persistent storage */
int mat_zero_persist(matrix_t *m, int rows, int cols)
{
    size_t bytes;
    int i, n;
    
    bytes = (size_t)rows * (size_t)cols * sizeof(apfc);
    m->data = (apfc *)mat_persist_alloc(bytes);
    if (!m->data) {
        m->rows = 0;
        m->cols = 0;
        return 0;
    }
    
    m->rows = rows;
    m->cols = cols;
    n = rows * cols;
    for (i = 0; i < n; i++) {
        apf_zero(&m->data[i].re);
        apf_zero(&m->data[i].im);
    }
    return 1;
}

/* ========== Initialization ========== */

/* Allocate matrix data from arena */
static int mat_alloc_data(matrix_t *m, int rows, int cols)
{
    size_t bytes = (size_t)rows * (size_t)cols * sizeof(apfc);
    m->rows = rows;
    m->cols = cols;
    m->data = (apfc *)mat_arena_alloc(bytes);
    return m->data != NULL;
}

void mat_init(matrix_t *m, int rows, int cols)
{
    mat_alloc_data(m, rows, cols);
}

void mat_zero(matrix_t *m, int rows, int cols)
{
    int i, n;
    /* Always allocate fresh from arena - prevents uninitialized pointer bugs */
    if (!mat_alloc_data(m, rows, cols)) return;
    n = rows * cols;
    for (i = 0; i < n; i++) {
        apf_zero(&m->data[i].re);
        apf_zero(&m->data[i].im);
    }
}

void mat_identity(matrix_t *m, int n)
{
    int i;
    mat_zero(m, n, n);
    if (!m->data) return;
    for (i = 0; i < n; i++) {
        apf_from_int(&MAT_AT(m, i, i).re, 1);
    }
}

void mat_ones(matrix_t *m, int rows, int cols)
{
    int i, n;
    if (!mat_alloc_data(m, rows, cols)) return;
    n = rows * cols;
    for (i = 0; i < n; i++) {
        apf_from_int(&m->data[i].re, 1);
        apf_zero(&m->data[i].im);
    }
}

void mat_copy(matrix_t *dst, const matrix_t *src)
{
    int i, n;
    if (!src || src->rows <= 0 || src->cols <= 0 || !src->data) {
        dst->rows = 0;
        dst->cols = 0;
        dst->data = NULL;
        return;
    }
    if (!mat_alloc_data(dst, src->rows, src->cols)) return;
    n = src->rows * src->cols;
    for (i = 0; i < n; i++) {
        dst->data[i] = src->data[i];
    }
}

/* Resize matrix - allocates from arena if needed */
int mat_resize(matrix_t *m, int rows, int cols)
{
    /* If data is NULL or size changed, allocate new space */
    if (!m->data || m->rows != rows || m->cols != cols) {
        return mat_alloc_data(m, rows, cols);
    }
    return 1;  /* Already correct size */
}

/* ========== Basic Info ========== */

int mat_is_square(const matrix_t *m)
{
    return m->rows == m->cols;
}

int mat_is_vector(const matrix_t *m)
{
    return m->rows == 1 || m->cols == 1;
}

int mat_same_size(const matrix_t *a, const matrix_t *b)
{
    return a->rows == b->rows && a->cols == b->cols;
}

int mat_length(const matrix_t *m)
{
    return m->rows > m->cols ? m->rows : m->cols;
}

/* ========== Matrix Multiplication ========== */

void mat_mul(matrix_t *r, const matrix_t *a, const matrix_t *b)
{
    int i, j, k;
    matrix_t tmp = {0, 0, NULL};
    static int mul_debug = -1;
    
    if (mul_debug < 0) {
        mul_debug = (getenv("MUL_DEBUG")) ? 1 : 0;
    }
    
    if (!a || !b || !a->data || !b->data) {
        if (mul_debug) printf("MUL: null input\n");
        r->rows = 0;
        r->cols = 0;
        r->data = NULL;
        return;
    }
    
    if (a->cols != b->rows) {
        printf("Error: inner dimensions must match for multiplication\n");
        mat_zero(r, 1, 1);
        return;
    }
    
    if (mul_debug) printf("MUL: %dx%d * %dx%d, arena remaining: %lu\n", 
        a->rows, a->cols, b->rows, b->cols, (unsigned long)mat_arena_remaining());
    
    mat_zero(&tmp, a->rows, b->cols);
    if (!tmp.data) {
        if (mul_debug) printf("MUL: alloc failed for %dx%d\n", a->rows, b->cols);
        /* Arena exhausted */
        r->rows = 0;
        r->cols = 0;
        r->data = NULL;
        return;
    }
    
    for (i = 0; i < a->rows; i++) {
        for (j = 0; j < b->cols; j++) {
            for (k = 0; k < a->cols; k++) {
                apfc prod, sum;
                apfc_mul(&prod, &MAT_AT(a, i, k), &MAT_AT(b, k, j));
                apfc_add(&sum, &MAT_AT(&tmp, i, j), &prod);
                MAT_AT(&tmp, i, j) = sum;
            }
        }
    }
    
    mat_copy(r, &tmp);
}

/* ========== Transpose ========== */

void mat_transpose(matrix_t *r, const matrix_t *a)
{
    int i, j;
    int rows_a = a->rows;
    int cols_a = a->cols;
    apfc *src_data = a->data;
    
    /* Allocate result directly (cols_a x rows_a) */
    mat_zero(r, cols_a, rows_a);
    if (!r->data) return;
    
    /* Copy with transposition - read from source, not from a anymore */
    for (i = 0; i < rows_a; i++) {
        for (j = 0; j < cols_a; j++) {
            MAT_AT(r, j, i) = src_data[i * cols_a + j];
        }
    }
}

/* Conjugate transpose (Hermitian transpose) */
void mat_conj_transpose(matrix_t *r, const matrix_t *a)
{
    int i, j;
    matrix_t tmp = {0, 0, NULL};
    
    /* Allocate tmp from arena */
    mat_zero(&tmp, a->cols, a->rows);
    if (!tmp.data) return;
    
    for (i = 0; i < a->rows; i++) {
        for (j = 0; j < a->cols; j++) {
            apfc_conj(&MAT_AT(&tmp, j, i), &MAT_AT(a, i, j));
        }
    }
    
    mat_copy(r, &tmp);
}

/* ========== Matrix Power ========== */

int mat_pow(matrix_t *r, const matrix_t *a, long n)
{
    matrix_t result = {0, 0, NULL}, base = {0, 0, NULL}, tmp = {0, 0, NULL};
    
    if (!mat_is_square(a)) {
        printf("Error: matrix power requires square matrix\n");
        return 0;
    }
    
    if (n < 0) {
        /* A^(-n) = inv(A)^n */
        if (!mat_inv(&base, a)) {
            return 0;
        }
        n = -n;
    } else {
        mat_copy(&base, a);
    }
    
    /* Start with identity */
    mat_identity(&result, a->rows);
    
    /* Binary exponentiation */
    while (n > 0) {
        if (n & 1) {
            mat_mul(&tmp, &result, &base);
            mat_copy(&result, &tmp);
        }
        mat_mul(&tmp, &base, &base);
        mat_copy(&base, &tmp);
        n >>= 1;
    }
    
    mat_copy(r, &result);
    return 1;
}

/* ========== Diagonal Operations ========== */

/* Extract diagonal from matrix into column vector */
void mat_diag_extract(matrix_t *r, const matrix_t *m)
{
    int i, n;
    
    n = m->rows < m->cols ? m->rows : m->cols;
    mat_zero(r, n, 1);
    
    for (i = 0; i < n; i++) {
        MAT_AT(r, i, 0) = MAT_AT(m, i, i);
    }
}

/* Create diagonal matrix from vector */
void mat_diag_create(matrix_t *r, const matrix_t *v)
{
    int i, n;
    
    n = v->rows * v->cols;
    if (n > MAT_MAX_ROWS) n = MAT_MAX_ROWS;
    
    mat_zero(r, n, n);
    
    for (i = 0; i < n; i++) {
        MAT_AT(r, i, i) = v->data[i];
    }
}

/* ========== Trace ========== */

void mat_trace(apfc *r, const matrix_t *m)
{
    int i, n;
    apfc sum;
    
    if (!mat_is_square(m)) {
        printf("Error: trace requires square matrix\n");
        apf_zero(&r->re);
        apf_zero(&r->im);
        return;
    }
    
    n = m->rows;
    apf_zero(&sum.re);
    apf_zero(&sum.im);
    
    for (i = 0; i < n; i++) {
        apfc tmp;
        apfc_add(&tmp, &sum, &MAT_AT(m, i, i));
        sum = tmp;
    }
    
    *r = sum;
}

/* ========== Determinant ========== */

/* Recursive determinant using cofactor expansion */
static void det_recursive(apfc *r, const matrix_t *m, int n, int *cols, int depth)
{
    int i, j, sign;
    apfc sum, term, minor_det, elem;
    int subcols[MAT_MAX_COLS];
    
    if (depth == n - 1) {
        /* Base case: find the remaining column */
        for (i = 0; i < n; i++) {
            if (cols[i] >= 0) {
                *r = MAT_AT(m, depth, cols[i]);
                return;
            }
        }
        apf_zero(&r->re);
        apf_zero(&r->im);
        return;
    }
    
    apf_zero(&sum.re);
    apf_zero(&sum.im);
    sign = 1;
    
    for (i = 0; i < n; i++) {
        if (cols[i] < 0) continue;
        
        /* Copy cols, marking this one as used */
        for (j = 0; j < n; j++) subcols[j] = cols[j];
        subcols[i] = -1;
        
        elem = MAT_AT(m, depth, cols[i]);
        
        /* Skip if element is zero */
        if (apf_is_zero(&elem.re) && apf_is_zero(&elem.im)) {
            sign = -sign;
            continue;
        }
        
        det_recursive(&minor_det, m, n, subcols, depth + 1);
        apfc_mul(&term, &elem, &minor_det);
        
        if (sign > 0) {
            apfc_add(&sum, &sum, &term);
        } else {
            apfc_sub(&sum, &sum, &term);
        }
        sign = -sign;
    }
    
    *r = sum;
}

void mat_det(apfc *r, const matrix_t *m)
{
    int cols[MAT_MAX_COLS];
    int i, n;
    
    if (!mat_is_square(m)) {
        printf("Error: determinant requires square matrix\n");
        apf_zero(&r->re);
        apf_zero(&r->im);
        return;
    }
    
    n = m->rows;
    
    /* Special cases */
    if (n == 1) {
        *r = MAT_AT(m, 0, 0);
        return;
    }
    
    if (n == 2) {
        apfc a, b, c, d, ad, bc;
        a = MAT_AT(m, 0, 0);
        b = MAT_AT(m, 0, 1);
        c = MAT_AT(m, 1, 0);
        d = MAT_AT(m, 1, 1);
        apfc_mul(&ad, &a, &d);
        apfc_mul(&bc, &b, &c);
        apfc_sub(r, &ad, &bc);
        return;
    }
    
    /* General case */
    for (i = 0; i < n; i++) cols[i] = i;
    det_recursive(r, m, n, cols, 0);
}

/* ========== Inverse ========== */

/* Gauss-Jordan elimination */
int mat_inv(matrix_t *r, const matrix_t *m)
{
    matrix_t aug = {0, 0, NULL};
    int i, j, k, n;
    int pivot_row;
    apfc pivot, factor, temp;
    apf abs_pivot, abs_test, zero_thresh;
    
    if (!mat_is_square(m)) {
        printf("Error: inverse requires square matrix\n");
        return 0;
    }
    
    n = m->rows;
    
    /* Create augmented matrix [M | I] stored as n x 2n */
    mat_zero(&aug, n, 2 * n);
    if (!aug.data) return 0;
    
    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            MAT_AT(&aug, i, j) = MAT_AT(m, i, j);
        }
        for (j = n; j < 2 * n; j++) {
            apf_zero(&MAT_AT(&aug, i, j).re);
            apf_zero(&MAT_AT(&aug, i, j).im);
            if (j - n == i) {
                apf_from_int(&MAT_AT(&aug, i, j).re, 1);
            }
        }
    }
    
    /* Small threshold for detecting zero pivot - create a proper tiny value */
    apf_from_int(&zero_thresh, 1);
    /* After apf_from_int(1), exp is about -127 (normalized MSB at bit 127).
     * To make it tiny (e.g., 2^-100), we need exp = -100 - 127 = -227 */
    zero_thresh.exp = -200;  /* Very small: about 2^(-200+(-127)) = 2^-327 or so */
    
    /* Forward elimination with partial pivoting */
    for (k = 0; k < n; k++) {
        /* Find pivot */
        pivot_row = k;
        apfc_abs(&abs_pivot, &MAT_AT(&aug, k, k));
        
        for (i = k + 1; i < n; i++) {
            apfc_abs(&abs_test, &MAT_AT(&aug, i, k));
            if (apf_lt(&abs_pivot, &abs_test)) {
                abs_pivot = abs_test;
                pivot_row = i;
            }
        }
        
        /* Check for singular matrix */
        if (apf_lt(&abs_pivot, &zero_thresh)) {
            printf("Error: matrix is singular\n");
            return 0;
        }
        
        /* Swap rows if needed */
        if (pivot_row != k) {
            for (j = 0; j < 2 * n; j++) {
                temp = MAT_AT(&aug, k, j);
                MAT_AT(&aug, k, j) = MAT_AT(&aug, pivot_row, j);
                MAT_AT(&aug, pivot_row, j) = temp;
            }
        }
        
        /* Scale pivot row */
        pivot = MAT_AT(&aug, k, k);
        for (j = 0; j < 2 * n; j++) {
            apfc_div(&MAT_AT(&aug, k, j), &MAT_AT(&aug, k, j), &pivot);
        }
        
        /* Eliminate column */
        for (i = 0; i < n; i++) {
            if (i == k) continue;
            factor = MAT_AT(&aug, i, k);
            for (j = 0; j < 2 * n; j++) {
                apfc prod, diff;
                apfc_mul(&prod, &factor, &MAT_AT(&aug, k, j));
                apfc_sub(&diff, &MAT_AT(&aug, i, j), &prod);
                MAT_AT(&aug, i, j) = diff;
            }
        }
    }
    
    /* Extract inverse from right half */
    mat_zero(r, n, n);
    if (!r->data) return 0;
    
    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            MAT_AT(r, i, j) = MAT_AT(&aug, i, j + n);
        }
    }
    
    return 1;
}

/* ========== Eigenvalues (2x2 only) ========== */

int mat_eig2(apfc *eig1, apfc *eig2, const matrix_t *m)
{
    apfc a, b, c, d;
    apfc tr, det, disc, sqrt_disc;
    apfc two, four;
    apfc tmp1, tmp2;
    
    if (m->rows != 2 || m->cols != 2) {
        printf("Error: eig2 only works for 2x2 matrices\n");
        return 0;
    }
    
    a = MAT_AT(m, 0, 0);
    b = MAT_AT(m, 0, 1);
    c = MAT_AT(m, 1, 0);
    d = MAT_AT(m, 1, 1);
    
    /* trace = a + d */
    apfc_add(&tr, &a, &d);
    
    /* det = ad - bc */
    apfc_mul(&tmp1, &a, &d);
    apfc_mul(&tmp2, &b, &c);
    apfc_sub(&det, &tmp1, &tmp2);
    
    /* discriminant = tr^2 - 4*det */
    apfc_mul(&tmp1, &tr, &tr);
    apf_from_int(&four.re, 4);
    apf_zero(&four.im);
    apfc_mul(&tmp2, &four, &det);
    apfc_sub(&disc, &tmp1, &tmp2);
    
    /* sqrt(discriminant) */
    apfc_sqrt(&sqrt_disc, &disc);
    
    /* eigenvalues = (tr ± sqrt(disc)) / 2 */
    apf_from_int(&two.re, 2);
    apf_zero(&two.im);
    
    apfc_add(&tmp1, &tr, &sqrt_disc);
    apfc_div(eig1, &tmp1, &two);
    
    apfc_sub(&tmp1, &tr, &sqrt_disc);
    apfc_div(eig2, &tmp1, &two);
    
    return 2;
}

/* ========== Norms ========== */

void mat_norm_frobenius(apfc *r, const matrix_t *m)
{
    int i;
    apf sum, abs_val, sq;
    
    apf_zero(&sum);
    
    for (i = 0; i < m->rows * m->cols; i++) {
        apfc_abs(&abs_val, &m->data[i]);
        apf_mul(&sq, &abs_val, &abs_val);
        apf_add(&sum, &sum, &sq);
    }
    
    apf_sqrt(&r->re, &sum);
    apf_zero(&r->im);
}

/* ========== Output ========== */

void mat_print(const matrix_t *m)
{
    int i, j;
    static char buf[128];
    int digits = display_digits ? display_digits : 10;  /* Default 10 for matrix display */
    apf max_val, threshold, abs_val, eps;
    int *col_widths;
    
    if (m->rows == 0 || m->cols == 0) {
        printf("[]\n");
        return;
    }
    
    /* Allocate column widths array */
    col_widths = (int *)malloc(m->cols * sizeof(int));
    if (!col_widths) {
        /* Fallback to unaligned printing */
        for (i = 0; i < m->rows; i++) {
            for (j = 0; j < m->cols; j++) {
                apfc_to_str(buf, sizeof(buf), &MAT_AT(m, i, j), digits);
                printf("%s", buf);
                if (j < m->cols - 1) printf("  ");
            }
            printf("\n");
        }
        return;
    }
    
    /* Initialize column widths to 0 */
    for (j = 0; j < m->cols; j++) {
        col_widths[j] = 0;
    }
    
    /* Find maximum absolute value in matrix for threshold calculation */
    apf_zero(&max_val);
    for (i = 0; i < m->rows; i++) {
        for (j = 0; j < m->cols; j++) {
            apfc_abs(&abs_val, &MAT_AT(m, i, j));
            if (apf_cmp(&abs_val, &max_val) > 0) {
                max_val = abs_val;
            }
        }
    }
    
    /* threshold = max_val * eps where eps = 10^-14 */
    apf_from_int(&eps, 10);
    {
        apf exp;
        apf_from_int(&exp, -14);
        apfx_pow(&eps, &eps, &exp);
    }
    apf_mul(&threshold, &eps, &max_val);
    
    /* First pass: compute column widths */
    for (i = 0; i < m->rows; i++) {
        for (j = 0; j < m->cols; j++) {
            apfc val = MAT_AT(m, i, j);
            int len;
            
            /* Clean up near-zero values relative to matrix magnitude */
            apfc_abs(&abs_val, &val);
            if (apf_cmp(&abs_val, &threshold) < 0 && !apf_is_zero(&max_val)) {
                apf_zero(&val.re);
                apf_zero(&val.im);
            }
            
            apfc_to_str(buf, sizeof(buf), &val, digits);
            len = strlen(buf);
            if (len > col_widths[j]) {
                col_widths[j] = len;
            }
        }
    }
    
    /* Second pass: print with alignment */
    for (i = 0; i < m->rows; i++) {
        for (j = 0; j < m->cols; j++) {
            apfc val = MAT_AT(m, i, j);
            int len, pad;
            
            /* Clean up near-zero values relative to matrix magnitude */
            apfc_abs(&abs_val, &val);
            if (apf_cmp(&abs_val, &threshold) < 0 && !apf_is_zero(&max_val)) {
                apf_zero(&val.re);
                apf_zero(&val.im);
            }
            
            apfc_to_str(buf, sizeof(buf), &val, digits);
            len = strlen(buf);
            
            /* Right-justify within column */
            pad = col_widths[j] - len;
            while (pad > 0) {
                printf(" ");
                pad--;
            }
            printf("%s", buf);
            if (j < m->cols - 1) printf("  ");
        }
        printf("\n");
    }
    
    free(col_widths);
}

void mat_to_str(char *buf, int bufsize, const matrix_t *m)
{
    int i, j, pos = 0;
    static char elem[128];
    int len;
    int digits = display_digits ? display_digits : 10;
    apf max_val, threshold, abs_val, eps;
    
    /* Find maximum absolute value in matrix for threshold calculation */
    apf_zero(&max_val);
    for (i = 0; i < m->rows; i++) {
        for (j = 0; j < m->cols; j++) {
            apfc_abs(&abs_val, &MAT_AT(m, i, j));
            if (apf_cmp(&abs_val, &max_val) > 0) {
                max_val = abs_val;
            }
        }
    }
    
    /* threshold = max_val * eps where eps = 10^-14 */
    apf_from_int(&eps, 10);
    {
        apf exp;
        apf_from_int(&exp, -14);
        apfx_pow(&eps, &eps, &exp);
    }
    apf_mul(&threshold, &eps, &max_val);
    
    buf[0] = '[';
    pos = 1;
    
    for (i = 0; i < m->rows && pos < bufsize - 10; i++) {
        for (j = 0; j < m->cols && pos < bufsize - 10; j++) {
            apfc val = MAT_AT(m, i, j);
            
            /* Clean up near-zero values relative to matrix magnitude */
            apfc_abs(&abs_val, &val);
            if (apf_cmp(&abs_val, &threshold) < 0 && !apf_is_zero(&max_val)) {
                apf_zero(&val.re);
                apf_zero(&val.im);
            }
            
            apfc_to_str(elem, sizeof(elem), &val, digits);
            len = strlen(elem);
            if (pos + len < bufsize - 2) {
                strcpy(buf + pos, elem);
                pos += len;
            }
            if (j < m->cols - 1 && pos < bufsize - 2) {
                buf[pos++] = ' ';
            }
        }
        if (i < m->rows - 1 && pos < bufsize - 3) {
            buf[pos++] = ';';
            buf[pos++] = ' ';
        }
    }
    
    if (pos < bufsize - 1) {
        buf[pos++] = ']';
    }
    buf[pos] = '\0';
}
