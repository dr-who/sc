/* matrix.h - Matrix support for scalc
 * MATLAB-style matrices with arbitrary precision complex elements
 * C89 compliant for Watcom C / DOS
 * 16-bit clean: int is 16-bit, long is 32-bit
 */
#ifndef MATRIX_H
#define MATRIX_H

#include "config.h"
#include "apfc.h"

/* Far memory macro for DOS - puts large arrays outside 64KB DGROUP */
#ifndef SC_FAR
  #ifdef __WATCOMC__
    #define SC_FAR __far
  #else
    #define SC_FAR
  #endif
#endif

/* Maximum matrix dimensions */
#ifndef MAT_MAX_ROWS
  #if defined(PLATFORM_DOS) || defined(SCALC_MEDIUM)
    #define MAT_MAX_ROWS 4
  #else
    #define MAT_MAX_ROWS 10
  #endif
#endif
#ifndef MAT_MAX_COLS
  #if defined(PLATFORM_DOS) || defined(SCALC_MEDIUM)
    #define MAT_MAX_COLS 4
  #else
    #define MAT_MAX_COLS 10
  #endif
#endif
#define MAT_MAX_ELEM (MAT_MAX_ROWS * MAT_MAX_COLS)

/* Matrix structure */
typedef struct {
    int rows;
    int cols;
    apfc data[MAT_MAX_ELEM];  /* Row-major storage */
} matrix_t;

/* Access element at (row, col) - 0-indexed */
#define MAT_AT(m, r, c) ((m)->data[(r) * (m)->cols + (c)])

/* ========== Initialization (matrix.c) ========== */
void mat_init(matrix_t *m, int rows, int cols);
void mat_zero(matrix_t *m, int rows, int cols);
void mat_ones(matrix_t *m, int rows, int cols);
void mat_identity(matrix_t *m, int n);
void mat_copy(matrix_t *dst, const matrix_t *src);

/* ========== Basic info (matrix.c) ========== */
int mat_is_square(const matrix_t *m);
int mat_is_vector(const matrix_t *m);
int mat_same_size(const matrix_t *a, const matrix_t *b);
int mat_length(const matrix_t *m);

/* ========== Element-wise operations (matrix.c, matrix_ops.c) ========== */
void mat_add(matrix_t *r, const matrix_t *a, const matrix_t *b);
void mat_sub(matrix_t *r, const matrix_t *a, const matrix_t *b);
void mat_neg(matrix_t *r, const matrix_t *a);
void mat_scale(matrix_t *r, const matrix_t *a, const apfc *s);
void mat_add_scalar(matrix_t *r, const matrix_t *a, const apfc *s);
void mat_elemwise_mul(matrix_t *r, const matrix_t *a, const matrix_t *b);
void mat_elemwise_div(matrix_t *r, const matrix_t *a, const matrix_t *b);
void mat_elemwise_pow(matrix_t *r, const matrix_t *a, const apfc *p);
void mat_elemwise_pow_mat(matrix_t *r, const matrix_t *a, const matrix_t *b);

/* ========== Matrix multiplication (matrix.c) ========== */
void mat_mul(matrix_t *r, const matrix_t *a, const matrix_t *b);

/* ========== Matrix power (matrix.c) ========== */
int mat_pow(matrix_t *r, const matrix_t *a, long n);

/* ========== Transpose (matrix.c) ========== */
void mat_transpose(matrix_t *r, const matrix_t *a);
void mat_conj_transpose(matrix_t *r, const matrix_t *a);

/* ========== Diagonal (matrix.c) ========== */
void mat_diag_extract(matrix_t *r, const matrix_t *m);
void mat_diag_create(matrix_t *r, const matrix_t *v);

/* ========== Determinant (matrix.c) ========== */
void mat_det(apfc *r, const matrix_t *m);

/* ========== Trace (matrix.c) ========== */
void mat_trace(apfc *r, const matrix_t *m);

/* ========== Inverse (matrix.c) ========== */
int mat_inv(matrix_t *r, const matrix_t *m);

/* ========== Linear algebra (matrix_linalg.c) ========== */
int mat_lu(matrix_t *L, matrix_t *U, int *perm, const matrix_t *m);
int mat_solve(matrix_t *x, const matrix_t *A, const matrix_t *b);
int mat_mldivide(matrix_t *r, const matrix_t *A, const matrix_t *B);
int mat_mrdivide(matrix_t *r, const matrix_t *A, const matrix_t *B);
int mat_svd2(apfc *s1, apfc *s2, const matrix_t *m);

/* ========== Eigenvalues (matrix.c) ========== */
int mat_eig2(apfc *eig1, apfc *eig2, const matrix_t *m);

/* ========== Norms (matrix.c, matrix_linalg.c) ========== */
void mat_norm_frobenius(apfc *r, const matrix_t *m);
void mat_norm_1(apfc *r, const matrix_t *m);
void mat_norm_inf(apfc *r, const matrix_t *m);

/* ========== Aggregations (matrix_linalg.c) ========== */
void mat_sum(apfc *r, const matrix_t *m);
void mat_mean(apfc *r, const matrix_t *m);
void mat_min(apfc *r, const matrix_t *m);
void mat_max(apfc *r, const matrix_t *m);
void mat_prod(apfc *r, const matrix_t *m);
void mat_std(apfc *r, const matrix_t *m);
void mat_var(apfc *r, const matrix_t *m);
void mat_median(apfc *r, const matrix_t *m);
void mat_sum_rows(matrix_t *r, const matrix_t *m);
void mat_sum_cols(matrix_t *r, const matrix_t *m);

/* ========== QR Decomposition (matrix_linalg.c) ========== */
int mat_qr(matrix_t *Q, matrix_t *R, const matrix_t *A);

/* ========== Eigenvalues (matrix_linalg.c) ========== */
int mat_eig(matrix_t *eigenvalues, const matrix_t *A);

/* ========== Concatenation (matrix_ops.c) ========== */
int mat_hcat(matrix_t *r, const matrix_t *a, const matrix_t *b);
int mat_vcat(matrix_t *r, const matrix_t *a, const matrix_t *b);

/* ========== Indexing/Slicing (matrix_ops.c) ========== */
void mat_get_row(matrix_t *r, const matrix_t *m, int row);
void mat_get_col(matrix_t *r, const matrix_t *m, int col);
void mat_set_row(matrix_t *m, int row, const matrix_t *v);
void mat_set_col(matrix_t *m, int col, const matrix_t *v);
void mat_submat(matrix_t *r, const matrix_t *m, int r1, int r2, int c1, int c2);

/* ========== Random (matrix_rand.c) ========== */
void mat_rand_seed(unsigned long seed);
void rand_auto_seed(void);

/* Scalar random */
void rand_uniform(apfc *r);
void rand_uniform_range(apfc *r, const apfc *low, const apfc *high);
void rand_normal(apfc *r);
void rand_int(apfc *r, long imax);
void rand_int_range(apfc *r, long imin, long imax);

/* Matrix random */
void mat_rand(matrix_t *r, int rows, int cols);
void mat_rand_range(matrix_t *r, int rows, int cols,
                    const apfc *low, const apfc *high);
void mat_randn(matrix_t *r, int rows, int cols);
void mat_randi(matrix_t *r, int rows, int cols, long imax);
void mat_randi_range(matrix_t *r, int rows, int cols, long imin, long imax);

/* ========== Output (matrix.c) ========== */
void mat_print(const matrix_t *m);
void mat_to_str(char *buf, int bufsize, const matrix_t *m);

#endif /* MATRIX_H */
