/* matrix.h - Matrix support for scalc
 * MATLAB-style matrices with arbitrary precision complex elements
 * C89 compliant for Watcom C / DOS
 * 16-bit clean: int is 16-bit, long is 32-bit
 *
 * Matrices use arena allocation - data is allocated from a static pool
 * that gets reset at the start of each REPL command. This allows large
 * matrices (like iris 150x4) without stack overflow.
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

/* ========== Matrix Arena ========== */
/* Default arena sizes - can be changed at runtime via mat_memory_init() */
#if defined(PLATFORM_DOS) || defined(SCALC_MEDIUM)
  #define MAT_ARENA_DEFAULT  (32 * 1024)       /* 32KB for DOS */
  #define MAT_PERSIST_DEFAULT (16 * 1024)      /* 16KB persistent */
#else
  #define MAT_ARENA_DEFAULT  (1024UL * 1024 * 1024)  /* 1GB for Linux */
  #define MAT_PERSIST_DEFAULT (1024UL * 1024 * 1024) /* 1GB persistent */
#endif

/* Initialize memory pools (call once at startup, 0 = use defaults) */
int mat_memory_init(size_t arena_bytes, size_t persist_bytes);
void mat_memory_free(void);

/* Arena functions */
void mat_arena_reset(void);              /* Reset at start of each command */
void *mat_arena_alloc(size_t bytes);     /* Allocate from arena */
size_t mat_arena_remaining(void);        /* Bytes remaining */
size_t mat_arena_size(void);             /* Total arena size */

void mat_persist_reset(void);            /* Clear all persistent storage */
void *mat_persist_alloc(size_t bytes);   /* Allocate persistent storage */
void mat_persist_free(void *ptr);        /* Mark block as free (simple) */
size_t mat_persist_remaining(void);      /* Bytes remaining */
size_t mat_persist_size(void);           /* Total persist size */

/* Maximum dimensions (for validation only, not static allocation) */
#if defined(PLATFORM_DOS) || defined(SCALC_MEDIUM)
  #define MAT_MAX_ROWS 64
  #define MAT_MAX_COLS 64
  #define MAT_MAX_ELEM 256   /* Max elements for stack-based temporaries */
#else
  #define MAT_MAX_ROWS 500000
  #define MAT_MAX_COLS 1000
  #define MAT_MAX_ELEM 4096  /* Max elements for stack-based temporaries */
#endif

/* Matrix structure - data allocated from arena */
typedef struct {
    int rows;
    int cols;
    apfc *data;  /* Pointer to arena-allocated storage */
} matrix_t;

/* Initialize matrix to empty state */
#define MAT_INIT_EMPTY(m) do { (m).rows = 0; (m).cols = 0; (m).data = NULL; } while(0)

/* Ensure matrix has space allocated - use before assigning elements */
#define MAT_ENSURE(m, r, c) do { \
    if (!(m)->data || (m)->rows != (r) || (m)->cols != (c)) { \
        mat_zero((m), (r), (c)); \
    } \
} while(0)

/* Access element at (row, col) - 0-indexed */
#define MAT_AT(m, r, c) ((m)->data[(r) * (m)->cols + (c)])

/* Copy matrix to persistent storage (for named variables) */
int mat_copy_persist(matrix_t *dst, const matrix_t *src);
int mat_zero_persist(matrix_t *m, int rows, int cols);

/* Resize/allocate matrix from arena - safe to call on existing matrix */
int mat_resize(matrix_t *m, int rows, int cols);

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

/* ========== Column-wise Statistics (MATLAB compatible) ========== */
void mat_mean_cols(matrix_t *r, const matrix_t *m);
void mat_std_cols(matrix_t *r, const matrix_t *m);
void mat_cov(matrix_t *r, const matrix_t *m);
void mat_corrcoef(matrix_t *r, const matrix_t *m);
void mat_crosstab(matrix_t *r, const matrix_t *a, const matrix_t *b);
void mat_randindex(apfc *r, const matrix_t *a, const matrix_t *b);

/* ========== QR Decomposition (matrix_linalg.c) ========== */
int mat_qr(matrix_t *Q, matrix_t *R, const matrix_t *A);

/* ========== Eigenvalues (matrix_linalg.c) ========== */
int mat_eig(matrix_t *eigenvalues, const matrix_t *A);

/* ========== Cholesky Decomposition (matrix_linalg.c) ========== */
int mat_chol(matrix_t *L, const matrix_t *A);

/* ========== Full SVD (matrix_linalg.c) ========== */
int mat_svd(matrix_t *U, matrix_t *S, matrix_t *V, const matrix_t *A);

/* ========== Null Space (matrix_linalg.c) ========== */
int mat_null(matrix_t *N, const matrix_t *A);

/* ========== Schur Decomposition (matrix_linalg.c) ========== */
int mat_schur(matrix_t *Q, matrix_t *T, const matrix_t *A);

/* ========== PCA and Clustering (matrix_linalg.c) ========== */
int mat_pca(matrix_t *coeff, matrix_t *score, matrix_t *latent, const matrix_t *X);
int mat_pca_reduce(matrix_t *result, const matrix_t *X, int k);
int mat_kmeans(matrix_t *idx, matrix_t *centroids, const matrix_t *X, int k);
int mat_pdist(matrix_t *D, const matrix_t *X);
void mat_silhouette(apfc *score, const matrix_t *X, const matrix_t *idx);

/* ========== Concatenation (matrix_ops.c) ========== */
int mat_hcat(matrix_t *r, const matrix_t *a, const matrix_t *b);
int mat_vcat(matrix_t *r, const matrix_t *a, const matrix_t *b);

/* ========== Indexing/Slicing (matrix_ops.c) ========== */
void mat_get_row(matrix_t *r, const matrix_t *m, int row);
void mat_get_col(matrix_t *r, const matrix_t *m, int col);
void mat_set_row(matrix_t *m, int row, const matrix_t *v);
void mat_set_col(matrix_t *m, int col, const matrix_t *v);
void mat_submat(matrix_t *r, const matrix_t *m, int r1, int r2, int c1, int c2);

/* Array manipulation */
void mat_fliplr(matrix_t *r, const matrix_t *m);
void mat_flipud(matrix_t *r, const matrix_t *m);
void mat_flip(matrix_t *r, const matrix_t *m);
void mat_rot90(matrix_t *r, const matrix_t *m);
void mat_triu(matrix_t *r, const matrix_t *m);
void mat_tril(matrix_t *r, const matrix_t *m);
void mat_cumsum(matrix_t *r, const matrix_t *m);
void mat_cumprod(matrix_t *r, const matrix_t *m);
void mat_diff(matrix_t *r, const matrix_t *m);
void mat_normalize(matrix_t *r, const matrix_t *m);
void mat_sort(matrix_t *r, const matrix_t *m);
void mat_squeeze(matrix_t *r, const matrix_t *m);
void mat_cummax(matrix_t *r, const matrix_t *m);
void mat_cummin(matrix_t *r, const matrix_t *m);
void mat_unique(matrix_t *r, const matrix_t *m);
void mat_find(matrix_t *r, const matrix_t *m);
void mat_sortrows(matrix_t *r, const matrix_t *m);
void mat_rescale(matrix_t *r, const matrix_t *m);

/* Statistics reduce functions */
void mat_range(apfc *r, const matrix_t *m);
void mat_rms(apfc *r, const matrix_t *m);
void mat_geomean(apfc *r, const matrix_t *m);
void mat_harmmean(apfc *r, const matrix_t *m);
void mat_sumsq(apfc *r, const matrix_t *m);
void mat_meansq(apfc *r, const matrix_t *m);

/* Boolean/query reduce functions */
void mat_any(apfc *r, const matrix_t *m);
void mat_all(apfc *r, const matrix_t *m);
void mat_isempty(apfc *r, const matrix_t *m);
void mat_numel(apfc *r, const matrix_t *m);
void mat_isscalar(apfc *r, const matrix_t *m);
void mat_isvector(apfc *r, const matrix_t *m);
void mat_isrow(apfc *r, const matrix_t *m);
void mat_iscolumn(apfc *r, const matrix_t *m);
void mat_ismatrix(apfc *r, const matrix_t *m);
void mat_issquare(apfc *r, const matrix_t *m);
void mat_skewness(apfc *r, const matrix_t *m);
void mat_kurtosis(apfc *r, const matrix_t *m);
void mat_nnz(apfc *r, const matrix_t *m);
void mat_mode(apfc *r, const matrix_t *m);
void mat_iqr(apfc *r, const matrix_t *m);

/* Cumulative functions (matrix -> matrix) */
void mat_cumsum(matrix_t *r, const matrix_t *m);
void mat_cumprod(matrix_t *r, const matrix_t *m);
void mat_cummax(matrix_t *r, const matrix_t *m);
void mat_cummin(matrix_t *r, const matrix_t *m);

/* Matrix property functions (reduce versions) */
void mat_numel_r(apfc *r, const matrix_t *m);
void mat_length_r(apfc *r, const matrix_t *m);
void mat_ndims_r(apfc *r, const matrix_t *m);
void mat_nrows_r(apfc *r, const matrix_t *m);
void mat_ncols_r(apfc *r, const matrix_t *m);
void mat_width_r(apfc *r, const matrix_t *m);
void mat_height_r(apfc *r, const matrix_t *m);

/* Higher moment statistics */
void mat_mad(apfc *r, const matrix_t *m);

/* Matrix transform functions */
void mat_center(matrix_t *r, const matrix_t *m);
void mat_gradient(matrix_t *r, const matrix_t *m);
void mat_diff2(matrix_t *r, const matrix_t *m);
void mat_detrend(matrix_t *r, const matrix_t *m);

/* ========== MATLAB Compatibility (matrix_linalg.c) ========== */
void mat_rref(matrix_t *r, const matrix_t *a);       /* Row reduced echelon form */
void mat_orth(matrix_t *r, const matrix_t *a);       /* Orthonormal basis */
void mat_poly(matrix_t *r, const matrix_t *a);       /* Characteristic polynomial */
void mat_roots(matrix_t *r, const matrix_t *coeffs); /* Polynomial roots */
void mat_expm(matrix_t *r, const matrix_t *a);       /* Matrix exponential */
void mat_logm(matrix_t *r, const matrix_t *a);       /* Matrix logarithm */
void mat_fft(matrix_t *r, const matrix_t *a);        /* Fast Fourier Transform */
void mat_ifft(matrix_t *r, const matrix_t *a);       /* Inverse FFT */
void mat_meshgrid(matrix_t *X, matrix_t *Y, const matrix_t *x, const matrix_t *y);

/* Data analysis */
void mat_histcounts(matrix_t *r, const matrix_t *data, int nbins);
void mat_xcorr(matrix_t *r, const matrix_t *x, const matrix_t *y);
void mat_xcov(matrix_t *r, const matrix_t *x, const matrix_t *y);
void mat_isoutlier(matrix_t *r, const matrix_t *data);
void mat_simpson(matrix_t *r, const matrix_t *y, double h);
void mat_hess(matrix_t *H, const matrix_t *A);
void mat_balance(matrix_t *B, const matrix_t *A);

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
