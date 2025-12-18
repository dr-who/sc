/* stubs_dos.c - Stub functions for DOS builds
 * Provides empty implementations for features not included in DOS build
 * to satisfy linker without including full implementations
 */

#include "config.h"

#if defined(SCALC_MEDIUM) || defined(SCALC_TINY) || defined(SCALC_MINIMAL)

#include "apf.h"
#include "apfc.h"
#include "matrix.h"

/* ========== Matrix stubs (only if matrix not compiled) ========== */
#ifndef HAVE_MATRIX

void mat_init(matrix_t *m, int rows, int cols) { m->rows = rows; m->cols = cols; }
void mat_zero(matrix_t *m, int rows, int cols) { m->rows = rows; m->cols = cols; }
void mat_ones(matrix_t *m, int rows, int cols) { (void)m; (void)rows; (void)cols; }
void mat_identity(matrix_t *m, int n) { (void)m; (void)n; }
void mat_copy(matrix_t *dst, const matrix_t *src) { (void)dst; (void)src; }
int mat_is_square(const matrix_t *m) { (void)m; return 0; }
int mat_is_vector(const matrix_t *m) { (void)m; return 0; }
int mat_same_size(const matrix_t *a, const matrix_t *b) { (void)a; (void)b; return 0; }
int mat_length(const matrix_t *m) { (void)m; return 0; }
void mat_add(matrix_t *r, const matrix_t *a, const matrix_t *b) { (void)r; (void)a; (void)b; }
void mat_sub(matrix_t *r, const matrix_t *a, const matrix_t *b) { (void)r; (void)a; (void)b; }
void mat_neg(matrix_t *r, const matrix_t *a) { (void)r; (void)a; }
void mat_scale(matrix_t *r, const matrix_t *a, const apfc *s) { (void)r; (void)a; (void)s; }
void mat_add_scalar(matrix_t *r, const matrix_t *a, const apfc *s) { (void)r; (void)a; (void)s; }
void mat_elemwise_mul(matrix_t *r, const matrix_t *a, const matrix_t *b) { (void)r; (void)a; (void)b; }
void mat_elemwise_div(matrix_t *r, const matrix_t *a, const matrix_t *b) { (void)r; (void)a; (void)b; }
void mat_elemwise_pow(matrix_t *r, const matrix_t *a, const apfc *p) { (void)r; (void)a; (void)p; }
void mat_elemwise_pow_mat(matrix_t *r, const matrix_t *a, const matrix_t *b) { (void)r; (void)a; (void)b; }
void mat_mul(matrix_t *r, const matrix_t *a, const matrix_t *b) { (void)r; (void)a; (void)b; }
int mat_pow(matrix_t *r, const matrix_t *a, long n) { (void)r; (void)a; (void)n; return 0; }
void mat_transpose(matrix_t *r, const matrix_t *a) { (void)r; (void)a; }
void mat_conj_transpose(matrix_t *r, const matrix_t *a) { (void)r; (void)a; }
void mat_print(const matrix_t *m) { (void)m; }
void mat_det(apfc *result, const matrix_t *m) { (void)result; (void)m; }
void mat_trace(apfc *result, const matrix_t *m) { (void)result; (void)m; }
int mat_inv(matrix_t *r, const matrix_t *a) { (void)r; (void)a; return 0; }
int mat_mldivide(matrix_t *r, const matrix_t *a, const matrix_t *b) { (void)r; (void)a; (void)b; return 0; }
int mat_mrdivide(matrix_t *r, const matrix_t *a, const matrix_t *b) { (void)r; (void)a; (void)b; return 0; }
void mat_norm_frobenius(apfc *result, const matrix_t *m) { (void)result; (void)m; }
int mat_eig(matrix_t *eigenvalues, const matrix_t *m) { (void)eigenvalues; (void)m; return 0; }
int mat_eig2(apfc *eig1, apfc *eig2, const matrix_t *m) { (void)eig1; (void)eig2; (void)m; return 0; }
void mat_sum(apfc *result, const matrix_t *m) { (void)result; (void)m; }
void mat_mean(apfc *result, const matrix_t *m) { (void)result; (void)m; }
void mat_min(apfc *result, const matrix_t *m) { (void)result; (void)m; }
void mat_max(apfc *result, const matrix_t *m) { (void)result; (void)m; }
void mat_median(apfc *result, const matrix_t *m) { (void)result; (void)m; }
void mat_std(apfc *result, const matrix_t *m) { (void)result; (void)m; }
void mat_diag_create(matrix_t *r, const matrix_t *v) { (void)r; (void)v; }
void mat_diag_extract(matrix_t *r, const matrix_t *m) { (void)r; (void)m; }

#endif /* HAVE_MATRIX */

/* ========== Matrix random stubs (matrix_rand.c not compiled for DOS) ========== */
/* These are needed even when HAVE_MATRIX is defined because matrix_rand.c is excluded */

void mat_rand(matrix_t *m, int rows, int cols) { (void)m; (void)rows; (void)cols; }
void mat_randn(matrix_t *m, int rows, int cols) { (void)m; (void)rows; (void)cols; }
void mat_randi_range(matrix_t *m, int rows, int cols, long lo, long hi) { (void)m; (void)rows; (void)cols; (void)lo; (void)hi; }

/* ========== Random stubs (rand is in matrix_rand.c, not compiled for DOS) ========== */

void rand_auto_seed(void) { }
void rand_uniform(apfc *r) { apfc_zero(r); }
void rand_normal(apfc *r) { apfc_zero(r); }
void rand_int(apfc *r, long imax) { (void)imax; apfc_zero(r); }
void rand_uniform_range(apfc *r, const apfc *low, const apfc *high) { (void)low; (void)high; apfc_zero(r); }
void rand_int_range(apfc *r, long imin, long imax) { (void)imin; (void)imax; apfc_zero(r); }

/* ========== Plot stubs (plot uses double, not available on DOS) ========== */

int do_plot(const char *expr, char var, apf *xmin, apf *xmax) { (void)expr; (void)var; (void)xmin; (void)xmax; return 0; }
int do_textplot(const char *expr, char var, apf *xmin, apf *xmax) { (void)expr; (void)var; (void)xmin; (void)xmax; return 0; }

#endif /* SCALC_MEDIUM || SCALC_TINY || SCALC_MINIMAL */
