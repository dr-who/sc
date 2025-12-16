/*
 * func_linalg.h - Linear Algebra function declarations
 * C89 compliant
 */
#ifndef FUNC_LINALG_H
#define FUNC_LINALG_H

#include "func_registry.h"
#include "matrix.h"

/* Category export */
extern const func_category_t func_linalg_category;

/* Function declarations */
int func_svd(matrix_t *U, matrix_t *S, matrix_t *V, const matrix_t *A);
int func_qr(matrix_t *Q, matrix_t *R, const matrix_t *A);
int func_lu(matrix_t *L, matrix_t *U, int *perm, const matrix_t *A);
int func_schur(matrix_t *Q, matrix_t *T, const matrix_t *A);
int func_eig(matrix_t *result, const matrix_t *A);
int func_chol(matrix_t *L, const matrix_t *A);
int func_det(apfc *result, const matrix_t *A);
int func_inv(matrix_t *result, const matrix_t *A);
int func_trace(apfc *result, const matrix_t *A);
int func_rank(int *result, const matrix_t *A);
int func_cond(apfc *result, const matrix_t *A);
int func_norm(apfc *result, const matrix_t *A);
int func_null(matrix_t *result, const matrix_t *A);
int func_linsolve(matrix_t *result, const matrix_t *A, const matrix_t *b);
int func_trans(matrix_t *result, const matrix_t *A);

#endif /* FUNC_LINALG_H */
