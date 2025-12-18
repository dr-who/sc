/*
 * decomp_funcs.c - Table-driven decomposition function dispatch
 * 
 * C89 compliant for Watcom C / DOS
 */
#include <stdio.h>
#include <string.h>
#include "decomp_funcs.h"
#include "matrix.h"

/* Hash table for decomposition functions - sized for ~80 entries */
#define DF_HASH_SIZE 61
static DecompFuncEntry *df_hash[DF_HASH_SIZE];

/* Static storage - room for 8x growth */
#define DF_MAX_ENTRIES 96
static DecompFuncEntry df_entries[DF_MAX_ENTRIES];
static int df_count = 0;

/* Case-insensitive string comparison */
static int df_str_eq(const char *a, const char *b) {
    while (*a && *b) {
        char ca = *a, cb = *b;
        if (ca >= 'A' && ca <= 'Z') ca += 32;
        if (cb >= 'A' && cb <= 'Z') cb += 32;
        if (ca != cb) return 0;
        a++; b++;
    }
    return *a == *b;
}

/* Register a decomposition function */
static void reg_decomp(const char *name, int min_out, int max_out, 
                       int min_in, int max_in) {
    DecompFuncEntry *e;
    if (df_count >= DF_MAX_ENTRIES) return;
    e = &df_entries[df_count++];
    e->name = name;
    e->min_outputs = min_out;
    e->max_outputs = max_out;
    e->min_inputs = min_in;
    e->max_inputs = max_in;
}

void decomp_funcs_init(void) {
    int i;
    for (i = 0; i < DF_HASH_SIZE; i++) df_hash[i] = NULL;
    df_count = 0;
    
    /* SVD: 1 output = singular values, 3 outputs = [U,S,V] */
    reg_decomp("svd", 1, 3, 1, 1);
    
    /* QR: 1 output = R, 2 outputs = [Q,R] */
    reg_decomp("qr", 1, 2, 1, 1);
    
    /* LU: 2 outputs = [L,U], 3 outputs = [L,U,P] */
    reg_decomp("lu", 2, 3, 1, 1);
    
    /* Cholesky: 1 output = L */
    reg_decomp("chol", 1, 1, 1, 1);
    reg_decomp("cholesky", 1, 1, 1, 1);
    
    /* Eigenvalues: 1 output = eigenvalues, 2 outputs = [V,D] */
    reg_decomp("eig", 1, 2, 1, 1);
    reg_decomp("eigenvalues", 1, 2, 1, 1);
    
    /* PCA: 1 output = scores, 3 outputs = [coeff,score,latent] */
    reg_decomp("pca", 1, 3, 1, 1);
    
    /* Schur: 1 output = T, 2 outputs = [Q,T] */
    reg_decomp("schur", 1, 2, 1, 1);
    
    /* K-means: 2 outputs = [idx, centroids] */
    reg_decomp("kmeans", 2, 2, 2, 2);
}

DecompFuncEntry *decomp_func_lookup(const char *name) {
    int i;
    for (i = 0; i < df_count; i++) {
        if (df_str_eq(df_entries[i].name, name)) {
            return &df_entries[i];
        }
    }
    return NULL;
}

/* Execute SVD decomposition */
static int exec_svd(int nout, DecompResult *result, const matrix_t *A) {
    matrix_t U = {0,0,NULL}, S = {0,0,NULL}, V = {0,0,NULL};
    
    if (!mat_svd(&U, &S, &V, A)) {
        return 0;
    }
    
    if (nout == 1) {
        /* Return singular values as column vector */
        int i, min_dim = (S.rows < S.cols) ? S.rows : S.cols;
        mat_zero(&result->outputs[0], min_dim, 1);
        for (i = 0; i < min_dim; i++) {
            MAT_AT(&result->outputs[0], i, 0) = MAT_AT(&S, i, i);
        }
        result->num_outputs = 1;
    } else {
        /* Return [U, S, V] */
        mat_copy(&result->outputs[0], &U);
        mat_copy(&result->outputs[1], &S);
        mat_copy(&result->outputs[2], &V);
        result->num_outputs = 3;
    }
    return 1;
}

/* Execute QR decomposition */
static int exec_qr(int nout, DecompResult *result, const matrix_t *A) {
    matrix_t Q = {0,0,NULL}, R = {0,0,NULL};
    
    if (!mat_qr(&Q, &R, A)) {
        return 0;
    }
    
    if (nout == 1) {
        /* Return R only */
        mat_copy(&result->outputs[0], &R);
        result->num_outputs = 1;
    } else {
        /* Return [Q, R] */
        mat_copy(&result->outputs[0], &Q);
        mat_copy(&result->outputs[1], &R);
        result->num_outputs = 2;
    }
    return 1;
}

/* Execute LU decomposition */
static int exec_lu(int nout, DecompResult *result, const matrix_t *A) {
    matrix_t L = {0,0,NULL}, U = {0,0,NULL};
    int perm[MAT_MAX_ROWS];
    int i, n;
    
    if (!mat_lu(&L, &U, perm, A)) {
        return 0;
    }
    
    mat_copy(&result->outputs[0], &L);
    mat_copy(&result->outputs[1], &U);
    
    if (nout >= 3) {
        /* Build permutation matrix P */
        n = A->rows;
        mat_zero(&result->outputs[2], n, n);
        for (i = 0; i < n; i++) {
            apf_from_int(&MAT_AT(&result->outputs[2], i, perm[i]).re, 1);
        }
        result->num_outputs = 3;
    } else {
        result->num_outputs = 2;
    }
    return 1;
}

/* Execute Cholesky decomposition */
static int exec_chol(int nout, DecompResult *result, const matrix_t *A) {
    (void)nout;  /* Always 1 output */
    if (!mat_chol(&result->outputs[0], A)) {
        return 0;
    }
    result->num_outputs = 1;
    return 1;
}

/* Execute eigenvalue decomposition */
static int exec_eig(int nout, DecompResult *result, const matrix_t *A) {
    if (nout == 1) {
        /* Return eigenvalues only */
        if (!mat_eig(&result->outputs[0], A)) {
            return 0;
        }
        result->num_outputs = 1;
    } else {
        /* Return [V, D] - eigenvectors and diagonal eigenvalues */
        /* Note: full eigenvector computation not implemented, fall back */
        if (!mat_eig(&result->outputs[0], A)) {
            return 0;
        }
        /* For now, just return eigenvalues */
        result->num_outputs = 1;
    }
    return 1;
}

/* Execute PCA */
static int exec_pca(int nout, DecompResult *result, const matrix_t *A) {
    matrix_t coeff = {0,0,NULL}, score = {0,0,NULL}, latent = {0,0,NULL};
    
    if (!mat_pca(&coeff, &score, &latent, A)) {
        return 0;
    }
    
    if (nout == 1) {
        /* Return scores only */
        mat_copy(&result->outputs[0], &score);
        result->num_outputs = 1;
    } else {
        /* Return [coeff, score, latent] */
        mat_copy(&result->outputs[0], &coeff);
        mat_copy(&result->outputs[1], &score);
        mat_copy(&result->outputs[2], &latent);
        result->num_outputs = 3;
    }
    return 1;
}

/* Execute Schur decomposition */
static int exec_schur(int nout, DecompResult *result, const matrix_t *A) {
    matrix_t Q = {0,0,NULL}, T = {0,0,NULL};
    
    if (!mat_schur(&Q, &T, A)) {
        return 0;
    }
    
    if (nout == 1) {
        /* Return T only */
        mat_copy(&result->outputs[0], &T);
        result->num_outputs = 1;
    } else {
        /* Return [Q, T] */
        mat_copy(&result->outputs[0], &Q);
        mat_copy(&result->outputs[1], &T);
        result->num_outputs = 2;
    }
    return 1;
}

/* Execute k-means clustering */
static int exec_kmeans(int nout, DecompResult *result, 
                       const matrix_t *X, const matrix_t *k_mat) {
    matrix_t idx = {0,0,NULL}, centroids = {0,0,NULL};
    int k;
    
    (void)nout;  /* Always 2 outputs */
    
    /* Extract k from k_mat (should be 1x1 scalar) */
    if (k_mat->rows == 1 && k_mat->cols == 1) {
        k = apf_to_long(&MAT_AT(k_mat, 0, 0).re);
    } else {
        k = 2;  /* Default */
    }
    
    if (!mat_kmeans(&idx, &centroids, X, k)) {
        return 0;
    }
    
    mat_copy(&result->outputs[0], &idx);
    mat_copy(&result->outputs[1], &centroids);
    result->num_outputs = 2;
    return 1;
}

int decomp_func_exec(const char *name, int nout, 
                     DecompResult *result, 
                     const matrix_t *input, 
                     const matrix_t *input2) {
    
    result->num_outputs = 0;
    
    if (df_str_eq(name, "svd")) {
        return exec_svd(nout, result, input);
    } else if (df_str_eq(name, "qr")) {
        return exec_qr(nout, result, input);
    } else if (df_str_eq(name, "lu")) {
        return exec_lu(nout, result, input);
    } else if (df_str_eq(name, "chol") || df_str_eq(name, "cholesky")) {
        return exec_chol(nout, result, input);
    } else if (df_str_eq(name, "eig") || df_str_eq(name, "eigenvalues")) {
        return exec_eig(nout, result, input);
    } else if (df_str_eq(name, "pca")) {
        return exec_pca(nout, result, input);
    } else if (df_str_eq(name, "schur")) {
        return exec_schur(nout, result, input);
    } else if (df_str_eq(name, "kmeans")) {
        if (!input2) {
            printf("Error: kmeans requires two arguments\n");
            return 0;
        }
        return exec_kmeans(nout, result, input, input2);
    }
    
    printf("Error: unknown decomposition function '%s'\n", name);
    return 0;
}
