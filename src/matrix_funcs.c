/*
 * matrix_funcs.c - Table-driven matrix function dispatch
 * 
 * C89 compliant for Watcom C / DOS
 */
#include <stdio.h>
#include <string.h>
#include "matrix_funcs.h"
#include "matrix.h"
#include "ml.h"
#include "saas.h"

/* Hash table for O(1) lookup - sized for ~700 entries */
#define MF_HASH_SIZE 251
static MatFuncEntry *mf_hash[MF_HASH_SIZE];

/* Static storage for entries - room for 8x growth */
#define MF_MAX_ENTRIES 768
static MatFuncEntry mf_entries[MF_MAX_ENTRIES];
static int mf_count = 0;

/* Simple string hash */
static unsigned mf_hash_name(const char *name) {
    unsigned h = 0;
    while (*name) {
        char c = *name++;
        if (c >= 'A' && c <= 'Z') c += 32;
        h = h * 31 + (unsigned char)c;
    }
    return h % MF_HASH_SIZE;
}

/* Case-insensitive compare */
static int mf_streq(const char *a, const char *b) {
    while (*a && *b) {
        char ca = *a++, cb = *b++;
        if (ca >= 'A' && ca <= 'Z') ca += 32;
        if (cb >= 'A' && cb <= 'Z') cb += 32;
        if (ca != cb) return 0;
    }
    return *a == *b;
}

/* Register a matrix->matrix function (no error return) */
static void reg_mat1(const char *name, mat_unary_fn fn) {
    unsigned h;
    MatFuncEntry *e;
    if (mf_count >= MF_MAX_ENTRIES) return;
    e = &mf_entries[mf_count++];
    e->name = name;
    e->type = MF_MAT_1;
    e->fn.m1 = fn;
    e->can_fail = 0;
    h = mf_hash_name(name);
    e->next = mf_hash[h];
    mf_hash[h] = e;
}

/* Register a matrix->matrix function (with error return) */
static void reg_mat1_err(const char *name, mat_unary_fn_err fn) {
    unsigned h;
    MatFuncEntry *e;
    if (mf_count >= MF_MAX_ENTRIES) return;
    e = &mf_entries[mf_count++];
    e->name = name;
    e->type = MF_MAT_1;
    e->fn.m1_err = fn;
    e->can_fail = 1;
    h = mf_hash_name(name);
    e->next = mf_hash[h];
    mf_hash[h] = e;
}

/* Register a matrix->scalar reduction function */
static void reg_reduce(const char *name, mat_reduce_fn fn) {
    unsigned h;
    MatFuncEntry *e;
    if (mf_count >= MF_MAX_ENTRIES) return;
    e = &mf_entries[mf_count++];
    e->name = name;
    e->type = MF_REDUCE_1;
    e->fn.red = fn;
    e->can_fail = 0;
    h = mf_hash_name(name);
    e->next = mf_hash[h];
    mf_hash[h] = e;
}

void matrix_funcs_init(void) {
    int i;
    for (i = 0; i < MF_HASH_SIZE; i++) mf_hash[i] = NULL;
    mf_count = 0;
    
    /* === MATRIX TRANSFORMS === */
    reg_mat1("transpose", mat_transpose);
    reg_mat1("trans", mat_transpose);
    reg_mat1("ctranspose", mat_conj_transpose);
    reg_mat1("neg", mat_neg);
    reg_mat1("cov", mat_cov);
    reg_mat1("sum_rows", mat_sum_rows);
    reg_mat1("sum_cols", mat_sum_cols);
    reg_mat1("mean_cols", mat_mean_cols);
    reg_mat1("std_cols", mat_std_cols);
    /* corrcoef has 2-arg variant, handled in parser */
    
    /* === LINEAR ALGEBRA (can fail) === */
    reg_mat1_err("inv", mat_inv);
    reg_mat1_err("inverse", mat_inv);
    reg_mat1_err("chol", mat_chol);
    reg_mat1_err("cholesky", mat_chol);
    reg_mat1_err("null", mat_null);
    reg_mat1_err("nullspace", mat_null);
    
    /* Array manipulation */
    reg_mat1("fliplr", mat_fliplr);
    reg_mat1("flipud", mat_flipud);
    reg_mat1("flip", mat_flip);
    reg_mat1("rot90", mat_rot90);
    reg_mat1("triu", mat_triu);
    reg_mat1("tril", mat_tril);
    reg_mat1("cumsum", mat_cumsum);
    reg_mat1("cumprod", mat_cumprod);
    reg_mat1("diff", mat_diff);
    reg_mat1("normalize", mat_normalize);
    reg_mat1("center", mat_center);
    reg_mat1("gradient", mat_gradient);
    reg_mat1("diff2", mat_diff2);
    reg_mat1("detrend", mat_detrend);
    reg_mat1("sort", mat_sort);
    reg_mat1("squeeze", mat_squeeze);
    reg_mat1("cummax", mat_cummax);
    reg_mat1("cummin", mat_cummin);
    reg_mat1("unique", mat_unique);
    reg_mat1("find", mat_find);
    reg_mat1("sortrows", mat_sortrows);
    reg_mat1("rescale", mat_rescale);
    reg_mat1_err("zscore", mat_zscore);
    reg_mat1_err("pdist", mat_pdist);
    reg_mat1_err("eig", mat_eig);
    reg_mat1_err("eigenvalues", mat_eig);
    
    /* === REDUCTIONS (matrix -> scalar) === */
    reg_reduce("sum", mat_sum);
    reg_reduce("prod", mat_prod);
    reg_reduce("mean", mat_mean);
    reg_reduce("det", mat_det);
    reg_reduce("trace", mat_trace);
    reg_reduce("tr", mat_trace);
    reg_reduce("min", mat_min);
    reg_reduce("max", mat_max);
    reg_reduce("std", mat_std);
    reg_reduce("sd", mat_std);
    reg_reduce("var", mat_var);
    reg_reduce("median", mat_median);
    reg_reduce("norm", mat_norm_frobenius);
    reg_reduce("norm1", mat_norm_1);
    reg_reduce("norminf", mat_norm_inf);
    reg_reduce("range", mat_range);
    reg_reduce("rms", mat_rms);
    reg_reduce("geomean", mat_geomean);
    reg_reduce("harmmean", mat_harmmean);
    reg_reduce("sumsq", mat_sumsq);
    reg_reduce("meansq", mat_meansq);
    reg_reduce("any", mat_any);
    reg_reduce("all", mat_all);
    reg_reduce("isempty", mat_isempty);
    reg_reduce("numel", mat_numel);
    reg_reduce("length", mat_length_r);
    reg_reduce("ndims", mat_ndims_r);
    reg_reduce("isscalar", mat_isscalar);
    reg_reduce("isvector", mat_isvector);
    reg_reduce("isrow", mat_isrow);
    reg_reduce("iscolumn", mat_iscolumn);
    reg_reduce("ismatrix", mat_ismatrix);
    reg_reduce("issquare", mat_issquare);
    reg_reduce("skewness", mat_skewness);
    reg_reduce("kurtosis", mat_kurtosis);
    
    /* Matrix properties */
    reg_reduce("numel", mat_numel_r);
    reg_reduce("length", mat_length_r);
    reg_reduce("ndims", mat_ndims_r);
    reg_reduce("rows", mat_nrows_r);
    reg_reduce("cols", mat_ncols_r);
    reg_reduce("width", mat_width_r);
    reg_reduce("height", mat_height_r);
    
    /* Higher moment statistics */
    reg_reduce("mad", mat_mad);
    reg_reduce("nnz", mat_nnz);
    reg_reduce("mode", mat_mode);
    reg_reduce("iqr", mat_iqr);
    
    /* Cumulative functions */
    reg_mat1("cumsum", mat_cumsum);
    reg_mat1("cumprod", mat_cumprod);
    reg_mat1("cummax", mat_cummax);
    reg_mat1("cummin", mat_cummin);
    
    /* MATLAB compatibility functions */
    reg_mat1("rref", mat_rref);
    reg_mat1("orth", mat_orth);
    reg_mat1("poly", mat_poly);
    reg_mat1("roots", mat_roots);
    reg_mat1("expm", mat_expm);
    reg_mat1("logm", mat_logm);
    reg_mat1("fft", mat_fft);
    reg_mat1("ifft", mat_ifft);
    reg_mat1("isoutlier", mat_isoutlier);
    reg_mat1("hess", mat_hess);
    reg_mat1("balance", mat_balance);
    
    /* === SaaS METRICS === */
    reg_mat1("mrr", mat_mrr);
    reg_mat1("arpu", mat_arpu);
    reg_mat1("mrrbridge", mat_mrrbridge);
    reg_mat1("customercount", mat_customercount);
    reg_mat1("newcustomers", mat_newcustomers);
    reg_mat1("churn", mat_churn);
    reg_mat1("reactivated", mat_reactivated);
    reg_mat1("churnrate", mat_churnrate);
    reg_mat1("nrr", mat_nrr);
    reg_mat1("grr", mat_grr);
    reg_mat1("retention", mat_retention);
    reg_mat1("tenure", mat_tenure);
    reg_mat1("concentration", mat_concentration);
    reg_mat1("revchurn", mat_revchurn);
    reg_mat1("netchurn", mat_netchurn);
    reg_mat1("quickratio", mat_quickratio);
    reg_mat1("ltv", mat_ltv);
}

MatFuncEntry *matrix_func_lookup(const char *name) {
    unsigned h = mf_hash_name(name);
    MatFuncEntry *e = mf_hash[h];
    while (e) {
        if (mf_streq(e->name, name)) return e;
        e = e->next;
    }
    return NULL;
}

/* Ensure arg is a matrix (convert scalar if needed) */
static void ensure_matrix(mf_value_t *v) {
    if (v->type == MF_VAL_SCALAR) {
        apfc saved = v->v.scalar;
        mat_zero(&v->v.matrix, 1, 1);
        MAT_AT(&v->v.matrix, 0, 0) = saved;
        v->type = MF_VAL_MATRIX;
    }
}

int matrix_func_eval1(MatFuncEntry *entry, mf_value_t *result, mf_value_t *arg) {
    ensure_matrix(arg);
    
    switch (entry->type) {
        case MF_MAT_1:
            result->type = MF_VAL_MATRIX;
            if (entry->can_fail) {
                if (!entry->fn.m1_err(&result->v.matrix, &arg->v.matrix)) {
                    return 0;
                }
            } else {
                entry->fn.m1(&result->v.matrix, &arg->v.matrix);
            }
            return 1;
            
        case MF_REDUCE_1:
            result->type = MF_VAL_SCALAR;
            entry->fn.red(&result->v.scalar, &arg->v.matrix);
            return 1;
            
        default:
            return 0;
    }
}
