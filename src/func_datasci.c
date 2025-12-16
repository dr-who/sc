/*
 * func_datasci.c - Data Science / Machine Learning functions
 * C89 compliant
 */
#include "func_registry.h"
#include "matrix.h"
#include <stdio.h>

/**
 * @func pca
 * @category Data Science
 * @syntax [coeff,score,latent] = pca(X)
 * @desc Principal Component Analysis. X is m-by-n (m observations, n variables).
 *       Returns coeff (principal components), score (projected data), 
 *       latent (explained variance).
 * @example X = [1,2;3,4;5,6;7,8] -> 1 2 / 3 4 / 5 6 / 7 8
 * @example [C,S,L] = pca(X) ->
 * @example S -> -4.24 0 / -1.41 0 / 1.41 0 / 4.24 0
 */
int func_pca(matrix_t *coeff, matrix_t *score, matrix_t *latent, const matrix_t *X)
{
    return mat_pca(coeff, score, latent, X);
}

/**
 * @func pcareduce
 * @category Data Science
 * @syntax pcareduce(X, k)
 * @desc Reduce data to k principal components for visualization.
 * @example X = [1,2,3;4,5,6;7,8,9;10,11,12] -> 1 2 3 / 4 5 6 / 7 8 9 / 10 11 12
 * @example pcareduce(X, 2) -> 2D projection
 */
int func_pcareduce(matrix_t *result, const matrix_t *X, int k)
{
    return mat_pca_reduce(result, X, k);
}

/**
 * @func kmeans
 * @category Data Science
 * @syntax [idx,C] = kmeans(X, k)
 * @desc K-means clustering. Partitions data into k clusters.
 *       Returns idx (cluster labels 0 to k-1) and C (centroids).
 * @example X = [1,1;1,2;2,1;8,8;8,9;9,8] -> 1 1 / 1 2 / 2 1 / 8 8 / 8 9 / 9 8
 * @example [idx,C] = kmeans(X, 2) ->
 * @example idx -> 0 / 0 / 0 / 1 / 1 / 1
 * @example C -> 1.33 1.33 / 8.33 8.33
 */
int func_kmeans(matrix_t *idx, matrix_t *centroids, const matrix_t *X, int k)
{
    return mat_kmeans(idx, centroids, X, k);
}

/**
 * @func silhouette
 * @category Data Science
 * @syntax silhouette(X, idx)
 * @desc Silhouette score for clustering quality. Range -1 to 1, higher is better.
 *       Measures how well-separated clusters are.
 * @example X = [1,1;1,2;8,8;8,9] -> 1 1 / 1 2 / 8 8 / 8 9
 * @example idx = [0;0;1;1] -> 0 / 0 / 1 / 1
 * @example silhouette(X, idx) -> 0.88
 */
int func_silhouette(apfc *result, const matrix_t *X, const matrix_t *idx)
{
    mat_silhouette(result, X, idx);
    return 1;
}

/**
 * @func pdist
 * @category Data Science
 * @syntax pdist(X)
 * @desc Pairwise Euclidean distances between rows.
 *       Returns condensed distance vector.
 * @example pdist([0,0;1,0;0,1]) -> 1 / 1 / 1.41
 */
int func_pdist(matrix_t *result, const matrix_t *X)
{
    return mat_pdist(result, X);
}

/**
 * @func zscore
 * @category Data Science
 * @syntax zscore(X)
 * @desc Standardize columns to mean=0, std=1.
 * @example zscore([1,2;3,4;5,6]) -> -1.22 -1.22 / 0 0 / 1.22 1.22
 */
int func_zscore(matrix_t *result, const matrix_t *X)
{
    return mat_zscore(result, X);
}

/**
 * @func normalize
 * @category Data Science
 * @syntax normalize(X)
 * @desc Normalize rows to unit length (L2 norm = 1).
 * @example normalize([3,4;1,0]) -> 0.6 0.8 / 1 0
 */
int func_normalize(matrix_t *result, const matrix_t *X)
{
    return mat_normalize(result, X);
}

/**
 * @func corrcoef
 * @category Data Science
 * @syntax corrcoef(X)
 * @desc Correlation coefficient matrix.
 * @example corrcoef([1,2,3;4,5,6]') -> correlation matrix
 */
int func_corrcoef(matrix_t *result, const matrix_t *X)
{
    return mat_corrcoef(result, X);
}

/**
 * @func cov
 * @category Data Science
 * @syntax cov(X)
 * @desc Covariance matrix.
 * @example cov([1,2;3,4;5,6]) -> 4 4 / 4 4
 */
int func_cov(matrix_t *result, const matrix_t *X)
{
    return mat_cov(result, X);
}

/*
 * Function registry for this category
 */
static const func_entry_t datasci_functions[] = {
    {"pca", "Data Science", "[coeff,score,latent] = pca(X)",
     "Principal Component Analysis",
     {"X = [1,2;3,4;5,6;7,8]", "[C,S,L] = pca(X)", "S", "L"},
     FUNC_MULTI_RET, (void*)func_pca},
    
    {"pcareduce", "Data Science", "pcareduce(X, k)",
     "Reduce data to k principal components",
     {"X = [1,2,3;4,5,6;7,8,9;10,11,12]", "pcareduce(X, 2)", NULL, NULL},
     FUNC_MATRIX2, (void*)func_pcareduce},
    
    {"kmeans", "Data Science", "[idx,C] = kmeans(X, k)",
     "K-means clustering",
     {"X = [1,1;1,2;2,1;8,8;8,9;9,8]", "[idx,C] = kmeans(X, 2)", "idx", "C"},
     FUNC_MULTI_RET, (void*)func_kmeans},
    
    {"silhouette", "Data Science", "silhouette(X, idx)",
     "Silhouette score for clustering quality (-1 to 1)",
     {"X = [1,1;1,2;8,8;8,9]; idx = [0;0;1;1]", "silhouette(X, idx)", NULL, NULL},
     FUNC_MATRIX2, (void*)func_silhouette},
    
    {"pdist", "Data Science", "pdist(X)",
     "Pairwise distances between rows",
     {"pdist([0,0;1,0;0,1])", NULL, NULL, NULL},
     FUNC_MATRIX, (void*)func_pdist},
    
    {"zscore", "Data Science", "zscore(X)",
     "Standardize columns to mean=0, std=1",
     {"zscore([1,2;3,4;5,6])", NULL, NULL, NULL},
     FUNC_MATRIX, (void*)func_zscore},
    
    {"normalize", "Data Science", "normalize(X)",
     "Normalize rows to unit length",
     {"normalize([3,4;1,0])", NULL, NULL, NULL},
     FUNC_MATRIX, (void*)func_normalize},
    
    {"corrcoef", "Data Science", "corrcoef(X)",
     "Correlation coefficient matrix",
     {"corrcoef([1,2,3;4,5,6]')", NULL, NULL, NULL},
     FUNC_MATRIX, (void*)func_corrcoef},
    
    {"cov", "Data Science", "cov(X)",
     "Covariance matrix",
     {"cov([1,2;3,4;5,6])", NULL, NULL, NULL},
     FUNC_MATRIX, (void*)func_cov},
    
    {NULL, NULL, NULL, NULL, {NULL, NULL, NULL, NULL}, 0, NULL}
};

/* Category export */
const func_category_t func_datasci_category = {
    "Data Science",
    datasci_functions,
    sizeof(datasci_functions) / sizeof(datasci_functions[0]) - 1
};
