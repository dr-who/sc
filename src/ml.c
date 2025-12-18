/*
 * ml.c - Machine Learning functions for scalc
 * C89 compliant, DOS compatible (16-bit clean)
 * 
 * Implements:
 * - zscore: Standardize data
 * - fitcknn: k-Nearest Neighbors classifier
 * - fitcsvm: Support Vector Machine (linear)
 * - fitctree: Decision Tree classifier
 * - fitcnb: Naive Bayes classifier (Gaussian)
 * - predict: Predict using trained model
 * - cvpartition: Cross-validation partitioning
 * - confusionmat: Confusion matrix
 * - accuracy: Classification accuracy
 */

#include <stdio.h>
#include <string.h>
#include "matrix.h"
#include "apf.h"
#include "apfx.h"

/* Model types */
#define ML_MODEL_NONE   0
#define ML_MODEL_KNN    1
#define ML_MODEL_SVM    2
#define ML_MODEL_TREE   3
#define ML_MODEL_NB     4

/* Limits - keep small for DOS */
#define MAX_MODELS 4
#define MAX_CLASSES 10
#define MAX_ML_FEATURES 20

/* Model storage - uses matrices from arena */
typedef struct {
    int type;
    int n_samples;
    int n_features;
    int n_classes;
    int k;  /* for KNN */
    
    /* Training data stored as matrices (uses arena) */
    matrix_t X_train;
    matrix_t Y_train;
    
    /* SVM weights and bias */
    matrix_t weights;
    matrix_t bias;
    
    /* Decision tree (simple stump) */
    int tree_feature;
    apf tree_threshold;
    int tree_left_class;
    int tree_right_class;
    
    /* Naive Bayes - store as matrices */
    matrix_t nb_means;   /* n_classes x n_features */
    matrix_t nb_vars;    /* n_classes x n_features */
    matrix_t nb_priors;  /* n_classes x 1 */
    
} ml_model_t;

static ml_model_t models[MAX_MODELS];
static int next_model_id = 0;

/* Initialize ML subsystem */
void ml_init(void)
{
    int i;
    for (i = 0; i < MAX_MODELS; i++) {
        models[i].type = ML_MODEL_NONE;
    }
    next_model_id = 0;
}

/* zscore: Standardize matrix columns to mean=0, std=1 */
int mat_zscore(matrix_t *r, const matrix_t *X)
{
    int i, j;
    int m, n;
    apf mean, std, sum, sumsq, diff, var, count;
    
    m = X->rows;
    n = X->cols;
    
    mat_zero(r, m, n);
    if (!r->data) return 0;
    
    apf_from_int(&count, m);
    
    for (j = 0; j < n; j++) {
        /* Compute mean */
        apf_zero(&sum);
        for (i = 0; i < m; i++) {
            apf_add(&sum, &sum, &MAT_AT(X, i, j).re);
        }
        apf_div(&mean, &sum, &count);
        
        /* Compute std */
        apf_zero(&sumsq);
        for (i = 0; i < m; i++) {
            apf_sub(&diff, &MAT_AT(X, i, j).re, &mean);
            apf_mul(&var, &diff, &diff);
            apf_add(&sumsq, &sumsq, &var);
        }
        apf_div(&var, &sumsq, &count);
        apf_sqrt(&std, &var);
        
        /* Standardize */
        for (i = 0; i < m; i++) {
            apf_sub(&diff, &MAT_AT(X, i, j).re, &mean);
            if (!apf_is_zero(&std)) {
                apf_div(&MAT_AT(r, i, j).re, &diff, &std);
            } else {
                apf_zero(&MAT_AT(r, i, j).re);
            }
            apf_zero(&MAT_AT(r, i, j).im);
        }
    }
    
    return 1;
}

/* Helper: Euclidean distance between two rows */
static void row_distance(apf *dist, const matrix_t *X1, int r1, 
                         const matrix_t *X2, int r2, int n_features)
{
    int j;
    apf sum, diff, sq;
    
    apf_zero(&sum);
    for (j = 0; j < n_features; j++) {
        apf_sub(&diff, &MAT_AT(X1, r1, j).re, &MAT_AT(X2, r2, j).re);
        apf_mul(&sq, &diff, &diff);
        apf_add(&sum, &sum, &sq);
    }
    apf_sqrt(dist, &sum);
}

/* fitcknn: Train k-Nearest Neighbors classifier */
int ml_fit_knn(int *model_id, const matrix_t *X, const matrix_t *Y, int k)
{
    int id, i;
    ml_model_t *m;
    int max_class;
    long label;
    
    if (next_model_id >= MAX_MODELS) {
        printf("Error: Maximum number of models reached\n");
        return 0;
    }
    
    id = next_model_id++;
    m = &models[id];
    
    m->type = ML_MODEL_KNN;
    m->n_samples = X->rows;
    m->n_features = X->cols;
    m->k = (k > 0) ? k : 5;
    
    /* Copy training data */
    mat_copy_persist(&m->X_train, X);
    mat_copy_persist(&m->Y_train, Y);
    
    /* Find number of classes */
    max_class = 0;
    for (i = 0; i < Y->rows; i++) {
        label = apf_to_long(&MAT_AT(Y, i, 0).re);
        if (label > max_class) max_class = (int)label;
    }
    m->n_classes = max_class;
    
    *model_id = id;
    printf("Trained k-NN (k=%d, %d samples, %d classes)\n", 
           m->k, m->n_samples, m->n_classes);
    return 1;
}

/* Predict using KNN - uses static arrays to avoid stack issues */
static int predict_knn(matrix_t *Ypred, const ml_model_t *m, const matrix_t *X)
{
    /* Static to avoid stack overflow on DOS */
    static apf distances[150];  /* Max 150 training samples */
    static int indices[150];
    int votes[MAX_CLASSES];
    int i, j, c, ki;
    int n_test;
    int n_train;
    apf tmp_d;
    int tmp_i;
    int min_idx;
    int max_votes, pred_class;
    long label;
    
    n_test = X->rows;
    n_train = m->n_samples;
    if (n_train > 150) n_train = 150;
    
    mat_zero(Ypred, n_test, 1);
    if (!Ypred->data) return 0;
    
    for (i = 0; i < n_test; i++) {
        /* Compute distances to all training points */
        for (j = 0; j < n_train; j++) {
            row_distance(&distances[j], X, i, &m->X_train, j, m->n_features);
            indices[j] = j;
        }
        
        /* Partial selection sort for k nearest */
        for (ki = 0; ki < m->k && ki < n_train; ki++) {
            min_idx = ki;
            for (j = ki + 1; j < n_train; j++) {
                if (apf_lt(&distances[j], &distances[min_idx])) {
                    min_idx = j;
                }
            }
            if (min_idx != ki) {
                tmp_d = distances[ki];
                tmp_i = indices[ki];
                distances[ki] = distances[min_idx];
                indices[ki] = indices[min_idx];
                distances[min_idx] = tmp_d;
                indices[min_idx] = tmp_i;
            }
        }
        
        /* Vote among k nearest */
        for (c = 0; c < MAX_CLASSES; c++) votes[c] = 0;
        for (ki = 0; ki < m->k && ki < n_train; ki++) {
            label = apf_to_long(&MAT_AT(&m->Y_train, indices[ki], 0).re);
            if (label > 0 && label <= MAX_CLASSES) {
                votes[label - 1]++;
            }
        }
        
        /* Find majority class */
        max_votes = 0;
        pred_class = 1;
        for (c = 0; c < m->n_classes; c++) {
            if (votes[c] > max_votes) {
                max_votes = votes[c];
                pred_class = c + 1;
            }
        }
        apf_from_int(&MAT_AT(Ypred, i, 0).re, pred_class);
        apf_zero(&MAT_AT(Ypred, i, 0).im);
    }
    
    return 1;
}

/* fitcsvm: Train linear SVM (one-vs-all) using simple SGD */
int ml_fit_svm(int *model_id, const matrix_t *X, const matrix_t *Y, int standardize)
{
    int id, i, j, c, iter;
    ml_model_t *m;
    int n, d;
    int max_class;
    long label;
    apf lr, margin, one;
    int yi;
    apf score, yi_apf, prod, update;
    apf hundred;
    
    (void)standardize;
    
    if (next_model_id >= MAX_MODELS) {
        printf("Error: Maximum number of models reached\n");
        return 0;
    }
    
    id = next_model_id++;
    m = &models[id];
    
    n = X->rows;
    d = X->cols;
    
    m->type = ML_MODEL_SVM;
    m->n_samples = n;
    m->n_features = d;
    
    /* Find number of classes */
    max_class = 0;
    for (i = 0; i < Y->rows; i++) {
        label = apf_to_long(&MAT_AT(Y, i, 0).re);
        if (label > max_class) max_class = (int)label;
    }
    m->n_classes = max_class;
    
    /* Initialize weights and bias */
    mat_zero_persist(&m->weights, max_class, d);
    mat_zero_persist(&m->bias, max_class, 1);
    
    /* lr = 0.01 = 1/100 */
    apf_from_int(&one, 1);
    apf_from_int(&hundred, 100);
    apf_div(&lr, &one, &hundred);
    
    /* Train one-vs-all classifiers */
    for (c = 0; c < max_class; c++) {
        /* SGD for linear SVM - 50 iterations */
        for (iter = 0; iter < 50; iter++) {
            for (i = 0; i < n; i++) {
                label = apf_to_long(&MAT_AT(Y, i, 0).re);
                yi = (label == c + 1) ? 1 : -1;
                
                /* Compute score = w.x + b */
                apf_copy(&score, &MAT_AT(&m->bias, c, 0).re);
                for (j = 0; j < d; j++) {
                    apf_mul(&prod, &MAT_AT(&m->weights, c, j).re, 
                            &MAT_AT(X, i, j).re);
                    apf_add(&score, &score, &prod);
                }
                
                /* Compute yi * score */
                apf_from_int(&yi_apf, yi);
                apf_mul(&margin, &yi_apf, &score);
                
                /* Hinge loss update: if yi * score < 1 */
                if (apf_lt(&margin, &one)) {
                    for (j = 0; j < d; j++) {
                        apf_mul(&prod, &yi_apf, &MAT_AT(X, i, j).re);
                        apf_mul(&update, &lr, &prod);
                        apf_add(&MAT_AT(&m->weights, c, j).re, 
                                &MAT_AT(&m->weights, c, j).re, &update);
                    }
                    apf_mul(&update, &lr, &yi_apf);
                    apf_add(&MAT_AT(&m->bias, c, 0).re, 
                            &MAT_AT(&m->bias, c, 0).re, &update);
                }
            }
        }
    }
    
    /* Store training data for reference */
    mat_copy_persist(&m->X_train, X);
    mat_copy_persist(&m->Y_train, Y);
    
    *model_id = id;
    printf("Trained linear SVM (%d samples, %d classes)\n", n, max_class);
    return 1;
}

/* Predict using SVM */
static int predict_svm(matrix_t *Ypred, const ml_model_t *m, const matrix_t *X)
{
    int i, j, c;
    int n_test;
    apf max_score, score, prod;
    int pred_class;
    
    n_test = X->rows;
    
    mat_zero(Ypred, n_test, 1);
    if (!Ypred->data) return 0;
    
    for (i = 0; i < n_test; i++) {
        apf_from_int(&max_score, -10000);
        pred_class = 1;
        
        for (c = 0; c < m->n_classes; c++) {
            apf_copy(&score, &MAT_AT(&m->bias, c, 0).re);
            
            for (j = 0; j < m->n_features; j++) {
                apf_mul(&prod, &MAT_AT(&m->weights, c, j).re, 
                        &MAT_AT(X, i, j).re);
                apf_add(&score, &score, &prod);
            }
            
            if (apf_lt(&max_score, &score)) {
                apf_copy(&max_score, &score);
                pred_class = c + 1;
            }
        }
        
        apf_from_int(&MAT_AT(Ypred, i, 0).re, pred_class);
        apf_zero(&MAT_AT(Ypred, i, 0).im);
    }
    
    return 1;
}

/* fitctree: Train decision tree (simple single split) */
int ml_fit_tree(int *model_id, const matrix_t *X, const matrix_t *Y)
{
    int id, i, j, c;
    ml_model_t *m;
    int n, d;
    int max_class;
    long label;
    apf best_gini, gini, threshold;
    int best_feature;
    apf sum, mean, count_apf;
    int left_counts[MAX_CLASSES];
    int right_counts[MAX_CLASSES];
    int left_total, right_total;
    apf gini_left, gini_right, p, p2, one_apf;
    apf left_n, right_n, total_n;
    int max_left, max_right;
    
    if (next_model_id >= MAX_MODELS) {
        printf("Error: Maximum number of models reached\n");
        return 0;
    }
    
    id = next_model_id++;
    m = &models[id];
    
    n = X->rows;
    d = X->cols;
    
    m->type = ML_MODEL_TREE;
    m->n_samples = n;
    m->n_features = d;
    
    /* Find number of classes */
    max_class = 0;
    for (i = 0; i < n; i++) {
        label = apf_to_long(&MAT_AT(Y, i, 0).re);
        if (label > max_class) max_class = (int)label;
    }
    m->n_classes = max_class;
    
    /* Find best single split */
    apf_from_int(&best_gini, 2);
    best_feature = 0;
    apf_zero(&threshold);
    
    for (j = 0; j < d; j++) {
        /* Try splitting at mean of feature j */
        apf_zero(&sum);
        for (i = 0; i < n; i++) {
            apf_add(&sum, &sum, &MAT_AT(X, i, j).re);
        }
        apf_from_int(&count_apf, n);
        apf_div(&mean, &sum, &count_apf);
        
        /* Count samples on each side */
        for (c = 0; c < MAX_CLASSES; c++) {
            left_counts[c] = 0;
            right_counts[c] = 0;
        }
        left_total = 0;
        right_total = 0;
        
        for (i = 0; i < n; i++) {
            label = apf_to_long(&MAT_AT(Y, i, 0).re) - 1;
            if (label >= 0 && label < MAX_CLASSES) {
                if (apf_lt(&MAT_AT(X, i, j).re, &mean)) {
                    left_counts[label]++;
                    left_total++;
                } else {
                    right_counts[label]++;
                    right_total++;
                }
            }
        }
        
        /* Compute Gini impurity */
        apf_from_int(&one_apf, 1);
        apf_from_int(&gini_left, 1);
        apf_from_int(&gini_right, 1);
        
        if (left_total > 0) {
            apf_from_int(&left_n, left_total);
            for (c = 0; c < max_class; c++) {
                apf_from_int(&p, left_counts[c]);
                apf_div(&p, &p, &left_n);
                apf_mul(&p2, &p, &p);
                apf_sub(&gini_left, &gini_left, &p2);
            }
        } else {
            apf_zero(&gini_left);
        }
        
        if (right_total > 0) {
            apf_from_int(&right_n, right_total);
            for (c = 0; c < max_class; c++) {
                apf_from_int(&p, right_counts[c]);
                apf_div(&p, &p, &right_n);
                apf_mul(&p2, &p, &p);
                apf_sub(&gini_right, &gini_right, &p2);
            }
        } else {
            apf_zero(&gini_right);
        }
        
        /* Weighted Gini */
        apf_from_int(&total_n, n);
        apf_from_int(&left_n, left_total);
        apf_from_int(&right_n, right_total);
        
        apf_mul(&gini_left, &gini_left, &left_n);
        apf_mul(&gini_right, &gini_right, &right_n);
        apf_add(&gini, &gini_left, &gini_right);
        apf_div(&gini, &gini, &total_n);
        
        if (apf_lt(&gini, &best_gini)) {
            apf_copy(&best_gini, &gini);
            best_feature = j;
            apf_copy(&threshold, &mean);
        }
    }
    
    /* Store tree */
    m->tree_feature = best_feature;
    apf_copy(&m->tree_threshold, &threshold);
    
    /* Determine majority class for each side */
    for (c = 0; c < MAX_CLASSES; c++) {
        left_counts[c] = 0;
        right_counts[c] = 0;
    }
    
    for (i = 0; i < n; i++) {
        label = apf_to_long(&MAT_AT(Y, i, 0).re) - 1;
        if (label >= 0 && label < MAX_CLASSES) {
            if (apf_lt(&MAT_AT(X, i, best_feature).re, &threshold)) {
                left_counts[label]++;
            } else {
                right_counts[label]++;
            }
        }
    }
    
    max_left = 0;
    max_right = 0;
    m->tree_left_class = 1;
    m->tree_right_class = 1;
    for (c = 0; c < max_class; c++) {
        if (left_counts[c] > max_left) {
            max_left = left_counts[c];
            m->tree_left_class = c + 1;
        }
        if (right_counts[c] > max_right) {
            max_right = right_counts[c];
            m->tree_right_class = c + 1;
        }
    }
    
    mat_copy_persist(&m->X_train, X);
    mat_copy_persist(&m->Y_train, Y);
    
    *model_id = id;
    {
        /* Print without using native double */
        char thresh_str[64];
        apf_to_str(thresh_str, sizeof(thresh_str), &threshold, 3);
        printf("Trained Decision Tree (split: feature %d < %s)\n", 
               best_feature + 1, thresh_str);
    }
    return 1;
}

/* Predict using decision tree */
static int predict_tree(matrix_t *Ypred, const ml_model_t *m, const matrix_t *X)
{
    int i;
    int n_test;
    int pred;
    
    n_test = X->rows;
    
    mat_zero(Ypred, n_test, 1);
    if (!Ypred->data) return 0;
    
    for (i = 0; i < n_test; i++) {
        if (apf_lt(&MAT_AT(X, i, m->tree_feature).re, &m->tree_threshold)) {
            pred = m->tree_left_class;
        } else {
            pred = m->tree_right_class;
        }
        apf_from_int(&MAT_AT(Ypred, i, 0).re, pred);
        apf_zero(&MAT_AT(Ypred, i, 0).im);
    }
    
    return 1;
}

/* fitcnb: Train Gaussian Naive Bayes */
int ml_fit_nb(int *model_id, const matrix_t *X, const matrix_t *Y)
{
    int id, i, j, c;
    ml_model_t *m;
    int n, d;
    int max_class;
    long label;
    int class_counts[MAX_CLASSES];
    apf count, total, eps;
    apf diff, sq;
    
    if (next_model_id >= MAX_MODELS) {
        printf("Error: Maximum number of models reached\n");
        return 0;
    }
    
    id = next_model_id++;
    m = &models[id];
    
    n = X->rows;
    d = X->cols;
    if (d > MAX_ML_FEATURES) d = MAX_ML_FEATURES;
    
    m->type = ML_MODEL_NB;
    m->n_samples = n;
    m->n_features = d;
    
    /* Find number of classes and count samples per class */
    max_class = 0;
    for (c = 0; c < MAX_CLASSES; c++) class_counts[c] = 0;
    
    for (i = 0; i < n; i++) {
        label = apf_to_long(&MAT_AT(Y, i, 0).re);
        if (label > max_class) max_class = (int)label;
        if (label > 0 && label <= MAX_CLASSES) {
            class_counts[label - 1]++;
        }
    }
    m->n_classes = max_class;
    
    /* Allocate NB parameter matrices */
    mat_zero_persist(&m->nb_means, max_class, d);
    mat_zero_persist(&m->nb_vars, max_class, d);
    mat_zero_persist(&m->nb_priors, max_class, 1);
    
    /* Compute means per class */
    for (i = 0; i < n; i++) {
        label = apf_to_long(&MAT_AT(Y, i, 0).re) - 1;
        if (label >= 0 && label < max_class) {
            for (j = 0; j < d; j++) {
                apf_add(&MAT_AT(&m->nb_means, (int)label, j).re, 
                        &MAT_AT(&m->nb_means, (int)label, j).re, 
                        &MAT_AT(X, i, j).re);
            }
        }
    }
    
    for (c = 0; c < max_class; c++) {
        if (class_counts[c] > 0) {
            apf_from_int(&count, class_counts[c]);
            for (j = 0; j < d; j++) {
                apf_div(&MAT_AT(&m->nb_means, c, j).re, 
                        &MAT_AT(&m->nb_means, c, j).re, &count);
            }
        }
    }
    
    /* Compute variances per class */
    for (i = 0; i < n; i++) {
        label = apf_to_long(&MAT_AT(Y, i, 0).re) - 1;
        if (label >= 0 && label < max_class) {
            for (j = 0; j < d; j++) {
                apf_sub(&diff, &MAT_AT(X, i, j).re, 
                        &MAT_AT(&m->nb_means, (int)label, j).re);
                apf_mul(&sq, &diff, &diff);
                apf_add(&MAT_AT(&m->nb_vars, (int)label, j).re, 
                        &MAT_AT(&m->nb_vars, (int)label, j).re, &sq);
            }
        }
    }
    
    /* eps = 1e-6 = 1/1000000 */
    {
        apf million;
        apf_from_int(&eps, 1);
        apf_from_int(&million, 1000000);
        apf_div(&eps, &eps, &million);
    }
    for (c = 0; c < max_class; c++) {
        if (class_counts[c] > 0) {
            apf_from_int(&count, class_counts[c]);
            for (j = 0; j < d; j++) {
                apf_div(&MAT_AT(&m->nb_vars, c, j).re, 
                        &MAT_AT(&m->nb_vars, c, j).re, &count);
                /* Add epsilon to avoid division by zero */
                apf_add(&MAT_AT(&m->nb_vars, c, j).re, 
                        &MAT_AT(&m->nb_vars, c, j).re, &eps);
            }
        }
    }
    
    /* Compute priors */
    apf_from_int(&total, n);
    for (c = 0; c < max_class; c++) {
        apf_from_int(&MAT_AT(&m->nb_priors, c, 0).re, class_counts[c]);
        apf_div(&MAT_AT(&m->nb_priors, c, 0).re, 
                &MAT_AT(&m->nb_priors, c, 0).re, &total);
    }
    
    mat_copy_persist(&m->X_train, X);
    mat_copy_persist(&m->Y_train, Y);
    
    *model_id = id;
    printf("Trained Naive Bayes (%d samples, %d classes)\n", n, max_class);
    return 1;
}

/* Predict using Naive Bayes */
static int predict_nb(matrix_t *Ypred, const ml_model_t *m, const matrix_t *X)
{
    int i, j, c;
    int n_test;
    apf max_log_prob, log_prob, log_prior;
    int pred_class;
    apf diff, sq, exp_term, log_term, pi2, half, var;
    apf two;
    
    n_test = X->rows;
    
    mat_zero(Ypred, n_test, 1);
    if (!Ypred->data) return 0;
    
    /* pi2 = 2*pi */
    apfx_pi(&pi2);
    apf_from_int(&two, 2);
    apf_mul(&pi2, &pi2, &two);
    
    /* half = -1/2 */
    apf_from_int(&half, -1);
    apf_div(&half, &half, &two);
    
    for (i = 0; i < n_test; i++) {
        apf_from_int(&max_log_prob, -10000);
        pred_class = 1;
        
        for (c = 0; c < m->n_classes; c++) {
            /* Start with log prior */
            if (!apf_is_zero(&MAT_AT(&m->nb_priors, c, 0).re)) {
                apfx_log(&log_prior, &MAT_AT(&m->nb_priors, c, 0).re);
            } else {
                apf_from_int(&log_prior, -100);
            }
            apf_copy(&log_prob, &log_prior);
            
            /* Add log likelihoods (Gaussian) */
            for (j = 0; j < m->n_features; j++) {
                apf_copy(&var, &MAT_AT(&m->nb_vars, c, j).re);
                
                apf_sub(&diff, &MAT_AT(X, i, j).re, 
                        &MAT_AT(&m->nb_means, c, j).re);
                apf_mul(&sq, &diff, &diff);
                
                /* -0.5 * (x-mu)^2 / var */
                if (!apf_is_zero(&var)) {
                    apf_div(&exp_term, &sq, &var);
                    apf_mul(&exp_term, &exp_term, &half);
                } else {
                    apf_zero(&exp_term);
                }
                
                /* -0.5 * log(2*pi*var) */
                apf_mul(&log_term, &pi2, &var);
                if (!apf_is_zero(&log_term)) {
                    apfx_log(&log_term, &log_term);
                    apf_mul(&log_term, &log_term, &half);
                } else {
                    apf_zero(&log_term);
                }
                
                apf_add(&log_prob, &log_prob, &exp_term);
                apf_add(&log_prob, &log_prob, &log_term);
            }
            
            if (apf_lt(&max_log_prob, &log_prob)) {
                apf_copy(&max_log_prob, &log_prob);
                pred_class = c + 1;
            }
        }
        
        apf_from_int(&MAT_AT(Ypred, i, 0).re, pred_class);
        apf_zero(&MAT_AT(Ypred, i, 0).im);
    }
    
    return 1;
}

/* predict: Unified prediction function */
int ml_predict(matrix_t *Ypred, int model_id, const matrix_t *X)
{
    ml_model_t *m;
    
    if (model_id < 0 || model_id >= MAX_MODELS) {
        printf("Error: Invalid model ID\n");
        return 0;
    }
    
    m = &models[model_id];
    
    switch (m->type) {
        case ML_MODEL_KNN:
            return predict_knn(Ypred, m, X);
        case ML_MODEL_SVM:
            return predict_svm(Ypred, m, X);
        case ML_MODEL_TREE:
            return predict_tree(Ypred, m, X);
        case ML_MODEL_NB:
            return predict_nb(Ypred, m, X);
        default:
            printf("Error: Unknown model type\n");
            return 0;
    }
}

/* cvpartition: Create train/test split indices 
 * holdout_pct is percentage (0-100) for test set */
int ml_cvpartition(matrix_t *train_idx, matrix_t *test_idx, 
                   int n, int holdout_pct)
{
    int i;
    int n_test, n_train;
    int test_count, train_count;
    
    /* n_test = n * holdout_pct / 100 */
    n_test = (n * holdout_pct) / 100;
    n_train = n - n_test;
    
    mat_zero(train_idx, n_train, 1);
    mat_zero(test_idx, n_test, 1);
    
    if (!train_idx->data || !test_idx->data) return 0;
    
    /* Simple split: last n_test go to test */
    train_count = 0;
    test_count = 0;
    for (i = 0; i < n; i++) {
        if (i < n_train) {
            apf_from_int(&MAT_AT(train_idx, train_count, 0).re, i + 1);
            apf_zero(&MAT_AT(train_idx, train_count, 0).im);
            train_count++;
        } else {
            apf_from_int(&MAT_AT(test_idx, test_count, 0).re, i + 1);
            apf_zero(&MAT_AT(test_idx, test_count, 0).im);
            test_count++;
        }
    }
    
    return 1;
}

/* confusionmat: Create confusion matrix */
int ml_confusionmat(matrix_t *C, const matrix_t *Ytrue, const matrix_t *Ypred)
{
    int i, n;
    int max_class;
    long c1, c2;
    int true_c, pred_c;
    apf one;
    
    n = Ytrue->rows;
    
    /* Find max class */
    max_class = 0;
    for (i = 0; i < n; i++) {
        c1 = apf_to_long(&MAT_AT(Ytrue, i, 0).re);
        c2 = apf_to_long(&MAT_AT(Ypred, i, 0).re);
        if (c1 > max_class) max_class = (int)c1;
        if (c2 > max_class) max_class = (int)c2;
    }
    
    mat_zero(C, max_class, max_class);
    if (!C->data) return 0;
    
    apf_from_int(&one, 1);
    for (i = 0; i < n; i++) {
        true_c = (int)apf_to_long(&MAT_AT(Ytrue, i, 0).re) - 1;
        pred_c = (int)apf_to_long(&MAT_AT(Ypred, i, 0).re) - 1;
        if (true_c >= 0 && true_c < max_class && 
            pred_c >= 0 && pred_c < max_class) {
            apf_add(&MAT_AT(C, true_c, pred_c).re, 
                    &MAT_AT(C, true_c, pred_c).re, &one);
        }
    }
    
    return 1;
}

/* accuracy: Compute classification accuracy */
void ml_accuracy(apfc *acc, const matrix_t *Ytrue, const matrix_t *Ypred)
{
    int i, n, correct;
    long t, p;
    apf total;
    
    n = Ytrue->rows;
    correct = 0;
    
    for (i = 0; i < n; i++) {
        t = apf_to_long(&MAT_AT(Ytrue, i, 0).re);
        p = apf_to_long(&MAT_AT(Ypred, i, 0).re);
        if (t == p) correct++;
    }
    
    apf_from_int(&acc->re, correct);
    apf_from_int(&total, n);
    apf_div(&acc->re, &acc->re, &total);
    apf_zero(&acc->im);
}

/* Get model info */
const char *ml_model_type_name(int model_id)
{
    if (model_id < 0 || model_id >= MAX_MODELS) return "invalid";
    
    switch (models[model_id].type) {
        case ML_MODEL_KNN: return "knn";
        case ML_MODEL_SVM: return "svm";
        case ML_MODEL_TREE: return "tree";
        case ML_MODEL_NB: return "nb";
        default: return "none";
    }
}

int ml_model_n_classes(int model_id)
{
    if (model_id < 0 || model_id >= MAX_MODELS) return 0;
    return models[model_id].n_classes;
}
