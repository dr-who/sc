/*
 * ml.h - Machine Learning functions for scalc
 * C89 compliant, DOS compatible
 */

#ifndef ML_H
#define ML_H

#include "matrix.h"
#include "apfc.h"

/* Initialize ML subsystem */
void ml_init(void);

/* Data preprocessing */
int mat_zscore(matrix_t *r, const matrix_t *X);

/* Training functions - return model ID through first parameter */
int ml_fit_knn(int *model_id, const matrix_t *X, const matrix_t *Y, int k);
int ml_fit_svm(int *model_id, const matrix_t *X, const matrix_t *Y, int standardize);
int ml_fit_tree(int *model_id, const matrix_t *X, const matrix_t *Y);
int ml_fit_nb(int *model_id, const matrix_t *X, const matrix_t *Y);

/* Prediction */
int ml_predict(matrix_t *Ypred, int model_id, const matrix_t *X);

/* Cross-validation */
int ml_cvpartition(matrix_t *train_idx, matrix_t *test_idx, 
                   int n, int holdout_pct);

/* Evaluation */
int ml_confusionmat(matrix_t *C, const matrix_t *Ytrue, const matrix_t *Ypred);
void ml_accuracy(apfc *acc, const matrix_t *Ytrue, const matrix_t *Ypred);

/* Model info */
const char *ml_model_type_name(int model_id);
int ml_model_n_classes(int model_id);

#endif /* ML_H */
