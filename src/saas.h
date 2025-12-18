/*
 * saas.h - SaaS metrics and time series analysis
 */
#ifndef SAAS_H
#define SAAS_H

#include "matrix.h"

/* Retime intervals */
#define RETIME_HOURLY    0
#define RETIME_DAILY     1
#define RETIME_WEEKLY    2
#define RETIME_MONTHLY   3
#define RETIME_QUARTERLY 4
#define RETIME_YEARLY    5

/* Revenue metrics */
void mat_mrr(matrix_t *r, const matrix_t *data);
void mat_arpu(matrix_t *r, const matrix_t *data);
void mat_mrrbridge(matrix_t *r, const matrix_t *data);

/* Customer metrics */
void mat_customercount(matrix_t *r, const matrix_t *data);
void mat_newcustomers(matrix_t *r, const matrix_t *data);
void mat_churn(matrix_t *r, const matrix_t *data);
void mat_reactivated(matrix_t *r, const matrix_t *data);
void mat_churnrate(matrix_t *r, const matrix_t *data);
void mat_topcustomers(matrix_t *r, const matrix_t *data, int n_top);
void mat_tenure(matrix_t *r, const matrix_t *data);

/* Retention metrics */
void mat_nrr(matrix_t *r, const matrix_t *data);
void mat_grr(matrix_t *r, const matrix_t *data);
void mat_retention(matrix_t *r, const matrix_t *data);

/* Retiming */
void mat_retime(matrix_t *r, const matrix_t *data, int interval);

#endif /* SAAS_H */
