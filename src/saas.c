/*
 * saas.c - SaaS metrics and time series analysis functions
 * 
 * Comprehensive metrics for subscription business analysis:
 * - Revenue: MRR, ARR, ARPU, MRR Bridge components
 * - Customers: Active, New, Churned, Reactivated
 * - Retention: NRR, GRR, Cohort tables
 * - Forecasting: Per-customer and aggregate
 * 
 * C89 compliant for DOS compatibility.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "matrix.h"
#include "apf.h"
#include "apfc.h"
#include "apf_native.h"

/* External datetime functions */
extern void unix_to_datetime(long ts, int *year, int *month, int *day, 
                             int *hour, int *min, int *sec);
extern long datetime_to_unix(int year, int month, int day, 
                             int hour, int min, int sec);

/* Safety limits */
#define MAX_MONTHS 10000
#define MAX_CUSTOMERS 1000000

/* ============================================================
 * HELPER FUNCTIONS
 * ============================================================ */

static int ts_to_month(long ts)
{
    int year, month, day, hour, min, sec;
    unix_to_datetime(ts, &year, &month, &day, &hour, &min, &sec);
    return year * 12 + (month - 1);
}

static long month_to_ts(int month_idx)
{
    int year = month_idx / 12;
    int month = (month_idx % 12) + 1;
    return datetime_to_unix(year, month, 1, 0, 0, 0);
}

/* Find data ranges - returns 0 on failure */
static int get_ranges(const matrix_t *data, int *min_m, int *max_m, 
                      int *max_c, int *n_months)
{
    int i, n = data->rows;
    int mi = 999999999, ma = -999999999, mc = 0;
    
    for (i = 0; i < n; i++) {
        long ts = apf_to_long(&MAT_AT(data, i, 0).re);
        int cid = (int)apf_to_long(&MAT_AT(data, i, 1).re);
        int m = ts_to_month(ts);
        if (m < mi) mi = m;
        if (m > ma) ma = m;
        if (cid > mc) mc = cid;
    }
    
    *min_m = mi;
    *max_m = ma;
    *max_c = mc;
    *n_months = ma - mi + 1;
    
    return (*n_months > 0 && *n_months <= MAX_MONTHS && 
            mc >= 0 && mc < MAX_CUSTOMERS);
}

/* Allocate and build per-customer-per-month revenue matrix */
static apf *build_revenue_matrix(const matrix_t *data, int min_m, int n_months,
                                  int max_c, int billing_col)
{
    int i, n = data->rows;
    size_t size = (size_t)n_months * (size_t)(max_c + 1);
    apf *rev = (apf *)calloc(size, sizeof(apf));
    
    if (!rev) return NULL;
    
    for (i = 0; i < (int)size; i++) apf_zero(&rev[i]);
    
    for (i = 0; i < n; i++) {
        long ts = apf_to_long(&MAT_AT(data, i, 0).re);
        int cid = (int)apf_to_long(&MAT_AT(data, i, 1).re);
        int m = ts_to_month(ts) - min_m;
        int idx = m * (max_c + 1) + cid;
        apf_add(&rev[idx], &rev[idx], &MAT_AT(data, i, billing_col).re);
    }
    
    return rev;
}

/* Build active customer bitmap */
static int *build_active_matrix(const matrix_t *data, int min_m, int n_months, int max_c)
{
    int i, n = data->rows;
    size_t size = (size_t)n_months * (size_t)(max_c + 1);
    int *active = (int *)calloc(size, sizeof(int));
    
    if (!active) return NULL;
    
    for (i = 0; i < n; i++) {
        long ts = apf_to_long(&MAT_AT(data, i, 0).re);
        int cid = (int)apf_to_long(&MAT_AT(data, i, 1).re);
        int m = ts_to_month(ts) - min_m;
        active[m * (max_c + 1) + cid] = 1;
    }
    
    return active;
}

/* Find first seen month for each customer */
static int *build_first_seen(const matrix_t *data, int max_c)
{
    int i, n = data->rows;
    int *first = (int *)malloc((size_t)(max_c + 1) * sizeof(int));
    
    if (!first) return NULL;
    
    for (i = 0; i <= max_c; i++) first[i] = 999999999;
    
    for (i = 0; i < n; i++) {
        long ts = apf_to_long(&MAT_AT(data, i, 0).re);
        int cid = (int)apf_to_long(&MAT_AT(data, i, 1).re);
        int m = ts_to_month(ts);
        if (m < first[cid]) first[cid] = m;
    }
    
    return first;
}

/* ============================================================
 * MRR: Monthly Recurring Revenue
 * Output: [month, mrr]
 * ============================================================ */
void mat_mrr(matrix_t *r, const matrix_t *data)
{
    int i, c, min_m, max_m, max_c, n_months;
    int billing_col = data->cols - 1;
    apf *rev;
    
    if (data->rows == 0 || data->cols < 3) {
        mat_zero(r, 0, 0);
        return;
    }
    
    if (!get_ranges(data, &min_m, &max_m, &max_c, &n_months)) {
        mat_zero(r, 0, 0);
        return;
    }
    
    rev = build_revenue_matrix(data, min_m, n_months, max_c, billing_col);
    if (!rev) {
        mat_zero(r, 0, 0);
        return;
    }
    
    mat_zero(r, n_months, 2);
    
    for (i = 0; i < n_months; i++) {
        apf sum;
        apf_zero(&sum);
        for (c = 0; c <= max_c; c++) {
            apf_add(&sum, &sum, &rev[i * (max_c + 1) + c]);
        }
        apf_from_int(&MAT_AT(r, i, 0).re, (int)month_to_ts(min_m + i));
        apf_copy(&MAT_AT(r, i, 1).re, &sum);
    }
    
    free(rev);
}

/* ============================================================
 * ACTIVE CUSTOMERS: Count per month
 * Output: [month, count]
 * ============================================================ */
void mat_customercount(matrix_t *r, const matrix_t *data)
{
    int i, c, min_m, max_m, max_c, n_months;
    int *active;
    
    if (data->rows == 0) {
        mat_zero(r, 0, 0);
        return;
    }
    
    if (!get_ranges(data, &min_m, &max_m, &max_c, &n_months)) {
        mat_zero(r, 0, 0);
        return;
    }
    
    active = build_active_matrix(data, min_m, n_months, max_c);
    if (!active) {
        mat_zero(r, 0, 0);
        return;
    }
    
    mat_zero(r, n_months, 2);
    
    for (i = 0; i < n_months; i++) {
        int count = 0;
        for (c = 0; c <= max_c; c++) {
            if (active[i * (max_c + 1) + c]) count++;
        }
        apf_from_int(&MAT_AT(r, i, 0).re, (int)month_to_ts(min_m + i));
        apf_from_int(&MAT_AT(r, i, 1).re, count);
    }
    
    free(active);
}

/* ============================================================
 * NEW CUSTOMERS: First appearance this month
 * Output: [month, new_count]
 * ============================================================ */
void mat_newcustomers(matrix_t *r, const matrix_t *data)
{
    int i, min_m, max_m, max_c, n_months;
    int *first, *new_count;
    
    if (data->rows == 0) {
        mat_zero(r, 0, 0);
        return;
    }
    
    if (!get_ranges(data, &min_m, &max_m, &max_c, &n_months)) {
        mat_zero(r, 0, 0);
        return;
    }
    
    first = build_first_seen(data, max_c);
    new_count = (int *)calloc((size_t)n_months, sizeof(int));
    
    if (!first || !new_count) {
        free(first); free(new_count);
        mat_zero(r, 0, 0);
        return;
    }
    
    for (i = 0; i <= max_c; i++) {
        if (first[i] != 999999999) {
            int m = first[i] - min_m;
            if (m >= 0 && m < n_months) new_count[m]++;
        }
    }
    
    mat_zero(r, n_months, 2);
    for (i = 0; i < n_months; i++) {
        apf_from_int(&MAT_AT(r, i, 0).re, (int)month_to_ts(min_m + i));
        apf_from_int(&MAT_AT(r, i, 1).re, new_count[i]);
    }
    
    free(first);
    free(new_count);
}

/* ============================================================
 * CHURNED CUSTOMERS: Active last month, not this month
 * Output: [month, churned_count]
 * ============================================================ */
void mat_churn(matrix_t *r, const matrix_t *data)
{
    int i, c, min_m, max_m, max_c, n_months;
    int *active, *churned;
    
    if (data->rows == 0) {
        mat_zero(r, 0, 0);
        return;
    }
    
    if (!get_ranges(data, &min_m, &max_m, &max_c, &n_months)) {
        mat_zero(r, 0, 0);
        return;
    }
    
    active = build_active_matrix(data, min_m, n_months, max_c);
    churned = (int *)calloc((size_t)n_months, sizeof(int));
    
    if (!active || !churned) {
        free(active); free(churned);
        mat_zero(r, 0, 0);
        return;
    }
    
    for (i = 1; i < n_months; i++) {
        for (c = 0; c <= max_c; c++) {
            int prev = active[(i-1) * (max_c + 1) + c];
            int curr = active[i * (max_c + 1) + c];
            if (prev && !curr) churned[i]++;
        }
    }
    
    mat_zero(r, n_months, 2);
    for (i = 0; i < n_months; i++) {
        apf_from_int(&MAT_AT(r, i, 0).re, (int)month_to_ts(min_m + i));
        apf_from_int(&MAT_AT(r, i, 1).re, churned[i]);
    }
    
    free(active);
    free(churned);
}

/* ============================================================
 * REACTIVATED: Not active last month, active this month, seen before
 * Output: [month, reactivated_count]
 * ============================================================ */
void mat_reactivated(matrix_t *r, const matrix_t *data)
{
    int i, c, min_m, max_m, max_c, n_months;
    int *active, *first, *reactivated;
    
    if (data->rows == 0) {
        mat_zero(r, 0, 0);
        return;
    }
    
    if (!get_ranges(data, &min_m, &max_m, &max_c, &n_months)) {
        mat_zero(r, 0, 0);
        return;
    }
    
    active = build_active_matrix(data, min_m, n_months, max_c);
    first = build_first_seen(data, max_c);
    reactivated = (int *)calloc((size_t)n_months, sizeof(int));
    
    if (!active || !first || !reactivated) {
        free(active); free(first); free(reactivated);
        mat_zero(r, 0, 0);
        return;
    }
    
    for (i = 1; i < n_months; i++) {
        for (c = 0; c <= max_c; c++) {
            int prev = active[(i-1) * (max_c + 1) + c];
            int curr = active[i * (max_c + 1) + c];
            int first_month = first[c] - min_m;
            /* Reactivated = not active last month, active now, and not new */
            if (!prev && curr && first_month < i) reactivated[i]++;
        }
    }
    
    mat_zero(r, n_months, 2);
    for (i = 0; i < n_months; i++) {
        apf_from_int(&MAT_AT(r, i, 0).re, (int)month_to_ts(min_m + i));
        apf_from_int(&MAT_AT(r, i, 1).re, reactivated[i]);
    }
    
    free(active);
    free(first);
    free(reactivated);
}

/* ============================================================
 * MRR BRIDGE: New, Expansion, Contraction, Churned MRR
 * Output: [month, new_mrr, expansion_mrr, contraction_mrr, churned_mrr, net_mrr]
 * ============================================================ */
void mat_mrrbridge(matrix_t *r, const matrix_t *data)
{
    int i, c, min_m, max_m, max_c, n_months;
    int billing_col = data->cols - 1;
    int *first;
    apf *rev;
    apf zero_val, new_mrr, exp_mrr, con_mrr, churn_mrr, net_mrr;
    apf prev_rev, curr_rev, diff;
    
    if (data->rows == 0 || data->cols < 3) {
        mat_zero(r, 0, 0);
        return;
    }
    
    if (!get_ranges(data, &min_m, &max_m, &max_c, &n_months)) {
        mat_zero(r, 0, 0);
        return;
    }
    
    rev = build_revenue_matrix(data, min_m, n_months, max_c, billing_col);
    first = build_first_seen(data, max_c);
    
    if (!rev || !first) {
        free(rev); free(first);
        mat_zero(r, 0, 0);
        return;
    }
    
    apf_zero(&zero_val);
    mat_zero(r, n_months, 6);
    
    for (i = 0; i < n_months; i++) {
        apf_zero(&new_mrr);
        apf_zero(&exp_mrr);
        apf_zero(&con_mrr);
        apf_zero(&churn_mrr);
        
        for (c = 0; c <= max_c; c++) {
            int idx = i * (max_c + 1) + c;
            int prev_idx = (i > 0) ? (i-1) * (max_c + 1) + c : -1;
            int first_month = first[c] - min_m;
            
            apf_copy(&curr_rev, &rev[idx]);
            if (prev_idx >= 0) {
                apf_copy(&prev_rev, &rev[prev_idx]);
            } else {
                apf_zero(&prev_rev);
            }
            
            if (first_month == i) {
                /* New customer this month */
                apf_add(&new_mrr, &new_mrr, &curr_rev);
            } else if (apf_cmp(&prev_rev, &zero_val) > 0) {
                /* Existing customer */
                if (apf_cmp(&curr_rev, &zero_val) == 0) {
                    /* Churned */
                    apf_add(&churn_mrr, &churn_mrr, &prev_rev);
                } else {
                    apf_sub(&diff, &curr_rev, &prev_rev);
                    if (apf_cmp(&diff, &zero_val) > 0) {
                        /* Expansion */
                        apf_add(&exp_mrr, &exp_mrr, &diff);
                    } else if (apf_cmp(&diff, &zero_val) < 0) {
                        /* Contraction */
                        apf_neg(&diff, &diff);
                        apf_add(&con_mrr, &con_mrr, &diff);
                    }
                }
            }
        }
        
        /* Net MRR = new + expansion - contraction - churn */
        apf_add(&net_mrr, &new_mrr, &exp_mrr);
        apf_sub(&net_mrr, &net_mrr, &con_mrr);
        apf_sub(&net_mrr, &net_mrr, &churn_mrr);
        
        apf_from_int(&MAT_AT(r, i, 0).re, (int)month_to_ts(min_m + i));
        apf_copy(&MAT_AT(r, i, 1).re, &new_mrr);
        apf_copy(&MAT_AT(r, i, 2).re, &exp_mrr);
        apf_copy(&MAT_AT(r, i, 3).re, &con_mrr);
        apf_copy(&MAT_AT(r, i, 4).re, &churn_mrr);
        apf_copy(&MAT_AT(r, i, 5).re, &net_mrr);
    }
    
    free(rev);
    free(first);
}

/* ============================================================
 * NRR: Net Revenue Retention (%)
 * (Revenue from last month's customers) / (Their revenue last month) * 100
 * Output: [month, nrr_percent]
 * ============================================================ */
void mat_nrr(matrix_t *r, const matrix_t *data)
{
    int i, c, min_m, max_m, max_c, n_months;
    int billing_col = data->cols - 1;
    apf *rev;
    apf prev_total, curr_total, nrr, hundred;
    
    if (data->rows == 0 || data->cols < 3) {
        mat_zero(r, 0, 0);
        return;
    }
    
    if (!get_ranges(data, &min_m, &max_m, &max_c, &n_months)) {
        mat_zero(r, 0, 0);
        return;
    }
    
    rev = build_revenue_matrix(data, min_m, n_months, max_c, billing_col);
    if (!rev) {
        mat_zero(r, 0, 0);
        return;
    }
    
    apf_from_int(&hundred, 100);
    mat_zero(r, n_months, 2);
    
    for (i = 0; i < n_months; i++) {
        apf_from_int(&MAT_AT(r, i, 0).re, (int)month_to_ts(min_m + i));
        
        if (i == 0) {
            apf_from_int(&MAT_AT(r, i, 1).re, 100);
            continue;
        }
        
        apf_zero(&prev_total);
        apf_zero(&curr_total);
        
        /* Only count customers who had revenue last month */
        for (c = 0; c <= max_c; c++) {
            int prev_idx = (i-1) * (max_c + 1) + c;
            int curr_idx = i * (max_c + 1) + c;
            
            if (apf_cmp_int(&rev[prev_idx], 0) > 0) {
                apf_add(&prev_total, &prev_total, &rev[prev_idx]);
                apf_add(&curr_total, &curr_total, &rev[curr_idx]);
            }
        }
        
        if (apf_cmp_int(&prev_total, 0) > 0) {
            apf_div(&nrr, &curr_total, &prev_total);
            apf_mul(&MAT_AT(r, i, 1).re, &nrr, &hundred);
        } else {
            apf_from_int(&MAT_AT(r, i, 1).re, 100);
        }
    }
    
    free(rev);
}

/* ============================================================
 * GRR: Gross Revenue Retention (%) - ignores expansion
 * min(current, previous) / previous * 100 for existing customers
 * Output: [month, grr_percent]
 * ============================================================ */
void mat_grr(matrix_t *r, const matrix_t *data)
{
    int i, c, min_m, max_m, max_c, n_months;
    int billing_col = data->cols - 1;
    apf *rev;
    apf prev_total, retained_total, grr, hundred, prev_rev, curr_rev, retained;
    
    if (data->rows == 0 || data->cols < 3) {
        mat_zero(r, 0, 0);
        return;
    }
    
    if (!get_ranges(data, &min_m, &max_m, &max_c, &n_months)) {
        mat_zero(r, 0, 0);
        return;
    }
    
    rev = build_revenue_matrix(data, min_m, n_months, max_c, billing_col);
    if (!rev) {
        mat_zero(r, 0, 0);
        return;
    }
    
    apf_from_int(&hundred, 100);
    mat_zero(r, n_months, 2);
    
    for (i = 0; i < n_months; i++) {
        apf_from_int(&MAT_AT(r, i, 0).re, (int)month_to_ts(min_m + i));
        
        if (i == 0) {
            apf_from_int(&MAT_AT(r, i, 1).re, 100);
            continue;
        }
        
        apf_zero(&prev_total);
        apf_zero(&retained_total);
        
        for (c = 0; c <= max_c; c++) {
            int prev_idx = (i-1) * (max_c + 1) + c;
            int curr_idx = i * (max_c + 1) + c;
            
            apf_copy(&prev_rev, &rev[prev_idx]);
            apf_copy(&curr_rev, &rev[curr_idx]);
            
            if (apf_cmp_int(&prev_rev, 0) > 0) {
                apf_add(&prev_total, &prev_total, &prev_rev);
                /* Retained = min(curr, prev) - no expansion credit */
                if (apf_cmp(&curr_rev, &prev_rev) < 0) {
                    apf_copy(&retained, &curr_rev);
                } else {
                    apf_copy(&retained, &prev_rev);
                }
                apf_add(&retained_total, &retained_total, &retained);
            }
        }
        
        if (apf_cmp_int(&prev_total, 0) > 0) {
            apf_div(&grr, &retained_total, &prev_total);
            apf_mul(&MAT_AT(r, i, 1).re, &grr, &hundred);
        } else {
            apf_from_int(&MAT_AT(r, i, 1).re, 100);
        }
    }
    
    free(rev);
}

/* ============================================================
 * ARPU: Average Revenue Per User
 * Output: [month, arpu]
 * ============================================================ */
void mat_arpu(matrix_t *r, const matrix_t *data)
{
    matrix_t mrr_data, cust_data;
    int i, n;
    
    mat_mrr(&mrr_data, data);
    mat_customercount(&cust_data, data);
    
    n = mrr_data.rows;
    if (n == 0) {
        mat_zero(r, 0, 0);
        return;
    }
    
    mat_zero(r, n, 2);
    
    for (i = 0; i < n; i++) {
        apf_copy(&MAT_AT(r, i, 0).re, &MAT_AT(&mrr_data, i, 0).re);
        if (apf_cmp_int(&MAT_AT(&cust_data, i, 1).re, 0) > 0) {
            apf_div(&MAT_AT(r, i, 1).re, &MAT_AT(&mrr_data, i, 1).re, 
                    &MAT_AT(&cust_data, i, 1).re);
        } else {
            apf_zero(&MAT_AT(r, i, 1).re);
        }
    }
}

/* ============================================================
 * CHURN RATE: Logo churn rate (%)
 * churned / previous_active * 100
 * Output: [month, churn_rate_percent]
 * ============================================================ */
void mat_churnrate(matrix_t *r, const matrix_t *data)
{
    matrix_t churn_data, cust_data;
    int i, n;
    apf hundred, rate;
    
    mat_churn(&churn_data, data);
    mat_customercount(&cust_data, data);
    
    n = churn_data.rows;
    if (n == 0) {
        mat_zero(r, 0, 0);
        return;
    }
    
    apf_from_int(&hundred, 100);
    mat_zero(r, n, 2);
    
    for (i = 0; i < n; i++) {
        apf_copy(&MAT_AT(r, i, 0).re, &MAT_AT(&churn_data, i, 0).re);
        if (i > 0 && apf_cmp_int(&MAT_AT(&cust_data, i-1, 1).re, 0) > 0) {
            apf_div(&rate, &MAT_AT(&churn_data, i, 1).re, 
                    &MAT_AT(&cust_data, i-1, 1).re);
            apf_mul(&MAT_AT(r, i, 1).re, &rate, &hundred);
        } else {
            apf_zero(&MAT_AT(r, i, 1).re);
        }
    }
}

/* ============================================================
 * RETENTION MATRIX: Cohort retention table
 * Row = cohort month, Col = months since start
 * Value = % of cohort still active
 * Output: NxN matrix
 * ============================================================ */
void mat_retention(matrix_t *r, const matrix_t *data)
{
    int i, m, c, min_m, max_m, max_c, n_months;
    int *active, *first, *cohort_size;
    apf hundred, pct;
    
    if (data->rows == 0) {
        mat_zero(r, 0, 0);
        return;
    }
    
    if (!get_ranges(data, &min_m, &max_m, &max_c, &n_months)) {
        mat_zero(r, 0, 0);
        return;
    }
    
    active = build_active_matrix(data, min_m, n_months, max_c);
    first = build_first_seen(data, max_c);
    cohort_size = (int *)calloc((size_t)n_months, sizeof(int));
    
    if (!active || !first || !cohort_size) {
        free(active); free(first); free(cohort_size);
        mat_zero(r, 0, 0);
        return;
    }
    
    /* Count cohort sizes */
    for (c = 0; c <= max_c; c++) {
        if (first[c] != 999999999) {
            int cohort = first[c] - min_m;
            if (cohort >= 0 && cohort < n_months) cohort_size[cohort]++;
        }
    }
    
    apf_from_int(&hundred, 100);
    mat_zero(r, n_months, n_months);
    
    for (i = 0; i < n_months; i++) {
        if (cohort_size[i] == 0) continue;
        
        for (m = 0; m < n_months - i; m++) {
            int month_idx = i + m;
            int retained = 0;
            
            for (c = 0; c <= max_c; c++) {
                if (first[c] - min_m == i && active[month_idx * (max_c + 1) + c]) {
                    retained++;
                }
            }
            
            apf_from_int(&pct, retained * 100);
            apf_from_int(&MAT_AT(r, i, m).re, cohort_size[i]);
            apf_div(&MAT_AT(r, i, m).re, &pct, &MAT_AT(r, i, m).re);
        }
    }
    
    free(active);
    free(first);
    free(cohort_size);
}

/* ============================================================
 * TOP CUSTOMERS: Top N by total revenue
 * Output: [customerid, total_revenue] sorted descending
 * ============================================================ */
void mat_topcustomers(matrix_t *r, const matrix_t *data, int n_top)
{
    int i, j, c, n, min_m, max_m, max_c, n_months;
    int billing_col = data->cols - 1;
    apf *rev, *totals;
    int *sorted_idx;
    
    int tmp_idx;
    
    if (data->rows == 0 || data->cols < 3 || n_top <= 0) {
        mat_zero(r, 0, 0);
        return;
    }
    
    if (!get_ranges(data, &min_m, &max_m, &max_c, &n_months)) {
        mat_zero(r, 0, 0);
        return;
    }
    
    rev = build_revenue_matrix(data, min_m, n_months, max_c, billing_col);
    totals = (apf *)malloc((size_t)(max_c + 1) * sizeof(apf));
    sorted_idx = (int *)malloc((size_t)(max_c + 1) * sizeof(int));
    
    if (!rev || !totals || !sorted_idx) {
        free(rev); free(totals); free(sorted_idx);
        mat_zero(r, 0, 0);
        return;
    }
    
    /* Sum revenue per customer */
    for (c = 0; c <= max_c; c++) {
        apf_zero(&totals[c]);
        sorted_idx[c] = c;
        for (i = 0; i < n_months; i++) {
            apf_add(&totals[c], &totals[c], &rev[i * (max_c + 1) + c]);
        }
    }
    
    /* Sort descending by total */
    n = max_c + 1;
    for (i = 0; i < n - 1; i++) {
        for (j = 0; j < n - i - 1; j++) {
            if (apf_cmp(&totals[sorted_idx[j]], &totals[sorted_idx[j+1]]) < 0) {
                tmp_idx = sorted_idx[j];
                sorted_idx[j] = sorted_idx[j+1];
                sorted_idx[j+1] = tmp_idx;
            }
        }
    }
    
    /* Output top N */
    if (n_top > n) n_top = n;
    
    /* Only include customers with revenue */
    for (i = 0; i < n_top; i++) {
        if (apf_cmp_int(&totals[sorted_idx[i]], 0) <= 0) {
            n_top = i;
            break;
        }
    }
    
    mat_zero(r, n_top, 2);
    for (i = 0; i < n_top; i++) {
        apf_from_int(&MAT_AT(r, i, 0).re, sorted_idx[i]);
        apf_copy(&MAT_AT(r, i, 1).re, &totals[sorted_idx[i]]);
    }
    
    free(rev);
    free(totals);
    free(sorted_idx);
}

/* ============================================================
 * CUSTOMER TENURE: Months active for each customer
 * Output: [customerid, tenure_months, first_month, last_month]
 * ============================================================ */
void mat_tenure(matrix_t *r, const matrix_t *data)
{
    int i, c, n, min_m, max_m, max_c, n_months;
    int *first, *last, out_count;
    
    if (data->rows == 0) {
        mat_zero(r, 0, 0);
        return;
    }
    
    if (!get_ranges(data, &min_m, &max_m, &max_c, &n_months)) {
        mat_zero(r, 0, 0);
        return;
    }
    
    first = build_first_seen(data, max_c);
    last = (int *)malloc((size_t)(max_c + 1) * sizeof(int));
    
    if (!first || !last) {
        free(first); free(last);
        mat_zero(r, 0, 0);
        return;
    }
    
    /* Initialize last seen */
    for (c = 0; c <= max_c; c++) last[c] = -1;
    
    /* Find last seen */
    n = data->rows;
    for (i = 0; i < n; i++) {
        long ts = apf_to_long(&MAT_AT(data, i, 0).re);
        int cid = (int)apf_to_long(&MAT_AT(data, i, 1).re);
        int m = ts_to_month(ts);
        if (m > last[cid]) last[cid] = m;
    }
    
    /* Count customers with activity */
    out_count = 0;
    for (c = 0; c <= max_c; c++) {
        if (first[c] != 999999999) out_count++;
    }
    
    mat_zero(r, out_count, 4);
    i = 0;
    for (c = 0; c <= max_c; c++) {
        if (first[c] != 999999999) {
            int tenure = last[c] - first[c] + 1;
            apf_from_int(&MAT_AT(r, i, 0).re, c);
            apf_from_int(&MAT_AT(r, i, 1).re, tenure);
            apf_from_int(&MAT_AT(r, i, 2).re, (int)month_to_ts(first[c]));
            apf_from_int(&MAT_AT(r, i, 3).re, (int)month_to_ts(last[c]));
            i++;
        }
    }
    
    free(first);
    free(last);
}

/* ============================================================
 * RETIME: Aggregate hourly data to monthly
 * Input: [timestamp, customerid, col3, col4, ...]
 * Output: [month_ts, customerid, sum(col3), sum(col4), ...]
 * ============================================================ */
void mat_retime(matrix_t *r, const matrix_t *data, int interval)
{
    int i, j, c, n, ncols, min_b, max_b, max_c, n_buckets;
    int *bucket;
    size_t alloc_size;
    apf *sums;
    int out_count, k;
    
    n = data->rows;
    ncols = data->cols;
    
    if (n == 0 || ncols < 2) {
        mat_zero(r, 0, 0);
        return;
    }
    
    bucket = (int *)malloc((size_t)n * sizeof(int));
    if (!bucket) {
        mat_zero(r, 0, 0);
        return;
    }
    
    /* Compute buckets and find ranges */
    min_b = 999999999;
    max_b = -999999999;
    max_c = 0;
    
    for (i = 0; i < n; i++) {
        long ts = apf_to_long(&MAT_AT(data, i, 0).re);
        int cid = (int)apf_to_long(&MAT_AT(data, i, 1).re);
        int b;
        
        switch (interval) {
            case 0: b = (int)(ts / 3600); break;      /* hourly */
            case 1: b = (int)(ts / 86400); break;     /* daily */
            case 2: b = (int)(ts / 604800); break;    /* weekly */
            case 3: b = ts_to_month(ts); break;       /* monthly */
            case 4: b = ts_to_month(ts) / 3; break;   /* quarterly */
            case 5: b = ts_to_month(ts) / 12; break;  /* yearly */
            default: b = ts_to_month(ts);
        }
        
        bucket[i] = b;
        if (b < min_b) min_b = b;
        if (b > max_b) max_b = b;
        if (cid > max_c) max_c = cid;
    }
    
    n_buckets = max_b - min_b + 1;
    
    if (n_buckets <= 0 || n_buckets > MAX_MONTHS || max_c > MAX_CUSTOMERS) {
        free(bucket);
        mat_zero(r, 0, 0);
        return;
    }
    
    /* Allocate sums array */
    alloc_size = (size_t)n_buckets * (size_t)(max_c + 1) * (size_t)ncols;
    sums = (apf *)calloc(alloc_size, sizeof(apf));
    
    if (!sums) {
        free(bucket);
        mat_zero(r, 0, 0);
        return;
    }
    
    for (i = 0; i < (int)alloc_size; i++) apf_zero(&sums[i]);
    
    /* Aggregate */
    for (i = 0; i < n; i++) {
        int b = bucket[i] - min_b;
        int cid = (int)apf_to_long(&MAT_AT(data, i, 1).re);
        int base = (b * (max_c + 1) + cid) * ncols;
        
        for (j = 0; j < ncols; j++) {
            apf_add(&sums[base + j], &sums[base + j], &MAT_AT(data, i, j).re);
        }
    }
    
    /* Count non-empty rows */
    out_count = 0;
    for (i = 0; i < n_buckets; i++) {
        for (c = 0; c <= max_c; c++) {
            int base = (i * (max_c + 1) + c) * ncols;
            /* Check if any column has data (check billing column) */
            if (apf_cmp_int(&sums[base + ncols - 1], 0) != 0) {
                out_count++;
            }
        }
    }
    
    mat_zero(r, out_count, ncols);
    
    k = 0;
    for (i = 0; i < n_buckets; i++) {
        for (c = 0; c <= max_c; c++) {
            int base = (i * (max_c + 1) + c) * ncols;
            if (apf_cmp_int(&sums[base + ncols - 1], 0) != 0) {
                long ts;
                
                switch (interval) {
                    case 0: ts = (long)(min_b + i) * 3600; break;
                    case 1: ts = (long)(min_b + i) * 86400; break;
                    case 2: ts = (long)(min_b + i) * 604800; break;
                    case 3: ts = month_to_ts(min_b + i); break;
                    case 4: ts = month_to_ts((min_b + i) * 3); break;
                    case 5: ts = month_to_ts((min_b + i) * 12); break;
                    default: ts = month_to_ts(min_b + i);
                }
                
                apf_from_int(&MAT_AT(r, k, 0).re, (int)ts);
                apf_from_int(&MAT_AT(r, k, 1).re, c);
                for (j = 2; j < ncols; j++) {
                    apf_copy(&MAT_AT(r, k, j).re, &sums[base + j]);
                }
                k++;
            }
        }
    }
    
    free(bucket);
    free(sums);
}

/* ========== Revenue Concentration ========== */

/* toprevenue - Get top N customers by revenue
 * Returns: [customerid, total_revenue, pct_of_total] sorted descending */
void mat_toprevenue(matrix_t *result, const matrix_t *data, int top_n)
{
    int i, j, n;
    long *cust_ids;
    double *cust_rev;
    int num_custs = 0;
    double total_rev = 0;
    
    if (!data || data->cols < 4) {
        mat_zero(result, 0, 0);
        return;
    }
    
    n = data->rows;
    cust_ids = (long *)malloc(n * sizeof(long));
    cust_rev = (double *)malloc(n * sizeof(double));
    if (!cust_ids || !cust_rev) {
        free(cust_ids); free(cust_rev);
        mat_zero(result, 0, 0);
        return;
    }
    
    /* Aggregate revenue by customer */
    for (i = 0; i < n; i++) {
        long cid = apf_to_long(&MAT_AT(data, i, 1).re);
        double rev = apf_to_double(&MAT_AT(data, i, 3).re);
        int found = -1;
        
        for (j = 0; j < num_custs; j++) {
            if (cust_ids[j] == cid) { found = j; break; }
        }
        
        if (found < 0) {
            found = num_custs++;
            cust_ids[found] = cid;
            cust_rev[found] = 0;
        }
        cust_rev[found] += rev;
        total_rev += rev;
    }
    
    /* Sort by revenue descending (simple bubble sort) */
    for (i = 0; i < num_custs - 1; i++) {
        for (j = i + 1; j < num_custs; j++) {
            if (cust_rev[j] > cust_rev[i]) {
                double tmp_rev = cust_rev[i];
                long tmp_id = cust_ids[i];
                cust_rev[i] = cust_rev[j];
                cust_ids[i] = cust_ids[j];
                cust_rev[j] = tmp_rev;
                cust_ids[j] = tmp_id;
            }
        }
    }
    
    /* Return top N */
    if (top_n <= 0 || top_n > num_custs) top_n = num_custs;
    mat_zero(result, top_n, 3);
    
    for (i = 0; i < top_n; i++) {
        apf_from_int(&MAT_AT(result, i, 0).re, (int)cust_ids[i]);
        apf_from_double(&MAT_AT(result, i, 1).re, cust_rev[i]);
        apf_from_double(&MAT_AT(result, i, 2).re, (cust_rev[i] / total_rev) * 100);
    }
    
    free(cust_ids);
    free(cust_rev);
}

/* concentration - Revenue concentration metrics
 * Returns: [top1%, top5%, top10%, top20%, herfindahl_index] */
void mat_concentration(matrix_t *result, const matrix_t *data)
{
    int i, j, n;
    long *cust_ids;
    double *cust_rev;
    int num_custs = 0;
    double total_rev = 0;
    double top1 = 0, top5 = 0, top10 = 0, top20 = 0;
    double hhi = 0;  /* Herfindahl-Hirschman Index */
    
    if (!data || data->cols < 4) {
        mat_zero(result, 1, 5);
        return;
    }
    
    n = data->rows;
    cust_ids = (long *)malloc(n * sizeof(long));
    cust_rev = (double *)malloc(n * sizeof(double));
    if (!cust_ids || !cust_rev) {
        free(cust_ids); free(cust_rev);
        mat_zero(result, 1, 5);
        return;
    }
    
    /* Aggregate revenue by customer */
    for (i = 0; i < n; i++) {
        long cid = apf_to_long(&MAT_AT(data, i, 1).re);
        double rev = apf_to_double(&MAT_AT(data, i, 3).re);
        int found = -1;
        
        for (j = 0; j < num_custs; j++) {
            if (cust_ids[j] == cid) { found = j; break; }
        }
        
        if (found < 0) {
            found = num_custs++;
            cust_ids[found] = cid;
            cust_rev[found] = 0;
        }
        cust_rev[found] += rev;
        total_rev += rev;
    }
    
    /* Sort by revenue descending */
    for (i = 0; i < num_custs - 1; i++) {
        for (j = i + 1; j < num_custs; j++) {
            if (cust_rev[j] > cust_rev[i]) {
                double tmp = cust_rev[i];
                cust_rev[i] = cust_rev[j];
                cust_rev[j] = tmp;
            }
        }
    }
    
    /* Calculate concentration */
    for (i = 0; i < num_custs; i++) {
        double pct = cust_rev[i] / total_rev;
        hhi += pct * pct;  /* Sum of squared market shares */
        
        if (i < 1) top1 += pct;
        if (i < 5) top5 += pct;
        if (i < 10) top10 += pct;
        if (i < (num_custs / 5 + 1)) top20 += pct;  /* Top 20% of customers */
    }
    
    mat_zero(result, 1, 5);
    apf_from_double(&MAT_AT(result, 0, 0).re, top1 * 100);
    apf_from_double(&MAT_AT(result, 0, 1).re, top5 * 100);
    apf_from_double(&MAT_AT(result, 0, 2).re, top10 * 100);
    apf_from_double(&MAT_AT(result, 0, 3).re, top20 * 100);
    apf_from_double(&MAT_AT(result, 0, 4).re, hhi * 10000);  /* HHI scaled 0-10000 */
    
    free(cust_ids);
    free(cust_rev);
}

/* revchurn - Revenue churn rate per month
 * Returns: [month, churned_mrr, start_mrr, rev_churn_rate%] */
void mat_revchurn(matrix_t *result, const matrix_t *data)
{
    matrix_t bridge;
    matrix_t mrr_data;
    int i, n;
    
    mat_mrrbridge(&bridge, data);
    mat_mrr(&mrr_data, data);
    
    n = bridge.rows;
    mat_zero(result, n, 4);
    
    for (i = 0; i < n; i++) {
        double month_ts = apf_to_double(&MAT_AT(&bridge, i, 0).re);
        double churned = apf_to_double(&MAT_AT(&bridge, i, 4).re);
        double prev_mrr = (i > 0) ? apf_to_double(&MAT_AT(&mrr_data, i-1, 1).re) : 0;
        double churn_rate = (prev_mrr > 0) ? (churned / prev_mrr) * 100 : 0;
        
        apf_from_double(&MAT_AT(result, i, 0).re, month_ts);
        apf_from_double(&MAT_AT(result, i, 1).re, churned);
        apf_from_double(&MAT_AT(result, i, 2).re, prev_mrr);
        apf_from_double(&MAT_AT(result, i, 3).re, churn_rate);
    }
}

/* netchurn - Net churn rate (including contraction, excluding expansion)
 * Returns: [month, net_churn_rate%] where net_churn = (churn + contraction - expansion) / start_mrr */
void mat_netchurn(matrix_t *result, const matrix_t *data)
{
    matrix_t bridge;
    matrix_t mrr_data;
    int i, n;
    
    mat_mrrbridge(&bridge, data);
    mat_mrr(&mrr_data, data);
    
    n = bridge.rows;
    mat_zero(result, n, 2);
    
    for (i = 0; i < n; i++) {
        double month_ts = apf_to_double(&MAT_AT(&bridge, i, 0).re);
        double expansion = apf_to_double(&MAT_AT(&bridge, i, 2).re);
        double contraction = apf_to_double(&MAT_AT(&bridge, i, 3).re);
        double churned = apf_to_double(&MAT_AT(&bridge, i, 4).re);
        double prev_mrr = (i > 0) ? apf_to_double(&MAT_AT(&mrr_data, i-1, 1).re) : 0;
        double net_churn = (prev_mrr > 0) ? ((churned + contraction - expansion) / prev_mrr) * 100 : 0;
        
        apf_from_double(&MAT_AT(result, i, 0).re, month_ts);
        apf_from_double(&MAT_AT(result, i, 1).re, net_churn);
    }
}

/* quickratio - SaaS Quick Ratio = (new + expansion) / (churn + contraction)
 * Higher is better. >4 is excellent, >2 is good, <1 is shrinking */
void mat_quickratio(matrix_t *result, const matrix_t *data)
{
    matrix_t bridge;
    int i, n;
    
    mat_mrrbridge(&bridge, data);
    n = bridge.rows;
    mat_zero(result, n, 2);
    
    for (i = 0; i < n; i++) {
        double month_ts = apf_to_double(&MAT_AT(&bridge, i, 0).re);
        double new_mrr = apf_to_double(&MAT_AT(&bridge, i, 1).re);
        double expansion = apf_to_double(&MAT_AT(&bridge, i, 2).re);
        double contraction = apf_to_double(&MAT_AT(&bridge, i, 3).re);
        double churned = apf_to_double(&MAT_AT(&bridge, i, 4).re);
        
        double growth = new_mrr + expansion;
        double loss = contraction + churned;
        double qr = (loss > 0) ? growth / loss : (growth > 0 ? 999 : 0);
        
        apf_from_double(&MAT_AT(result, i, 0).re, month_ts);
        apf_from_double(&MAT_AT(result, i, 1).re, qr);
    }
}

/* ltv - Customer Lifetime Value per customer
 * Returns: [customerid, total_revenue, tenure_months, monthly_avg, estimated_ltv] */
void mat_ltv(matrix_t *result, const matrix_t *data)
{
    matrix_t ten;
    int i, j, n;
    long *cust_ids;
    double *cust_rev;
    int num_custs = 0;
    double avg_churn = 0.05;  /* Assume 5% if can't calculate */
    
    if (!data || data->cols < 4) {
        mat_zero(result, 0, 0);
        return;
    }
    
    n = data->rows;
    cust_ids = (long *)malloc(n * sizeof(long));
    cust_rev = (double *)malloc(n * sizeof(double));
    if (!cust_ids || !cust_rev) {
        free(cust_ids); free(cust_rev);
        mat_zero(result, 0, 0);
        return;
    }
    
    /* Get tenure data */
    mat_tenure(&ten, data);
    
    /* Aggregate revenue by customer */
    for (i = 0; i < n; i++) {
        long cid = apf_to_long(&MAT_AT(data, i, 1).re);
        double rev = apf_to_double(&MAT_AT(data, i, 3).re);
        int found = -1;
        
        for (j = 0; j < num_custs; j++) {
            if (cust_ids[j] == cid) { found = j; break; }
        }
        
        if (found < 0) {
            found = num_custs++;
            cust_ids[found] = cid;
            cust_rev[found] = 0;
        }
        cust_rev[found] += rev;
    }
    
    /* Build result with LTV estimates */
    mat_zero(result, num_custs, 5);
    
    for (i = 0; i < num_custs; i++) {
        int tenure = 1;
        double monthly_avg, ltv;
        
        /* Find tenure for this customer */
        for (j = 0; j < ten.rows; j++) {
            if ((int)apf_to_double(&MAT_AT(&ten, j, 0).re) == (int)cust_ids[i]) {
                tenure = (int)apf_to_double(&MAT_AT(&ten, j, 1).re);
                if (tenure < 1) tenure = 1;
                break;
            }
        }
        
        monthly_avg = cust_rev[i] / tenure;
        ltv = monthly_avg / avg_churn;  /* Simple LTV = ARPU / churn */
        
        apf_from_int(&MAT_AT(result, i, 0).re, (int)cust_ids[i]);
        apf_from_double(&MAT_AT(result, i, 1).re, cust_rev[i]);
        apf_from_int(&MAT_AT(result, i, 2).re, tenure);
        apf_from_double(&MAT_AT(result, i, 3).re, monthly_avg);
        apf_from_double(&MAT_AT(result, i, 4).re, ltv);
    }
    
    free(cust_ids);
    free(cust_rev);
}

/* ========== Data Generation ========== */

/* Generate synthetic hourly billing data
 * Creates datasets/hourlybilling.csv (or datasets\hourlybilling.csv on Windows) */
void saas_generate_billing_data(int target_rows)
{
    FILE *fp;
    int i, month, hour;
    int rows_written = 0;
    long start_ts = 1704067200;  /* Jan 1, 2024 */
    int hours_per_month = 730;
    int num_months = 24;
    int total_customers = 113;
    int hours_per_row;
    
    /* Customer state */
    struct {
        int id;
        int segment;  /* 0=enterprise, 1=midmarket, 2=smb, 3=startup */
        double base_size;
        double base_rate;
        int start_month;
        double churn_prob;
        double expand_prob;
        int is_active;
        double current_size;
    } customers[120];
    
    /* Seed random */
    srand((unsigned int)time(NULL));
    
#ifdef _WIN32
    fp = fopen("datasets\\hourlybilling.csv", "w");
#else
    fp = fopen("datasets/hourlybilling.csv", "w");
#endif
    
    if (!fp) {
        /* Try creating directory first */
        int dummy;
#ifdef _WIN32
        dummy = system("mkdir datasets 2>nul");
        (void)dummy;
        fp = fopen("datasets\\hourlybilling.csv", "w");
#else
        dummy = system("mkdir -p datasets 2>/dev/null");
        (void)dummy;
        fp = fopen("datasets/hourlybilling.csv", "w");
#endif
    }
    
    if (!fp) {
        printf("Error: Cannot create hourlybilling.csv\n");
        return;
    }
    
    /* Initialize customers */
    for (i = 0; i < total_customers; i++) {
        customers[i].id = i + 1;
        
        if (i < 8) {  /* Enterprise */
            customers[i].segment = 0;
            customers[i].base_size = 50 + (rand() % 150);
            customers[i].base_rate = 0.08 + (rand() % 40) / 1000.0;
            customers[i].start_month = 1 + (rand() % 3);
            customers[i].churn_prob = 0.01;
            customers[i].expand_prob = 0.15;
        } else if (i < 33) {  /* Mid-market */
            customers[i].segment = 1;
            customers[i].base_size = 10 + (rand() % 40);
            customers[i].base_rate = 0.10 + (rand() % 50) / 1000.0;
            customers[i].start_month = 1 + (rand() % 6);
            customers[i].churn_prob = 0.03;
            customers[i].expand_prob = 0.12;
        } else if (i < 93) {  /* SMB */
            customers[i].segment = 2;
            customers[i].base_size = 1 + (rand() % 9);
            customers[i].base_rate = 0.12 + (rand() % 80) / 1000.0;
            customers[i].start_month = 1 + (rand() % 12);
            customers[i].churn_prob = 0.06;
            customers[i].expand_prob = 0.08;
        } else {  /* Startup */
            customers[i].segment = 3;
            customers[i].base_size = 0.5 + (rand() % 50) / 10.0;
            customers[i].base_rate = 0.15 + (rand() % 100) / 1000.0;
            customers[i].start_month = 1 + (rand() % 18);
            customers[i].churn_prob = 0.08;
            customers[i].expand_prob = 0.25;
        }
        
        customers[i].is_active = 0;
        customers[i].current_size = customers[i].base_size;
    }
    
    /* Calculate skip rate */
    hours_per_row = (num_months * hours_per_month * 70) / target_rows;
    if (hours_per_row < 1) hours_per_row = 1;
    
    /* Write header */
    fprintf(fp, "timestamp,customerid,sizetb,billing\n");
    
    /* Generate data */
    for (month = 1; month <= num_months && rows_written < target_rows; month++) {
        long month_ts = start_ts + (long)(month - 1) * hours_per_month * 3600;
        
        /* Update customer states */
        for (i = 0; i < total_customers; i++) {
            /* Activate new customers */
            if (!customers[i].is_active && month >= customers[i].start_month) {
                customers[i].is_active = 1;
                customers[i].current_size = customers[i].base_size * (0.8 + (rand() % 40) / 100.0);
            }
            
            /* Check for churn */
            if (customers[i].is_active && (rand() % 100) < (int)(customers[i].churn_prob * 100)) {
                customers[i].is_active = 0;
            }
            
            /* Check for expansion/contraction */
            if (customers[i].is_active) {
                if ((rand() % 100) < (int)(customers[i].expand_prob * 100)) {
                    customers[i].current_size *= 1.1 + (rand() % 40) / 100.0;
                } else if ((rand() % 100) < 10) {
                    customers[i].current_size *= 0.7 + (rand() % 20) / 100.0;
                }
                /* Monthly noise */
                customers[i].current_size *= 0.95 + (rand() % 10) / 100.0;
            }
            
            /* Reactivation chance */
            if (!customers[i].is_active && month > customers[i].start_month && (rand() % 100) < 5) {
                customers[i].is_active = 1;
                customers[i].current_size = customers[i].base_size * 0.6;
            }
        }
        
        /* Generate hourly records */
        for (hour = 0; hour < hours_per_month && rows_written < target_rows; hour += hours_per_row) {
            long ts = month_ts + (long)hour * 3600;
            
            for (i = 0; i < total_customers && rows_written < target_rows; i++) {
                if (customers[i].is_active) {
                    double size = customers[i].current_size * (0.9 + (rand() % 20) / 100.0);
                    double billing = size * customers[i].base_rate;
                    
                    fprintf(fp, "%ld,%d,%.4f,%.2f\n", ts, customers[i].id, size, billing);
                    rows_written++;
                }
            }
        }
    }
    
    fclose(fp);
    printf("Generated datasets%chourlybilling.csv with %d rows\n", 
#ifdef _WIN32
           '\\',
#else
           '/',
#endif
           rows_written);
}
